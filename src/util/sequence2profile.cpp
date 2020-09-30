// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte


#include "CSProfile.h"
#include "MathUtil.h"
#include "DBReader.h"
#include "Parameters.h"
#include "DBWriter.h"

#include <string>
#include <PSSMMasker.h>


#ifdef OPENMP
#include <omp.h>
#endif

int sequence2profile(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, 0.0);

    DBReader<unsigned int> sequenceDb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    sequenceDb.open(DBReader<unsigned int>::NOSORT);
    DBWriter resultDbw(par.db2.c_str(), par.db2Index.c_str(), par.threads,  par.compressed, Parameters::DBTYPE_HMM_PROFILE);
    resultDbw.open();
    Debug::Progress progress( sequenceDb.getSize());

#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, sequenceDb.getDbtype(), &subMat, 0, false, false);
        CSProfile ps(par.maxSeqLen);
        char * data = new char[sequenceDb.getMaxSeqLen() * Sequence::PROFILE_READIN_SIZE];
        ProbabilityMatrix probMatrix(subMat);
        PSSMMasker masker(sequenceDb.getMaxSeqLen(), probMatrix, subMat);

#pragma omp for schedule(static)
        for (size_t id = 0; id < sequenceDb.getSize(); id++) {
            int thread_idx = 0;
            progress.updateProgress();
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char *seqData     = sequenceDb.getData(id, thread_idx);
            unsigned int queryKey = sequenceDb.getDbKey(id);
            unsigned int seqLen = sequenceDb.getSeqLen(id);

            seq.mapSequence(id, queryKey, seqData, seqLen);
            float * prob =  ps.computeProfile(&seq, par.neff, par.tau);
//            if (par.maskProfile == true) {
//                masker.mask(seq, pssmRes);
//            }
//            pssmRes.toBuffer(centerSequence, subMat, result);
            size_t idx = 0;
            for (int i = 0; i < seq.L; i++) {
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE;aa++)
                {
                    data[idx++] = prob[i*Sequence::PROFILE_READIN_SIZE + aa];
                    //std::cout<<"\t"<<(int)data[idx-1];
                }
                //std::cout<<std::endl;
                data[idx++] = static_cast<unsigned char>(seq.numSequence[i]); // query
                data[idx++] = static_cast<unsigned char>(seq.numSequence[i]); // consensus
                data[idx++] = MathUtil::convertNeffToChar(prob[i*Sequence::PROFILE_READIN_SIZE + Sequence::PROFILE_AA_SIZE]);
            }
            resultDbw.writeData(data, idx, queryKey, thread_idx);
        }
        delete [] data;
    }
    sequenceDb.close();
    resultDbw.close();
    return 0;
}
