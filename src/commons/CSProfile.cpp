//
// Created by mad on 11/29/16.
//

#include "CSProfile.h"
#include "K4000.crf.h"


ContextLibrary::ContextLibrary(){
    std::string lib((const char *)K4000_crf,  K4000_crf_len);
    read(lib);
}
ContextLibrary *ContextLibrary::contextLibrary = NULL;

ContextLibrary::~ContextLibrary() {
    for (size_t k = 0; k < libSize; ++k) {
        for(size_t i = 0; i < wlen_; i++){
            free(context_weights[k][i]);
        }
        delete [] context_weights[k];
        free(pc_weights[k]);
        free(pc[k]);
    }
    delete [] context_weights;
    delete [] pc_weights;
    delete [] pc;
}

void ContextLibrary::read(std::string &libStr){
    std::stringstream in(libStr);
    // Parse and check header information
    if (!reader.StreamStartsWith(in, "CRF")){
        Debug(Debug::WARNING) << "Stream does not start with class id 'CRF'!\n";
        EXIT(EXIT_FAILURE);
    }
    std::string line = reader.getline(in);

    libSize = reader.ReadInt(line.c_str(), "SIZE", "Unable to parse context library 'SIZE'!");
    line = reader.getline(in);
    wlen_ = reader.ReadInt(line.c_str(), "LENG", "Unable to parse context library 'LENG'!");
    center = (wlen_ - 1) / 2;
    // Read context profiles
    //profiles_.Resize(size);
    size_t k;
    context_weights = new float**[libSize];
    pc_weights = new float*[libSize];
    pc = new float*[libSize];
    for (k = 0; k < libSize && !in.eof(); ++k){
        context_weights[k]           = new float*[wlen_];
        for(size_t i = 0; i < wlen_; i++) { // 24 instead of 20 for memory alignment
            context_weights[k][i] = (float *)mem_align(16, 24 * sizeof(float));
        }
        pc_weights[k] = (float *)mem_align(16, 24 * sizeof(float));
        pc[k] = (float *)mem_align(16, 24 * sizeof(float));
        readContextProfile(in, reader, context_weights[k], pc_weights[k], pc[k]);
        //            for(size_t i = 0; i < wlen_; i++) {
        //                std::cout << i<< "\t";
        //                for(size_t a = 0; a < 20; a++){
        //                   printf("%.3f\t",pc_weights[k][a]);
        //                }
        //                std::cout << std::endl;
        //            }
    }
    if (k != libSize){
        Debug(Debug::WARNING) << "Serialized context library should have "
                              << libSize << " profiles but actually has " << k << "!\n";
        EXIT(EXIT_FAILURE);
    }
}

void ContextLibrary::readContextProfile(std::stringstream &in, LibraryReader &reader,
                                        float ** context_weight, float * pc_weight, float * pc) {
    // Parse and check header information
    if (!reader.StreamStartsWith(in, "CrfState")){
        Debug(Debug::ERROR) << "Stream does not start with class id 'CrfState'!\n";
        EXIT(EXIT_FAILURE);
    }
    std::string line = reader.getline(in);
    if (strstr(line.c_str(), "NAME")) {
        names.push_back(reader.ReadString(line.c_str(), "NAME", "Unable to parse CRF state 'NAME'!"));
        line = reader.getline(in);
    }
    bias_weight.push_back(reader.ReadDouble(line.c_str(), "BIAS", "Unable to parse CRF state 'BIAS'!"));
    line = reader.getline(in);
    size_t len = reader.ReadInt(line.c_str(), "LENG", "Unable to parse CRF state 'LENG'!");
    line = reader.getline(in);
    size_t nalph = reader.ReadInt(line.c_str(), "ALPH", "Unable to parse CRF state 'ALPH'!");
    if (nalph != 20){
        Debug(Debug::ERROR) << "Alphabet size of serialized CRF state should be 20 "
                "but is acutally" << nalph << "!\n";
        EXIT(EXIT_FAILURE);
    }
    // If everything went fine we can resize our data memmbers

    // Read context weights and pseudocount weights
    size_t i = 0;
    line = reader.getline(in); // skip alphabet description line
    line = reader.getline(in);
//    float minFlt = 0.0;
//    float maxFlt = 0.0;

    while (line[0] != '/' && line[1] != '/') {
        char * ptr = (char *)line.c_str();
        if (line[0] != 'P' && line[1] != 'C') {
            i = strtol(ptr, NULL, 10) - 1;
            // TODO: include ANY char in serialization
            for (size_t a = 0; a < 20; ++a) {
                ptr += Util::skipNoneWhitespace(ptr);
                ptr += Util::skipWhitespace(ptr);
                int mm_a = ProfileStates::hh2mmseqsAAorder(a);
                context_weight[i][mm_a] = static_cast<double>(strtol(ptr, NULL, 10)) / kScale;
//                minFlt = std::min(minFlt, context_weight[i][mm_a]);
//                maxFlt = std::max(maxFlt, context_weight[i][mm_a]);
            }
//            std::cout << minFlt << "\t" << maxFlt << std::endl;

            context_weight[i][20] = 0.0;
        } else {
            for (size_t a = 0; a < 20; ++a){
                ptr += Util::skipNoneWhitespace(ptr);
                ptr += Util::skipWhitespace(ptr);
                int mm_a = ProfileStates::hh2mmseqsAAorder(a);
                pc_weight[mm_a] = static_cast<double>(strtol(ptr,NULL,10)) / kScale;
            }
        }
        context_weight[i][20] = 0.0;
        pc_weight[20] = 0.0;
        line = reader.getline(in);
    }
//    std::cout << minFlt << "\t" << maxFlt << std::endl;

    if (i != len - 1){
        Debug(Debug::ERROR) << "CRF state should have " << len <<" columns but actually has "<< i+1 << "!\n";
        EXIT(EXIT_FAILURE);
    }
    // Calculate maximum of pseudocount weights
    double max = -DBL_MAX;
    double mean = 0.0;
    for (size_t a = 0; a < 20; ++a) {
        mean += pc_weight[a];
        if (pc_weight[a] > max) max = pc_weight[a];
    }
    mean /= 20.0;

    // Rescale pseudocount weights and calculate their sum in lin-space
    long double sum = 0.0;
    for (size_t a = 0; a < 20; ++a)
        sum += exp(pc_weight[a] - max);

    // Update emission pseudocount vector
    double tmp = max + log(sum);
    for (size_t a = 0; a < 20; ++a) {
        pc[a] = DBL_MIN + exp(pc_weight[a] - tmp);
        // state.pc_weights[a] -= mean; // Not necessary if pc_weights are centered on central context weights
    }
}

inline float CSProfile::computeContextScore(float ** context_weights,
                                            const unsigned char * seq, const int L,
                                            size_t idx, size_t center) {
    const size_t beg = std::max(0, static_cast<int>(idx - center));
    const size_t end = std::min(static_cast<size_t >(L), idx + center + 1);
    const size_t diff = end - beg;
    float score = 0.0;
    float score1, score2, score3, score4;
    size_t j = beg - idx + center;
    switch(diff){
        case 1:
            score += context_weights[j][seq[beg]];
            break;
        case 2:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score = score1 + score2;
            break;
        case 3:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score = score1 + score2 + score3;
            break;
        case 4:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score = score1 + score2 + score3 + score4;
            break;
        case 5:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score = score1 + score2 + score3 + score4;
            break;
        case 6:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score = score1 + score2 + score3 + score4;
            break;
        case 7:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score = score1 + score2 + score3 + score4;
            break;
        case 8:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score = score1 + score2 + score3 + score4;
            break;
        case 9:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score = score1 + score2 + score3 + score4;
            break;
        case 10:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score2 += context_weights[j+9][seq[beg+9]];
            score = score1 + score2 + score3 + score4;
            break;
        case 11:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score2 += context_weights[j+9][seq[beg+9]];
            score3 += context_weights[j+10][seq[beg+10]];
            score = score1 + score2 + score3 + score4;
            break;
        case 12:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score2 += context_weights[j+9][seq[beg+9]];
            score3 += context_weights[j+10][seq[beg+10]];
            score4 += context_weights[j+11][seq[beg+11]];
            score = score1 + score2 + score3 + score4;
            break;
        case 13:
            //simdi_loadu()
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score2 += context_weights[j+9][seq[beg+9]];
            score3 += context_weights[j+10][seq[beg+10]];
            score4 += context_weights[j+11][seq[beg+11]];
            score1 += context_weights[j+12][seq[beg+12]];
            score = score1 + score2 + score3 + score4;
            break;
    }
//    for(size_t i = beg, j = beg - idx + center; i < end; ++i, ++j)
//        score += context_weights[j][seq[i]];
    return score;
}






PSSMCalculator::Profile CSProfile::computeProfile(BaseMatrix * subMat, Sequence * seq, float neff, float pTau){
    //std::cout << "Adding pseudocounts ...\n";
    const int center = ctxLib->center;
    // Calculate posterior probability ppi[k] of state k given sequence window
    // around position 'i'
    std::fill(maximums, maximums + seq->L, -FLT_MAX);
    for (size_t k = 0; k < ctxLib->libSize; ++k) {
        float* ppi = &pp[k * seq->L];
        float bias = ctxLib->bias_weight[k];
        for (int i = 0; i < seq->L; i++) {
            float contextScore = computeContextScore(ctxLib->context_weights[k], seq->numSequence, seq->L, i, center);
            ppi[i] = bias + contextScore;
            if (ppi[i] > maximums[i]){
                maximums[i] = ppi[i];  // needed for log-sum-exp trick
            }
        }
    }
    // Log-sum-exp trick begins here
    std::fill(sums, sums + seq->L, 0.0);
    const float scalingLog2 = 1/log(2.0);
    for (size_t k = 0; k < ctxLib->libSize; ++k){
        float* ppi = &pp[k * seq->L];
        for (int i = 0; i < seq->L; i++) {
            // exp(x) = 2^(1/log(2) * x)
            // http://www.wolframalpha.com/input/?i=exp(x)+%3D++2%5E(y+*+x)+solve+for+y
            sums[i] += MathUtil::fpow2((ppi[i] - maximums[i]) * scalingLog2);
        }
    }
    for (int i = 0; i < seq->L; i++) {
        maximums[i] = maximums[i] + log(sums[i]);
    }
    for (size_t k = 0; k < ctxLib->libSize; ++k) {
        float *ppi = &pp[k * seq->L];
        for (int i = 0; i < seq->L; i++) {
            ppi[i] = MathUtil::fpow2((ppi[i] - maximums[i]) * scalingLog2);
        }
    }
    // Calculate pseudocount vector P(a|X_i)
    std::fill(profile, profile + (seq->L * Sequence::PROFILE_AA_SIZE), 0.0);
    calculatePseudoeCounts(seq, pp, profile, ctxLib);
    for (int i = 0; i < seq->L; i++) {
        MathUtil::NormalizeTo1(&profile[i * Sequence::PROFILE_AA_SIZE], Sequence::PROFILE_AA_SIZE);
        //for(size_t a = 0; a < Sequence::PROFILE_AA_SIZE; ++a){
        //    printf("%f\t",profile[i * Sequence::PROFILE_AA_SIZE + a]);
        //}
        //printf("\n");
    }

    //AdmixTo(seq, p, admix);
    double tau = pTau; //TODO
    double t = 1 - tau;
    for (int i = 0; i < seq->L; ++i) {
        for (int a = 0; a < 20; ++a) {
            profile[i*Sequence::PROFILE_AA_SIZE + a] *= tau;
        }
        profile[i * Sequence::PROFILE_AA_SIZE + seq->numSequence[i]] += t;
        // set neff
        Neff_M[i] = neff;
        consensusSequence[i] = seq->numSequence[i];
    }
    for (int i = 0; i < seq->L; ++i) {
        MathUtil::NormalizeTo1(&profile[i*Sequence::PROFILE_AA_SIZE], Sequence::PROFILE_AA_SIZE);
    }
    PSSMCalculator::computeLogPSSM(subMat, pssm, profile, 8.0,  seq->L, 0.0);

    return PSSMCalculator::Profile(pssm, profile, Neff_M, consensusSequence);
}
