//
// Created by mad on 11/29/16.
//

#ifndef MMSEQS_CSBLASTPROFILES_H
#define MMSEQS_CSBLASTPROFILES_H

#include <stdlib.h>
#include <simd/simd.h>
#include "LibraryReader.h"
#include "Debug.h"
#include "Util.h"
#include "ProfileStates.h"
#include "Sequence.h"
class ContextLibrary{

public:

    size_t wlen_;                            // size of context window.
    size_t center;                           // center of window
    size_t libSize;                          // size of library (e.g. 4000 K4000.crf)
    LibraryReader reader;
    std::vector<std::string> names;
    std::vector<float>  bias_weight;
    // keep ContextLibrary probs
    float *** context_weights; // k * states * aa
    float ** pc_weights;
    float ** pc;

    static ContextLibrary* getContextLibraryInstance()
    {
#pragma omp critical
        {
            if (contextLibrary==NULL) {
                contextLibrary = new ContextLibrary();
            }
        }
        return contextLibrary;
    }

private:
    static ContextLibrary * contextLibrary;

    ContextLibrary();
    ~ContextLibrary();

    // read library file
    void read(std::string &libStr);

    // read context profiles from library
    void readContextProfile(std::stringstream &in, LibraryReader &reader,
                            float ** context_weight, float * pc_weight, float * pc);
};

class CSProfile {
    ContextLibrary * ctxLib;
    float * profile;
    float * pp;
    float * maximums;
    float * sums;
public:

    CSProfile(size_t maxSeqLen) {
        ctxLib = ContextLibrary::getContextLibraryInstance();
        this->profile = (float * )mem_align(16, Sequence::PROFILE_AA_SIZE * maxSeqLen * sizeof(float));
        this->pp =  (float * ) mem_align(16, 4000 * maxSeqLen * sizeof(float)); //TODO how can I avoid th 4000?
        this->maximums = new float[maxSeqLen];
        this->sums =  new float[maxSeqLen];
    };

    ~CSProfile(){
        free(profile);
        free(pp);
        delete [] maximums;
        delete [] sums;
    }

    float computeContextScore(float ** context_weights,
                              const unsigned char * seq, const int L,
                              size_t idx, size_t center);

    void calculatePseudoeCounts(Sequence * seq, float * pp, float * profile,ContextLibrary * ctxLib){
        for (size_t k = 0; k < ctxLib->libSize; ++k){
            float * ppi = &pp[k * seq->L];
            float * ctxLib_pc = ctxLib->pc[k];
            for (int i = 0; i < seq->L; i++) {
                float *pc = &profile[i * Sequence::PROFILE_AA_SIZE];
                __m128 simd_ppi = _mm_set_ps1(ppi[i]);
                for (size_t a = 0; a < Sequence::PROFILE_AA_SIZE; a += 4) {
                    //pc[a] += ppi[i] * ctxLib_pc[a];
                    __m128 ctxLib_pc_a = _mm_load_ps(&ctxLib_pc[a]);
                    __m128 pc_a = _mm_load_ps(&pc[a]);
                    __m128 pc_res = _mm_add_ps(pc_a, _mm_mul_ps(ctxLib_pc_a, simd_ppi));
                    _mm_store_ps(&pc[a], pc_res);
                }
            }
        }
    }

    float * computeProfile(Sequence * seq, float neff, float tau);

};


#endif //MMSEQS_CSBLASTPROFILES_H
