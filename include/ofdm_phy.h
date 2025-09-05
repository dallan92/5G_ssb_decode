#ifndef OFDM_PHY_H
#define OFDM_PHY_H
#include "../include/dsp.h"
#include <complex.h>

void ofdmDemod(float complex *ofdmSymbs, struct fftParam ft,
               const unsigned int cpLen, const unsigned int numSymbs,
               float complex *symbs);
float nVarEstnull(float complex *symb, unsigned int nfft, unsigned int numAvg,
                  unsigned int used);
float cpFreqSynch(float complex *x, int length, int nfft, int cpLen);
#endif
