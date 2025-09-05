#ifndef FEC_H
#define FEC_H
#include <complex.h>

unsigned char *getFrozen(const unsigned int N, const unsigned int K);
void polarSCD(float *llrs, unsigned char *frozen, unsigned int N,
              unsigned char *d);
void qpskLLR(float complex *qpsk_symbs, int length, float nVar, float *llrs);
void csiScale(float *llrs, const float *csi, const unsigned int llrLen);
unsigned char crcCheck(unsigned char *bits, unsigned char *poly,
                       unsigned int polLen, unsigned int bitLen);
#endif
