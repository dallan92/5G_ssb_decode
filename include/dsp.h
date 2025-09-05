#ifndef DSP_H
#define DSP_H
#include <complex.h>
#include <fftw3.h>

struct fftParam {
  float complex *in;
  float complex *out;
  fftwf_plan p;
  unsigned int nfft;
};

void fft1(float complex *data, unsigned int nfft, fftwf_plan plan,
          float complex *in, float complex *out);
void fftInit(struct fftParam *ft, const unsigned int nfft,
             const unsigned int dir);
void fftDestroy(struct fftParam *ft);
void fftshift(float complex *data, const unsigned int nfft);
void fftNorm(float complex *x, unsigned int nfft);
void compFIR(float complex *y, float complex *x, float complex *h, int y_len,
             int h_len);
void reFIR(float *y, float *x, float *h, int y_len, int h_len);
void compFreqShift(float complex *x, float complex *y, float normFreq, int n);
float complex compCrossCorr(float complex *x, float complex *y, int length);
void fastConv(float complex *x, int xLen, float complex *h, unsigned int bs,
              struct fftParam *ft, struct fftParam *ift, float complex *y);
float *movingAvg(float *y, int sigLen, int ml);
float complex *linInterp(float complex *y, int yLen, int *x, int *q);
// float complex *polyDec2hb(float complex *in, const unsigned int sigLen);
#endif
