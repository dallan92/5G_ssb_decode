#ifndef HELPERS_H
#define HELPERS_H
#include <complex.h>

#define PI 3.141592653589793

void compFliplr(float complex *d_orig, float complex *flip, int length);
void reFliplr(float *d_orig, float *d_flip, int length);
void intFliplr(int *d_orig, int *d_flip, int length);
void max(float *arr, int length, float *mx, int *indx);
void find(int *arr, int *indices, int length, int K, int dir);
void maxThresh(float *arr, float *thresh, int length, float *max, int *indx);
void readTestSig(float complex *x, int length);
void sig2file(void *x, int type, int length, char *path);
float complex compAvg(float complex *x, int length);
float reAvg(float *x, int length);
void getReal(float *y, float complex *x, int length);
void getImag(float *y, float complex *x, int length);
float absre(float x);
float sign(float x);
float min(float x, float y);
void mag2(float complex *x, float *y, int length);
float abs2(float complex x);
float var(float complex *x, unsigned int len);
void dec2bin(unsigned int x, unsigned int n, unsigned char *bin);
#endif
