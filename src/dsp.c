#include "../include/dsp.h"
#include "../include/helpers.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Function to perform FFT */
void fft1(float complex *data, unsigned int nfft, fftwf_plan plan,
          float complex *in, float complex *out) {
  memcpy(in, data, sizeof(float complex) * nfft);
  fftwf_execute_dft(plan, in, out);
  memcpy(data, out, sizeof(float complex) * nfft);
}

/* Function to initialise fftParam object */
void fftInit(struct fftParam *ft, const unsigned int nfft,
             const unsigned int dir) {
  ft->in = (float complex *)malloc(sizeof(float complex) * nfft);
  ft->out = (float complex *)malloc(sizeof(float complex) * nfft);
  if (dir == 0) {
    ft->p =
        fftwf_plan_dft_1d(nfft, ft->in, ft->out, FFTW_FORWARD, FFTW_ESTIMATE);
  } else {
    ft->p =
        fftwf_plan_dft_1d(nfft, ft->in, ft->out, FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  ft->nfft = nfft;
}

/* Function to destroy fftParam object */
void fftDestroy(struct fftParam *ft) {
  fftwf_destroy_plan(ft->p);
  free(ft->in);
  free(ft->out);
}

/* Function to normalise the FFT */
void fftNorm(float complex *x, unsigned int nfft) {
  for (int i = 0; i < nfft; i++) {
    x[i] = x[i] / (float)nfft;
  }
}

/* Function to perform fftshift */
void fftshift(float complex *data, const unsigned int nfft) {

  float complex last;

  for (int i = 0; i < nfft / 2; i++) {

    /* Last value of array */
    last = data[nfft - 1];

    /* shift elements of array 1 to right */
    for (int j = nfft - 1; j > 0; j--)
      data[j] = data[j - 1];

    /* Update 1st element of array to previous last */
    data[0] = last;
  }
}

/* Complex FIR filter */
void compFIR(float complex *y, float complex *x, float complex *h, int y_len,
             int h_len) {
  float complex temp;
  int l;
  for (int i = 0; i < y_len; i++) {
    temp = 0 + 0 * I;
    l = h_len - 1;
    for (int j = 0; j < h_len; j++) {
      temp += x[i + j] * h[l];
      l--;
    }
    y[i] = temp;
  }
}

/* Fast convolution using Overlap Add method  */
void fastConv(float complex *x, int xLen, float complex *h, unsigned int bs,
              struct fftParam *ft, struct fftParam *ift, float complex *y) {

  /*  Perform FFT of filter coefficients as they don't change */
  unsigned int hLen = ft->nfft - (bs - 1);
  float complex *hZP = (float complex *)calloc(ft->nfft, sizeof(float complex));
  for (int i = 0; i < hLen; i++)
    hZP[i] = h[i];
  fft1(hZP, ft->nfft, ft->p, ft->in, ft->out);

  /* Array to hold corelation */
  float complex *corr =
      (float complex *)calloc(ft->nfft, sizeof(float complex));
  float complex *blockZP =
      (float complex *)calloc(ft->nfft, sizeof(float complex));

  /* Perform FFT convolution using overlap add method */
  for (int i = 0; i < xLen; i += bs) {

    /* Re-initialise FFT block */
    for (int z = 0; z < ft->nfft; z++) {
      blockZP[z] = 0.0 + 0.0 * I;
    }

    /* Take block of samples from input signal */
    for (int j = 0; j < bs && (i + j) < xLen; j++)
      blockZP[j] = x[i + j];

    /* Perform FFT */
    fft1(blockZP, ft->nfft, ft->p, ft->in, ft->out);

    /* Perform element wise multiplication */
    for (int k = 0; k < ft->nfft; k++)
      corr[k] = blockZP[k] * hZP[k];

    /* Perfrom IFFT */
    fft1(corr, ft->nfft, ift->p, ift->in, ift->out);

    /* Oerlap and Add */
    if ((i + bs) > xLen) {
      unsigned int zp = bs - (xLen % bs);
      for (int l = 0; l < ft->nfft - zp; l++)
        y[i + l] += corr[l] / (float)ft->nfft;
    } else {
      for (int l = 0; l < ft->nfft; l++)
        y[i + l] += corr[l] / (float)ft->nfft;
    }
  }

  /* Free memory */
  free(hZP);
  free(corr);
  free(blockZP);
}

/* Real FIR filter */
void reFIR(float *y, float *x, float *h, int y_len, int h_len) {
  float temp;
  int l;
  for (int i = 0; i < y_len; i++) {
    temp = 0;
    l = h_len - 1;
    for (int j = 0; j < h_len; j++) {
      temp += x[i + j] * h[l];
      l--;
    }
    y[i] = temp;
  }
}

/* Perfrom cross correlation of complex sequences */
float complex compCrossCorr(float complex *x, float complex *y, int length) {
  /* Perform cross correlation */
  float complex *corr = (float complex *)malloc(sizeof(float complex) * length);
  for (int i = 0; i < length; i++)
    corr[i] = x[i] * conjf(y[i]);

  /* Compute average */
  float complex avg = compAvg(corr, length);

  /* Free memory */
  free(corr);

  return avg;
}

/* Function to perform complex frequency shift */
void compFreqShift(float complex *x, float complex *y, float normFreq, int n) {
  for (int i = 0; i < n; i++)
    y[i] =
        x[i] * (cosf(2 * PI * normFreq * i) + sinf(2 * PI * normFreq * i) * I);
}

/* Moving Average filter */
float *movingAvg(float *x, int sigLen, int ml) {

  /* Zero pad by ml-1 samples so output is same length as input */
  int zpLen = sigLen + ml - 1;
  float *xZP = (float *)calloc(zpLen, sizeof(float));
  for (int i = 0; i < sigLen; i++)
    xZP[i + ml - 1] = x[i];

  /* Get sum of 1st ml samples */
  float initSum = reAvg(xZP, ml) * (float)ml;
  float temp = initSum;

  /* Create y and assign 1st output */
  float *y = (float *)malloc(sigLen * sizeof(float));
  y[0] = temp;

  for (int i = 1; i < sigLen; i++) {
    /* Add newest sample + subtract oldest sample from sum */
    temp = temp + xZP[i + ml - 1] - x[i - 1];
    y[i] = temp;
  }

  free(xZP);
  return y;
}

/* Function to perform compalex linear interpolation between two points (Y1, X1)
 * and (Y2, X2) at points q */
static void interp(float complex y1, float complex y2, int x1, int x2, int num,
                   int step, float complex *yInterp) {
  int X = x1 + step;
  yInterp[0] = y1;

  // Gradient
  float complex m = (y2 - y1) / (x2 - x1);

  // Compute interpolated values
  int dx;
  for (int i = 1; i < num; i++) {
    dx = X - x1;
    yInterp[i] = y1 + m * (float)dx;
    X += step;
  }
}

/* Function to perform linear interpolation of a complex sequence at points x
 * evaluated at query points q  */
float complex *linInterp(float complex *y, int yLen, int *x, int *q) {

  // No. of interpolations between adjacent points
  unsigned int nInterp = yLen - 1;

  // Step along x-axis of interpolated sequence
  unsigned int step = q[1] - q[0];

  // No. of interpolated values between adjacent samples of orignal sequence
  unsigned int nInterpVals = ((x[1] - x[0]) / step) - 1;

  // Length of interpolated sequence
  unsigned int seqLen = yLen + nInterpVals * nInterp;

  // Create vectors to hold interpolated + final sequences
  float complex *yInt =
      (float complex *)malloc((seqLen - 1) * sizeof(float complex));
  float complex *yFinal =
      (float complex *)malloc((seqLen) * sizeof(float complex));

  // Perform linear interpolation
  int j = 0;
  int k = 0;
  const unsigned int num = nInterpVals + 1;
  float complex yInterp[num];

  for (int i = 0; i < nInterp; i++) {
    interp(y[j], y[j + 1], x[j], x[j + 1], num, step, yInterp);
    for (int l = 0; l < num; l++) {
      yInt[l + k] = yInterp[l];
    }
    j += 1;
    k += num;
  }

  // Construct final sequence
  yFinal[seqLen - 1] = y[yLen - 1]; // copy last sample
  for (int i = 0; i < (seqLen - 1); i++) {
    yFinal[i] = yInt[i];
  }

  // Free memory
  free(yInt);

  return yFinal;
}
