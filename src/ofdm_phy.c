#include "../include/ofdm_phy.h"
#include "../include/dsp.h"
#include "../include/helpers.h"
#include <math.h>
#include <stdlib.h>

/* ----- OFDM PHY Utilities ----- */

/* Performs OFDM demodulation of multiple OFDM symbols
 Assumes that input begins with 1st sample of the CP of the 1st symbol */
void ofdmDemod(float complex *ofdmSymbs, struct fftParam ft,
               const unsigned int cpLen, const unsigned int numSymbs,
               float complex *symbs) {

  /*  Peform FFT on each OFDM symbol */
  int l = 0;
  int k = cpLen;

  for (int i = 0; i < numSymbs; i++) {
    for (int j = 0; j < ft.nfft; j++) {
      symbs[l + j] = ofdmSymbs[k + j];
    }
    fft1(&symbs[l], ft.nfft, ft.p, ft.in, ft.out);
    fftshift(&symbs[l], ft.nfft);
    fftNorm(&symbs[l], ft.nfft);
    l += ft.nfft;
    k += ft.nfft + cpLen;
  }
}

/* Function to estimate noise variance using null sub-carriers of OFDM signal.
   Assumes DC is at centre of FFT (two sided spectrum)
*/
float nVarEstnull(float complex *symb, unsigned int nfft, unsigned int numAvg,
                  unsigned int used) {

  /* Estimate using numAvg samples from null sub-carriers
     on lower edge of carrier. */
  unsigned int low = (nfft - used) / 2;
  unsigned int start;
  float nVar = 0.0;
  if (numAvg > low) {
    /* Not enough samples available for desired averaging length.
       Use all available samples. */
    start = 0;
    nVar = var(&symb[start], low);
  } else {
    start = (low - 1) - (numAvg - 1);
    nVar = var(&symb[start], numAvg);
  }
  return nVar;
}

/* Function to estimate fractional frequency offset using Cyclic Prefix (CP)
 * algorithm. */
float cpFreqSynch(float complex *x, int length, int nfft, int cpLen) {
  /* Create signal delayed by lag= nfft */
  float complex *x_lag = (float complex *)calloc(length, sizeof(float complex));
  for (int i = 0; i < (length - nfft); i++)
    x_lag[i + nfft] = x[i];

  /* Perform auto-correlation at lag = nfft */
  float complex *x_corr =
      (float complex *)calloc(length, sizeof(float complex));
  for (int i = 0; i < length; i++)
    x_corr[i] = x[i] * conjf(x_lag[i]);

  /* Pass through MA filter of length = cpLen */
  int zpLen = length + cpLen - 1;
  float complex *y = (float complex *)calloc(zpLen, sizeof(float complex));
  for (int i = 0; i < length; i++)
    y[i + cpLen - 1] = x_corr[i];

  float complex h[cpLen]; // MA filter weights
  for (int i = 0; i < cpLen; i++)
    h[i] = 1.0 + 0.0 * I;

  float complex *y_filt =
      (float complex *)calloc(length, sizeof(float complex));
  compFIR(y_filt, y, h, length, cpLen);

  /* Average correlation result over 4 SSB symbols to improve SNR */
  float complex *y_avg =
      (float complex *)calloc(nfft + cpLen, sizeof(float complex));
  float complex val[4] = {0.0 + 0.0 * I};
  int k;
  for (int i = 0; i < (nfft + cpLen); i++) {
    k = 0;
    for (int j = 0; j < 4; j++) {
      val[j] = y_filt[i + k];
      k += nfft + cpLen;
    }
    // Average samples in val array
    y_avg[i] = compAvg(val, 4);
  }

  /* Find mag2 of averaged output */
  float *y_avg_abs2 = (float *)calloc(nfft + cpLen, sizeof(float));
  mag2(y_avg, y_avg_abs2, nfft + cpLen);
  float mx;
  int indx;
  max(y_avg_abs2, nfft + cpLen, &mx, &indx);

  /* Find normalised frequency offset */
  float f_err = (cargf(y_avg[indx])) / (2 * PI * nfft);

  /* Free allocated memory */
  free(x_lag);
  free(x_corr);
  free(y);
  free(y_filt);
  free(y_avg);
  free(y_avg_abs2);

  return f_err;
}
