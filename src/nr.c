#include "../include/nr.h"
#include "../include/dsp.h"
#include "../include/fec.h"
#include "../include/helpers.h"
#include "../include/ofdm_phy.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ---------- 5G NR Utilities ---------- */

/* Pseudo random sequence generation - Clause 5.2 of TS 38.211 */
static void nrPRS(float *prs, uint32_t x2init, uint32_t length) {
  uint32_t Nc = 1600;
  uint32_t x1, x2;       /* LFSR registers */
  uint32_t x1_fb, x2_fb; /* Feedback values */

  /* Initialise registers */
  x1 = 0x1;
  x2 = x2init;

  for (uint32_t n = 0; n < Nc + length; n++) {
    if (n >= Nc)
      /* Only take last length values */
      prs[n - Nc] = (float)((x1 ^ x2) & 0x1);

    /* Update x1 */
    x1_fb = (x1 >> 3 ^ x1) & 0x1;
    x1 = (x1 >> 1) | (x1_fb << 30);

    /* Update x2 */
    x2_fb = (x2 >> 3 ^ x2 >> 2 ^ x2 >> 1 ^ x2) & 0x1;
    x2 = (x2 >> 1) | (x2_fb << 30);
  }
}

/* Frequency Domain PSS sequence */
static void nrPSS(float *d_pss, int n_id_2) {
  /*m-sequence */
  int x[] = {0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0,
             0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0,
             0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
             0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
             1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1,
             0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1};

  /* Compute d_pss */
  for (int n = 0; n < 127; n++)
    d_pss[n] = 1 - 2 * x[(n + 43 * n_id_2) % 127];
}

/* Frequency Domaain SSS sequence */
static void nrSSS(float complex *d_sss, int n_id_1, int n_id_2) {
  /* 1st m-sequence */
  int x0[] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1,
              1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1,
              0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1,
              0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0,
              1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1,
              0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0};

  /* 2nd m-sequence */
  int x1[] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1,
              0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1,
              0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0,
              0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,
              0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0,
              0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1};

  /* Cyclic shifts */
  int m0 = 15 * (int)floor((float)n_id_1 / 112.0) + 5 * n_id_2;
  int m1 = n_id_1 % 112;

  /* Generate SSS sequence */
  for (int n = 0; n < 127; n++)
    d_sss[n] =
        ((1 - 2 * x0[(n + m0) % 127]) * (1 - 2 * x1[(n + m1) % 127])) + 0 * I;
}

/* Time Domain PSS sequence */
static void nrPSSTime(float complex *pssTime, float *pssFreq,
                      struct fftParam *ft) {

  /* Map to REs */
  unsigned int scStart = (ft->nfft / 2 - 1 - 63);
  for (int i = 0; i < 127; i++)
    pssTime[i + scStart] = (float complex)pssFreq[i];

  /* Perform IFFT */
  fftshift(pssTime, ft->nfft);
  fft1(pssTime, ft->nfft, ft->p, ft->in, ft->out);

  /* Normalise IFFT */
  for (int i = 0; i < ft->nfft; i++)
    pssTime[i] = pssTime[i] / (float)ft->nfft;
}

/* ---------- PSS / SSS detection ---------- */

/* Perform PSS detection to find SSB and detremine n_id_2 part of PCI. */
static void pssCorr(float complex *x, const unsigned int sigLen,
                    struct ssbConf *s, int *ssbStart, int *n_id_2) {

  int convLen = sigLen + s->nfft - 1;
  float complex *y = (float complex *)calloc(convLen, sizeof(float complex));
  float *y_abs2 = (float *)calloc(convLen, sizeof(float));
  float complex *h = (float complex *)calloc(s->nfft, sizeof(float complex));

  float mx = 0.0f;
  float mx_curr = mx;
  int init_est_ssb_start = 0, init_est_n_id_2 = 0, indx_max = 0;

  float d_pss[127];

  // Initialise FFT plans
  struct fftParam conv_ft = {0};
  struct fftParam conv_ift = {0};
  struct fftParam pss_ift = {0};

  unsigned int convLenSeg = 2 * s->nfft;
  unsigned int blockSize = convLenSeg - s->nfft + 1;

  fftInit(&conv_ft, convLenSeg, 0);
  fftInit(&conv_ift, convLenSeg, 1);
  fftInit(&pss_ift, s->nfft, 1);

  /* For each value of n_id_2 */
  for (int i = 0; i < 3; i++) {

    // Generate time domain PSS
    nrPSS(d_pss, i);
    nrPSSTime(h, d_pss, &pss_ift);

    // Cross correlation
    for (int j = 0; j < convLen; j++)
      y[j] = 0.0 + 0.0 * I;
    fastConv(x, sigLen, h, blockSize, &conv_ft, &conv_ift, y);

    // Compute y_abs2
    mag2(y, y_abs2, convLen);

    // Find max of |abs(y)|^2
    max(y_abs2, convLen, &mx, &indx_max);

    // Update correlation result if new maximum peak is found
    if (mx > mx_curr) {
      mx_curr = mx;
      init_est_ssb_start = indx_max - (s->nfft + s->cpLen -
                                       1); // Peak occurs at end of PSS symbol
      init_est_n_id_2 = i;
    }
  }

  /* If capture begins during PSS and init_est_ssb_strat is -ve, force
     to 0 in order to avoid seg fault */
  if (init_est_ssb_start < 0) {
    printf("WARNING: Signal capture may have started during PSS so MIB "
           "decoding will likely "
           "fail.\n");
    printf("The physical cell group is: %d.\n", init_est_n_id_2);
    init_est_ssb_start = 0;
  }

  /* If ssb extends beynd end of capture, adjust estimated start to avoid a
    seg fault */
  unsigned int ssbLen = 2192;
  unsigned int last = init_est_ssb_start + ssbLen;
  if (last > sigLen) {
    printf("WARNING: SSB may extend beyond edge of captured signal so MIB "
           "decoding will likely fail.\n");
    printf("The physical cell group is: %d.\n", init_est_n_id_2);
    init_est_ssb_start -= (last - sigLen);
  }
  /* Return ssb_start and n_id_2 */
  *ssbStart = init_est_ssb_start;
  *n_id_2 = init_est_n_id_2;

  // Free allocated memory
  free(y);
  free(y_abs2);
  free(h);
  fftDestroy(&conv_ft);
  fftDestroy(&conv_ift);
  fftDestroy(&pss_ift);
}

/* Perform SSS correlations to find n_id_1 part of PCI */
static unsigned int sssCorr(float complex *rgSSB, unsigned int n_id_2) {

  // Perform SSS correlations to estimate n_id_1
  float sssCorr[336] = {0.0};
  float complex d_sss[127] = {0.0 + 0.0 * I};
  float complex avg;
  const int index = 296;

  for (int i = 0; i < 336; i++) {
    nrSSS(d_sss, i, n_id_2); // Generate local SSS sequence
    avg = compCrossCorr(&rgSSB[index], d_sss, 127);
    sssCorr[i] = abs2(avg);
  }

  // Index of max of SSS correlation = n_id_1
  float mx;
  int n_id_1;
  max(sssCorr, 336, &mx, &n_id_1);

  // Return estimated n_id_1
  return n_id_1;
}

/* Find transmitted PBCH DMRS sequence and i_bar_ssb for de-srambling of MIB */
static void dmrsSearch(float complex *rgSSB, int n_id_cell, float complex *dmrs,
                       int *i_bar_ssb) {
  uint32_t length = 288;
  float prs[length];
  float complex *dmrsLocal =
      (float complex *)malloc(sizeof(float complex) * (length / 2));
  float complex *c =
      (float complex *)malloc(sizeof(float complex) * (length / 2));
  float *cabs = (float *)malloc(sizeof(float) * (length / 2));
  float c_abs_avg[8];
  int v = n_id_cell % 4;

  /* DMRS indices */
  int dmrsInd[] = {
      0 + v,   4 + v,   8 + v,   12 + v,  16 + v,  20 + v,  24 + v,  28 + v,
      32 + v,  36 + v,  40 + v,  44 + v,  48 + v,  52 + v,  56 + v,  60 + v,
      64 + v,  68 + v,  72 + v,  76 + v,  80 + v,  84 + v,  88 + v,  92 + v,
      96 + v,  100 + v, 104 + v, 108 + v, 112 + v, 116 + v, 120 + v, 124 + v,
      128 + v, 132 + v, 136 + v, 140 + v, 144 + v, 148 + v, 152 + v, 156 + v,
      160 + v, 164 + v, 168 + v, 172 + v, 176 + v, 180 + v, 184 + v, 188 + v,
      192 + v, 196 + v, 200 + v, 204 + v, 208 + v, 212 + v, 216 + v, 220 + v,
      224 + v, 228 + v, 232 + v, 236 + v, 240 + v, 244 + v, 248 + v, 252 + v,
      256 + v, 260 + v, 264 + v, 268 + v, 272 + v, 276 + v, 280 + v, 284 + v,
      432 + v, 436 + v, 440 + v, 444 + v, 448 + v, 452 + v, 456 + v, 460 + v,
      464 + v, 468 + v, 472 + v, 476 + v, 480 + v, 484 + v, 488 + v, 492 + v,
      496 + v, 500 + v, 504 + v, 508 + v, 512 + v, 516 + v, 520 + v, 524 + v,
      528 + v, 532 + v, 536 + v, 540 + v, 544 + v, 548 + v, 552 + v, 556 + v,
      560 + v, 564 + v, 568 + v, 572 + v, 576 + v, 580 + v, 584 + v, 588 + v,
      592 + v, 596 + v, 600 + v, 604 + v, 608 + v, 612 + v, 616 + v, 620 + v,
      624 + v, 628 + v, 632 + v, 636 + v, 640 + v, 644 + v, 648 + v, 652 + v,
      656 + v, 660 + v, 664 + v, 668 + v, 672 + v, 676 + v, 680 + v, 684 + v,
      688 + v, 692 + v, 696 + v, 700 + v, 704 + v, 708 + v, 712 + v, 716 + v};

  /* Find most likely DMRS sequence based on which leads to maximum value of
   * channel power */
  int l = 0;

  for (int i = 0; i < 8; i++) {

    /* Initialise LFSR for PRS sequence according to Clause 7.4.1.4 in TS38.211
     */
    uint32_t x2init = (2048 * (i + 1) * (floor(n_id_cell / 4) + 1)) +
                      (64 * (i + 1)) + (n_id_cell % 4);
    nrPRS(prs, x2init, length);

    /* Map to complex symbols */
    l = 0;
    for (int j = 0; j < length / 2; j++) {
      dmrsLocal[j] = (1 / sqrt(2)) * (1 - (2 * prs[l])) +
                     I * (1 / sqrt(2)) * (1 - (2 * prs[l + 1]));
      l += 2;
    }

    /* Perform cross correlation of received + local sequences */
    for (int k = 0; k < 144; k++) {
      c[k] = rgSSB[dmrsInd[k]] / dmrsLocal[k];
    }

    /* Compute power of cross - correlation */
    for (int k = 0; k < 144; k++) {
      cabs[k] = abs2(c[k]);
    }
    c_abs_avg[i] = reAvg(cabs, 144);
  }

  /* Find maximum value in array to determine i_bar_ssb */
  float maxVal;
  max(c_abs_avg, 8, &maxVal, i_bar_ssb);

  /* Regenerate most likely DMRS sequence based on i_bar_ssb & n_id_cell */
  uint32_t x2init = (2048 * (*i_bar_ssb + 1) * (floor(n_id_cell / 4) + 1)) +
                    (64 * (*i_bar_ssb + 1)) + (n_id_cell % 4);
  nrPRS(prs, x2init, length);

  /* Map to complex symbols */
  l = 0;
  for (int i = 0; i < length / 2; i++) {
    dmrs[i] = (1 / sqrt(2)) * (1 - (2 * prs[l])) +
              I * (1 / sqrt(2)) * (1 - (2 * prs[l + 1]));
    l += 2;
  }

  /* Free memory */
  free(dmrsLocal);
  free(c);
  free(cabs);
}

/* ---------- Physical Broadcast Channel (PBCH) ---------- */

/* Generate scrambling sequence for PBCH */
static void pbchPRS(int n_id_cell, int v, float *seq) {
  /* Generate PRS */
  int m_bit = 864;
  int len = m_bit + v * 864;
  float cs[len];
  nrPRS(cs, n_id_cell, len);

  /* Extract last m_bit bits from generated sequence
     and BPSK modulate */
  int j = len - m_bit;
  for (int i = 0; i < m_bit; i++) {
    seq[i] = 1 - 2 * cs[j];
    j += 1;
  }
}

/* Extract PBCH resource grid */
static void pbchRGExtract(float complex *rxPBCH, int nfft,
                          float complex *rgPBCH) {
  int start = (nfft - 240) / 2;
  int k = start;
  int l = 0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 240; j++)
      rgPBCH[j + l] = rxPBCH[j + k];
    l += 240;
    k += nfft;
  }
}

/* Channel estimation for PBCH - assume channel is static across all SSB symbols
 */
static void pbchChanEst(float complex *pbchSymbs, float complex *dmrs, int v,
                        float complex *H, float *rsrp) {

  /* DMRS Indices */
  int dmrsInd[] = {
      0 + v,   4 + v,   8 + v,   12 + v,  16 + v,  20 + v,  24 + v,  28 + v,
      32 + v,  36 + v,  40 + v,  44 + v,  48 + v,  52 + v,  56 + v,  60 + v,
      64 + v,  68 + v,  72 + v,  76 + v,  80 + v,  84 + v,  88 + v,  92 + v,
      96 + v,  100 + v, 104 + v, 108 + v, 112 + v, 116 + v, 120 + v, 124 + v,
      128 + v, 132 + v, 136 + v, 140 + v, 144 + v, 148 + v, 152 + v, 156 + v,
      160 + v, 164 + v, 168 + v, 172 + v, 176 + v, 180 + v, 184 + v, 188 + v,
      192 + v, 196 + v, 200 + v, 204 + v, 208 + v, 212 + v, 216 + v, 220 + v,
      224 + v, 228 + v, 232 + v, 236 + v, 240 + v, 244 + v, 248 + v, 252 + v,
      256 + v, 260 + v, 264 + v, 268 + v, 272 + v, 276 + v, 280 + v, 284 + v,
      432 + v, 436 + v, 440 + v, 444 + v, 448 + v, 452 + v, 456 + v, 460 + v,
      464 + v, 468 + v, 472 + v, 476 + v, 480 + v, 484 + v, 488 + v, 492 + v,
      496 + v, 500 + v, 504 + v, 508 + v, 512 + v, 516 + v, 520 + v, 524 + v,
      528 + v, 532 + v, 536 + v, 540 + v, 544 + v, 548 + v, 552 + v, 556 + v,
      560 + v, 564 + v, 568 + v, 572 + v, 576 + v, 580 + v, 584 + v, 588 + v,
      592 + v, 596 + v, 600 + v, 604 + v, 608 + v, 612 + v, 616 + v, 620 + v,
      624 + v, 628 + v, 632 + v, 636 + v, 640 + v, 644 + v, 648 + v, 652 + v,
      656 + v, 660 + v, 664 + v, 668 + v, 672 + v, 676 + v, 680 + v, 684 + v,
      688 + v, 692 + v, 696 + v, 700 + v, 704 + v, 708 + v, 712 + v, 716 + v};

  /* Get channel estimates for DMRS symbols */
  const unsigned int ndmrs = 144;
  float complex *H_dmrs =
      (float complex *)malloc(sizeof(float complex) * ndmrs);
  float *dmrs_abs2 = (float *)malloc(sizeof(float) * ndmrs);
  for (int i = 0; i < ndmrs; i++) {
    H_dmrs[i] = pbchSymbs[dmrsInd[i]] / dmrs[i];
    dmrs_abs2[i] = abs2(pbchSymbs[dmrsInd[i]]);
  }

  /* Compute the Reference Signal Received Power (RSRP) in dBm */
  float dmrs_pwr = reAvg(dmrs_abs2, ndmrs);
  *rsrp = 10 * log10(dmrs_pwr) -
          60; // Needs actual calibration but just using rx gain for now.

  /* Interpolate in frequency domain to get channel estimates for data positions
   */
  const unsigned int interplen = 237;
  const unsigned int ndmrsPerSymb = 60;
  float complex *H_interp_1 =
      (float complex *)malloc(sizeof(float complex) * interplen);
  float complex *H_interp_3 =
      (float complex *)malloc(sizeof(float complex) * interplen);
  int x[ndmrsPerSymb];
  int q[interplen];
  int j = 1;
  int k = 1;
  for (int i = 0; i < ndmrsPerSymb; i++) {
    x[i] = j;
    j += 4;
  }
  for (int i = 0; i < interplen; i++) {
    q[i] = k;
    k += 1;
  }
  const unsigned int start = 84;
  H_interp_1 = linInterp(H_dmrs, ndmrsPerSymb, x, q);
  H_interp_3 = linInterp(&H_dmrs[start], ndmrsPerSymb, x, q);

  /* Final channel estimates for each symbol */
  const unsigned int numREssbPerSymb = 240;
  float complex *H_final_1 =
      (float complex *)malloc(sizeof(float complex) * numREssbPerSymb);
  float complex *H_final_3 =
      (float complex *)malloc(sizeof(float complex) * numREssbPerSymb);

  /* Populate final channel estimate across entire resource grid including
   * missing sub-carriers with no channel estimate */
  if (v == 0) {
    memcpy(&H_final_1[0], &H_interp_1[0], sizeof(float complex) * interplen);
    memcpy(&H_final_3[0], &H_interp_3[0], sizeof(float complex) * interplen);
    int indOuter[] = {237, 238, 239};
    int indNear[] = {236, 236, 236};
    for (int i = 0; i < 3; i++) {
      H_final_1[indOuter[i]] = H_interp_1[indNear[i]];
      H_final_3[indOuter[i]] = H_interp_3[indNear[i]];
    }
  } else if (v == 1) {
    memcpy(&H_final_1[1], &H_interp_1[0], sizeof(float complex) * interplen);
    memcpy(&H_final_3[1], &H_interp_3[0], sizeof(float complex) * interplen);
    int indOuter[] = {0, 238, 239};
    int indNear[] = {0, 236, 236};
    for (int i = 0; i < 3; i++) {
      H_final_1[indOuter[i]] = H_interp_1[indNear[i]];
      H_final_3[indOuter[i]] = H_interp_3[indNear[i]];
    }
  } else if (v == 2) {
    memcpy(&H_final_1[2], &H_interp_1[0], sizeof(float complex) * interplen);
    memcpy(&H_final_3[2], &H_interp_3[0], sizeof(float complex) * interplen);
    int indOuter[] = {0, 1, 239};
    int indNear[] = {0, 0, 236};
    for (int i = 0; i < 3; i++) {
      H_final_1[indOuter[i]] = H_interp_1[indNear[i]];
      H_final_3[indOuter[i]] = H_interp_3[indNear[i]];
    }
  } else if (v == 3) {
    memcpy(&H_final_1[3], &H_interp_1[0], sizeof(float complex) * interplen);
    memcpy(&H_final_3[3], &H_interp_3[0], sizeof(float complex) * interplen);
    int indOuter[] = {0, 1, 2};
    int indNear[] = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
      H_final_1[indOuter[i]] = H_interp_1[indNear[i]];
      H_final_3[indOuter[i]] = H_interp_3[indNear[i]];
    }
  }

  /* Create channel estimate for symbol 2 */
  const unsigned int interpLenSymb2 = 45;
  const unsigned int ndmrsUpLow = 12;
  float complex *H_interp_2_lower =
      (float complex *)malloc(sizeof(float complex) * interpLenSymb2);
  float complex *H_interp_2_upper =
      (float complex *)malloc(sizeof(float complex) * interpLenSymb2);

  int x2[ndmrsUpLow];
  int q2[interpLenSymb2];
  j = 1;
  k = 1;
  for (int i = 0; i < ndmrsUpLow; i++) {
    x2[i] = j;
    j += 4;
  }
  for (int i = 0; i < interpLenSymb2; i++) {
    q2[i] = k;
    k += 1;
  }
  H_interp_2_lower = linInterp(&H_dmrs[ndmrsPerSymb], ndmrsUpLow, x2, q2);
  H_interp_2_upper =
      linInterp(&H_dmrs[ndmrsPerSymb + ndmrsUpLow], ndmrsUpLow, x2, q2);

  float complex *H_final_2_lower =
      (float complex *)malloc(sizeof(float complex) * (interpLenSymb2 + 3));
  float complex *H_final_2_upper =
      (float complex *)malloc(sizeof(float complex) * (interpLenSymb2 + 3));

  // Fill in missing outer sub-carriers
  if (v == 0) {
    memcpy(&H_final_2_lower[0], &H_interp_2_lower[0],
           sizeof(float complex) * interpLenSymb2);
    memcpy(&H_final_2_upper[0], &H_interp_2_upper[0],
           sizeof(float complex) * interpLenSymb2);
    int indOuter[] = {45, 46, 47};
    int indNear[] = {44, 44, 44};
    for (int i = 0; i < 3; i++) {
      H_final_2_lower[indOuter[i]] = H_interp_2_lower[indNear[i]];
      H_final_2_upper[indOuter[i]] = H_interp_2_upper[indNear[i]];
    }
  } else if (v == 1) {
    memcpy(&H_final_2_lower[1], &H_interp_2_lower[0],
           sizeof(float complex) * interpLenSymb2);
    memcpy(&H_final_2_upper[1], &H_interp_2_upper[0],
           sizeof(float complex) * interpLenSymb2);
    int indOuter[] = {0, 46, 47};
    int indNear[] = {0, 44, 44};
    for (int i = 0; i < 3; i++) {
      H_final_2_lower[indOuter[i]] = H_interp_2_lower[indNear[i]];
      H_final_2_upper[indOuter[i]] = H_interp_2_upper[indNear[i]];
    }
  } else if (v == 2) {
    memcpy(&H_final_2_lower[2], &H_interp_2_lower[0],
           sizeof(float complex) * interpLenSymb2);
    memcpy(&H_final_2_upper[2], &H_interp_2_upper[0],
           sizeof(float complex) * interpLenSymb2);
    int indOuter[] = {0, 1, 47};
    int indNear[] = {0, 0, 44};
    for (int i = 0; i < 3; i++) {
      H_final_2_lower[indOuter[i]] = H_interp_2_lower[indNear[i]];
      H_final_2_upper[indOuter[i]] = H_interp_2_upper[indNear[i]];
    }
  } else if (v == 3) {
    memcpy(&H_final_2_lower[3], &H_interp_2_lower[0],
           sizeof(float complex) * interpLenSymb2);
    memcpy(&H_final_2_upper[3], &H_interp_2_upper[0],
           sizeof(float complex) * interpLenSymb2);
    int indOuter[] = {0, 1, 2};
    int indNear[] = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
      H_final_2_lower[indOuter[i]] = H_interp_2_lower[indNear[i]];
      H_final_2_upper[indOuter[i]] = H_interp_2_upper[indNear[i]];
    }
  }

  /* Create final estimate for channel 2 */
  float complex *H_final_2_full =
      (float complex *)calloc(sizeof(float complex), numREssbPerSymb);
  memcpy(&H_final_2_full[0], &H_final_2_lower[0],
         sizeof(float complex) * (interpLenSymb2 + 3));
  memcpy(&H_final_2_full[192], &H_final_2_upper[0],
         sizeof(float complex) * (interpLenSymb2 + 3));

  /* Create final channel estimate for entire resource grid  */
  memcpy(H, H_final_1, sizeof(float complex) * numREssbPerSymb);
  memcpy(&H[numREssbPerSymb], H_final_2_full,
         sizeof(float complex) * numREssbPerSymb);
  memcpy(&H[2 * numREssbPerSymb], H_final_3,
         sizeof(float complex) * numREssbPerSymb);

  /* Free memory resources */
  free(H_dmrs);
  free(dmrs_abs2);
  free(H_interp_1);
  free(H_final_1);
  free(H_interp_3);
  free(H_final_3);
  free(H_interp_2_lower);
  free(H_interp_2_upper);
  free(H_final_2_lower);
  free(H_final_2_upper);
  free(H_final_2_full);
}

/* Equalise PBCH symbols */
static void pbchEqualise(float complex *pbchSymb, float complex *H, int v,
                         float complex *pbchEq, float *csi, float nVar) {
  /* PBCH symbol indices */
  const unsigned int numREpbch = 432;
  int pbchInd[numREpbch];

  if (v == 0) {
    int pbchInd_0[] = {
        1,   2,   3,   5,   6,   7,   9,   10,  11,  13,  14,  15,  17,  18,
        19,  21,  22,  23,  25,  26,  27,  29,  30,  31,  33,  34,  35,  37,
        38,  39,  41,  42,  43,  45,  46,  47,  49,  50,  51,  53,  54,  55,
        57,  58,  59,  61,  62,  63,  65,  66,  67,  69,  70,  71,  73,  74,
        75,  77,  78,  79,  81,  82,  83,  85,  86,  87,  89,  90,  91,  93,
        94,  95,  97,  98,  99,  101, 102, 103, 105, 106, 107, 109, 110, 111,
        113, 114, 115, 117, 118, 119, 121, 122, 123, 125, 126, 127, 129, 130,
        131, 133, 134, 135, 137, 138, 139, 141, 142, 143, 145, 146, 147, 149,
        150, 151, 153, 154, 155, 157, 158, 159, 161, 162, 163, 165, 166, 167,
        169, 170, 171, 173, 174, 175, 177, 178, 179, 181, 182, 183, 185, 186,
        187, 189, 190, 191, 193, 194, 195, 197, 198, 199, 201, 202, 203, 205,
        206, 207, 209, 210, 211, 213, 214, 215, 217, 218, 219, 221, 222, 223,
        225, 226, 227, 229, 230, 231, 233, 234, 235, 237, 238, 239, 241, 242,
        243, 245, 246, 247, 249, 250, 251, 253, 254, 255, 257, 258, 259, 261,
        262, 263, 265, 266, 267, 269, 270, 271, 273, 274, 275, 277, 278, 279,
        281, 282, 283, 285, 286, 287, 433, 434, 435, 437, 438, 439, 441, 442,
        443, 445, 446, 447, 449, 450, 451, 453, 454, 455, 457, 458, 459, 461,
        462, 463, 465, 466, 467, 469, 470, 471, 473, 474, 475, 477, 478, 479,
        481, 482, 483, 485, 486, 487, 489, 490, 491, 493, 494, 495, 497, 498,
        499, 501, 502, 503, 505, 506, 507, 509, 510, 511, 513, 514, 515, 517,
        518, 519, 521, 522, 523, 525, 526, 527, 529, 530, 531, 533, 534, 535,
        537, 538, 539, 541, 542, 543, 545, 546, 547, 549, 550, 551, 553, 554,
        555, 557, 558, 559, 561, 562, 563, 565, 566, 567, 569, 570, 571, 573,
        574, 575, 577, 578, 579, 581, 582, 583, 585, 586, 587, 589, 590, 591,
        593, 594, 595, 597, 598, 599, 601, 602, 603, 605, 606, 607, 609, 610,
        611, 613, 614, 615, 617, 618, 619, 621, 622, 623, 625, 626, 627, 629,
        630, 631, 633, 634, 635, 637, 638, 639, 641, 642, 643, 645, 646, 647,
        649, 650, 651, 653, 654, 655, 657, 658, 659, 661, 662, 663, 665, 666,
        667, 669, 670, 671, 673, 674, 675, 677, 678, 679, 681, 682, 683, 685,
        686, 687, 689, 690, 691, 693, 694, 695, 697, 698, 699, 701, 702, 703,
        705, 706, 707, 709, 710, 711, 713, 714, 715, 717, 718, 719};

    memcpy(pbchInd, pbchInd_0, sizeof(int) * numREpbch);

  } else if (v == 1) {

    int pbchInd_1[] = {
        0,   2,   3,   4,   6,   7,   8,   10,  11,  12,  14,  15,  16,  18,
        19,  20,  22,  23,  24,  26,  27,  28,  30,  31,  32,  34,  35,  36,
        38,  39,  40,  42,  43,  44,  46,  47,  48,  50,  51,  52,  54,  55,
        56,  58,  59,  60,  62,  63,  64,  66,  67,  68,  70,  71,  72,  74,
        75,  76,  78,  79,  80,  82,  83,  84,  86,  87,  88,  90,  91,  92,
        94,  95,  96,  98,  99,  100, 102, 103, 104, 106, 107, 108, 110, 111,
        112, 114, 115, 116, 118, 119, 120, 122, 123, 124, 126, 127, 128, 130,
        131, 132, 134, 135, 136, 138, 139, 140, 142, 143, 144, 146, 147, 148,
        150, 151, 152, 154, 155, 156, 158, 159, 160, 162, 163, 164, 166, 167,
        168, 170, 171, 172, 174, 175, 176, 178, 179, 180, 182, 183, 184, 186,
        187, 188, 190, 191, 192, 194, 195, 196, 198, 199, 200, 202, 203, 204,
        206, 207, 208, 210, 211, 212, 214, 215, 216, 218, 219, 220, 222, 223,
        224, 226, 227, 228, 230, 231, 232, 234, 235, 236, 238, 239, 240, 242,
        243, 244, 246, 247, 248, 250, 251, 252, 254, 255, 256, 258, 259, 260,
        262, 263, 264, 266, 267, 268, 270, 271, 272, 274, 275, 276, 278, 279,
        280, 282, 283, 284, 286, 287, 432, 434, 435, 436, 438, 439, 440, 442,
        443, 444, 446, 447, 448, 450, 451, 452, 454, 455, 456, 458, 459, 460,
        462, 463, 464, 466, 467, 468, 470, 471, 472, 474, 475, 476, 478, 479,
        480, 482, 483, 484, 486, 487, 488, 490, 491, 492, 494, 495, 496, 498,
        499, 500, 502, 503, 504, 506, 507, 508, 510, 511, 512, 514, 515, 516,
        518, 519, 520, 522, 523, 524, 526, 527, 528, 530, 531, 532, 534, 535,
        536, 538, 539, 540, 542, 543, 544, 546, 547, 548, 550, 551, 552, 554,
        555, 556, 558, 559, 560, 562, 563, 564, 566, 567, 568, 570, 571, 572,
        574, 575, 576, 578, 579, 580, 582, 583, 584, 586, 587, 588, 590, 591,
        592, 594, 595, 596, 598, 599, 600, 602, 603, 604, 606, 607, 608, 610,
        611, 612, 614, 615, 616, 618, 619, 620, 622, 623, 624, 626, 627, 628,
        630, 631, 632, 634, 635, 636, 638, 639, 640, 642, 643, 644, 646, 647,
        648, 650, 651, 652, 654, 655, 656, 658, 659, 660, 662, 663, 664, 666,
        667, 668, 670, 671, 672, 674, 675, 676, 678, 679, 680, 682, 683, 684,
        686, 687, 688, 690, 691, 692, 694, 695, 696, 698, 699, 700, 702, 703,
        704, 706, 707, 708, 710, 711, 712, 714, 715, 716, 718, 719};

    memcpy(pbchInd, pbchInd_1, sizeof(int) * numREpbch);

  } else if (v == 2) {

    int pbchInd_2[] = {
        0,   1,   3,   4,   5,   7,   8,   9,   11,  12,  13,  15,  16,  17,
        19,  20,  21,  23,  24,  25,  27,  28,  29,  31,  32,  33,  35,  36,
        37,  39,  40,  41,  43,  44,  45,  47,  48,  49,  51,  52,  53,  55,
        56,  57,  59,  60,  61,  63,  64,  65,  67,  68,  69,  71,  72,  73,
        75,  76,  77,  79,  80,  81,  83,  84,  85,  87,  88,  89,  91,  92,
        93,  95,  96,  97,  99,  100, 101, 103, 104, 105, 107, 108, 109, 111,
        112, 113, 115, 116, 117, 119, 120, 121, 123, 124, 125, 127, 128, 129,
        131, 132, 133, 135, 136, 137, 139, 140, 141, 143, 144, 145, 147, 148,
        149, 151, 152, 153, 155, 156, 157, 159, 160, 161, 163, 164, 165, 167,
        168, 169, 171, 172, 173, 175, 176, 177, 179, 180, 181, 183, 184, 185,
        187, 188, 189, 191, 192, 193, 195, 196, 197, 199, 200, 201, 203, 204,
        205, 207, 208, 209, 211, 212, 213, 215, 216, 217, 219, 220, 221, 223,
        224, 225, 227, 228, 229, 231, 232, 233, 235, 236, 237, 239, 240, 241,
        243, 244, 245, 247, 248, 249, 251, 252, 253, 255, 256, 257, 259, 260,
        261, 263, 264, 265, 267, 268, 269, 271, 272, 273, 275, 276, 277, 279,
        280, 281, 283, 284, 285, 287, 432, 433, 435, 436, 437, 439, 440, 441,
        443, 444, 445, 447, 448, 449, 451, 452, 453, 455, 456, 457, 459, 460,
        461, 463, 464, 465, 467, 468, 469, 471, 472, 473, 475, 476, 477, 479,
        480, 481, 483, 484, 485, 487, 488, 489, 491, 492, 493, 495, 496, 497,
        499, 500, 501, 503, 504, 505, 507, 508, 509, 511, 512, 513, 515, 516,
        517, 519, 520, 521, 523, 524, 525, 527, 528, 529, 531, 532, 533, 535,
        536, 537, 539, 540, 541, 543, 544, 545, 547, 548, 549, 551, 552, 553,
        555, 556, 557, 559, 560, 561, 563, 564, 565, 567, 568, 569, 571, 572,
        573, 575, 576, 577, 579, 580, 581, 583, 584, 585, 587, 588, 589, 591,
        592, 593, 595, 596, 597, 599, 600, 601, 603, 604, 605, 607, 608, 609,
        611, 612, 613, 615, 616, 617, 619, 620, 621, 623, 624, 625, 627, 628,
        629, 631, 632, 633, 635, 636, 637, 639, 640, 641, 643, 644, 645, 647,
        648, 649, 651, 652, 653, 655, 656, 657, 659, 660, 661, 663, 664, 665,
        667, 668, 669, 671, 672, 673, 675, 676, 677, 679, 680, 681, 683, 684,
        685, 687, 688, 689, 691, 692, 693, 695, 696, 697, 699, 700, 701, 703,
        704, 705, 707, 708, 709, 711, 712, 713, 715, 716, 717, 719};

    memcpy(pbchInd, pbchInd_2, sizeof(int) * numREpbch);

  } else if (v == 3) {

    int pbchInd_3[] = {
        0,   1,   2,   4,   5,   6,   8,   9,   10,  12,  13,  14,  16,  17,
        18,  20,  21,  22,  24,  25,  26,  28,  29,  30,  32,  33,  34,  36,
        37,  38,  40,  41,  42,  44,  45,  46,  48,  49,  50,  52,  53,  54,
        56,  57,  58,  60,  61,  62,  64,  65,  66,  68,  69,  70,  72,  73,
        74,  76,  77,  78,  80,  81,  82,  84,  85,  86,  88,  89,  90,  92,
        93,  94,  96,  97,  98,  100, 101, 102, 104, 105, 106, 108, 109, 110,
        112, 113, 114, 116, 117, 118, 120, 121, 122, 124, 125, 126, 128, 129,
        130, 132, 133, 134, 136, 137, 138, 140, 141, 142, 144, 145, 146, 148,
        149, 150, 152, 153, 154, 156, 157, 158, 160, 161, 162, 164, 165, 166,
        168, 169, 170, 172, 173, 174, 176, 177, 178, 180, 181, 182, 184, 185,
        186, 188, 189, 190, 192, 193, 194, 196, 197, 198, 200, 201, 202, 204,
        205, 206, 208, 209, 210, 212, 213, 214, 216, 217, 218, 220, 221, 222,
        224, 225, 226, 228, 229, 230, 232, 233, 234, 236, 237, 238, 240, 241,
        242, 244, 245, 246, 248, 249, 250, 252, 253, 254, 256, 257, 258, 260,
        261, 262, 264, 265, 266, 268, 269, 270, 272, 273, 274, 276, 277, 278,
        280, 281, 282, 284, 285, 286, 432, 433, 434, 436, 437, 438, 440, 441,
        442, 444, 445, 446, 448, 449, 450, 452, 453, 454, 456, 457, 458, 460,
        461, 462, 464, 465, 466, 468, 469, 470, 472, 473, 474, 476, 477, 478,
        480, 481, 482, 484, 485, 486, 488, 489, 490, 492, 493, 494, 496, 497,
        498, 500, 501, 502, 504, 505, 506, 508, 509, 510, 512, 513, 514, 516,
        517, 518, 520, 521, 522, 524, 525, 526, 528, 529, 530, 532, 533, 534,
        536, 537, 538, 540, 541, 542, 544, 545, 546, 548, 549, 550, 552, 553,
        554, 556, 557, 558, 560, 561, 562, 564, 565, 566, 568, 569, 570, 572,
        573, 574, 576, 577, 578, 580, 581, 582, 584, 585, 586, 588, 589, 590,
        592, 593, 594, 596, 597, 598, 600, 601, 602, 604, 605, 606, 608, 609,
        610, 612, 613, 614, 616, 617, 618, 620, 621, 622, 624, 625, 626, 628,
        629, 630, 632, 633, 634, 636, 637, 638, 640, 641, 642, 644, 645, 646,
        648, 649, 650, 652, 653, 654, 656, 657, 658, 660, 661, 662, 664, 665,
        666, 668, 669, 670, 672, 673, 674, 676, 677, 678, 680, 681, 682, 684,
        685, 686, 688, 689, 690, 692, 693, 694, 696, 697, 698, 700, 701, 702,
        704, 705, 706, 708, 709, 710, 712, 713, 714, 716, 717, 718};

    memcpy(pbchInd, pbchInd_3, sizeof(int) * numREpbch);
  }

  /* Zero Forcing Equaliser or MMSE equaliser */
  const unsigned int zf = 0;
  if (zf == 1) {
    for (int i = 0; i < numREpbch; i++) {
      // ZF equaliser
      pbchEq[i] = pbchSymb[pbchInd[i]] / H[pbchInd[i]];
      csi[i] = abs2(H[pbchInd[i]]);
    }
  } else {
    for (int i = 0; i < numREpbch; i++) {
      // MMSE equaliser
      pbchEq[i] = (pbchSymb[pbchInd[i]] * conjf(H[pbchInd[i]])) /
                  (abs2(H[pbchInd[i]]) + nVar);
      csi[i] = abs2(H[pbchInd[i]]);
    }
  }
}

/* PBCH rate recovery */
static void pbchrateRecover(float *e, float *d) {

  /* Estimate the interleaved bits, y */
  float y[512];

  /* Average LLRs for repeated bits */
  float eAvg[352];
  for (int i = 0; i < 352; i++)
    eAvg[i] = (e[i] + e[i + 512]);

  /* Populate final LLRs */
  for (int i = 0; i < 512; i++) {
    if (i < 352) {
      y[i] = eAvg[i];
    } else {
      y[i] = e[i];
    }
  }

  /* Perform sub-block de-interleaving to obtain codewoard d */
  int P[] = {0,  1,  2,  4,  3,  5,  6,  7,  8,  16, 9,  17, 10, 18, 11, 19,
             12, 20, 13, 21, 14, 22, 15, 23, 24, 25, 26, 28, 27, 29, 30, 31};

  int J[512];
  int i;
  for (int n = 0; n < 512; n++) {
    i = (int)floor(32 * n / 512);
    J[n] = P[i] * 16 + n % 16;
    d[J[n]] = y[n];
  }
}

/* De-scramble PBCH */
static void pbchDescramble(float *llrs, unsigned int n_id_cell,
                           unsigned int i_bar_ssb, unsigned int llrLen) {

  /* Generate scrambling sequence for PBCH */
  float seq[llrLen];
  pbchPRS(n_id_cell, i_bar_ssb, seq);

  /* Descramble PBCH */
  for (int i = 0; i < llrLen; i++) {
    llrs[i] = llrs[i] * seq[i];
  }
}

/* Function to de-interleave polar code in accordance with Section 5.3.1.1 of
 * TS 38.212 */
static void pcDeIntrlv(const unsigned char *cbar, unsigned int Il,
                       const unsigned int K, unsigned char *c) {

  int pi_il_max[] = {
      0,   2,   4,   7,   9,   14,  19,  20,  24,  25,  26,  28,  31,  34,  42,
      45,  49,  50,  51,  53,  54,  56,  58,  59,  61,  62,  65,  66,  67,  69,
      70,  71,  72,  76,  77,  81,  82,  83,  87,  88,  89,  91,  93,  95,  98,
      101, 104, 106, 108, 110, 111, 113, 115, 118, 119, 120, 122, 123, 126, 127,
      129, 132, 134, 138, 139, 140, 1,   3,   5,   8,   10,  15,  21,  27,  29,
      32,  35,  43,  46,  52,  55,  57,  60,  63,  68,  73,  78,  84,  90,  92,
      94,  96,  99,  102, 105, 107, 109, 112, 114, 116, 121, 124, 128, 130, 133,
      135, 141, 6,   11,  16,  22,  30,  33,  36,  44,  47,  64,  74,  79,  85,
      97,  100, 103, 117, 125, 131, 136, 142, 12,  17,  23,  37,  48,  75,  80,
      86,  137, 143, 13,  18,  38,  144, 39,  145, 40,  146, 41,  147, 148, 149,
      150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163};

  /* Generate interleaving pattern */
  unsigned int P[K];
  unsigned int k_il_max = 164;

  if (Il == 0) {
    for (int k = 0; k < K; k++)
      P[k] = k;
  } else {
    int k = 0;
    for (int m = 0; m < k_il_max; m++) {
      if (pi_il_max[m] >= (k_il_max - K)) {
        P[k] = pi_il_max[m] - (k_il_max - K);
        k += 1;
      }
    }
  }

  /* De-interleave input sequence */
  for (int k = 0; k < K; k++)
    c[P[k]] = cbar[k];
}

static void mibDescramble(unsigned char *adash, const unsigned int n_id_cell,
                          unsigned char *a) {

  // Interleaving pattern.
  unsigned int G[] = {16, 23, 18, 17, 8,  30, 10, 6,  24, 7,  0,
                      5,  3,  2,  1,  4,  9,  11, 12, 13, 14, 15,
                      19, 20, 21, 22, 25, 26, 27, 28, 29, 31};

  /* Determine v using 3rd and 2nd LSB of SFN.
  The 3rd LSB of the SFN is the 8th bit in the SFN.
  The SFN is a 10-bit value (0-1023).*/
  unsigned int j_sfn = 7;
  unsigned int v = 2 * adash[G[j_sfn]] + adash[G[j_sfn + 1]];

  // Initialise parameters.
  unsigned int j = 0;
  unsigned int A = 32;
  unsigned int M = A - 3;

  // Generate scrambling sequence.
  const uint32_t length = v * M + M;
  float prs[length];
  nrPRS(prs, n_id_cell, length);

  // Get de-scrambling sequence
  unsigned int j_hrf = 10;
  unsigned char s[A];
  for (int i = 0; i < A; i++) {
    if (i == G[j_hrf] || i == G[j_sfn] || i == G[j_sfn + 1]) {
      s[i] = 0;
    } else {
      s[i] = (unsigned char)prs[j + v * M];
      j += 1;
    }
  }

  // De-scramble
  for (int i = 0; i < A; i++) {
    a[i] = (adash[i] + s[i]) % 2;
  }
}

static void mibDeIntrlv(unsigned char *a, unsigned char *abar) {

  // Initialise params
  unsigned int Abar = 24;    // Length of MIB payload
  unsigned int A = Abar + 8; // Total length of payload.
  unsigned int j_sfn = 0;
  unsigned int j_hrf = 10;
  unsigned int j_ssb = 11;
  unsigned int j_other = 14;

  // Interleaving pattern
  unsigned int G[] = {16, 23, 18, 17, 8,  30, 10, 6,  24, 7,  0,
                      5,  3,  2,  1,  4,  9,  11, 12, 13, 14, 15,
                      19, 20, 21, 22, 25, 26, 27, 28, 29, 31};

  // Perform de-interleaving
  for (int i = 0; i < A; i++) {
    if ((i > 0 && i <= 6) || (i >= 24 && i <= 27)) {
      abar[i] = a[G[j_sfn]];
      j_sfn += 1;
    } else if (i == 28) {
      abar[i] = a[G[j_hrf]];
    } else if (i >= 29 && i <= 31) {
      abar[i] = a[G[j_ssb]];
      j_ssb += 1;
    } else {
      abar[i] = a[G[j_other]];
      j_other += 1;
    }
  }
}

static void parseMIB(unsigned char *abar, struct ssbParam *s) {

  // System Frame Number (SFN)
  unsigned int sfn_ind[] = {1, 2, 3, 4, 5, 6, 24, 25, 26, 27};
  s->sfn = 0;
  unsigned int j = 9;
  for (int i = 0; i < 10; i++) {
    s->sfn += abar[sfn_ind[i]] * pow(2, j);
    j -= 1;
  }

  // Common sub-carrier spacing.
  s->scsCommon = (abar[7] == 0) ? 15 : 30;

  // k_ssb
  unsigned int kssb_ind[] = {8, 9, 10, 11};
  j = 3;
  for (int i = 0; i < 4; i++) {
    s->k_ssb += abar[kssb_ind[i]] * pow(2, j);
    j -= 1;
  }

  // DMRS Type A Position
  s->dmrsTypeApos = (abar[12] == 0) ? 2 : 3;

  // Coreset Zero config
  unsigned int cz_ind[] = {13, 14, 15, 16};
  j = 3;
  for (int i = 0; i < 4; i++) {
    s->coresetZero += abar[cz_ind[i]] * pow(2, j);
    j -= 1;
  }

  // Search Space Zero
  unsigned int ssz_ind[] = {17, 18, 19, 20};
  j = 3;
  for (int i = 0; i < 4; i++) {
    s->searchSpaceZero += abar[ssz_ind[i]] * pow(2, j);
    j -= 1;
  }

  // Cell barred
  s->cellBarred = abar[21];

  // Intra Freq Reselection
  s->intraFreqReselection = abar[22];
}

void ssbConfig(struct ssbConf *s, unsigned int band, float fs) {

  /* Determine SSB config (scs, pattern etc) based
   on Table 5.4.3.3-1 of TS 38.104
   TODO: Support for bands beynd n77 / 78 */
  if (band == 77 || band == 78) {
    s->subCarrierSpacing = 30e3;
    s->pattern = 'C';
    s->nfft = fs / s->subCarrierSpacing;
    int num = 1;
    float kappa = fs / 30.72e6;
    s->cpLen = 144 * kappa * pow(2, -num);
    s->symbLen = s->nfft + s->cpLen;
    s->fs = fs;
  }
}

void pbchDecode(float complex *rxSig, const unsigned int sigLen,
                struct ssbParam *s, struct ssbConf *ssb_conf) {

  // Initialise parameters
  const unsigned int ssbLen = 4 * ssb_conf->symbLen;
  const unsigned int pbchRGsize = 720;
  const unsigned int pbchUsedSC = 240;
  const unsigned int npbchRE = 432;

  // FFT Init
  struct fftParam ofdm_ft = {0};
  fftInit(&ofdm_ft, ssb_conf->nfft, 0);

  // Estimate ssbStart and n_id_2
  int ssbStart, n_id_2;
  pssCorr(rxSig, sigLen, ssb_conf, &ssbStart, &n_id_2);
  s->ssbStart = ssbStart;

  // Extract SSB from received signal
  float complex *rxSSB =
      (float complex *)malloc(sizeof(float complex) * ssbLen);
  for (int i = 0; i < ssbLen; i++) {
    rxSSB[i] = rxSig[ssbStart + i];
  }

  // Estimate and correct fractional frequency offset using CP method
  float frac_err = cpFreqSynch(rxSSB, ssbLen, ssb_conf->nfft, ssb_conf->cpLen);

  // Correct fractional frequency offset
  compFreqShift(rxSSB, rxSSB, -frac_err, ssbLen);
  float f_offs = frac_err * ssb_conf->fs;
  s->freqOffs = f_offs / 1e3;

  // Demodulate only PBCH and SSS symbols, i.e. last 3 in block.
  float complex *pbchSymbs =
      (float complex *)malloc(sizeof(float complex) * (3 * ssb_conf->nfft));
  ofdmDemod(&rxSSB[ssb_conf->symbLen], ofdm_ft, ssb_conf->cpLen, 3, pbchSymbs);

  // Estimate noise variance using null sub-carriers
  const unsigned int numAvg = 100;
  float nVar = nVarEstnull(pbchSymbs, ofdm_ft.nfft, numAvg, pbchUsedSC);

  // Extract PBCH resource grid
  float complex *rgPBCH =
      (float complex *)malloc(sizeof(float complex) * pbchRGsize);
  pbchRGExtract(pbchSymbs, ssb_conf->nfft, rgPBCH);

  // Estimate n_id_1
  unsigned int n_id_1 = sssCorr(rgPBCH, n_id_2);

  // Compute final Physical Cell ID (PCI)
  s->n_id_cell = 3 * n_id_1 + n_id_2;

  // Perform DMRS search for channel estimation and i_bar_ssb needed for
  //   descrambling of MIB
  float complex *dmrs = (float complex *)malloc(sizeof(float complex) * 144);
  int i_bar_ssb;
  dmrsSearch(rgPBCH, s->n_id_cell, dmrs, &i_bar_ssb);
  for (int i = 0; i < 8; i++) {
    s->ssb_pos[i] = (i == i_bar_ssb) ? '1' : '0';
  }
  s->ibarSSB = i_bar_ssb;

  // Perform channel estimation for PBCH
  float complex *H = (float complex *)calloc(pbchRGsize, sizeof(float complex));
  int v = s->n_id_cell % 4;
  pbchChanEst(rgPBCH, dmrs, v, H, &s->rsrp);

  // Extract data sub-carriers and perform MMSE / ZF equalisation
  float complex *pbchEq =
      (float complex *)malloc(sizeof(float complex) * npbchRE);
  float *csi = (float *)malloc(sizeof(float) * npbchRE);
  pbchEqualise(rgPBCH, H, v, pbchEq, csi, nVar);

  // Calculate the QPSK LLRs
  const unsigned int llrLen = 2 * npbchRE;
  float llrs[llrLen];
  qpskLLR(pbchEq, llrLen, nVar, llrs);

  // De-scramble PBCH
  pbchDescramble(llrs, s->n_id_cell, i_bar_ssb, llrLen);

  // Apply CSI scaling
  csiScale(llrs, csi, llrLen);

  // PBCH Rate Recovery
  unsigned int N = 512;
  float *d = (float *)malloc(sizeof(float) * N);
  pbchrateRecover(llrs, d);

  // Get frozen bit positions
  const unsigned int K = 56;
  unsigned char *frozen = getFrozen(N, K);

  // Perform polar SCD
  unsigned char *cbar = (unsigned char *)malloc(sizeof(unsigned char) * K);
  polarSCD(d, frozen, N, cbar);

  // De-interleave the decoded message
  unsigned char *c = (unsigned char *)malloc(sizeof(unsigned char) * K);
  pcDeIntrlv(cbar, 1, K, c);

  // Perform CRC check using CRC 24C.
  unsigned char poly[] = {1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1,
                          0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1};
  unsigned int polLen = 25;
  s->crcRes = crcCheck(c, poly, polLen, K);

  // De-scramble MIB.
  unsigned char adash[32];
  unsigned char a[32];
  for (int i = 0; i < 32; i++) {
    adash[i] = c[i];
  }
  mibDescramble(adash, s->n_id_cell, a);

  // De-interleave MIB.
  unsigned char abar[32] = {0};
  mibDeIntrlv(a, abar);

  // Parse MIB
  parseMIB(abar, s);

  // Free memory
  free(rxSSB);
  free(pbchSymbs);
  free(rgPBCH);
  free(dmrs);
  free(H);
  free(csi);
  free(pbchEq);
  free(frozen);
  free(d);
  free(c);
  free(cbar);
  fftDestroy(&ofdm_ft);
  fftwf_cleanup();
}
