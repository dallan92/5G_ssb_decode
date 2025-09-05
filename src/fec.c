#include "../include/fec.h"
#include "../include/helpers.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ---------- Forward Error Correction (FEC)  ---------- */

/* Successive Cancellation Decoder (SCD) for Polar codes */

/* f function for SCD decoding */
static void f(float *llrs, float *L, unsigned int n) {

  for (int i = 0; i < n; i++) {
    float sl = sign(llrs[i]);
    float sr = sign(llrs[i + n]);
    float la = absre(llrs[i]);
    float ra = absre(llrs[i + n]);
    float min_lr = min(la, ra);
    L[i] = sl * sr * min_lr;
  }
}

/* g function for SCD decoding */
static void g(float *llrs, float *R, unsigned int n, unsigned char **u,
              unsigned int s, unsigned int strt) {

  for (int i = 0; i < n; i++)
    R[i] = llrs[i + n] + (1 - 2 * (float)u[s][strt + i]) * llrs[i];
}

/* Function to determine frozen bit positions in polar coding */
unsigned char *getFrozen(const unsigned int N, const unsigned int K) {

  /* Polarised channel indices in ascending order of reliability */
  unsigned int q_0_n_max[] = {
      0,    1,    2,    4,   8,    16,   32,   3,    5,    64,   9,    6,
      17,   10,   18,   128, 12,   33,   65,   20,   256,  34,   24,   36,
      7,    129,  66,   512, 11,   40,   68,   130,  19,   13,   48,   14,
      72,   257,  21,   132, 35,   258,  26,   513,  80,   37,   25,   22,
      136,  260,  264,  38,  514,  96,   67,   41,   144,  28,   69,   42,
      516,  49,   74,   272, 160,  520,  288,  528,  192,  544,  70,   44,
      131,  81,   50,   73,  15,   320,  133,  52,   23,   134,  384,  76,
      137,  82,   56,   27,  97,   39,   259,  84,   138,  145,  261,  29,
      43,   98,   515,  88,  140,  30,   146,  71,   262,  265,  161,  576,
      45,   100,  640,  51,  148,  46,   75,   266,  273,  517,  104,  162,
      53,   193,  152,  77,  164,  768,  268,  274,  518,  54,   83,   57,
      521,  112,  135,  78,  289,  194,  85,   276,  522,  58,   168,  139,
      99,   86,   60,   280, 89,   290,  529,  524,  196,  141,  101,  147,
      176,  142,  530,  321, 31,   200,  90,   545,  292,  322,  532,  263,
      149,  102,  105,  304, 296,  163,  92,   47,   267,  385,  546,  324,
      208,  386,  150,  153, 165,  106,  55,   328,  536,  577,  548,  113,
      154,  79,   269,  108, 578,  224,  166,  519,  552,  195,  270,  641,
      523,  275,  580,  291, 59,   169,  560,  114,  277,  156,  87,   197,
      116,  170,  61,   531, 525,  642,  281,  278,  526,  177,  293,  388,
      91,   584,  769,  198, 172,  120,  201,  336,  62,   282,  143,  103,
      178,  294,  93,   644, 202,  592,  323,  392,  297,  770,  107,  180,
      151,  209,  284,  648, 94,   204,  298,  400,  608,  352,  325,  533,
      155,  210,  305,  547, 300,  109,  184,  534,  537,  115,  167,  225,
      326,  306,  772,  157, 656,  329,  110,  117,  212,  171,  776,  330,
      226,  549,  538,  387, 308,  216,  416,  271,  279,  158,  337,  550,
      672,  118,  332,  579, 540,  389,  173,  121,  553,  199,  784,  179,
      228,  338,  312,  704, 390,  174,  554,  581,  393,  283,  122,  448,
      353,  561,  203,  63,  340,  394,  527,  582,  556,  181,  295,  285,
      232,  124,  205,  182, 643,  562,  286,  585,  299,  354,  211,  401,
      185,  396,  344,  586, 645,  593,  535,  240,  206,  95,   327,  564,
      800,  402,  356,  307, 301,  417,  213,  568,  832,  588,  186,  646,
      404,  227,  896,  594, 418,  302,  649,  771,  360,  539,  111,  331,
      214,  309,  188,  449, 217,  408,  609,  596,  551,  650,  229,  159,
      420,  310,  541,  773, 610,  657,  333,  119,  600,  339,  218,  368,
      652,  230,  391,  313, 450,  542,  334,  233,  555,  774,  175,  123,
      658,  612,  341,  777, 220,  314,  424,  395,  673,  583,  355,  287,
      183,  234,  125,  557, 660,  616,  342,  316,  241,  778,  563,  345,
      452,  397,  403,  207, 674,  558,  785,  432,  357,  187,  236,  664,
      624,  587,  780,  705, 126,  242,  565,  398,  346,  456,  358,  405,
      303,  569,  244,  595, 189,  566,  676,  361,  706,  589,  215,  786,
      647,  348,  419,  406, 464,  680,  801,  362,  590,  409,  570,  788,
      597,  572,  219,  311, 708,  598,  601,  651,  421,  792,  802,  611,
      602,  410,  231,  688, 653,  248,  369,  190,  364,  654,  659,  335,
      480,  315,  221,  370, 613,  422,  425,  451,  614,  543,  235,  412,
      343,  372,  775,  317, 222,  426,  453,  237,  559,  833,  804,  712,
      834,  661,  808,  779, 617,  604,  433,  720,  816,  836,  347,  897,
      243,  662,  454,  318, 675,  618,  898,  781,  376,  428,  665,  736,
      567,  840,  625,  238, 359,  457,  399,  787,  591,  678,  434,  677,
      349,  245,  458,  666, 620,  363,  127,  191,  782,  407,  436,  626,
      571,  465,  681,  246, 707,  350,  599,  668,  790,  460,  249,  682,
      573,  411,  803,  789, 709,  365,  440,  628,  689,  374,  423,  466,
      793,  250,  371,  481, 574,  413,  603,  366,  468,  655,  900,  805,
      615,  684,  710,  429, 794,  252,  373,  605,  848,  690,  713,  632,
      482,  806,  427,  904, 414,  223,  663,  692,  835,  619,  472,  455,
      796,  809,  714,  721, 837,  716,  864,  810,  606,  912,  722,  696,
      377,  435,  817,  319, 621,  812,  484,  430,  838,  667,  488,  239,
      378,  459,  622,  627, 437,  380,  818,  461,  496,  669,  679,  724,
      841,  629,  351,  467, 438,  737,  251,  462,  442,  441,  469,  247,
      683,  842,  738,  899, 670,  783,  849,  820,  728,  928,  791,  367,
      901,  630,  685,  844, 633,  711,  253,  691,  824,  902,  686,  740,
      850,  375,  444,  470, 483,  415,  485,  905,  795,  473,  634,  744,
      852,  960,  865,  693, 797,  906,  715,  807,  474,  636,  694,  254,
      717,  575,  913,  798, 811,  379,  697,  431,  607,  489,  866,  723,
      486,  908,  718,  813, 476,  856,  839,  725,  698,  914,  752,  868,
      819,  814,  439,  929, 490,  623,  671,  739,  916,  463,  843,  381,
      497,  930,  821,  726, 961,  872,  492,  631,  729,  700,  443,  741,
      845,  920,  382,  822, 851,  730,  498,  880,  742,  445,  471,  635,
      932,  687,  903,  825, 500,  846,  745,  826,  732,  446,  962,  936,
      475,  853,  867,  637, 907,  487,  695,  746,  828,  753,  854,  857,
      504,  799,  255,  964, 909,  719,  477,  915,  638,  748,  944,  869,
      491,  699,  754,  858, 478,  968,  383,  910,  815,  976,  870,  917,
      727,  493,  873,  701, 931,  756,  860,  499,  731,  823,  922,  874,
      918,  502,  933,  743, 760,  881,  494,  702,  921,  501,  876,  847,
      992,  447,  733,  827, 934,  882,  937,  963,  747,  505,  855,  924,
      734,  829,  965,  938, 884,  506,  749,  945,  966,  755,  859,  940,
      830,  911,  871,  639, 888,  479,  946,  750,  969,  508,  861,  757,
      970,  919,  875,  862, 758,  948,  977,  923,  972,  761,  877,  952,
      495,  703,  935,  978, 883,  762,  503,  925,  878,  735,  993,  885,
      939,  994,  980,  926, 764,  941,  967,  886,  831,  947,  507,  889,
      984,  751,  942,  996, 971,  890,  509,  949,  973,  1000, 892,  950,
      863,  759,  1008, 510, 979,  953,  763,  974,  954,  879,  981,  982,
      927,  995,  765,  956, 887,  985,  997,  986,  943,  891,  998,  766,
      511,  988,  1001, 951, 1002, 893,  975,  894,  1009, 955,  1004, 1010,
      957,  983,  958,  987, 1012, 999,  1016, 767,  989,  1003, 990,  1005,
      959,  1011, 1013, 895, 1006, 1014, 1017, 1018, 991,  1020, 1007, 1015,
      1019, 1021, 1022, 1023};

  /* Extract N polarised channel indices in ascending order of reliability */
  int k = 0;
  unsigned int q_0_n[N];
  for (int i = 0; i < 1024; i++) {
    if (q_0_n_max[i] < N) {
      q_0_n[k] = q_0_n_max[i];
      k += 1;
    }
    if (k == N) {
      break;
    }
  }

  /* Get frozen bit indices */
  unsigned char *frozen = (unsigned char *)calloc(N, sizeof(unsigned char));
  for (int i = 0; i < (N - K); i++) {
    frozen[q_0_n[i]] = 1; // Frozen bits
  }

  /* return frozen bits */
  return frozen;
}

void polarSCD(float *llrs, unsigned char *frozen, unsigned int N,
              unsigned char *d) {

  /* Input checking */
  if ((N % 2) != 0) {
    fprintf(stderr, "N must be a power of 2.\n");
  }

  /* Matrix holding bit start indices for g operations performed at each stage
     g is performed N/2 times at each stage to decode entire codeword. */
  const unsigned int n = log2(N);
  unsigned int **gbitInds = (unsigned int **)malloc(n * sizeof(unsigned int *));
  for (int i = 0; i < n; i++)
    gbitInds[i] = (unsigned int *)calloc(N / 2, sizeof(unsigned int));
  unsigned int rep = 1;
  unsigned int step = 2;
  unsigned int gind = 0;
  unsigned int ind = 0;

  for (int i = n - 1; i >= 1; i--) {
    gind = 0;
    ind = 0;
    for (int j = 0; j < (N / (2 * rep)); j++) {
      for (int k = 0; k < rep; k++) {
        gbitInds[i][ind + k] = gind;
      }
      ind += rep;
      gind += step;
    }
    // Update step and rep
    step *= 2;
    rep *= 2;
  }

  /* Array holding index of latest g operation performed at each stage */
  unsigned int *g_ptr = (unsigned int *)calloc(n, sizeof(unsigned int));

  /* Matrix holding re-encoded bits at each stage */
  unsigned char **uhat = (unsigned char **)malloc(n * sizeof(unsigned char *));
  for (int i = 0; i < n; i++)
    uhat[i] = (unsigned char *)calloc(N, sizeof(unsigned char));

  unsigned char *c =
      (unsigned char *)calloc(N, sizeof(unsigned char)); // estimated cw
  float *ll = (float *)calloc(N / 2, sizeof(float));     // llr updates
  unsigned char *fgOP =
      (unsigned char *)calloc(n, sizeof(unsigned char)); // f/g operations
  unsigned int num = N - 1;
  unsigned int half;
  unsigned int start;

  // Main loop
  for (int i = 0; i < N; i++) {

    // f/g operaions are N-1 to 0 in n bit representation.
    dec2bin(num, n, fgOP);

    /* If bit is frozen, do not perform decoding process.
       But increment g index if a g was supposed to be perfomed
       Then jump to the next bit.
    */
    if (frozen[i] == 1) {
      for (int s = 0; s < n; s++) {
        if (fgOP[s] == 0) {
          g_ptr[s] += 1;
        }
      }
      num -= 1;
      continue;
    }

    // Start at top of decoding tree each time.
    half = N / 2;

    // For each stage.
    for (int s = 0; s < n; s++) {

      if (fgOP[s] == 1) {
        // Compute f function.
        if (s == 0) {
          f(llrs, ll, half); // stage 1
        } else if (s > 0 && s < n - 1) {
          f(ll, ll, half); // middle stages
        } else {
          f(ll, ll, half); // last stage

          // Make bit decision at last stage.
          if (ll[0] > 0) {
            c[i] = 0;
          } else {
            c[i] = 1;
          }
        }
      } else {
        // Compute g fuction
        if (s == 0) {
          // 1st stage
          start = gbitInds[0][g_ptr[0]];
          g(llrs, ll, half, uhat, 0, start);
        } else if (s > 0 && s < n - 1) {
          // Middle stages
          start = gbitInds[s][g_ptr[s]];
          g(ll, ll, half, uhat, s, start);
        } else {
          // Last stage
          start = gbitInds[n - 1][g_ptr[n - 1]];
          g(ll, ll, half, uhat, n - 1, start);

          // Make bit decision at last stage.
          if (ll[0] > 0) {
            c[i] = 0;
          } else {
            c[i] = 1;
          }
        }
        // Increment g pointer for stage s
        g_ptr[s] += 1;
      }
      half = half / 2;
    }

    // Update u_hat for next path descent
    for (int j = 0; j < N; j++)
      uhat[n - 1][j] = c[j];

    unsigned int vl = 2;
    unsigned int indx;
    for (int s = n - 1; s >= 1; s--) {
      indx = 0;
      for (int j = 0; j < (N / vl); j++) {
        for (int l = 0; l < (vl / 2); l++) {
          // uhat = [u_left + u_right, u_right]
          uhat[s - 1][l + indx] =
              (uhat[s][l + indx] + uhat[s][l + indx + (vl / 2)]) % 2;
          uhat[s - 1][l + indx + (vl / 2)] = uhat[s][l + indx + (vl / 2)];
        }
        indx += vl;
      }
      vl *= 2;
    }
    num = num - 1;
  }

  // Update estiated bit vector.
  int cnt = 0;
  for (int i = 0; i < N; i++) {
    if (frozen[i] == 0) {
      d[cnt] = c[i];
      cnt += 1;
    }
  }

  // Free memory resources
  for (int i = 0; i < n; i++) {
    free(gbitInds[i]);
    free(uhat[i]);
  }
  free(gbitInds);
  free(uhat);
  free(g_ptr);
  free(c);
  free(ll);
  free(fgOP);
}

/* QPSK LLR */
void qpskLLR(float complex *qpsk_symbs, int length, float nVar, float *llrs) {
  int j = 0;
  for (int i = 0; i < (length / 2); i++) {
    llrs[j] = (2 * crealf(qpsk_symbs[i])) / nVar;
    llrs[j + 1] = (2 * cimagf(qpsk_symbs[i])) / nVar;
    j += 2;
  }
}

/* Function to perform CSI scaling */
void csiScale(float *llrs, const float *csi, const unsigned int llrLen) {

  float *csiRep = (float *)malloc(sizeof(float) * llrLen);
  int j = 0;
  for (int i = 0; i < (llrLen / 2); i++) {
    csiRep[j] = csi[i];
    csiRep[j + 1] = csi[i];
    j += 2;
  }

  for (int i = 0; i < llrLen; i++) {
    llrs[i] = llrs[i] * csiRep[i];
  }
  free(csiRep);
}

unsigned char crcCheck(unsigned char *bits, unsigned char *poly,
                       unsigned int polLen, unsigned int bitLen) {

  // Initialise shift register
  const unsigned int shiftLen = polLen - 1;
  unsigned char *shiftreg =
      (unsigned char *)calloc(shiftLen, sizeof(unsigned char));

  // Main loop
  for (int i = 0; i < bitLen; i++) {
    // XOR input bit with leftmost bit of shift register
    unsigned char feedback = (bits[i] + shiftreg[0]) % 2;

    // Shift the register to the left
    for (int j = 0; j < (shiftLen - 1); j++) {
      shiftreg[j] = shiftreg[j + 1];
    }

    // Set last register bit to 0
    shiftreg[shiftLen - 1] = 0;

    // If feedback is 1, XOR the register with the polynomial (excluding highest
    // degree term)
    if (feedback == 1) {
      for (int k = 0; k < shiftLen; k++) {
        shiftreg[k] = (shiftreg[k] + poly[k + 1]) % 2;
      }
    }
  }

  // If shift register contains all zeros, CRC passed.
  // Otherwise, error(s) have been found.
  unsigned char e = 0;
  for (int i = 0; i < shiftLen; i++) {
    if (shiftreg[i] != 0) {
      // If any bit is non-zero, error(s) have been found in received message.
      e = 1;
      break;
    }
  }

  // Free memory
  free(shiftreg);

  return e;
}
