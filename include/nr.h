#ifndef NR_H
#define NR_H
#include "dsp.h"
#include <complex.h>
#include <stdint.h>

// Struct to contain the SSB parameters
struct ssbParam {
  float freqOffs;
  unsigned int ssbStart;
  unsigned int n_id_cell;
  char ssb_pos[8];
  unsigned int ibarSSB;
  unsigned char crcRes;
  unsigned int sfn;
  unsigned int scsCommon;
  unsigned int k_ssb;
  unsigned int dmrsTypeApos;
  unsigned int coresetZero;
  unsigned int searchSpaceZero;
  unsigned int cellBarred;
  unsigned int intraFreqReselection;
  float rsrp;
};

struct ssbConf {
  unsigned int nfft;
  unsigned int cpLen;
  unsigned int symbLen;
  float subCarrierSpacing;
  char pattern;
  float fs;
};

void pbchDecode(float complex *rxSig, const unsigned int sigLen,
                struct ssbParam *s, struct ssbConf *ssb_conf);
void ssbConfig(struct ssbConf *s, unsigned int band, float fs);
#endif
