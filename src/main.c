#include "../include/dsp.h"
#include "../include/helpers.h"
#include "../include/nr.h"
#include "../include/usrp_radio.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

int main() {

  unsigned int a = 1;
  struct ssbParam s = {0};
  struct ssbConf ssb_conf = {0};
  float complex *rxSig = NULL;

  if (a == 0) {
    float fs = 15.36e6;
    float dur = 20e-3; // 2 frames.
    float fc = 4080e6;
    float gain = 60.0;
    unsigned int sigLen = (unsigned int)(dur / (1 / fs));
    rxSig = (float complex *)malloc(sizeof(float complex) * sigLen);
    int cap = usrpIQcapture(rxSig, fs, fc, gain, sigLen);
    if (cap == 1) {
      printf("Failed to capture 5G signal from USRP!\n");
      return EXIT_FAILURE;
    }
    printf("Successfuly captured 5G signal from USRP!\n");

    /* Determine SSB config (scs, pattern, fft size) based
    on Table 5.4.3.3-1 of TS 38.104 */
    unsigned int band = 77;
    ssbConfig(&ssb_conf, band, fs);

    // Perform SSB detection + MIB decoding.
    pbchDecode(rxSig, sigLen, &s, &ssb_conf);
  } else {

    float fs = 15.36e6;
    unsigned int sigLen = 307200;
    float complex *rxSig =
        (float complex *)malloc(sizeof(float complex) * sigLen);
    readTestSig(rxSig, sigLen);

    /* Determine SSB config (scs, pattern, fft size) based
        on Table 5.4.3.3-1 of TS 38.104 */
    unsigned int band = 77;
    ssbConfig(&ssb_conf, band, fs);

    // Perform SSB detection + MIB decoding.
    pbchDecode(rxSig, sigLen, &s, &ssb_conf);
  }

  // Print results
  if (s.crcRes == 0) {
    float timeOffs = (s.ssbStart * (1 / ssb_conf.fs)) * 1000;
    printf("SSB found at sample index of %d or %0.2f ms\n", s.ssbStart,
           timeOffs);
    printf("Corrected a Frequency offset of %0.3fkHz\n", s.freqOffs);
    printf("The Physical Cell ID = %d\n", s.n_id_cell);
    printf("The SSB position bitmap = \"%s\"\n", s.ssb_pos);
    printf("The RSRP is %f dBm\n", s.rsrp);
    printf("The MIB was decoded successfully...\n");
    printf(" SFN = %d\n", s.sfn);
    printf(" Common sub-carrier spacing = %dkHz\n", s.scsCommon);
    printf(" kSSB = %d\n", s.k_ssb);
    printf(" dmrsTypeApos = %d\n", s.dmrsTypeApos);
    printf(" coresetZero = %d\n", s.coresetZero);
    printf(" searchSpaceZero = %d\n", s.searchSpaceZero);
    if (s.cellBarred == 0) {
      printf(" The cell is barred!\n");
    }
    if (s.intraFreqReselection == 0) {
      printf(" Intra Frequency cell reselection is allowed.\n");
    } else {
      printf(" Intra Frequency cell reselection is not allowed.\n");
    }
  } else {
    printf("No cell detected!\n");
  }

  // Free memory
  free(rxSig);

  return EXIT_SUCCESS;
}
