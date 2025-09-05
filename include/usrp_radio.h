#ifndef USRP_RADIO_H
#define USRP_RADIO_H
#include <complex.h>

int usrpIQcapture(float complex *rxSig, const float fs, const float fc,
                  const float gain, const unsigned int length);

#endif
