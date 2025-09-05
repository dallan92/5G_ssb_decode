#include "dsp.h"
#include "helpers.h"
#include "nr.h"
#include <complex.h>
#include <immintrin.h> /* AVX2 intrinsics */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {

  // read in llrs
  float llrs[8];
  readReTestSig(llrs, 8);

  unsigned char frozen[] = {
      1, 1, 1, 0, 1, 0, 0, 0,
  };
  unsigned char d[4];
  polarSCD(llrs, frozen, 8, d);

  for (int i = 0; i < 4; i++)
    printf("%d\n", d[i]);

  return 0;
}
