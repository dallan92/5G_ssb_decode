#!/bin/bash
gcc -o cellSearch -O2 -Wall ../src/main.c ../src/nr.c ../src/dsp.c ../src/helpers.c ../src/usrp_radio.c ../src/ofdm_phy.c ../src/fec.c -lm -lfftw3f -luhd
./cellSearch
