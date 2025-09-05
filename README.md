# 5G_ssb_decode
This software detects and decodes the SSB in 5G NR. Can be tested with captured file or USRP. Verified using Amarisoft gNB.   

## Dependencies 
libfftw3 (sudo apt install ibfftw3-dev)
UHD (sudo apt-get install libuhd-dev uhd-host)

## Build + run code
cd build </br>
./run.sh

## Notes
valgrind --leak-check=full --show-leak-kinds=all ./cellSearch </br> 


