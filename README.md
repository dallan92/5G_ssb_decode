# 5G_ssb_decode
This software detects and decodes the SSB in 5G NR. Can be tested with captured file or USRP. Verified using Amarisoft gNB.   

## Dependencies 
libfftw3 (sudo apt install ibfftw3-dev) </br>
UHD (sudo apt-get install libuhd-dev uhd-host)

## Build code 
cd build </br>
./build.sh

## Using the software
The software is configured using a simple YAML file. You can choose whether the signal source is the USRP or a file. If it is a USRP, you need to specify centre frequency, sample rate, rx gain and capture duration. If it is a file, you only need to specify the sample rate and length of the capture. </br>

The program is launched  by running ./cellSearch path_to_config_file path_to_test_signal if you are reading the signal from a file. If you are using the USRP, it is simply launched with ./cellSearch path_to_config_file.
