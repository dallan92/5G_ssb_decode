# 5G_ssb_decode
This software detects and decodes the SSB in 5G NR. Can be tested with captured file or USRP. Verified using Amarisoft gNB.   

## Dependencies 
libfftw3 (sudo apt install ibfftw3-dev) </br>
UHD (sudo apt-get install libuhd-dev uhd-host)

## Build code 
cd build </br>
./build.sh

## Using the software
The software is configured using a simple YAML file. You can choose whether the signal source is a USRP or file. If it is a USRP you need to specify </br>
centre frequency, sample rate, rx gain and capture duration. If it is a file, you need to specify the sample rate and length of the capture. The program 
is launched with by running ./cellSearch path_to_config_file where path_to_config_file is the location of the config file. 

## Notes
valgrind --leak-check=full --show-leak-kinds=all ./cellSearch </br> 


