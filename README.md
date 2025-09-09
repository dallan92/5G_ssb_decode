# 5G_ssb_decode
This software detects and decodes the SSB in 5G NR. Can be tested with captured file or USRP. Verified using Amarisoft gNB.   

## Dependencies 
libfftw3 (sudo apt install ibfftw3-dev) </br>
UHD (sudo apt-get install libuhd-dev uhd-host)

## Build + run code
cd build </br>
./run.sh

## Using the sofwtare
This software implements all steps to detect and decode the Synchronisation Signal Block (SSB) in 5G NR. The SSB consists of the Primary Synchronisation Signal (PSS), </br>
Secondary Synchronisation Signal (SSS) and Physical Broadcast Channel (PBCH). The SSB is "always on" and is usually broadcast every 10 or 20ms. It is used for initial cell search, </br>
and cell signal measurements for processes such as handover, cell-reselection and uplink power control. </br> 

## Notes
valgrind --leak-check=full --show-leak-kinds=all ./cellSearch </br> 


