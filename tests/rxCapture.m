% Load captured waveform
clear; close all; 
file = '~/cellSearcher/captures/IQRxSignals-PCI-1.csv';
rxSig = readsigfromcsv(file); 

% Write to file for input to cell seracher code 
fid = fopen('~/cellSearcher/tests/testSignal.txt','wt');
for i = 1:numel(rxSig)
    for j = 1:2 
        if (j == 1)
           fprintf(fid,'%f\n',real(rxSig(i)));
        else 
           fprintf(fid,'%f\n',imag(rxSig(i)));
        end 
    end
end

