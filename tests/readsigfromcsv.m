function rxSig = readsigfromcsv(file)

sig_iq = csvread(file); 
rxSig = zeros(size(sig_iq,1),1); 

for i = 1:size(sig_iq,1)
    rxSig(i) = sig_iq(i,1) + 1i*sig_iq(i,2);   
end

end