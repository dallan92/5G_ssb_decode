function fid = writeTestSig(vec)
fid = fopen('testSignal.txt', 'w');
fprintf(fid, '%f\n', vec);
fclose(fid);
end