function DA_WriteFFT(outfile,f,psi,A)
%   writeFile(outfile,t,z)
%       
%       Writes a complex number in file "outfile" 
%       a three column file
%

fid  = fopen(outfile,'w');
lfft = length(f);
for nn=1:lfft;
    fprintf(fid, '%30.20f %30.20f %30.20f \n',[f(nn); psi(nn); A(nn)]);
end 

return    
