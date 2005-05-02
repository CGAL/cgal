function v = taucs_vec_write_binary(A,filename)
% v = taucs_vec_read_binary(A,filename)
% 
% read a vector v from a binary
% file in TAUCS format.
% v must conform in dimension and
% data type (real/complex) to A

[nrows ncols] = size(A);

f = fopen(filename,'rb');

if isreal(A)
 disp('reading double vector');
else
 disp('reading double complex vector');
end;

if isreal(A)
 v = fread(f,ncols,'double');
else
 RI = zeros(2*ncols,1);
 RI = fread(f,2*ncols,'double');
 v = RI(1:2:2*ncols-1) + i*RI(2:2:2*ncols);
end;

fclose(f);



