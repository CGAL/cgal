function taucs_vec_write_binary(A,v,filename)
% taucs_vec_write_binary(A,v,filename)
% 
% writes a vector v a binary
% file in TAUCS format.
% v must conform in dimension and
% data type (real/complex) to A

[nrows ncols] = size(A);

f = fopen(filename,'wb');

if isreal(A)
 disp('writing double vector');
else
 disp('writing double complex vector');
end;

if isreal(A)
 fwrite(f,v,'double');
else
 RI = zeros(2*ncols,1);
 RI(1:2:2*ncols-1) = real(v);
 RI(2:2:2*ncols  ) = imag(v);
 fwrite(f,RI,'double');
end;

fclose(f);


