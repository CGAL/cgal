function taucs_ccs_write_binary(A,filename)
% taucs_ccs_write_binary(A,filename)
% 
% writes a sparse matrix A into a binary
% file in TAUCS format

[I J V] = find(A);
[nrows ncols] = size(A);
nnz = length(V);

f = fopen(filename,'wb');

fwrite(f,nrows,'int');
fwrite(f,ncols,'int');

if isreal(V)
 flags = 2048; % double
 disp('writing double matrix');
else
 flags = 8192; % dcomplex
 disp('writing double complex matrix');
end;

fwrite(f,flags,'int');

colptr = zeros(ncols,1);

for ip=2:length(J)
  if J(ip-1) ~= J(ip) 
    colptr(J(ip),1) = ip-1;
  end;
end;
colptr(ncols+1,1) = length(J);

fwrite(f,colptr,'int');

fwrite(f,I-1,'int');

if isreal(V)
 fwrite(f,V,'double');
else
 RI = zeros(2*length(V),1);
 RI(1:2:2*nnz-1) = real(V);
 RI(2:2:2*nnz  ) = imag(V);
 fwrite(f,RI,'double');
end;

fclose(f);



  



