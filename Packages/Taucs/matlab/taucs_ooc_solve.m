function x = taucs_ccs_solve(A,b)

taucs_ccs_write_binary(A,'/tmp/matlab-taucs-A.bin');
taucs_vec_write_binary(A,b,'/tmp/matlab-taucs-b.bin');

disp('calling taucs to factor and solve');

unix('../bin/linux/ooc_factor_solve /tmp/matlab-taucs')

disp('taucs returned');

x = taucs_vec_read_binary(A,'/tmp/matlab-taucs-x.bin');

% !/bin/rm /tmp/matlab-taucs-*

relative_residual_norm = norm(A*x-b)/norm(b)
