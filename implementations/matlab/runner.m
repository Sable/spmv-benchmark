function runner(dim,density,normal_stdev,iterations)
% Example: runner(5000,2000,0.01,100);
% rand_csr
csr_num_rows     = dim;
csr_num_cols     = dim;
csr_density_perc = density/10000;
csr_nz_per_row   = dim * density / 1000000;
csr_num_nonzeros = round(csr_nz_per_row * dim);
csr_stddev       = normal_stdev * csr_nz_per_row;
csr_Ap           = ones(1,csr_num_rows + 1);
csr_Aj           = zeros(1,csr_num_nonzeros);

high_bound       = min(csr_num_cols, 2 * csr_nz_per_row);
used_cols        = zeros(1, csr_num_cols);
used_cols        = blanks(csr_num_cols);
update_interval  = round(csr_num_rows / 10);
if update_interval == 0
    update_interval = csr_num_rows;
end
seed = uint32(10000);
[ kn, fn, wn ] = r4_nor_setup ( );

for i = 1:csr_num_rows
    if rem(i, update_interval) == 0
        disp(sprintf('\t%d of %d (%5.1f%%) Rows Generated. Continuing...\n',i,csr_num_rows,i/csr_num_rows*100));
    end
    [ value, seed ] = r4_nor(seed, kn, fn, wn);
    nnz_ith_row_double  = value;
    nnz_ith_row_double  = nnz_ith_row_double * csr_stddev + csr_nz_per_row;
    if nnz_ith_row_double < 0
        nnz_ith_row = 0;
    elseif nnz_ith_row_double > high_bound
        nnz_ith_row = high_bound;
    else
        nnz_ith_row = round(nnz_ith_row_double);
    end
    csr_Ap(i+1) = csr_Ap(i) + nnz_ith_row;
    if csr_Ap(i+1) > csr_num_nonzeros + 1
        csr_Ap = ones(1, csr_Ap(i+1));
    end
    used_cols(:) = 0;
    for j = 1:nnz_ith_row
        rand_col = abs(gen_rand(0,csr_num_cols-1)) + 1; %rand from [1,csr_num_cols]
        while used_cols(rand_col)
          rand_col = abs(gen_rand(0,csr_num_cols-1)) + 1; %rand from [1,csr_num_cols]
        end
        csr_Aj(csr_Ap(i) + j - 1) = rand_col;
        used_cols(rand_col) = 1;
    end
    sta = csr_Ap(i);
    csr_Aj(sta:nnz_ith_row+sta-1) = sort(csr_Aj(sta:nnz_ith_row+sta-1));
end
nz_error = abs(csr_num_nonzeros - csr_Ap(csr_num_rows + 1)) / csr_num_nonzeros;
if nz_error >= 0.05
    error('WARNING: Actual NNZ differs from Theoretical NNZ by %5.2f%%!\n',nz_error*100);
end
csr_num_nonzeros = csr_Ap(csr_num_rows + 1) - 1;
csr_density_perc = csr_num_nonzeros * 100 / csr_num_cols / csr_num_rows;
csr_density_ppm  = round(csr_density_perc * 10000);
csr_Ax           = zeros(1,csr_num_nonzeros);
for i = 1:csr_num_nonzeros
    csr_Ax(i) = single(1.0 - 2.0 * commonRandomJS());
    while csr_Ax(i) == 0
        csr_Ax(i) = single(1.0 - 2.0 * commonRandomJS());
    end
end
% the end of rand_csr
vec = zeros(1,dim);
for i = 1:dim
  vec(i) = single(commonRandomJS());
end

tic
for i = 1:iterations
    res = spmv_core(dim,csr_num_rows,csr_Ap,csr_Ax,csr_Aj,vec);
end
elapsedTime = toc;

csr_Aj(1:csr_num_nonzeros) = csr_Aj(1:csr_num_nonzeros) - 1;
csr_Ap = csr_Ap - 1;

fprintf(1, '{ \"status\": %d, \"options\": \"-n %d -d %d -s %f\", \"time\": %f, \"output\": {\"row_ptr\": %d, \"col\": %d, \"val\": %d, \"x\": %d, \"y\": %d} }\n', 1, dim, density, normal_stdev, elapsedTime, floor(fletcherSum(csr_Ap)), floor(fletcherSum(csr_Aj)), floor(fletcherSum(csr_Ax)), floor(fletcherSum(vec)), floor(fletcherSum(res)));
end
