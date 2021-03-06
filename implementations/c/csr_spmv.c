/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014, Erick Lavoie, Faiz Khan, Sujay Kathrotia, Vincent
 * Foley-Bourgon, Laurie Hendren
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef SERIAL
#define SERIAL
#endif

#include "sparse_formats.h"
#include "./common/common.h"
#include "./common/common_rand.h"
#include <getopt.h>
#include <stdlib.h>

func_ret_t create_vector_from_random(float **vp, int size) {
  float *v;
  int i;

  srand(time(NULL));

  v = (float *) malloc(size*sizeof(float));
  if( v == NULL){
    return RET_FAILURE;
  }

  for(i = 0; i < size; ++i){
    v[i] = common_randJS();
  }

  *vp = v;
  return RET_SUCCESS;
}
/**
 * Sparse Matrix-Vector Multiply
 *
 * Multiplies csr matrix by vector x, adds vector y, and stores output in vector out
 */
void spmv_csr_cpu(const csr_matrix* csr,const float* x,const float* y,float* out)
{
    unsigned int row,row_start,row_end,jj;
    float sum = 0.0;
    for(row=0; row < csr->num_rows; row++)
    {
        sum = y[row];
        row_start = csr->Ap[row];
        row_end   = csr->Ap[row+1];

        for (jj = row_start; jj < row_end; jj++){
            sum += csr->Ax[jj] * x[csr->Aj[jj]];
        }
        out[row] = sum;
    }
}

static struct option long_options[] = {
    /* name, has_arg, flag, val */
    {"stddev", 1, NULL, 's'},
    {"density", 1, NULL, 'd'},
    {"size", 1, NULL, 'n'},
    {"iterations", 1, NULL, 'i'},
    {0,0,0,0}
};

int main(int argc, char *argv[]){
    int opt, option_index=0;
    unsigned int dim=1024, density=5000;
    double normal_stdev=0.01;
    unsigned int seed = 10000;
    float *v;
    stopwatch sw;
    unsigned int iterations = 1;
    int i;

    while ((opt = getopt_long(argc, argv, "s:d:n:i:", long_options, &option_index)) != -1){
        switch(opt){
        case 's':
            normal_stdev = atof(optarg);
            break;
        case 'd':
            density =  atoi(optarg);
            break;
        case 'n':
            dim  = atoi(optarg);
            break ;
        case 'i':
            iterations = atoi(optarg);
            break ;
        default:
            fprintf(stderr, "Usage: %s [-s stddev] [-d density] [-n dimension]", argv[0]);
            break;
        }
    }

    float *sum = calloc(dim, sizeof(float));
    float *result = calloc(dim, sizeof(float));
    //memset(sum, 0.0, sizeof(sum));
    //memset(result, 0.0, sizeof(result));

    csr_matrix sm =  rand_csr(dim, density, normal_stdev, &seed, stderr);
    create_vector_from_random(&v, dim);

    stopwatch_start(&sw);
    for(i=0; i< iterations; ++i) spmv_csr_cpu(&sm,v,sum, result);
    stopwatch_stop(&sw);

    int Ajlen = (double)dim*dim*density/1000000.0;

    printf("{ \"status\": %d, \"options\": \"-n %d -d %d -s %f\", \"time\": %f, \"output\": {\"row_ptr\": %d, \"col\": %d, \"val\": %d, \"x\": %d, \"y\": %d} }\n", 1, dim, density, normal_stdev, get_interval_by_sec(&sw), fletcher_sum_1d_array_unsigned_int(sm.Ap, dim+1), fletcher_sum_1d_array_unsigned_int(sm.Aj, Ajlen), fletcher_sum_1d_array_float(sm.Ax, sm.num_nonzeros), fletcher_sum_1d_array_float(v, dim), fletcher_sum_1d_array_float(result, dim));

    free(sum);
    free(result);
    free(v);

}
