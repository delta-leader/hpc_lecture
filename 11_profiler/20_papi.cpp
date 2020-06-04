#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <immintrin.h>
#include <sys/time.h>
#include <papi.h>

void avx_matmult(float ** A, float ** B, float ** C, int N) {
#pragma omp parallel for
  for (int i=0; i<N; i++) {
    for (int k=0; k<N; k++) {
      __m256 Aik = _mm256_broadcast_ss(&A[i][k]);
      for (int j=0; j<N; j+=8) {
        __m256 Cij = _mm256_loadu_ps(&C[i][j]);
        __m256 Bkj = _mm256_loadu_ps(&B[k][j]);
        Cij = _mm256_fmadd_ps(Aik, Bkj, Cij);
        _mm256_storeu_ps(&C[i][j], Cij);
      }
    }
  }
}

void check_result(float ** A, float ** B, float ** C, int N) {
#pragma omp parallel for
  for (int i=0; i<N; i++) {
    for (int k=0; k<N; k++) {
      for (int j=0; j<N; j++) {
        C[i][j] -= A[i][k] * B[k][j];
      }
    }
  }
}

int main(int argc, char **argv) {
  int N = atoi(argv[1]);
  float ** A = new float * [N];
  float ** B = new float * [N];
  float ** C = new float * [N];
  for (int i=0; i<N; i++) {
    A[i] = new float [N];
    B[i] = new float [N];
    C[i] = new float [N];
    for (int j=0; j<N; j++) {
      A[i][j] = drand48();
      B[i][j] = drand48();
      C[i][j] = 0;
    }
  }
  int events[2] = {PAPI_L2_DCM, PAPI_L2_DCH};
  struct timeval tic, toc;
  gettimeofday(&tic, NULL);
  PAPI_start_counters(events, 2);
  avx_matmult(A, B, C, N);
  long long values[2];
  PAPI_read_counters(values, 2);
  gettimeofday(&toc, NULL);
  double time = toc.tv_sec-tic.tv_sec+(toc.tv_usec-tic.tv_usec)*1e-6;
  printf("N=%d: %lf s (%lf GFlops)\n",N,time,2.*N*N*N/time/1e9);
  gettimeofday(&tic, NULL);
  check_result(A, B, C, N);
  gettimeofday(&toc, NULL);
  time = toc.tv_sec-tic.tv_sec+(toc.tv_usec-tic.tv_usec)*1e-6;
  printf("N=%d: %lf s (%lf GFlops)\n",N,time,2.*N*N*N/time/1e9);
  float err = 0;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      err += fabs(C[i][j]);
    }
  }
  printf("error: %f\n",err/N/N);
  printf("L2 cache misses:  %lld\n", values[0]);
  printf("L2 cache hits  :  %lld\n", values[1]);
  for (int i=0; i<N; i++) {
    delete[] A[i];
    delete[] B[i];
    delete[] C[i];
  }
  delete[] A;
  delete[] B;
  delete[] C;
}
