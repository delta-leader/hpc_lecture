#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  // load the result vectors
  __m256 fxvec = _mm256_load_ps(fx);
  __m256 fyvec = _mm256_load_ps(fy);
  for(int i=0; i<N; i++) {
    // load data
    __m256 xvec = _mm256_load_ps(x);
    __m256 yvec = _mm256_load_ps(y);
    __m256 xsubvec = _mm256_set1_ps(x[i]);
    __m256 ysubvec = _mm256_set1_ps(y[i]);

    // create mask
    __m256 mask = _mm256_cmp_ps(xvec, xsubvec, _CMP_NEQ_OQ);
    __m256 zerovec = _mm256_setzero_ps();

    // use mask to avoid div 0 after subtraction
    xsubvec = _mm256_blendv_ps(zerovec, xsubvec, mask);
    ysubvec = _mm256_blendv_ps(zerovec, ysubvec, mask);
    __m256 rxvec = _mm256_sub_ps(xvec, xsubvec);
    __m256 ryvec = _mm256_sub_ps(yvec, ysubvec);
    
    __m256 rvec = _mm256_add_ps(_mm256_mul_ps(rxvec, rxvec), _mm256_mul_ps(ryvec, ryvec));
    rvec = _mm256_sqrt_ps(rvec);
    rvec = _mm256_div_ps(_mm256_set1_ps(1), rvec);
    //rvec = _mm256_rsqrt_ps(rvec);
    __m256 rmulvec = _mm256_mul_ps(rvec, _mm256_mul_ps(rvec, rvec));
    
   // load m
    __m256 mvec = _mm256_set1_ps(m[i]);
    mvec = _mm256_blendv_ps(zerovec, mvec, mask);
    rxvec = _mm256_mul_ps(rxvec, _mm256_mul_ps(mvec, rmulvec));
    ryvec = _mm256_mul_ps(ryvec, _mm256_mul_ps(mvec, rmulvec));
    fxvec = _mm256_sub_ps(fxvec, rxvec);
    fyvec = _mm256_sub_ps(fyvec, ryvec);
}
  // store the result vectors
  _mm256_store_ps(fx, fxvec);
  _mm256_store_ps(fy, fyvec);

  for(int i=0; i<N; i++) 
    printf("%d %g %g\n",i,fx[i],fy[i]);
}
