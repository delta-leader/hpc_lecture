#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
  #pragma omp parallel for
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  std::vector<int> bucket(range);
  #pragma omp parallel for 
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  #pragma omp parallel for
  for (int i=0; i<n; i++) {
    #pragma omp atomic update
    bucket[key[i]]++;
  }

  // parallel prefix sum, probably large overhead for small ranges, but a nice exercise
  std::vector<int> tmp_bucket(range);
  #pragma omp parallel
  for (int i=1; i<range; i<<=1){
    #pragma omp for
    for (int j=0; j<range; ++j)
      tmp_bucket[j] = bucket[j];
    #pragma omp for
    for (int j=i; j<range; ++j)
      bucket[j] += tmp_bucket[j-i];
  }

  // test output
  //printf("Summed bucket counts: ");
  //for (int i=0; i<range; ++i)
    //printf("%*d ", 2, bucket[i]);
  //printf("\n");

  #pragma omp parallel for
  for (int i=0; i<range; i++) {
    int j = (i==0 ? 0 : bucket[i-1]);
    int cnt = bucket[i] - j;
    // alternative for getting the bucket counts without previously summing
    //for (int k=0; k<i; ++k)
       //j += bucket[k];
    for (; cnt>0; cnt--) {
      key[j++] = i;
    }
    //printf("Thread %d, i:%d, j:%d\n", omp_get_thread_num(), i, j);
  }

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
