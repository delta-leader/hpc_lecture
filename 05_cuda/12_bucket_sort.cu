#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void bucket_sort(int *key, int* bucket){
  atomicAdd(&bucket[key[threadIdx.x]], 1);
  __syncthreads();
  int buck_val = 0;
  for (int i = threadIdx.x; i >= bucket[buck_val]; i-=bucket[buck_val++]);
  key[threadIdx.x] = buck_val;
}

int main() {
  int n = 50;
  int range = 5;
  int *key, *bucket;
  cudaMallocManaged(&key, n*sizeof(int));
  cudaMallocManaged(&bucket, range*sizeof(int));
  // initialize to 0
  cudaMemset(bucket, 0, range);

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  //only launching 1 block because n is small, could be extended but would need cooperative_groups for grid synchronization
  bucket_sort<<<1, n>>>(key, bucket);
  cudaDeviceSynchronize();
  

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  /*for (int i=0; i<range; i++) {
    printf("%d ",bucket[i]);
  }
  printf("\n");*/
  
  cudaFree(key);
  cudaFree(bucket);
}
