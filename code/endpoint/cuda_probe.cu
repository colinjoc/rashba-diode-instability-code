#include <cuda_runtime.h>
#include <cstdio>
int main(){int n=-1;cudaError_t e=cudaGetDeviceCount(&n);printf("{\"cuda_error\":%d,\"message\":\"%s\",\"device_count\":%d}\n",(int)e,cudaGetErrorString(e),n);return e==cudaSuccess?0:2;}
