#include <cstddef>
extern "C" void compensated(const double* a, size_t n, double* sum, double* abs_sum) {
 double s=0,c=0,ab=0,ac=0;
 for(size_t i=0;i<n;i++){double y=a[i]-c,t=s+y;c=(t-s)-y;s=t;
 double v=a[i]<0?-a[i]:a[i];y=v-ac;t=ab+y;ac=(t-ab)-y;ab=t;}
 *sum=s;*abs_sum=ab;
}
