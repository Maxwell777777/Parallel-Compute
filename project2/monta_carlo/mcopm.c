#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#define MAX_ITERS 100000000
double Rand(double L, double R)  
{  
    return L + (R - L) * rand() * 1.0 / RAND_MAX;  
}  
int main(int argc, char* argv[])
{
   long long i,count=0;
   long long niter = atoi(argv[1]);
   int numthreads = atoi(argv[2]);
   double x,y,pi,start_time,end_time;
	omp_set_num_threads(numthreads);
   time_t t;
   srand((unsigned) time(&t));
	start_time=omp_get_wtime();
   #pragma omp parallel for private(x,y) reduction(+:count)
   for(i=0;i<niter;i++)
   {
    x = Rand(-1, 1);  
	y = Rand(-1, 1);  
	
      if((x*x+y*y)<=1)
         count++;
   } 
   pi = count*4.0/niter;
   end_time=omp_get_wtime();
   printf("Pi:%0.8f\n",pi);
   printf("time%0.8f s\n",(end_time-start_time));
   return 0;
}
