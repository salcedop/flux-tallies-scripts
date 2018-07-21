#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#define NBIN 100
void test_vec();
char* concat(const char *, const char *);

int main()
{
 
 test_vec();
 return 0;
}

char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
       strcpy(result, s1);

       strcat(result, s2);
      return result;
}

void test_vec()
{
     clock_t start, end;
     long NPP=30000000;
     double cpu_time_used;
     long double flux_vec[NBIN] = {0};
     long double xs_vec[NBIN] = {0};
     long double results[NBIN] = {0};
     int counter = 0;
     char *path = "/home/salcedop/group_xs?";
     char *tot_string;
     char string[20];
     char data[5] = "data";
     int factor;
     double start_omp = omp_get_wtime();
     start = clock();
     
     for ( long i = 0; i < NPP; i++)
     {
        factor = i / NPP * NBIN ;
        if (factor > NBIN) {
          factor = NBIN-1;
         }
        
       counter = counter + 1;
       flux_vec[factor] = counter * (2 * i); 
       xs_vec[factor] = counter / (8*i +2);

     }
     
     for (long j = 0; j < NPP; j++)
     {  
         
        factor = j / NPP * NBIN ;
        if (factor > NBIN) {
          factor = NBIN-1;
         }
    
        //sprintf(string, "%d",j);
        
        /*tot_string = "data_";
        tot_string[5] = string;*/
        //tot_string = concat(data,string);
        //FILE *fp = fopen(tot_string,"w");
        results[factor] = flux_vec[factor] * xs_vec[factor];
        /*fprintf(fp,"%f\n",results[j]);*/
        //fclose(fp);
        //remove(tot_string);
        
     }
     double end_omp = omp_get_wtime();
     end = clock();
     double omp_time = (end_omp - start_omp);
     
     cpu_time_used = (end - start) / CLOCKS_PER_SEC;

     printf("loop takes %f seconds to execute \n", cpu_time_used);
     printf("loop takes %f seconds to execute(omp) \n", omp_time);
}
      
