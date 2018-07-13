#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>
#include <string.h>
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
     long NPP = 30000000; 
     clock_t start, end;
     double cpu_time_used;
     long double flux_vec[NBIN] = {0};
     long double xs_vec[NBIN] = {0};
     long double results[NBIN] = {0};
     int counter = 0;
     char *path = "/home/salcedop/group_xs?";
     char string[20];
     FILE *fp;
     char *tot_string;
     char data[5] = "data";
     /**int chdir(const char *path);
     chdir(path);*/
     start = clock();
     
     for ( long i = 0; i < NPP; i++)
     {
        int factor = i / NPP * NBIN ;
        if (factor > NBIN) {
          factor = NBIN-1;
         }
      
       counter = counter + 1;
       flux_vec[factor] = counter * (2 * i); 
       xs_vec[factor] = counter / (8*i +2);

     }
      omp_set_num_threads(omp_get_num_procs()); 
     //omp_set_num_threads(omp_get_num_procs()); 
     #pragma omp parallel default (none) \
     shared(data,NPP,flux_vec,xs_vec,results) \
     //private(fp,string,tot_string) 
     {
     #pragma omp for schedule(dynamic,10)
     for (long  j = 0; j < NPP; j++)
     {  
       
        int factor = j / NPP * NBIN ;
        if (factor > NBIN) {
          factor = NBIN-1;
         }
        
        //printf(" %d", factor); 
        //sprintf(string, "%d",j);
        /*tot_string = "data_";
        tot_string[5] = string;*/
        //char tot_string[30] = {'d','a','t','a','_',string[0],'.','d','a','t','\0'}
        //tot_string = concat(data,string);
        /*char tot_string[20] = j;*/
        //fp = fopen(tot_string,"w");
        results[factor] = flux_vec[factor] * xs_vec[factor];
        /*fprintf(fp,"%f\n",results[j]);*/
        //fclose(fp);
        //remove(tot_string);
        
     }
     }
     end = clock();

     cpu_time_used = (end - start) / CLOCKS_PER_SEC;

     printf("loop takes %f seconds to execute \n", cpu_time_used);
}
      
