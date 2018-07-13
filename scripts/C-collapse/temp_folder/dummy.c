#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#define NPP 30

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
     double cpu_time_used;
     double flux_vec[NPP] = {0};
     double xs_vec[NPP] = {0};
     double results[NPP] = {0};
     int counter = 0;
     char *path = "/home/salcedop/group_xs?";
     char *tot_string;
     char string[20];
     char data[5] = "data";
     start = clock();
     
     for ( int i = 0; i < NPP; i++)
     {
       counter = counter + 1;
       flux_vec[i] = cos(counter * (2 * i)); 
       xs_vec[i] = cos(counter / (8*i +2));

     }
 
     for (int j = 0; j < NPP; j++)
     {  
         
        sprintf(string, "%d",j);
        
        /*tot_string = "data_";
        tot_string[5] = string;*/
        tot_string = concat(data,string);
        FILE *fp = fopen(tot_string,"w");
        results[j] = flux_vec[j] * xs_vec[j];
        /*fprintf(fp,"%f\n",results[j]);*/
        fclose(fp);
        //remove(tot_string);
        
     }

     end = clock();

     cpu_time_used = (end - start) / CLOCKS_PER_SEC;

     printf("loop takes %f seconds to execute \n", cpu_time_used);
}
      
