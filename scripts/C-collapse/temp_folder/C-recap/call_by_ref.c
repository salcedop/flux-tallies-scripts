#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>




void main(){
  int a = 10;
  int b = 20;

  printf("before swap:\t %d\t %d",a,b);

  swap(a,b);

  printf("after swap:\t %d\t %d",a,b);
  
}
void swap(int x, int y){
  int temp = x;
  x = y
  y = temp;
}
