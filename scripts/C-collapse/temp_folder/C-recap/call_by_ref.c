#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

void swap(int*,int*);

void main(){
  int a = 10;
  int b = 20;
  printf("before swap:\t %d\t %d\n",a,b);
  swap(&a,&b);
  printf("after swap:\t %d\t %d\n",a,b);
}

void swap(int *x, int *y){
  int temp = *x;
  *x = *y;
  *y = temp;
}
