#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

int add(int, int);
int multiply(int, int, int);
int multiply2(int, int);

void main(){
  int r1, r2, r3,r4,r5,a,b,c;
  int (*ptr)(int, int);
  a = 10;
  b = 20;
  c = 30;
  r1 = add(a,b);
  r2 = multiply(a,b,c);
  printf("r1: %d\nr2:%d\n", r1, r2);
  ptr = &add;
  r3 = ptr(a,b);
  r4 = multiply2(b,c);
  ptr = &multiply2;
  r5 = ptr(b,c);
  printf("using ptr-func, this is r1: %d\n",r3);
  printf("r4: %d\n",r4);
  printf("using ptr-func, this is r4: %d\n",r5);
}
int add(int x, int y){
  int z = x + y;
  return z;
}
int multiply(int x, int y, int z){
  int w = x*y*z;
  return w;
}
int multiply2(int x, int y){
  int w = x*y;
  return w;
}
