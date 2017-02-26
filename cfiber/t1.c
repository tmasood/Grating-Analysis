#include<stdio.h>
#include<math.h>
#include "f2c.h"

int main(int argc, char *argv[])
{
int i;
long double x = 123456789;
for (i = 0; i < 10; i++)
{
x *= x;
printf("x = %g\n", (double)x);
}
return 0;
}
