#include <stdio.h>
#include <stdlib.h>

int main(){
  double a = 3.14159265358979;

  FILE *fp;
  fp = fopen("test.txt", "w+");
  fprintf(fp, "Printing double a (f)    = %f\n", a);
  fprintf(fp, "Printing double a (.17f) = %.17f\n", a);
  fclose(fp);

  return 0;

  //result:
  //Printing double a (f)    = 3.141593
  //Printing double a (.17f) = 3.14159265358979001
}

