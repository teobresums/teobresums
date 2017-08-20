#include "Array.h"

void initArray(Double_Array *a, int initialSize) {
  a->array = (double *)malloc(initialSize * sizeof(double));
  a->used = 0;
  a->size = initialSize;
}

void insertArray(Double_Array *a, double element) {
  // a->used is the number of used entries, because a->array[a->used++] updates a->used only *after* the array has been accessed.
  // Therefore a->used can go up to a->size 
  if (a->used == a->size) {
    a->size *= 2;
    a->array = (double *)realloc(a->array, a->size * sizeof(double));
    if (a->array == NULL){
      fprintf(stderr, "Fail in realloc memory in insertArray.");
    }
  }
  a->array[a->used++] = element;
}

void freeArray(Double_Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

