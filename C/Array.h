/**
 * Copyright (C) 2017 Alessandro Nagar, Gregorio Carullo, Ka Wa Tsang, Philipp Fleig, Sebastiano Bernuzzi, Walter Del Pozzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#ifndef _Array_h
#define _Array_h

#include <stdlib.h>
#include <stdio.h>

typedef struct {
  double *array;
  int used;
  int size;
} Double_Array;

void initArray(Double_Array *a, int initialSize);
void insertArray(Double_Array *a, double element);
void freeArray(Double_Array *a);


//void initArray(Double_Array *a, int initialSize) {
//  a->array = (double *)malloc(initialSize * sizeof(double));
//  a->used = 0;
//  a->size = initialSize;
//}
//
//void insertArray(Double_Array *a, double element) {
//  // a->used is the number of used entries, because a->array[a->used++] updates a->used only *after* the array has been accessed.
//  // Therefore a->used can go up to a->size 
//  if (a->used == a->size) {
//    a->size += 10;
//    a->array = (double *)realloc(a->array, a->size * sizeof(double));
//    if (a->array == NULL){
//      fprintf(stderr, "Fail in realloc memory in insertArray.");
//    }
//  }
//  a->array[a->used++] = element;
//}
//
//void freeArray(Double_Array *a) {
//  free(a->array);
//  a->array = NULL;
//  a->used = a->size = 0;
//}

#endif /* _Array_h */

// Example Usage :
// Array a;
// int i;
//
// initArray(&a, 5);  // initially 5 elements
// for (i = 0; i < 100; i++)
//   insertArray(&a, i);  // automatically resizes as necessary
//   printf("%d\n", a.array[9]);  // print 10th element
//   printf("%d\n", a.used);  // print number of elements
//   freeArray(&a);
