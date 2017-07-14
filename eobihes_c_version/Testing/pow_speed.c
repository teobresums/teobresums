#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <time.h>


int main(){
    
//    Testing pow speed
//    
//        double num = 3.373;
//    
//        clock_t start1, end1;
//        start1 = clock();
//        double result_1 = pow(num, 8);
//        end1 = clock();
//        printf("Time required for execution 1: %f seconds.\n\n", (double)(end1-start1)/CLOCKS_PER_SEC);
//    
//        clock_t start2, end2;
//        start2 = clock();
//        double result_2 = num*num*num*num*num*num*num*num;
//        end2 = clock();
//        printf("Time required for execution 2: %f seconds.\n\n", (double)(end2-start2)/
//        CLOCKS_PER_SEC);
    
//    Output:
//    Time required for execution 1: 0.000008 seconds.
//        
//    Time required for execution 2: 0.000001 seconds.

    
//  Cube root vs pow
    double num = 3.373;

    clock_t start1, end1;
    start1 = clock();
    double result_1 = pow(num, 2./3.);
    end1 = clock();
    printf("Result1: %f\nTime required for execution 1: %f seconds.\n\n", result_1, (double)(end1-start1)/CLOCKS_PER_SEC);
    
    clock_t start2, end2;
    start2 = clock();
    double result_2 = cbrt(num*num);
    end2 = clock();
    printf("Result2: %f\nTime required for execution 2: %f seconds.\n\n", result_2, (double)(end2-start2)/CLOCKS_PER_SEC);
    
//    Output:
//    Time required for execution 1: 0.000008 seconds.
//
//    Time required for execution 2: 0.000002 seconds.

}

