//
//  PearsonCorrelationParallel+Serial.c
//  C
//
//  Created by Krzysztof Nalborski on 20/10/2016.
//  Copyright © 2016 Krzysztof Nalborski. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


//Global variables

int number_processes;

int process_rank;

const int n= 5000000;

clock_t begin_serial, end_serial, begin_parallel, end_parallel;

double time_spent_parallel;

double time_spent_serial;


int main(void) {


//Methods
void ParallelPCC ();
    
void SerialPCC();


//Initializing MPI
    
MPI_Init(NULL, NULL);
    
MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    
MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    
if (process_rank == 0) {

SerialPCC();
}
    
if (number_processes > 1) {

ParallelPCC();
}

MPI_Finalize();
    
}

//Serial Method


void SerialPCC () {

double *a= malloc(n * sizeof(a));
    
double *b= malloc(n * sizeof(b));

int i; //declares the "ith" elements of the arrays

double sum_a= 0; //variable to hold the result of sum of array a

double sum_b= 0; //variable to hold the result of sum of array b

double mean_a= 0; //variable to hold the result of mean- array a

double mean_b= 0; //variable to hold the result of mean- array b

double sd_asum= 0; //variable to hold the result of the sum of (a1-mean_a)^2 + (a2-mean_a)^2...

double SD_a= 0; //variable to hold SD for array a

double sd_bsum= 0; //variable to hold the result of the sum of (b1-mean_b)^2 + (b2-mean_b)^2...

double SD_b= 0; //variable to hold SD for array b

double sum_r= 0; //sum of (a[1]-mean_a) * (b[1]-mean_b)...+ (a[n]-mean_a) * (b[n]-mean_b)

double r= 0; //variable to hold the result of Pearson Correlation Coefficient

clock_t begin_serial = clock();

//Initializing Arrays

for(i= 0; i < n; i++) {

a[i]= sin(i);
b[i]= sin(i +2);
}
    
//Computing sum of array a and b

for(i= 0; i < n; i++) {
    
    sum_a= sum_a + a[i];
    sum_b= sum_b + b[i];
    
}

//Computing mean

mean_a= sum_a/ n;
mean_b= (sum_b/ n);


//Computing standard deviation of array a and b

for(i= 0; i < n; i++) {
    
    
sd_asum= sd_asum + (a[i] - mean_a) * (a[i] - mean_a);
    
sd_bsum= sd_bsum + (b[i] - mean_b) * (b[i] - mean_b);

}

SD_a= sqrt(sd_asum/n);
SD_b= sqrt(sd_bsum/n);

//Computing Pearson Correlation between the two arrays

for(i= 0; i < n; i++) {
    
sum_r= sum_r + (a[i] - mean_a) * (b[i] - mean_b);
    
}

r= (sum_r / n) / (SD_a * SD_b);

clock_t end_serial = clock();

time_spent_serial = (double)(end_serial - begin_serial) / CLOCKS_PER_SEC;

printf("\n\n##########     Serial Algorithm Results:     ##########\n\n");

printf("\n     mean_a= %lf    mean_b= %lf\n", mean_a, mean_b);

printf("\n     SD_a= %lf    SD_b= %lf\n", SD_a, SD_b);

printf("\n     Pearson Correlation Coefficient = %lf\n", r);

printf("\n     Execution time is %lf s\n", time_spent_serial);

/*printf("\n#################################     Warning:     ################################\n \nIn order to initialize Parallel Algorithm number of processes HAVE TO be at least 2.\n \n(assumption made: n is evenly divisible by the number of processes).\n \n#################################################################################\n");*/
}

//Parallel Method

void ParallelPCC () {

//Methods used
int Mean(double a[], double b[],int local_n, int local_a, int local_b,
             double *localmean_a, double *localmean_b);
    
int Standard_Deviation (double a[], double b[],int local_n, int local_a, int local_b, double totalmean_a, double totalmean_b, double *tempSD_a, double *tempSD_b);
    
int Pearson_Correlation_Coefficient (double a[], double b[], int local_n, int local_a, int local_b, double totalmean_a, double totalmean_b, double *temp_r);

//Defining variables
    
        
double *a= malloc(n * sizeof(a));

double *b= malloc(n * sizeof(b));
        
int local_a= 0; //1st element of the subset (Trapezodial Rule)
        
int local_b= 0; //last element of the subset (Trapezodial Rule)
        
int local_n= 0; //number of elements in the subset
        
double localmean_a= 0; //local mean for subset a
        
double localmean_b= 0; //local mean for subset b
        
double totalmean_a= 0; //total mean aggregated from all local means- array a
        
double totalmean_b= 0; ////total mean aggregated from all local means- array b
        
double tempSD_a= 0; //numerator in the SD equation provided- array a
        
double tempSD_b= 0; //numerator in the SD equation provided- array b
        
double SD_a= 0; //variable to hold SD of a
        
double SD_b= 0; //variable to hold SD of a
        
double temp_r= 0; //numerator in the equation provided
    
double r= 0; //variable to hold the result of Pearson Correlation Coefficient


clock_t begin_parallel= clock();

//Block Partitioning (Trapezoidal Rule)
        
local_n= n/number_processes;
        

//Determining local subsets of the arrays
        
local_a= process_rank * local_n;
local_b= local_a+ local_n-1;
    

Mean(a, b, local_n, local_a, local_b, &localmean_a, &localmean_b);
        
if (process_rank != 0) {
            
MPI_Send(&localmean_a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            
MPI_Send(&localmean_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
}
        
else {
            
totalmean_a += localmean_a;
totalmean_b += localmean_b;
            

for (int q = 1; q < number_processes; q++) {
                
                
MPI_Recv(&localmean_a, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
MPI_Recv(&localmean_b, 1, MPI_DOUBLE, q, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
totalmean_a += localmean_a;
totalmean_b += localmean_b;
                
}
            
totalmean_a= totalmean_a/number_processes;
totalmean_b= totalmean_b/number_processes;
            
}

//Broadcasting the means to all the processes
    
MPI_Bcast(&totalmean_a, 1 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&totalmean_b, 1 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
Standard_Deviation(a, b, local_n, local_a, local_b, totalmean_a, totalmean_b, &tempSD_a, &tempSD_b);
        
if (process_rank != 0) {
            
            
MPI_Send(&tempSD_a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            
MPI_Send(&tempSD_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
}
        
else {
            
SD_a += tempSD_a;
SD_b += tempSD_b;
            
            
for (int q = 1; q < number_processes; q++) {
                
MPI_Recv(&tempSD_a, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
MPI_Recv(&tempSD_b, 1, MPI_DOUBLE, q, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
SD_a += tempSD_a;
SD_b += tempSD_b;
                
}
            
SD_a=sqrt(SD_a/n);
SD_b= sqrt(SD_b/n);
            
}

Pearson_Correlation_Coefficient(a, b, local_n, local_a, local_b, totalmean_a, totalmean_b, &temp_r);
        
if (process_rank != 0) {
            
MPI_Send(&temp_r, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            
}

else {
            
r += temp_r;
            
for (int q = 1; q < number_processes; q++) {
                

MPI_Recv(&temp_r, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
r += temp_r;
                
}

clock_t end_parallel = clock();
            
time_spent_parallel = (double)(end_parallel - begin_parallel) / CLOCKS_PER_SEC;

printf("\n\n##########     Parallel Algorithm Results:     ##########\n\n");
            
printf("\n     Mean a= %lf      Mean b= %lf\n", totalmean_a, totalmean_b);
            
printf("\n     SD a= %lf        SD b= %lf\n", SD_a, SD_b);
            
r= (r / n) / (SD_a * SD_b);
        
printf("\n     Pearson Correltion Coefficient = %lf\n", r);
            
printf("\n     Execution time: %lf s \n", time_spent_parallel);
            
printf("\n     Speedup Achieved =  %lf\n", time_spent_serial /time_spent_parallel);

printf("\n     Percentage Speedup Achieved = %0.2lf%%\n\n", (1-(time_spent_parallel/time_spent_serial))*100);
// it denotes a percentage decrease of time when running the serial algorithm vs. running the prarllel algorithm
}
    
}

int Mean(double a[], double b[],int local_n, int local_a, int local_b,
         double *localmean_a, double *localmean_b) {

int i;
        
double sum_a= 0; //local sum_a
        
double sum_b= 0; //local sum_a
        

for(i= local_a; i <= local_b; i++) {
a[i]=sin(i);
b[i]=sin(i+2);
            
sum_a= sum_a + a[i];
sum_b= sum_b + b[i];
            
}
    
*localmean_a= sum_a/local_n;
*localmean_b = sum_b/local_n;
        
return 0;

}
    

int Standard_Deviation (double a[], double b[],int local_n, int local_a, int local_b, double totalmean_a, double totalmean_b, double *tempSD_a, double *tempSD_b) {

int i;
        
for(i= local_a; i <= local_b; i++) {
            
*tempSD_a += (a[i] - totalmean_a) * (a[i] - totalmean_a);
*tempSD_b += (b[i] - totalmean_b) * (b[i] - totalmean_b);

}
        
return 0;

}


int Pearson_Correlation_Coefficient (double a[], double b[], int local_n, int local_a, int local_b, double totalmean_a, double totalmean_b, double *temp_r) {

int i;
        
for(i= local_a; i <= local_b; i++) {
            
*temp_r += (a[i] - totalmean_a) * (b[i] - totalmean_b);

}

return 0;

}
    

