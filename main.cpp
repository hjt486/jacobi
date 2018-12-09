/******************************************************************************
 * Texas A&M University
 * 
 * ECEN 751 Fall 2018 Final Project
 * 
 * December 9, 2018
 * 
 * Author: Jianhao Chen & Jiatai Han
 * 
 * Program      :   Linear Solver using Jacobi Method with Pthreads
 * 
 * Objective    :   Jacobi method to solve AX = b matrix system of linear equations. 
 * 
 * Input        :   .mtx file (MatrixMarket file format)
 *                  Number of Threads (Maximum 8)
 * 
 * Output       :   The solution of  Ax = b and
 *                  Other running information
 * 
 * 
 * E-mail       :   jiataihan@tamu.edu
 *                  chenjh@tamu.edu
 * 
 * Other reference:
 * mmio.c and mmio.h are from MatrixMarket (https://math.nist.gov/MatrixMarket/) to load .mtx file
 ***************************************************************************/

#include <stdio.h>
#include <pthread.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h> 
#include "mmio.c"
#include "mmio.h"
#include "read_mat.c"
#include "read_arr.c"

#define  MAX_ITERATIONS 1000000 // Define maximum iteraions while not able to converge
#define MAXTHREADS 8

double   Distance(double *X_Old, double *X_New, int matrix_size);
pthread_mutex_t  mutex1 = PTHREAD_MUTEX_INITIALIZER;

double   **Matrix_A, *RHS, *RHS_Vector;
double   *X_New, *X_Old, *Bloc_X, rno,sum;
int matrix_size;
int NumThreads,ithread;
int height;

int Number; 
void jacobi(int);

int main(int argc, char **argv)
{
    double tolerance  = 1.0E-10;
    int NoofRows, NoofCols,CLASS_SIZE,THREADS;
    double rowsum;
    double  sum;
    int irow, icol, index, Iteration,iteration,ret_count;
    double time_start, time_end,memoryused;
    struct timeval tv;
    struct timezone tz;
    char * FILENAME1;
    char * FILENAME2;
    FILE *fp;
    
    pthread_attr_t pta;
    pthread_t *threads;
    
    memoryused =0.0;
    
    printf("\n---------------------------------------------------------------------------");
    printf("\n Texas A&M University");
    printf("\n ECEN 751 Fall 2018 Final Project");
    printf("\n---------------------------------------------------------------------------");
    printf("\n Linear Solver using Jacobi Method with Pthreads\n ");
    printf("\n Jianhao Chen & Jiatai Han");
    printf("\n---------------------------------------------------------------------------\n");
    
    if( argc != 4 )
    {
        printf("\t\t Wrong Arguments\n ");
        printf("\t\t Syntax : exec <Matrix filename in .mtx)> <RHS filename in .mtx)> <Threads>\n");
        exit(-1);
    }
    else {
        FILENAME1 = argv[1];
        FILENAME2 = argv[2];
        THREADS = atoi(argv[3]);
    }
    
    if (THREADS > MAXTHREADS )
    {
        printf("\n Maximum threads is 8!\n");
        printf("Aborting ...\n");
        return(0); 
    }
    Matrix_A = read_mat(FILENAME1, &matrix_size);
    RHS = read_arr(FILENAME2, &matrix_size);
    NumThreads = THREADS;
    
    NoofRows = matrix_size; 
    NoofCols = matrix_size;
    
    // Allocate memory
    RHS_Vector = (double *) malloc(matrix_size * sizeof(double));
    
    for (irow = 0; irow < matrix_size; irow++)
    {
        RHS_Vector[irow] = RHS[irow];
    }
    
    memoryused+=(NoofRows * NoofCols * sizeof(double));
    memoryused+=(NoofRows * sizeof(double));  
    
    printf("\n");
    
    if (NoofRows != NoofCols) 
    {
        printf(" Input matrix must be square! \n");
        exit(-1);
    }
    
    // Allocate memory
    X_New = (double *) malloc(matrix_size * sizeof(double));
    memoryused+=(NoofRows * sizeof(double));
    X_Old = (double *) malloc(matrix_size * sizeof(double));
    memoryused+=(NoofRows * sizeof(double));
    Bloc_X = (double *) malloc(matrix_size * sizeof(double));
    memoryused+=(NoofRows * sizeof(double));
    
    // Get start tune
    gettimeofday(&tv, &tz);
    time_start= (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
    
    // Initailize X[i] = B[i]
    for (irow = 0; irow < matrix_size; irow++)
    {
        Bloc_X[irow] = RHS_Vector[irow];
        X_New[irow] =  RHS_Vector[irow];
    }

    // Allocate memory for threads
    threads = (pthread_t *) malloc(sizeof(pthread_t) * NumThreads);  
    
    // Initialize thread attribute
    pthread_attr_init(&pta);
    
    
    //Height is the rows that each thread should handle
    height = (int) (matrix_size-0.5)/NumThreads + 1;
    
    // Iteration begins
    Iteration = 0;
    printf(" Iterating for soluion ..... \n");
    do
    {
        for(index = 0; index < matrix_size; index++) 
            X_Old[index] = X_New[index];
        
        for(ithread=0; ithread <NumThreads;ithread++)
        {
            // Creating threads and do jacobi 
            ret_count = pthread_create(&threads[ithread], &pta,(void *(*) (void *))jacobi, (void *)(long)ithread);
            if(ret_count)
            {
                printf("\n ERROR : Return code from pthread_create() is %d ",ret_count);
                exit(-1);
            }
            
        }
        
        Iteration++;
        
        for (ithread=0; ithread<NumThreads; ithread++)
        {
            ret_count=pthread_join(threads[ithread], NULL); 
            if(ret_count)
            {
                printf("\n ERROR : Return code from pthread_join() is %d ",ret_count);
                exit(-1);
            }
        }
        
        ret_count=pthread_attr_destroy(&pta);
        
        if(ret_count)
        {
            printf("\n ERROR : Return code from pthread_attr_destroy() is %d ",ret_count);
            exit(-1);
        }
        // If reaches maximum iteraions or converges
    }while ((Iteration < MAX_ITERATIONS) && (Distance(X_Old, X_New, matrix_size) >= tolerance));

    
    // Calculating total running time
    gettimeofday(&tv, &tz);
    time_end= (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
    
    // Print out the result
    printf(" Solution :\n");
    for (irow = 0; irow < matrix_size; irow++)
    {
        printf(" X[%d] = %.10f\n", irow, X_Old[irow]);
    }

    //Show statistics
    printf("\n Matrix Size                  :  %d",matrix_size);
    printf("\n Threads                      :  %d",NumThreads);
    printf("\n Rows for each thread         :  %d\n",height);
    printf("\n Total Number Of Iterations   :  %d",Iteration);
    printf("\n Memory Consumed              :  %lf MB",(memoryused/(1024*1024)));
    printf("\n Time in  Seconds (T)         :  %lf",(time_end - time_start));
    printf("\n..........................................................................\n");
    
    // Freeing memory
    free(X_New);
    free(X_Old);
    free(Matrix_A);
    free(RHS_Vector);
    free(Bloc_X);
    free(threads);
    return 0;
}

//Calculate distance for convergence
double Distance(double *X_Old, double *X_New, int matrix_size)
{
    int             index;
    double          Sum;
    
    Sum = 0.0;
    
    for (index = 0; index < matrix_size; index++)
        Sum += (X_New[index] - X_Old[index]) * (X_New[index] - X_Old[index]);
    return (Sum);
}

void jacobi(int thread)
{
    // Calcuate rows that this thread should handle
    int firstRow, lastRow;
    firstRow = thread * height;
    lastRow = firstRow + height-1;
    lastRow = ((matrix_size-1) > lastRow) ? lastRow : (matrix_size-1);
    
    // Jacobi iteration begins
    int i,j;
    for(i = firstRow; i <= lastRow; i++)
    {
        Bloc_X[i] = RHS_Vector[i];
        //printf("\nRHS: %f\n", RHS_Vector[i]);
        j = 0;
        while (j < i) 
        {
            Bloc_X[i] -= X_Old[j] * Matrix_A[i][j];
            j += 1;
        }
        j = i+1;
        while (j < matrix_size) 
            
        {
            Bloc_X[i] -= X_Old[j] * Matrix_A[i][j];
            j += 1;
        }
        
        Bloc_X[i] = Bloc_X[i] / Matrix_A[i][i];
    }
    
    for(i = firstRow; i <= lastRow; i++)
        
    { 
        X_New[i] = Bloc_X[i];
    }
}  
