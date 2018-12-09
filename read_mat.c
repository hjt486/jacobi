/******************************************************************************
 * This file is used to read matrix from .mtx file
******************************************************************************/

double **read_mat(char *filename, int *size){
	int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    int i, j, *I, *J;
    double *val;

    if ((f = fopen(filename, "r")) == NULL) 
    {
        printf(" No file can be read\n");
        exit(1);
    }
    else
    {
        printf(" File \"%s\" readed\n", filename);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf(" Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf(" Sorry, this application does not support ");
        printf(" Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);

    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;
        J[i]--;
    }

    if (f !=stdin) fclose(f);

    double **A;
    A = (double **) malloc(M * sizeof(double *));
    for (i = 0; i < M; i += 1){
    	A[i] = (double *) malloc(M * sizeof(double));
    	for(j = 0; j < M; j += 1){
    		A[i][j] = 0;
		}
    }
    
    for (i=0; i<nz; i++){

    	A[I[i]][J[i]] = val[i];
    }
    *size = M;
    free(I);
    free(J);
    free(val);
    return A;
}
