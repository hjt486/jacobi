/******************************************************************************
 * This file is used to read right-hand-side array from .mtx file
******************************************************************************/

double *read_arr(char *filename, int *size){ 
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N;
    int i, j;
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


    if (!mm_is_real(matcode) && !mm_is_array(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }


    if ((ret_code = mm_read_mtx_array_size(f, &M, &N)) !=0)
        exit(1);


    val = (double *) malloc(M * sizeof(double));

    for (int i=0; i<M; i++)
    {
        fscanf(f, " %lg\n", &val[i]);
    }
    
    if (f !=stdin) fclose(f);

//     for (int i=0; i<M; i++)
//         fprintf(stdout, "%20.19g\n", val[i]);
    return val;
}
