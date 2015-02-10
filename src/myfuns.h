/* functions prototypes for helper functions defined in myfuns.c */

double *dmalloc(unsigned long n);
int getCount(double z, double *dat, int n);
int getSum(int *x, int n);
int *imalloc(unsigned long n);
void mcopy(double *x, double *copy, int n, int m);
void imcopy(int *x, int *copy, int n, int m);
double **pdmalloc(unsigned long n);
void printIntMatrix(int *x, int nrow, int ncol);
void printMatrix(double *y, int nrow, int ncol);

