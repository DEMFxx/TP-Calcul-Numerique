/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  
  // Fill AB with the tridiagonal Poisson operator
  
  // pointeur -> variable
  int n = *la; // Dimension de AB
  int labv = *lab; // Nombre de colonnes
  int kv_local = *kv; // Nombre de superdiagonales
  int kl = 1; // Matrice tridiagonale


  // Remplir AB pour qu'il soit nulle
  for(int j=0; j<n; j++){
    for(int i=0; i<labv; i++){

      AB[i+labv*j] = 0.0;

    }
  }

  double h = 1.0 / (double)(n+1); // valeur de h
  double den = 1.0 / (h*h); // denominateur

  // Remplir AB
  for(int j=0; j<n; j++){

    // diagonal
    AB[kv_local + kl + j * labv] = 2.0 * den;

    // subdiagonal
    if(j < n-1){

      AB[(kv_local+kl+1)+j*labv] = -1.0 * den;

    }

    // superdiagonal
    if(j > 0){

      AB[(kv_local+kl-1)+j*labv] = -1.0 * den;

    }

  }

}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){

  // Only the main diagonal should have 1, all other entries are 0

  // pointeur -> variable
  int n = *la; // Dimension de AB
  int labv = *lab; // Nombre de colonnes
  int kv_local = *kv; // Nombre de superdiagonaux

  // Remplir AB pour qu'il soit nulle
  for(int j=0; j<n; j++){
    for(int i=0; i<labv; i++){

      AB[j*labv+i] = 0.0;

    }
  }

  // Remplir la diagonale de AB
  for(int j=0; j<n; j++){

    int row_diag = kv_local;
    int idx_diag = row_diag + j*labv;
    AB[idx_diag] = 1.0;

  }

}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  
  // Compute RHS vector
  
  int n = *la;
  double h = 1.0 / (double)(n+1);

  for(int i=0; i<n; i++){
    RHS[i] = 0.0;
  }

  // Conditions aux limites
  RHS[0] += (*BC0) / (h*h);
  RHS[n-1] += (*BC1) / (h*h);

}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){

  // Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D

  int n = *la;

  for(int i=0; i<n; i++){

    EX_SOL[i] = (*BC0) + X[i] * ((*BC1)-(*BC0));

  }

}  

void set_grid_points_1D(double* x, int* la){

  // Generate uniformly spaced grid points in [0,1]

  int n = *la;
  double h = 1.0 / (double)(n+1);

  for(int i=0; i<n; i++){

    x[i] = (double)(i+1) * h;

  }

}

double relative_forward_error(double* x, double* y, int* la){
  
  // Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  
  int n = *la;
  double num = 0.0;
  double den = 0.0;

  for(int i=0; i<n; i++){

    num += (x[i]-y[i]) * (x[i]-y[i]);
    den += x[i] * x[i];

  }

  return sqrt(num) / sqrt(den);

}

int indexABCol(int i, int j, int *lab){
  
  // Return the correct index formula for column-major band storage
  
  int labv = *lab;
  int ku = (labv-1) / 2;

  int row = ku + i - j;

  if(row<0 || row>=labv) return -1;

  return row + j * labv;

}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  
  // Implement specialized LU factorization for tridiagonal matrices
  
  int N = *n;
  int labv = *lab;

  // Pour qu'on ne traite que les matrices tridiagonales
  if(*kl!=1 || *ku!=1){

    *info = -1;
    return *info;

  }

  int row_diag = *kl + *ku; // position de la diagonale
  int row_sup = row_diag - 1; // superdiagonal
  int row_sub = row_diag + 1; // subdiagonal

  *info = 0;

  for(int i=0; i<N; i++){

    ipiv[i] = i+1;

  }

  for(int i=0; i<N-1; i++){

    double pivot = AB[row_diag+i*labv];

    // Pour que le pivot ne soit pas nul
    if(pivot == 0.0){

      *info = i+1;
      return *info;

    }

    double mult = AB[row_sub+i*labv] / pivot;
    AB[row_sub+i*labv] = mult;

    AB[row_diag+(i+1)*labv] -= mult * AB[row_sup+(i+1)*labv];

  }

  if(AB[row_diag+(N-1)*labv] == 0.0) *info = N;

  return *info;

}
