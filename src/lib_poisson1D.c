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
  int kv_local = *kv; // Nombre de superdiagonaux

  // Remplir AB pour qu'il soit nulle
  for(int j=0; j<n; j++){
    for(int i=0; i<labv; i++){

      AB[j*labv+i] = 0.0;

    }
  }

  double h = 1.0 / (double)(n+1); // valeur de h
  double den = 1.0 / (h*h); // denominateur

  // Remplir AB
  for(int j=0; j<n; j++){

    // diagonale
    int idx_diag = kv_local + j * labv;
    AB[idx_diag] = 2.0 * inv_h2;

    // subdiagonale
    if(j+1 < n){

      int row_sub = kv_local + 1;
      int idx_sub = row_sub + j * labv;
      AB[idx_sub] = -1.0 * inv_h2;

    }

    // superdiagonale
    if(j-1 >= 0){

      int row_sup = kv_local - 1;
      int idx_sup = row_sup + j * labv;
      AB[idx_sup] = -1.0 * inv_h2;

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
  


}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  return 0.0;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  return *info;
}
