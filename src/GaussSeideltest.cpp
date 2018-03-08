// Program that tests the Gauss-Seidel method for solving system
// of equations of the form A*x = b. The matrix needs to NOT have zeroes
// on the diagonal. Convergence is only guaranteed if the matrix is either
//diagonally dominant, or symmetric and positive definite.

#include <cstdio>
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

const int Nmax = 4; // dimensions of the matrix
const int MAX_ITER = 100; // maximum number of iterations allowed. It iterates (MAX_ITER - 1) times over the 0th iteration.
const double res_threshold = 1.0E-7; // threshold of the residual to stop the iteration

// prints an array of N x N
void printarray(double arg[Nmax][Nmax]){
	for(int i=0; i<Nmax; i++){
		for(int j=0; j<Nmax; j++){
			cout << arg[i][j] << " ";
		}
		cout << endl;
	}
    cout << endl;
}

// prints a vector of dimension N
void printvector(double arg[Nmax]){
	for(int i=0; i<Nmax; i++){
		cout << arg[i] << " ";
	}
    cout << endl;
}

// prints the x vector on the correspondent iteration
void printxvector(double arg[Nmax][(MAX_ITER+1)], int iter){
	for(int i=0; i<Nmax; i++){
		cout << arg[i][iter] << " ";
	}
    cout << endl;
}

// initializes the x vector, in this case to all zeroes
void initvector(double arg[Nmax][(MAX_ITER+1)]){
	for(int i=0; i<Nmax; i++){
		arg[i][0] = 0.0;
	}
}

int main()
{
    int iter = 0; // iteration number
    double A_array[Nmax][Nmax] = {{10, -1,2,0},{-1,11,-1,3},{2,-1,10,-1},{0,3,-1,8}}; // Array that stores the coefficients of the system of equations
    double b_vector[Nmax] = {6,25,-11,15}; // Array that stores the vector b on the RHS
    double x_vector[Nmax][(MAX_ITER+1)];
    double residuals[Nmax];
    double max_residual = 1.0;
    double first_iter=1/1.0E-7;

    initvector(x_vector);
	printarray(A_array);
	printvector(b_vector);
    cout << endl;

    while ((max_residual > res_threshold) && (iter<MAX_ITER)) {

    	// Restarts maximum residual to zero at each iteration
    	max_residual = 0.0;

    	for (int i=0; i<Nmax; i++){
        	double s1=0.0;
        	double s2=0.0;
        	double sum_res_Ax=0.0;
        	residuals[i] = 0.0;

    		for (int j=0; j<Nmax; j++){
    			//This is the sum of A[i][j]*x[j] from j=0..N that it is used
    			// later to calculate the residual
    			sum_res_Ax = sum_res_Ax + A_array[i][j]*x_vector[j][iter];

    			// S1 and S2 are the sums that go in the Gauss-Seidel solution
    			// x[i]^(k+1) = ( b[i] - s1^(k+1) - s2^(k) )/A[i][i]
    			if (j<i){
    				s1=s1+A_array[i][j]*x_vector[j][iter+1];
    			}

    			if (j>i){
    				s2=s2+A_array[i][j]*x_vector[j][iter];
    			}
    		}
    		// Gauss-Seidel method of solution
        	x_vector[i][iter+1] = (b_vector[i] - s1 - s2 )/A_array[i][i];

        	// We store the residuals of x1 ... xN in this vector. It restarts for every iteration
        	residuals[i] = abs(sum_res_Ax - b_vector[i]);

        	// This variable stores the maximum residual of the vector residuals[i], by comparing each element
        	// with itself. If it is higher, it stores the value.
			max_residual = max(residuals[i],max_residual);
    	}
    	if (iter == 0){
    		// We save the first iteration value of the residual
    		//to use as the normalization factor for the the next residuals.
    		first_iter = max_residual;
    	}
    	// Value of the residual at each iteration normalized to the value
    	// of the first residual
    	max_residual = max_residual / first_iter;

    	// We print the iteration number and the values of the solution in each iteration
    	cout << "Iteration " << iter<<": ";
    	printxvector(x_vector,iter);
    	cout << endl;
    	iter++;
      	if (iter == MAX_ITER) {
      		cout << "You have reached the maximum number of iterations before converging to the set threshold" << endl;
      	}
    }

    cout << "The end" << endl;

    return 0;
}
