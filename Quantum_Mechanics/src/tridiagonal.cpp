#include "tridiagonal.hpp"

TridiagonalMatrix::TridiagonalMatrix()
{
	this->A = sp_mat();
}

TridiagonalMatrix::TridiagonalMatrix(vec LD, vec MD, vec UD)
{
	int N = MD.n_elem;
	sp_mat A = zeros<arma::sp_mat>(N, N);
	A.diag(0) = MD;
	A.diag(1) = UD;
	A.diag(-1) = LD;
	
	this->A = A;
}

void TridiagonalMatrix::solveEVP(EVPsol& sol, int eig_number)
{
	cout << "Solving EVP..." << endl;

	vec eigvals;
	mat eigvecs;

	eigs_sym(eigvals, eigvecs, this->A, eig_number, "sm");
	//eigvec *= -1;

	cout << "Done with eigenvalues and eigenvectors." << endl;
	sol.eigenvals = eigvals;
	sol.eigenvecs = eigvecs;
	sol.eigenvecs.insert_rows(0, 1); // Insert a row of zeros at the top
	sol.eigenvecs.insert_rows(sol.eigenvecs.n_rows, 1);

}

