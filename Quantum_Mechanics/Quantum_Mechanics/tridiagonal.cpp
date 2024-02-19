#include "tridiagonal.hpp"

TridiagonalMatrix::TridiagonalMatrix()
{
	this->A = sp_mat();
	cout << "Empty Matrix" << endl;
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

EVPsol TridiagonalMatrix::solveEVP()
{
	cout << "Solving EVP..." << endl;

	vec eigvals;
	mat eigvecs;

	eigs_sym(eigvals, eigvecs, this->A, this->k, "sm");
	//eigvec *= -1;

	cout << "Done with eigenvalues and eigenvectors." << endl;
	EVPsol sol;
	sol.eigenvals = eigvals;
	sol.eigenvecs = eigvecs;
	return sol;
}

