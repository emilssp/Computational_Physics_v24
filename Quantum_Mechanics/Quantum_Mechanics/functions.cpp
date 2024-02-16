#pragma once
#include "functions.hpp"

mat tridiagonal_matrix(vec LD, vec MD, vec UD){

	int N = MD.n_elem;
	mat A = zeros<arma::mat>(N, N);
	A.diag(0) = MD;
	A.diag(1) = UD;
	A.diag(-1) = LD;

	return A;
}

void solveEVP_to_txt(mat A, std::string path){
	
	vec eigval;
	mat eigvec;

	eig_sym(eigval, eigvec, A);

	eigvec *= -1;

	cout << "############################### Done with eigenvalues and vectors #####################################" << endl;

	// Write eigenvalues and eigenvectors to text files
	eigval.save(path + "/eigenvalues.txt", raw_ascii);
	eigvec.save(path + "/eigenvectors.txt", raw_ascii);

	cout << "########################################### Raw data exported ##########################################" << endl;


}

mat load_eigvec_from_txt(std::string path) {
	
	std::string eigenvectors_file = path + "/eigenvectors.txt";

	// Armadillo vectors and matrix to store eigenvalues and eigenvectors
	mat eigvecs;

	// Load eigenvectors
	if (!eigvecs.load(eigenvectors_file)) {
		std::cerr << "Error loading file: " << eigenvectors_file << std::endl;
		return mat();
	}
	return eigvecs;
}

vec load_eigval_from_txt(std::string path) {

	vec eigenvalues;
	std::string eigenvalues_file = path + "/eigenvalues.txt";

	// Armadillo vectors and matrix to store eigenvalues and eigenvectors
	vec eigvals;

	// Load eigenvalues
	if (!eigvals.load(eigenvalues_file)) {
		std::cerr << "Error loading file: " << eigenvalues_file << std::endl;
		return vec();
	}

	return eigvals;
}