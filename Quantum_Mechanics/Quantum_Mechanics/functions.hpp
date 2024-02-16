#pragma once
#include <armadillo>
#include "params.cpp"
#include <string>

using namespace arma;

mat tridiagonal_matrix(vec LD, vec MD, vec UD); //takes 3 vectors lower, middle, upper diagonal returns tridiagonal matrix
void solveEVP_to_txt(mat A, std::string path);	//Takes Matrix A and solve EVP. Writes the eigenvalues and vectors to a txt
mat load_eigvec_from_txt(std::string path);	//Loads the eigenvectors from the txt file
vec load_eigval_from_txt(std::string path);	//Loads the eigenvalues from the txt file
