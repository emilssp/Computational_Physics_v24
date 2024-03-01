#pragma once
#include "functions.hpp"


using namespace arma;
using namespace std;

//Further plans override the arithmetic operators for the tridiagonal matrices

struct TridiagonalMatrix {
	sp_mat A;
	TridiagonalMatrix();
	TridiagonalMatrix(vec LD, vec MD, vec UD);
	void solveEVP(EVPsol& sol, int eig_number);
};

struct ComplexTridiagonalMatrix {
	sp_cx_mat A;
	ComplexTridiagonalMatrix();
	ComplexTridiagonalMatrix(TridiagonalMatrix M);
	ComplexTridiagonalMatrix(cx_vec LD, cx_vec MD, cx_vec UD);
};