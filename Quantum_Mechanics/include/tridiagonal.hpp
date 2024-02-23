#pragma once
#include "functions.hpp"
#include <iostream>

using namespace arma;
using namespace std;


struct TridiagonalMatrix {
	sp_mat A;
	TridiagonalMatrix();
	TridiagonalMatrix(vec LD, vec MD, vec UD);
	void solveEVP(EVPsol& sol, int eig_number);
};