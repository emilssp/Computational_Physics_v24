#pragma once
#include "functions.hpp"
#include <iostream>

using namespace arma;
using namespace std;


class TridiagonalMatrix {
	sp_mat A;

	const int k = 200;

public:
	TridiagonalMatrix();
	TridiagonalMatrix(vec LD, vec MD, vec UD);
	void solveEVP(EVPsol& sol);
};