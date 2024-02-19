#pragma once
#include "functions.hpp"
#include <iostream>

using namespace arma;
using namespace std;


class TridiagonalMatrix {
	sp_mat A;
	const int k = 20;

public:
	TridiagonalMatrix();
	TridiagonalMatrix(vec LD, vec MD, vec UD);
	EVPsol solveEVP();
};