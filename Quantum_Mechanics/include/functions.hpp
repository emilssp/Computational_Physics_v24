#pragma once
#include <armadillo>
#include <iostream>
#include "params.cpp"
#include <string>

using namespace std;
using namespace arma;

struct EVPsol {

	vec eigenvals;
	mat eigenvecs;
};

void eigToFile(EVPsol sol, string path, string name);
EVPsol eigFromFile(string path);	//Loads the eigenvectors from the txt file
