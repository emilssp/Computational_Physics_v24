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

void eigToFile(EVPsol sol, std::string path);
EVPsol eigFromFile(std::string path);	//Loads the eigenvectors from the txt file
