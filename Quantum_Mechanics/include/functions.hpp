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
EVPsol eigFromFile(string path, string name);	//Loads the eigenvectors from the txt file
void potToFile(vec V, string path, string name);
vec potFromFile(string path, string name);
