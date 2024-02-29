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

struct ComplexEVPsol {

	cx_vec eigenvals;
	cx_mat eigenvecs;
};


void eigToFile(EVPsol sol, string path, string name);
EVPsol eigFromFile(string path, string name);	//Loads the eigenvectors from the txt file
void potToFile(vec V, string path, string name);
vec potFromFile(string path, string name);

double newtonRaphson(function<double(double, double)> f, double initial_guess, double V0, double tolerance = 1e-10, int max_iterations = 5000);
double derivative(function<double(double, double)> f, double x, double V0, const double h = 1e-12);
double f(double lambda, double V0);
