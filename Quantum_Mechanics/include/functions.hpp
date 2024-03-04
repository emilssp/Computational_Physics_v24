#pragma once
#include <armadillo>
#include <iostream>
#include "params.cpp"
#include <string>
//#include <cmath>


using namespace std;
using namespace arma;

struct EVPsol {

	vec eigenvals;
	mat eigenvecs;
};

void eigToFile(EVPsol sol, string path, string name);
EVPsol eigFromFile(string path, string name);	//Loads the eigenvectors from the file
void potToFile(vec V, string path, string name);
vec potFromFile(string path, string name);

cx_vec H_psi(cx_vec psi, double e0, double tau, double w, double t);
double f(double lambda, double V0);

double newtonRaphson(function<double(double, double)> f, double initial_guess, double V0, double tolerance = 1e-10, int max_iterations = 5000);
double derivative(function<double(double, double)> f, double x, double V0, const double h = 1e-12);
cx_vec extendedSimpsonsRule(vec init, double e0, double tau, double w, double start, double stop, int n);
cx_vec solveF(vec init, double e0, double tau, double w, double start, double stop, int n);
cx_vec thomasAlgorithm(cx_mat A, cx_vec y);

double tau_func(vec g, vec e, sp_mat H);