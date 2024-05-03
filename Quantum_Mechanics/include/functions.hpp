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

void eigToFile(EVPsol sol, string path, string name); // Write eigenvalues to file
EVPsol eigFromFile(string path, string name);	// Loads the eigenvectors from the file
void potToFile(vec V, string path, string name);// Write potential to file
vec potFromFile(string path, string name);// Load potential from eigenvalues to file

cx_vec H_psi(cx_vec psi, double e0, double tau, double w, double t); //multiplies the hamiltonian with the wavefunction
double f(double lambda, double V0); //support function

uvec nonzeroColsMatrix(cx_mat A);

double newtonRaphson(function<double(double, double)> f, double initial_guess, double V0, double tolerance = 1e-10, int max_iterations = 5000);
double derivative(function<double(double, double)> f, double x, double V0, const double h = 1e-12);
cx_vec trapz(cx_mat cx_init, double eps0, double tau, double w, double h);

cx_vec solveF(cx_mat init, double eps0, double tau, double w, double h);
cx_vec thomasAlgorithm(cx_mat A, cx_vec y);
