#pragma once
#include "functions.hpp"

void eigToFile(EVPsol sol, string path, string name)
{
	cout << "Exporting data..." << endl;

	// Write eigenvalues and eigenvectors to text files
	sol.eigenvals.save(path + "/eigenvalues_" + name + ".csv", csv_ascii);
	sol.eigenvecs.save(path + "/eigenvectors_" + name + ".csv", csv_ascii);

	cout << "Data exported to " << path << endl;
}

EVPsol eigFromFile(string path, string name)
{
	cout << "Loading eigenvalues..." << endl;
	EVPsol sol;
	vec eigvals;
	std::string eigenvalues_file = path + "/eigenvalues_" + name + ".csv";

	// Load eigenvalues
	if (!eigvals.load(eigenvalues_file)) {
		cerr << "Error loading eigenvalues file: " << eigenvalues_file << endl;
		return sol;
	}
	sol.eigenvals = eigvals;
	
	cout << "Loading eigenvectors..." << endl;

	std::string eigenvectors_file = path + "/eigenvectors_" + name + ".csv";

	// Armadillo vectors and matrix to store eigenvalues and eigenvectors
	mat eigvecs;

	// Load eigenvectors
	if (!eigvecs.load(eigenvectors_file)) {
		cerr << "Error loading eigenvectors file: " << eigenvectors_file << endl;
		return sol;
	}
	sol.eigenvecs = eigvecs;

	cout << "Finished loading data." << endl;
	return sol;
}

vec potFromFile(string path, string name) {
	cout << "Loading potential..." << endl;
	vec pot;

	std::string potential_file = path + "/potential_" + name + ".csv";

	// Load eigenvalues
	if (!pot.load(potential_file)) {
		cerr << "Error loading potential file: " << potential_file << endl;
		return pot;
	}
	return pot;
}

void potToFile(vec V, string path, string name) {
	V.save(path + "/potential_" + name + ".csv", csv_ascii);
}

double f(double lambda, double V0)
{
	if (lambda >= V0) return INFINITY;

	double k = sqrt(lambda);
	double kappa = sqrt(V0 - lambda);
	double term1 = (kappa * sin(k / 3) + k * cos(k / 3));
	double term2 = (kappa * sin(k / 3) - k * cos(k / 3));

	return exp(k / 3) * term1 * term1 - exp(-k / 3) * term2 * term2;
}

double derivative(function<double(double, double)> f, double x, double V0, const double h){
	return f(x + h, V0) - f(x - h, V0) / (2 * h);
}

double newtonRaphson(function<double(double, double)> f, double initial_guess, double V0, double tolerance, int max_iterations)
{
	double x = initial_guess;
	for (int i = 0; i < max_iterations; ++i) {
		
		double f_val = f(x, V0);
		double df_val = derivative(f,x,V0);

		if (abs(df_val) < tolerance) {
			cerr << "Derivative too small." << endl;
			return x;
		}

		double x_next = x - f_val / df_val;

		if (abs(x_next - x) < tolerance) {
			return x_next;
		}

		x = x_next;
	}
	cerr << "Max iterations reached without convergence." << endl;
	return x;
}