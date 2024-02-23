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