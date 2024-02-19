#pragma once
#include "functions.hpp"

void eigToFile(EVPsol sol, std::string path)
{
	cout << "Exporting data..." << endl;

	// Write eigenvalues and eigenvectors to text files
	sol.eigenvals.save(path + "/eigenvalues.txt", raw_ascii);
	sol.eigenvecs.save(path + "/eigenvectors.txt", raw_ascii);

	cout << "Data exported to " << path << endl;
}

EVPsol eigFromFile(std::string path)
{
	cout << "Loading eigenvalues..." << endl;
	EVPsol sol;
	vec eigvals;
	std::string eigenvalues_file = path + "/eigenvalues.txt";

	// Load eigenvalues
	if (!eigvals.load(eigenvalues_file)) {
		std::cerr << "Error loading eigenvalues file: " << eigenvalues_file << std::endl;
		return sol;
	}
	sol.eigenvals = eigvals;
	
	cout << "Loading eigenvectors..." << endl;

	std::string eigenvectors_file = path + "/eigenvectors.txt";

	// Armadillo vectors and matrix to store eigenvalues and eigenvectors
	mat eigvecs;

	// Load eigenvectors
	if (!eigvecs.load(eigenvectors_file)) {
		std::cerr << "Error loading eigenvectors file: " << eigenvectors_file << std::endl;
		return sol;
	}
	sol.eigenvecs = eigvecs;

	cout << "Finished loading data." << endl;
	return sol;
}
