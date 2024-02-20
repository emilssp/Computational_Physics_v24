// Quantum_Mechanics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <thread>

#include <armadillo>
#include "functions.hpp"
#include "tridiagonal.hpp"
#include "WaveFunction.hpp"

using namespace std;
using namespace arma;

int main()
{
	cout << "Armadillo version: " << arma_version::as_string() << endl;


	cout<<"######################################## PROGRAM STARTS HERE ################################################"<<endl;
	arma_rng::set_seed_random();
	vec x = linspace(0.0, L, SPACESTEPS+2);

	cout <<"dx = "<< dx << endl;
	
	//##############################################################
	// //Solve an eigenvalue problem. 
	// //Uncomment this block to solve and write to file.
	// 
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	TridiagonalMatrix A = { LD,MD,UD };
	// 
	EVPsol sol;
	thread th (&TridiagonalMatrix::solveEVP, &A, ref(sol));
	th.join();
	eigToFile(sol, RAW_PATH);


	//EVPsol sol = eigFromFile(RAW_PATH);
	//vec init = sol.eigenvecs.col(0);
	//WaveFunction psi{ sol, init };
	//psi.saveToFile(RAW_PATH, "wavefunction_boring");



	// // Builds the wavefunction
	// // Uncomment to build a wavefunction with
	// // psi_init = psi_1 + psi_2 / sqrt(2) 
	//
	//EVPsol sol = eigFromFile(RAW_PATH); 
	//vec init = sqrt(0.5) * (sol.eigenvecs.col(0) + sol.eigenvecs.col(1));
	//WaveFunction psi{ sol, init};
	//psi.saveToFile(RAW_PATH, "wavefunction1");
	//WaveFunction psi{ sol, init};
	//psi.saveToFile(RAW_PATH, "wavefunction1");


	//EVPsol sol = eigFromFile(RAW_PATH);
	vec init = zeros(SPACESTEPS + 2);
	init((SPACESTEPS + 2) / 2) = 1;
	WaveFunction psi{ sol, init };
	psi.saveToFile(RAW_PATH, "wavefunction_delta");


	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
