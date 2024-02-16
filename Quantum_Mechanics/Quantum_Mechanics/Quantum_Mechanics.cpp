// Quantum_Mechanics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "functions.hpp"

using namespace std;
using namespace arma;

int main()
{
	cout<<"############################### PROGRAM STARTS HERE ################################################"<<endl;
	arma_rng::set_seed_random();
	vec x = linspace(0.0, L, SPACESTEPS+2);

	vec UD = zeros(SPACESTEPS - 1) - 1.00 / (dx * dx);
	vec MD = zeros(SPACESTEPS) + 2.00 / (dx * dx);
	vec LD = zeros(SPACESTEPS - 1) - 1.00 / (dx * dx);

	mat A = tridiagonal_matrix(LD, MD, UD);

	//solveEVP_to_txt(A, "../data/raw");

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
