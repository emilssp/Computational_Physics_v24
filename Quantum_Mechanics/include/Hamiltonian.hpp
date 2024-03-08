#pragma once
#include "Tridiagonal.hpp"
#include "functions.hpp"

class Hamiltonian{
	TridiagonalMatrix H;
	EVPsol sol;
	vec V;
public:
	Hamiltonian();

	Hamiltonian(EVPsol sol);
	Hamiltonian(string path, string name);

	Hamiltonian(vec V_x);
	Hamiltonian(vec V_x, EVPsol sol);

	TridiagonalMatrix getH() { return this->H; }
	void solveEVP(int k = 250);
	EVPsol getSol();
	vec getV() { return this->V; };
	void toFile(string path, string name);
	//double tunnelingAmp();
	double tunnelingAmp(vec g, vec e);

	vec operator*(const vec& x) const;
	mat operator*(const mat& X) const;

	friend ostream& operator<< (ostream& os, const Hamiltonian& H);
};
