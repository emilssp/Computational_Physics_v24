#include "Hamiltonian.hpp"
#include <thread>

Hamiltonian::Hamiltonian()
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	this->H = TridiagonalMatrix(LD, MD, UD);

	solveEVP();
}
Hamiltonian::Hamiltonian(EVPsol sol)
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	this->H = TridiagonalMatrix(LD, MD, UD);
	this->sol = sol;
}

Hamiltonian::Hamiltonian(string path)
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	this->H = TridiagonalMatrix(LD, MD, UD);
	this->sol = eigFromFile(path);
}

Hamiltonian::Hamiltonian(vec V_x)
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	sp_mat V = sp_mat(diagmat(V_x));
	this->H = TridiagonalMatrix(LD, MD, UD);
	this->H.A = this->H.A + V;

	solveEVP();
}

Hamiltonian::Hamiltonian(vec V_x, EVPsol sol) 
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	sp_mat V = sp_mat(diagmat(V_x));
	
	this->H = TridiagonalMatrix(LD, MD, UD);
	this->H.A = this->H.A + V;

	this->sol = sol;
}

Hamiltonian::Hamiltonian(vec V_x, string path)
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	sp_mat V = sp_mat(diagmat(V_x));

	this->H = TridiagonalMatrix(LD, MD, UD);
	this->H.A = this->H.A + V;

	this->sol = eigFromFile(path);
}


void Hamiltonian::solveEVP()
{
	this->H.solveEVP(this->sol, k);
}

EVPsol Hamiltonian::getSol() 
{
	return this->sol;
}

void Hamiltonian::toFile(string path, string name)
{
	eigToFile(this->sol, path, name);
}

ostream& operator<<(ostream& os, const Hamiltonian& H)
{
	os << H.H.A;
	return os;
}
