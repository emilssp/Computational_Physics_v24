#include "Hamiltonian.hpp"
#include <thread>

Hamiltonian::Hamiltonian()
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	this->V = zeros(SPACESTEPS);
	this->H = TridiagonalMatrix(LD, MD, UD);

}
Hamiltonian::Hamiltonian(EVPsol sol)
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	this->V = zeros(SPACESTEPS);
	this->H = TridiagonalMatrix(LD, MD, UD);
	this->sol = sol;
}

Hamiltonian::Hamiltonian(string path, string name)
{

	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	this->V = potFromFile(path, name);
	this->H = TridiagonalMatrix(LD, MD, UD);
	this->sol = eigFromFile(path, name);

}

Hamiltonian::Hamiltonian(vec V_in)
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	//sp_mat V = sp_mat(diagmat(V_in.subvec(1,V_in.n_elem-2)));
	sp_mat V = sp_mat(diagmat(V_in));

	this->H = TridiagonalMatrix(LD, MD, UD);
	this->H.A = this->H.A + V;
	this->V = V_in;

}

Hamiltonian::Hamiltonian(vec V_in, EVPsol sol) 
{
	vec UD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));
	vec MD = zeros(SPACESTEPS) + (2.00 / (dx * dx));
	vec LD = zeros(SPACESTEPS - 1) - (1.00 / (dx * dx));

	sp_mat V = sp_mat(diagmat(V_in));
	
	this->V = V_in;
	this->H = TridiagonalMatrix(LD, MD, UD);
	this->H.A = this->H.A + V;

	this->sol = sol;
}

void Hamiltonian::solveEVP(int k)
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

	potToFile(this->V, path, name);
}

ostream& operator<<(ostream& os, const Hamiltonian& H)
{
	os << H.H.A;
	return os;
}

vec Hamiltonian::operator*(const vec& x) const {
	return this->H.A * x;
}
mat Hamiltonian::operator*(const mat& X) const {
	return this->H.A * X;
}


double Hamiltonian::tunnelingAmp(vec g, vec e)
{
	g = g.subvec(1, g.n_elem-2);
	e = e.subvec(1, e.n_elem-2);
	return cdot(g, H * e);
}