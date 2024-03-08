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


/*double Hamiltonian::tunnelingAmp()
{
	vec g = this->sol.eigenvecs.col(1);

	vec e = this->sol.eigenvecs.col(1);

	double lambda = this->sol.eigenvals(1);

	double res = g.at(0) * lambda * e.at(0) + g.at(g.n_elem - 1) * lambda * e.at(e.n_elem-1);
	// Sum for terms with coefficient 4
	for (int i = 1; i < SPACESTEPS; i += 2) {
		res = 4.00 * g.at(i) * lambda * e.at(i);
	}

	// Sum for terms with coefficient 2
	for (int i = 2; i < SPACESTEPS - 1; i += 2) {
	 	res = 2.00 * g.at(i) * lambda * e.at(i);
	}

	return (dx / 3) * res;
}*/

double Hamiltonian::tunnelingAmp(vec g, vec e)
{
	return cdot(g,this->H.A*e);
}
