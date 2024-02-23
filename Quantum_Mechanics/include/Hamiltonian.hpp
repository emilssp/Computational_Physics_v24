#include "tridiagonal.hpp"

class Hamiltonian{
	TridiagonalMatrix H;
	EVPsol sol;
	void solveEVP();
	int k = 250;
	vec V;
public:

	Hamiltonian();

	Hamiltonian(EVPsol sol);
	Hamiltonian(string path, string name);

	Hamiltonian(vec V_x);
	Hamiltonian(vec V_x, EVPsol sol);

	EVPsol getSol();
	void toFile(string path, string name);

	friend ostream& operator<< (ostream& os, const Hamiltonian& H);
};