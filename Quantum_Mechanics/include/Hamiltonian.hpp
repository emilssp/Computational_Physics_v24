#include "tridiagonal.hpp"

class Hamiltonian{
	TridiagonalMatrix H;
	EVPsol sol;
	void solveEVP();
	int k = 250;
public:

	Hamiltonian();
	Hamiltonian(EVPsol sol);
	Hamiltonian(string path);

	Hamiltonian(vec V_x);
	Hamiltonian(vec V_x, EVPsol sol);
	Hamiltonian(vec V_x, string path);

	EVPsol getSol();
	void toFile(string path, string name);

	friend ostream& operator<< (ostream& os, const Hamiltonian& H);
};