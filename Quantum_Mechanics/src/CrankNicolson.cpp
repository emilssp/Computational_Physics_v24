#include "CrankNicolson.hpp"

fwEuler::fwEuler(cx_vec init, Hamiltonian H, const double delta_t, const double delta_x, const unsigned int timesteps, const unsigned int spacesteps)
{
	ComplexTridiagonalMatrix A = H.getH();
	this->time = linspace(0, timesteps * delta_t, timesteps);

	this->res = cx_mat().zeros(spacesteps, timesteps);
	this->res.col(0) = init;

	A.A = speye<sp_cx_mat>(spacesteps, spacesteps) + complex<double>(0,dt)*A.A;
	for (int i = 1; i < time.n_elem; i++) 
	{
		
	}
}
