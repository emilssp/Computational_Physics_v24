#pragma once

#include "Hamiltonian.hpp"

struct fwEuler {
	cx_mat res;
	vec time;
	fwEuler(cx_vec init, Hamiltonian H, const double delta_t, const double delta_x, const unsigned int timesteps, const unsigned int spacesteps);

};