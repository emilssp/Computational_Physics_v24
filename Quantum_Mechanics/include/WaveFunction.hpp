#pragma once
#include "functions.hpp"
#include <complex>

class WaveFunction {
	
	cx_mat wave;

	vec alpha_n;
	cx_mat exponential;
	
	vec time;

	void alphas(mat eigenvecs, vec init);
	void exponent(vec time, vec eigenvals);

public:
	WaveFunction(string path, vec init);
	WaveFunction(EVPsol sol, vec init);
	void saveToFile(string path, string filename);
};