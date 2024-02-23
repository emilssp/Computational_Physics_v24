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
	WaveFunction(string path, string name, vec init, double time);
	WaveFunction(EVPsol sol, vec init, double time);
	void saveToFile(string path, string filename);
};