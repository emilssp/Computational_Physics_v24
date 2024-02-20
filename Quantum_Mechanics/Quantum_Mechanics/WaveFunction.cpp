#pragma once
#include "WaveFunction.hpp"
#include <chrono>

void WaveFunction::exponent(vec time, vec eigenvals)
{
	cout << "Calculating exponentials... \n";
	for (uword i = 0; i < eigenvals.n_rows; i++) {
		for (uword t_idx = 0; t_idx < time.n_rows; t_idx++) {
			std::complex<double> exponent(0.0, -eigenvals(i) * time(t_idx));
			this->exponential(t_idx, i) = std::exp(exponent);
		}
	}
	cout << "Calculating exponentials finished." << endl;
}

void WaveFunction::alphas(mat eigenvecs, vec init)
{
	cout << "Calculating coefficients... \n";
	for (uword i = 0; i < eigenvecs.n_cols; i++) {
		this->alpha_n.row(i) = cdot(eigenvecs.col(i), init);
	}
	cout << "Calculating coefficients finished." << endl;
}

WaveFunction::WaveFunction(string path, vec init) {

	EVPsol sol = eigFromFile("../data/raw");
	int N = sol.eigenvecs.n_cols;

	this->time = linspace(0, END_TIME, TIMESTEPS);
	this->alpha_n = zeros(N);
	this->exponential = zeros<cx_mat>(TIMESTEPS, N);
	this->wave = zeros<cx_mat>(sol.eigenvecs.n_rows, TIMESTEPS);

	thread th1(&WaveFunction::alphas, this, sol.eigenvecs, init);
	thread th2(&WaveFunction::exponent, this, this->time, sol.eigenvals);
	th1.join();
	th2.join();

	for (uword t_idx = 0; t_idx < TIMESTEPS; t_idx++) {
		cx_rowvec exp_vector = this->exponential.row(t_idx);
		for (uword n = 0; n < N; n++) {
			wave.col(t_idx) += alpha_n(n) * exp_vector(n) * sol.eigenvecs.col(n);
		}
	}

}

WaveFunction::WaveFunction(EVPsol sol, vec init) {
	
	int N = sol.eigenvecs.n_cols;

	this->time = linspace(0, END_TIME, TIMESTEPS);
	this->alpha_n = zeros(N);
	this->exponential = zeros<cx_mat>(TIMESTEPS, N);
	this->wave = zeros<cx_mat>(sol.eigenvecs.n_rows, TIMESTEPS);

	thread th1(&WaveFunction::alphas, this, sol.eigenvecs, init);
	thread th2(&WaveFunction::exponent, this, this->time, sol.eigenvals);
	th1.join();
	th2.join();
	cout << "Setting together wavefunctions..." << endl;
	auto start = std::chrono::high_resolution_clock::now();

	#pragma omp parallel for //parallelize the for loop 
	for (uword t_idx = 0; t_idx < TIMESTEPS; t_idx++) {
		cx_rowvec exp_vector = this->exponential.row(t_idx);

		//Use SIMD (Single Instruction, Multiple Data) 
		//instructions to process multiple data points with a single instruction.
		#pragma omp simd 
		for (uword n = 0; n < N; n++) {
			cx_double pref = alpha_n(n) * exp_vector(n);
			wave.col(t_idx) += pref * sol.eigenvecs.col(n);
		}
	
		wave.col(t_idx) = normalise(wave.col(t_idx));

	}
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
	cout << "Wavefunction finished." << endl;
}

void WaveFunction::saveToFile(string path, string filename)
{
	cout << "Exporting data..." << endl;

	// Write eigenvalues and eigenvectors to text files
	this->wave.save(path + "/"+filename+".csv", csv_ascii);

	cout << "Data exported to " << path << endl;
}
