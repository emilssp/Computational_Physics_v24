// Quantum_Mechanics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <thread>
#include <fstream>


#include <armadillo>

#include "functions.hpp"
#include "Tridiagonal.hpp"
#include "TwoLevelHamiltonian.hpp"

#include "Hamiltonian.hpp"
#include "WaveFunction.hpp"
#include "CrankNicolson.hpp"


using namespace std;
using namespace arma;

int main()
{
	cout << "Armadillo version: " << arma_version::as_string() << endl;


	cout<<"######################################## PROGRAM STARTS HERE ################################################"<<endl;
	arma_rng::set_seed_random();
	vec x = linspace(0.0, L, SPACESTEPS+2);
	vec x_ = x.subvec(1, x.n_rows - 2);
	cout << "dx = " << dx << endl;
	cout << "dt = " << dt << endl;



/*
//##############################################################
//####					Particle in box			   			####
//##############################################################
//   (Uncomment to build a free particle in a box Hamiltonian)
//	Solve an eigenvalue problem. 

	Hamiltonian H;
	H.solveEVP()
	H.toFile(RAW_PATH, "infwell");
*/


/*
//##############################################################
//####		Particle in box wave function different IC		####
//##############################################################

	//Builds the wavefunction
	//Uncomment to build a wavefunction with
 	//psi_init = psi_1
	//psi_init = psi_1 + psi_2 / sqrt(2)
	//psi_init = delta(x-L/2)


	Hamiltonian H{ RAW_PATH };
	
 	vec init = H.getSol().eigenvecs.col(0);
	WaveFunction psi{H.getSol(), init, END_TIME};
	psi.saveToFile(RAW_PATH, "wavefunction_boring");
	WaveFunction psi{ H.getSol(), init};
	psi.saveToFile(RAW_PATH, "wavefunction_boring");
	
	vec init = sqrt(0.5) * (H.getSol().eigenvecs.col(0) + H.getSol().eigenvecs.col(1));
	WaveFunction psi{H.getSol(), init,END_TIME};
	psi.saveToFile(RAW_PATH, "wavefunction1");
	WaveFunction psi{ H.getSol(), init};
	psi.saveToFile(RAW_PATH, "wavefunction1");

	vec init = zeros(SPACESTEPS + 2);
	init.subvec(1, init.n_rows-2) = arma::exp(-(x.subvec(1, init.n_rows - 2) - 0.5 * L) % (x.subvec(1, init.n_rows - 2) - 0.5 * L) / dx);
	init = normalise(init);
	WaveFunction psi{ H.getSol(), init, END_TIME};
	psi.saveToFile(RAW_PATH, "wavefunction_delta");

*/

/*
//##############################################################
//####			Particle in box	with a Barrier				####
//##############################################################
	const double V0_ = 22;

	ofstream outfile_clear(RAW_PATH + "/roots.txt", ios::trunc);
	outfile_clear.close();

	ofstream outfile(RAW_PATH + "/roots.txt", ios::app);
	
	for (double V0 = V0_; V0 < V0_+3; V0++)
	{
		vec V = ones(SPACESTEPS) * V0;

		uvec indices = find(x < (L / 3) || x >(2 * L / 3));
		V.elem(indices).zeros();
		Hamiltonian H{ V };
		H.solveEVP();
		//H.toFile(RAW_PATH, "barrier1000");

	//	vec init = (H.getSol().eigenvecs.col(0) + H.getSol().eigenvecs.col(1))/sqrt(2);
	//	double cutoffTime = 10 * pi / (H.getSol().eigenvals(1) - H.getSol().eigenvals(0));
	//	cout << "Cutoff time t' = " << cutoffTime << "\n";
	//	cout << "lambda_1 = " << H.getSol().eigenvals(0) << "\n";
	//	cout << "lambda_2 = " << H.getSol().eigenvals(1) << endl;


	//	WaveFunction psi{ H.getSol(), init, cutoffTime };
	//	psi.saveToFile(RAW_PATH, "barrier1000_10T");

		vec initial_guesses = H.getSol().eigenvals.elem(find(H.getSol().eigenvals <= V0));
		vec roots = zeros(initial_guesses.n_elem);
		for (int i = 0; i < initial_guesses.n_elem; i++) {
			double guess = initial_guesses(i);
			roots(i) = newtonRaphson(f, guess, V0);
		}

		outfile << fixed << setprecision(12);

		outfile << "v0 = " << V0 << setw(13) << "Eigenvalue" << setw(20) << "Roots" << "\n";
		for (arma::uword i = 0; i < roots.n_elem; ++i) {
			outfile <<setw(18)<<"    " << setw(20) << initial_guesses(i) << setw(20) << roots(i) << "\n";
		}

	}
	*/

/*
//##############################################################
//####	Particle in box	with a Barrier, Crank-Nicolson		####
//##############################################################

	const double V0 = 300;
	vec V = ones(SPACESTEPS) * V0;
	uvec indices = find(x_ < (L / 3) || x_ >(2 * L / 3));
	V.elem(indices).zeros();
	Hamiltonian H(V);
	H.solveEVP();

	cx_vec init = cx_vec(H.getSol().eigenvecs.col(0),zeros(H.getSol().eigenvecs.col(0).n_elem));
	CrankNicolson wave(init, H, dt, dx, TIMESTEPS, SPACESTEPS);
	wave.saveToFile(RAW_PATH, "test2");

	vec init1 = H.getSol().eigenvecs.col(0);
	cx_vec cx_init1 = cx_vec(init1, zeros(init1.n_elem));
	CrankNicolson wave_1 (cx_init1, H, dt, dx, TIMESTEPS, SPACESTEPS);
	wave_1.saveToFile(RAW_PATH, "psi_1");

	vec init2 = sqrt(0.5) * (H.getSol().eigenvecs.col(0) + H.getSol().eigenvecs.col(1));
	cx_vec cx_init2 = cx_vec(init2, zeros(init2.n_elem));
	CrankNicolson wave_12 (cx_init2, H, dt, dx, TIMESTEPS, SPACESTEPS);
	wave_12.saveToFile(RAW_PATH, "psi_12");

	vec init3 = zeros(x.n_rows);
	init3.subvec(1, x.n_rows - 2) = arma::exp(-(x.subvec(1, x.n_rows - 2) - 0.2 * L) % (x.subvec(1, x.n_rows - 2) - 0.2 * L) / dx);
	init3 = normalise(init3);
	cx_vec cx_init3 = cx_vec(init3, zeros(init3.n_elem));
	CrankNicolson wave_delta(cx_init3, H, dt, dx, TIMESTEPS, SPACESTEPS);
	wave_delta.saveToFile(RAW_PATH, "psi_delta");
	*/

/*

//##############################################################
//####		Test for different detuning potentials			####
//##############################################################
	const double V0 = 100;
	vec Vr_vec = linspace(-100, 100, 21);
	mat res = zeros<mat>(2,21);
	mat vec1 = zeros<mat>(SPACESTEPS + 2, 21);
	mat vec2 = zeros<mat>(SPACESTEPS + 2, 21);

	for (int i = 0; i < Vr_vec.n_elem; i++){
		cout << i+1 <<". Vr = " << Vr_vec[i]<<"\n";
		vec V = ones(SPACESTEPS) * V0;
		uvec indices_left = find(x_ < (L / 3));
		V.elem(indices_left).zeros();
		uvec indices_right = find(x_ > (2 * L / 3));

		V.elem(indices_right).ones();
		V.elem(indices_right) = V.elem(indices_right) * Vr_vec(i);
		Hamiltonian H(V);

		H.solveEVP(15);
		
		res.col(i) = H.getSol().eigenvals.subvec(0,1);
		vec1.col(i) = H.getSol().eigenvecs.col(0);
		vec2.col(i) = H.getSol().eigenvecs.col(1);

	}
	res.save(RAW_PATH + "/Eigenvalues_different_Vr.csv", csv_ascii);
	vec1.save(RAW_PATH + "/vec1_Vr.csv", csv_ascii);
	vec2.save(RAW_PATH + "/vec2_Vr.csv", csv_ascii);
*/
/*

*/
	int N = 101;
	const double V0 = 100;
	//double Vr = 0;
	vec g;
	vec e;
	vec tau = zeros(N);
	mat res = zeros<mat>(2, N);
	mat potential = zeros<mat>(SPACESTEPS, N);
	mat vec1 = zeros<mat>(SPACESTEPS + 2, N);
	mat vec2 = zeros<mat>(SPACESTEPS + 2, N);

	vec Vr_vec = linspace(-50, 50, N);
	for (int i = 0; i < N; i++) {
		cout << i + 1 << ". Vr = " << Vr_vec[i] << "\n";
		vec V = ones(SPACESTEPS) * V0;
		uvec indices_left = find(x_ < (L / 3));
		V.elem(indices_left).zeros();
		uvec indices_right = find(x_ > (2 * L / 3));

		V.elem(indices_right).ones();
		V.elem(indices_right) = V.elem(indices_right) * Vr_vec(i);
		Hamiltonian H(V);
		H.solveEVP(20);
		g = H.getSol().eigenvecs.col(0).subvec(1, SPACESTEPS);
		e = H.getSol().eigenvecs.col(1).subvec(1, SPACESTEPS);
		res.col(i) = H.getSol().eigenvals.subvec(0, 1);
		vec1.col(i) = H.getSol().eigenvecs.col(0);
		vec2.col(i) = H.getSol().eigenvecs.col(1);
		potential.col(i) = V;
		tau(i) = H.tunnelingAmp(g,e);
	}

	tau.save(RAW_PATH + "/tau.csv", csv_ascii);
	res.save(RAW_PATH + "/Eigenvalues_different_Vr.csv", csv_ascii);
	vec1.save(RAW_PATH + "/vec1_Vr.csv", csv_ascii);
	vec2.save(RAW_PATH + "/vec2_Vr.csv", csv_ascii);
	potential.save(RAW_PATH + "/potential_mat.csv", csv_ascii);
	return 0;
}
