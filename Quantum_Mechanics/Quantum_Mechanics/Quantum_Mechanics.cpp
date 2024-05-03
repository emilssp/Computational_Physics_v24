// Quantum_Mechanics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <thread>
#include <fstream>
#include <cmath>

#include <armadillo>

#include "functions.hpp"
#include "Tridiagonal.hpp"

#include "Hamiltonian.hpp"
#include "WaveFunction.hpp"
#include "CrankNicolson.hpp"


using namespace std;
using namespace arma;

int main()
{
	cout << "Armadillo version: " << arma_version::as_string() << endl;
	cout << "C++ version: " << __cplusplus << std::endl;

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
	H.solveEVP();
	H.toFile(RAW_PATH, "infwell");

	eigToFile(H.getSol(),RAW_PATH, "infwell");


//##############################################################
//####		Particle in box wave function different IC		####
//##############################################################

	//Builds the wavefunction
	//Uncomment to build a wavefunction with
 	//psi_init = psi_1
	//psi_init = psi_1 + psi_2 / sqrt(2)
	//psi_init = delta(x-L/2)
	
 	vec init1 = H.getSol().eigenvecs.col(0);
	WaveFunction psi1{H.getSol(), init1, END_TIME};
	psi1.saveToFile(RAW_PATH, "boring");

	vec init2 = sqrt(0.5) * (H.getSol().eigenvecs.col(0) + H.getSol().eigenvecs.col(1));
	WaveFunction psi2{H.getSol(), init2,END_TIME};
	psi2.saveToFile(RAW_PATH, "1");

	vec init3 = zeros(SPACESTEPS + 2);
	init3.subvec(1, init3.n_rows-2) = arma::exp(-(x.subvec(1, init3.n_rows - 2) - 0.5 * L) % (x.subvec(1, init3.n_rows - 2) - 0.5 * L) / dx);
	init3 = normalise(init3);
	WaveFunction psi3{ H.getSol(), init3, END_TIME};
	psi3.saveToFile(RAW_PATH, "delta");

*/
/*
//##############################################################
//####			Particle in box	with a Barrier				####
//##############################################################
	const double V0 = 500;
	vec V = ones(SPACESTEPS) * V0;
	uvec indices = find(x_ < (L / 3) || x_ >(2 * L / 3));
	
	V.elem(indices).zeros();
	Hamiltonian H{ V };
	H.solveEVP();
	double cutoffTime = 1 * pi / (H.getSol().eigenvals(1) - H.getSol().eigenvals(0)); //time to tunnel fully to the other side of the potential
	cout << cutoffTime << endl;


	vec init1 = H.getSol().eigenvecs.col(0);
	WaveFunction psi1{ H.getSol(), init1, cutoffTime };
	psi1.saveToFile(RAW_PATH, "boring_10T");
	vec init2 = sqrt(0.5) * (H.getSol().eigenvecs.col(0) + H.getSol().eigenvecs.col(1));
	WaveFunction psi2{ H.getSol(), init2,cutoffTime };
	psi2.saveToFile(RAW_PATH, "12_10T");
	vec init3 = zeros(SPACESTEPS + 2);
	init3.subvec(1, init3.n_rows - 2) = arma::exp(-(x.subvec(1, init3.n_rows - 2) - 0.8 * L) % (x.subvec(1, init3.n_rows - 2) - 0.8 * L) / dx);
	init3 = normalise(init3);
	WaveFunction psi3{ H.getSol(), init3, cutoffTime};
	psi3.saveToFile(RAW_PATH, "delta_10T");



	//	vec init = (H.getSol().eigenvecs.col(0) + H.getSol().eigenvecs.col(1))/sqrt(2);
	//	double cutoffTime = 10 * pi / (H.getSol().eigenvals(1) - H.getSol().eigenvals(0));
	//	cout << "Cutoff time t' = " << cutoffTime << "\n";
	//	cout << "lambda_1 = " << H.getSol().eigenvals(0) << "\n";
	//	cout << "lambda_2 = " << H.getSol().eigenvals(1) << endl;


	//	WaveFunction psi{ H.getSol(), init, cutoffTime };
	//	psi.saveToFile(RAW_PATH, "barrier1000_10T");


*/

//##############################################################
//####						Root finding					####
//##############################################################
/*
	const int max_iter = 10000;
	ofstream outfile_clear(RAW_PATH + "/roots.txt", ios::trunc);
	outfile_clear.close();
	ofstream outfile(RAW_PATH + "/roots.txt", ios::app);
	vec initial_guesses = H.getSol().eigenvals.elem(find(H.getSol().eigenvals <= V0));

	vec temp = initial_guesses;
	vec roots = zeros(initial_guesses.n_elem);
	int k = 0;
	while (k < max_iter) {
		for (int i = 0; i < initial_guesses.n_elem; i++) {
			double guess = temp(i);
			roots(i) = temp(i);
			temp(i) = newtonRaphson(f, guess, V0);	
		}
		if (arma::all((roots - temp) < 1e-16) &&arma::all((roots - temp) > 0)) {
			break;
		}
		k++;
	}
	outfile << fixed << setprecision(12);

	outfile << "v0 = " << V0 << setw(13) << "Eigenvalue" << setw(20) << "Roots" << "\n";
	for (arma::uword i = 0; i < roots.n_elem; ++i) {
		outfile <<setw(18)<<"    " << setw(20) << initial_guesses(i) << setw(20) << roots(i) << "\n";
	}

	
*/
/*
//##############################################################
//####	Particle in box	with a Barrier, Crank-Nicolson		####
//##############################################################

	const double V0 = 500;
	vec V = ones(SPACESTEPS) * V0;
	uvec indices = find(x_ < (L / 3) || x_ >(2 * L / 3));
	V.elem(indices).zeros();
	Hamiltonian H(V);
	H.solveEVP();
	double cutoffTime = 20 * pi / (H.getSol().eigenvals(1) - H.getSol().eigenvals(0));
	cout << cutoffTime << endl;
	vec init1 = H.getSol().eigenvecs.col(0);
	cx_vec cx_init1 = cx_vec(init1, zeros(init1.n_elem));
	CrankNicolson wave_1 (cx_init1, H, dt, dx, TIMESTEPS, SPACESTEPS);
	wave_1.saveToFile(RAW_PATH, "psi_1");
	WaveFunction psi2{ H.getSol(), init1,END_TIME };
	psi2.saveToFile(RAW_PATH, "1_compare");

	vec init3 = zeros(SPACESTEPS + 2);
	init3.subvec(1, init3.n_rows - 2) = arma::exp(-(x.subvec(1, init3.n_rows - 2) - 0.8 * L) % (x.subvec(1, init3.n_rows - 2) - 0.8 * L) / dx);
	init3 = normalise(init3);;
	cx_vec cx_init3 = cx_vec(init3, zeros(init3.n_elem));
	CrankNicolson wave_delta(cx_init3, H, dt, dx, TIMESTEPS, SPACESTEPS);
	wave_delta.saveToFile(RAW_PATH, "psi_delta");
	WaveFunction psi3{ H.getSol(), init3,END_TIME };
	psi3.saveToFile(RAW_PATH, "delta_compare");

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

//##############################################################
//####		Test for different detuning potentials			####
//##############################################################
	const double V0 = 100;
	vec V = ones(SPACESTEPS) * V0;
	vec Vr = linspace(-100, 100, 201);
	uvec indices_left = find(x_ < (L / 3));
	V.elem(indices_left).zeros();
	uvec indices_right = find(x_ > (2 * L / 3));
	vec tAmp = zeros(Vr.n_elem+1);

	V.elem(indices_right).zeros();

	Hamiltonian H(V);
	H.solveEVP(20);
	
	vec g0 = H.getSol().eigenvecs.col(0).subvec(1, H.getH().A.n_rows);

	vec e0 = H.getSol().eigenvecs.col(1).subvec(1, H.getH().A.n_rows);

	int i = 0;
	for (auto v : Vr) {

		V.elem(indices_right).ones();

		V.elem(indices_right) = V.elem(indices_right) * v;
		Hamiltonian H(V);
		tAmp[i] = H.tunnelingAmp(g0,e0);
		i++;
	}
	tAmp.save(RAW_PATH + "/tau3.csv", csv_ascii);
	*/

	/*
//##############################################################
//####					Rabi oscillations					####
//##############################################################
	const double V0 = 100;
	int n = 1001;
	vec V = ones(SPACESTEPS) * V0;
	uvec indices_left = find(x_ < (L / 3));
	V.elem(indices_left).zeros();
	uvec indices_right = find(x_ > (2 * L / 3));
	V.elem(indices_right).zeros();

	Hamiltonian H(V);
	H.solveEVP(20);
	vec init = ones(2);	
	init(1) = 0;
	vec e0 = ones(2);
	e0(0) = 0;

	cx_vec cx_init = cx_vec(init, zeros(2));
	cx_vec cx_e0 = cx_vec(e0, zeros(2));
	
	double eps0 = H.getSol().eigenvals(1) - H.getSol().eigenvals(0);
	double tau = 0.02 * eps0;//H.getSol().eigenvals(0);
	double w = eps0;
	double h = END_TIME * 1 / n;

	vec t = linspace(0, 100, n);

	cx_mat res = zeros<cx_mat>(2, n);
	res.col(0) = cx_init;
	cout << res.col(0);

	for (int i = 1; i<n; i++)
	{
		res.col(i) = solveF(res.cols(0,i-1), eps0, tau, w, h);
	}
	cout << "Time: " << n * h << endl;
	vec p_t = zeros(n);
	vec pg = zeros(n);
	for (int i = 0; i < n; i++) {
		
		p_t(i) = pow(abs(cdot(cx_e0, res.col(i))), 2);
	}
	vec pt_an = sin(t*tau/(2));
	res.save(RAW_PATH+"/testfile.csv", csv_ascii);
	p_t.save(RAW_PATH + "/pt.csv", csv_ascii);
	pt_an.save(RAW_PATH + "/pt_an.csv", csv_ascii);
	*/
	
	return 0;
}

