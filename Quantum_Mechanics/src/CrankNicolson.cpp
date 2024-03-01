#include "CrankNicolson.hpp"

fwEuler::fwEuler(cx_vec init, Hamiltonian H, const double delta_t, const double delta_x, const unsigned int timesteps, const unsigned int spacesteps)
{
	ComplexTridiagonalMatrix A = H.getH();
	this->time = linspace(0, timesteps * delta_t, timesteps);

	this->res = cx_mat().zeros(spacesteps+2, timesteps);
	this->res.col(0) = init;
    int N = res.col(0).n_elem;
	A.A = speye<sp_cx_mat>(spacesteps, spacesteps) + complex<double>(0,dt)*A.A;
    cout << "Solving SE..." << endl;
	for (int i = 1; i < time.n_elem; i++) 
	{
        cout << i << endl;
		res.col(i).subvec(1, N - 2) = A.A * res.col(i - 1).subvec(1, N - 2);
	}
    cout << "Done with solving SE.";
}

void fwEuler::saveToFile(string path, string name)
{
    cout << "Exporting data..." << endl;
    ofstream outFile(path + "/it_schemes/" + name + ".csv");

    if (!outFile) {
        cerr << "Error opening file for writing.\n";
        return;
    }

    // Write the time vector as the first row
    for (size_t i = 0; i < time.n_elem; ++i) {
        outFile << time(i);
        if (i < time.n_elem - 1) outFile << ", "; // CSV format
    }
    outFile << "\n\n"; // Empty row (two newlines)

    // Write the complex matrix res
    res.save(outFile, arma::csv_ascii);

    outFile.close();
    cout << "Data exported to " << path << "/it_schemes/" << name << ".csv" << endl;

}

CrankNicolson::CrankNicolson(cx_vec init, Hamiltonian H_in, const double delta_t, const double delta_x, const unsigned int timesteps, const unsigned int spacesteps)
{
    ComplexTridiagonalMatrix H = H_in.getH();
    this->time = linspace(0, timesteps * delta_t, timesteps);
    this->res = cx_mat().zeros(spacesteps + 2, timesteps);
    this->res.col(0) = init;
    int N = res.col(0).n_elem;
    cx_mat I = eye<cx_mat>(spacesteps, spacesteps);
    cx_mat A = I + complex<double>(0, 0.5*dt) * H.A;
    cx_mat B = I - complex<double>(0, 0.5*dt) * H.A;
    cx_vec fac;


    cout << "Check if operator is unitary: A\\B = I: " << (approx_equal(solve(conj(A),B), I,"absdiff", 1e-12) ? "True" : "False") << endl;
    cout << "Solving Schrodinger's Equation..." << endl;
    auto start = chrono::high_resolution_clock::now();

    #pragma omp parallel for //parallelize the for loop 
    for (int i = 1; i < time.n_elem; i++)
    {
        
        fac = B * res.col(i-1).subvec(1, N - 2);
        res.col(i).subvec(1, N - 2) = thomasAlgorithm(A, fac);
        //res.col(i).subvec(1, N - 2) = solve(A, fac);
    }
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Elapsed time: " << elapsed.count() << " s\n";

    cout << "Done with solving Schrodingers Equation.";
}

void CrankNicolson::saveToFile(string path, string name)
{
    cout << "Exporting data..." << endl;
    ofstream outFile(path + "/it_schemes/" + name + ".csv");

    if (!outFile) {
        cerr << "Error opening file for writing.\n";
        return;
    }

    // Write the time vector as the first row
    for (size_t i = 0; i < time.n_elem; ++i) {
        outFile << time(i);
        if (i < time.n_elem - 1) outFile << ", "; // CSV format
    }
    outFile << "\n\n"; // Empty row (two newlines)

    // Write the complex matrix res
    res.save(outFile, arma::csv_ascii);

    outFile.close();
    cout << "Data exported to " << path << "/it_schemes/" << name << ".csv" << endl;

}
