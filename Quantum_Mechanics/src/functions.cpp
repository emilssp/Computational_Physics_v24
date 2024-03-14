#include "functions.hpp"

void eigToFile(EVPsol sol, string path, string name)
{
	cout << "Exporting data..." << endl;

	// Write eigenvalues and eigenvectors to text files
	sol.eigenvals.save(path + "/eigenvalues_" + name + ".csv", csv_ascii);
	sol.eigenvecs.save(path + "/eigenvectors_" + name + ".csv", csv_ascii);

	cout << "Data exported to " << path << endl;
}

EVPsol eigFromFile(string path, string name)
{
	cout << "Loading eigenvalues..." << endl;
	EVPsol sol;
	vec eigvals;
	std::string eigenvalues_file = path + "/eigenvalues_" + name + ".csv";

	// Load eigenvalues
	if (!eigvals.load(eigenvalues_file)) {
		cerr << "Error loading eigenvalues file: " << eigenvalues_file << endl;
		return sol;
	}
	sol.eigenvals = eigvals;
	
	cout << "Loading eigenvectors..." << endl;

	std::string eigenvectors_file = path + "/eigenvectors_" + name + ".csv";

	// Armadillo vectors and matrix to store eigenvalues and eigenvectors
	mat eigvecs;

	// Load eigenvectors
	if (!eigvecs.load(eigenvectors_file)) {
		cerr << "Error loading eigenvectors file: " << eigenvectors_file << endl;
		return sol;
	}
	sol.eigenvecs = eigvecs;

	cout << "Finished loading data." << endl;
	return sol;
}

vec potFromFile(string path, string name) {
	cout << "Loading potential..." << endl;
	vec pot;

	std::string potential_file = path + "/potential_" + name + ".csv";

	// Load eigenvalues
	if (!pot.load(potential_file)) {
		cerr << "Error loading potential file: " << potential_file << endl;
		return pot;
	}
	return pot;
}

void potToFile(vec V, string path, string name) {
	V.save(path + "/potential_" + name + ".csv", csv_ascii);
}

double f(double lambda, double V0)
{
	if (lambda >= V0) return INFINITY;

	double k = sqrt(lambda);
	double kappa = sqrt(V0 - lambda);
	double term1 = (kappa * sin(k / 3) + k * cos(k / 3));
	double term2 = (kappa * sin(k / 3) - k * cos(k / 3));

	return exp(k / 3) * term1 * term1 - exp(-k / 3) * term2 * term2;
}

uvec nonzeroColsMatrix(cx_mat A)
{
	uvec nonzero;
	for (size_t i = 0; i < A.n_cols; ++i) {
		if (arma::accu(A.col(i) != 0) > 0) {
			nonzero.insert_rows(nonzero.n_elem, 1);
			nonzero(nonzero.n_elem - 1) = i;
		}
	}
	return nonzero;
}

cx_vec thomasAlgorithm(cx_mat A, cx_vec d)
{
	/*
	Idea for the alghorithm 
	https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	*/

	uword N = d.size();

	cx_vec a = A.diag(-1);
	cx_vec b = A.diag(0);
	cx_vec c = A.diag(1);
	cx_double w = (0.0,0.0);

	for (uword i = 1; i < N; i++) {
		w = a[i] / b[i - 1];
		b[i] = b[i] - w * c[i - 1];
		d[i] = d[i] - w * d[i - 1];
	}
	d[N - 1] = d[N - 1]/b[N-1];
	// This is the reverse sweep, used to update the solution vector f                                                                                                                                                 
	for (int i = N - 2; i>=0; i-- ) {
		d[i] = (d[i] - c[i] * d[i + 1])/b[i];
	}
	return d;
}

double derivative(function<double(double, double)> f, double x, double V0, const double h){
	return f(x + h, V0) - f(x - h, V0) / (2 * h);
}

double newtonRaphson(function<double(double, double)> f, double initial_guess, double V0, double tolerance, int max_iterations)
{
	double x = initial_guess;
	for (int i = 0; i < max_iterations; ++i) {
		
		double f_val = f(x, V0);
		double df_val = derivative(f,x,V0);

		if (abs(df_val) < tolerance) {
			cerr << "Derivative too small." << endl;
			return x;
		}

		double x_next = x - f_val / df_val;

		if (abs(x_next - x) < tolerance) {
			return x_next;
		}

		x = x_next;
	}
	cerr << "Max iterations reached without convergence." << endl;
	return x;
}

cx_vec H_psi(cx_vec psi, double eps0, double tau, double w, double t) {
	psi(0) = complex<double>(0, 1) * exp(complex<double>(0, -eps0 * t )) * tau * sin(w * t ) * psi(1);
	psi(1) = complex<double>(0, 1) * exp(complex<double>(0,  eps0 * t )) * tau * sin(w * t ) * psi(0);
	return psi;
}

cx_vec trapz(cx_mat cx_init, double eps0, double tau, double w, double h)
{
	cx_vec temp = zeros<cx_vec>(2);
	cx_vec sum = zeros<cx_vec>(2);
	sum += H_psi(cx_init.col(0), eps0, tau, w, 0);
	for (int i = 1; i < cx_init.n_cols - 1; i++) {
		temp = cx_init.col(i);
		sum += 2 * H_psi(temp, eps0, tau, w, i * h);
	}

	return (h / 2) * sum;
}

cx_vec solveF(cx_mat init, double eps0, double tau, double w, double h)
{	
	cx_vec intgrl = zeros<cx_vec>(2);
	intgrl = trapz(init, eps0, tau, w, h);
	int n = init.n_cols;
	cx_mat H = zeros<cx_mat>(2, 2);

	H(0, 1) = exp(complex<double>(0, -eps0 * n * h )) * tau * sin(w * n * h);
	H(1, 0) = exp(complex<double>(0, eps0 * n * h )) * tau * sin(w * n  * h);

	cx_mat A = (eye<cx_mat>(2, 2) + complex<double>(0,h/2) * H);
	cx_vec b = init.col(0) -  intgrl;
	cx_vec sol = zeros<cx_vec>(2);
	sol = solve(A, b);

	return sol;
}
