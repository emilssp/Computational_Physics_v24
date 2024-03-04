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

double tau_func(vec g, vec e, sp_mat H)
{
	mat A(H);
	double temp = 0.0;

	for (int i = 0; i < A.n_rows; i++) {
		for (int j = 0; j < A.n_cols; j++) {
			temp += g(i) * H(i, j) * e(j)*dx;
		}
	}
	return temp;
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

cx_vec H_psi(cx_vec psi, double e0, double tau, double w, double t) {
	psi(1) *= complex<double>(0, 1) * exp(complex<double>(0, -e0 * t)) * tau * sin(w * t);
	psi(0) *= complex<double>(0, 1) * exp(complex<double>(0, e0 * t)) * tau * sin(w * t);
	return psi;
}

cx_vec extendedSimpsonsRule(vec init, double e0, double tau, double w, double start, double stop, int n) {

	if (n % 2 != 0) {
		std::cerr << "n must be even for Simpson's Rule." << std::endl;
		return cx_vec();
	}
	cx_vec cx_init (init,zeros(2));

	double h = (stop - start) / n;
	cx_vec sum = H_psi(cx_init, e0, tau, w, 0);

	// Sum for terms with coefficient 4
	for (int i = 1; i < n; i += 2) {
		sum += 4 * H_psi(cx_init, e0, tau, w, start + i * h);
	}

	// Sum for terms with coefficient 2
	for (int i = 2; i < n - 1; i += 2) {
		sum += 2 * H_psi(cx_init, e0, tau, w, start + i * h);
	}

	return (h / 3) * sum;
}

cx_vec solveF(vec init, double e0, double tau, double w, double start, double stop, int n)
{
	cx_vec intgrl = extendedSimpsonsRule(init, e0, tau, w, start, stop, n);
	double h = (stop - start) / n;
	

	cx_mat H = zeros<cx_mat>(2, 2);
	
	H(1, 0) = complex<double>(0, 1)* exp(complex<double>(0, -e0 * n * h))* tau* sin(w * n * h);
	H(1, 0) = complex<double>(0, 1) * exp(complex<double>(0, e0 * n * h)) * tau * sin(w * n * h);

	cx_mat A = eye<cx_mat>(2, 2) - H;
	cx_vec b = cx_vec(init, zeros(2)) - intgrl;
	return solve(A,b);
	//return cx_vec();
}