#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <string>

#include <chrono>

#include "inverse_distances.h"

using namespace std;


int main(int argc, char *argv[]) {

	string filename = "";
	int nx;			                // количество ячеек по Ox;
    int ny;			                // количество ячеек по Oy;
	double beta;
	int cpu_threads = CPU_MAX_THREADS;

	if (argc < 5) {
		cout << "Error: Too few parameters." << endl;
		exit(1);
	}
	else {
		filename = argv[1];
		nx = stoi(argv[2]);
		ny = stoi(argv[3]);
		beta = stod(argv[4]);
	}
	if(argc > 5)
		cpu_threads = stoi(argv[5]);

	ifstream data_file(filename);
	double _x, _y, _z;
	vector<double> X, Y, Z;
	while ((data_file >> _x)) {
		data_file >> _y;
		data_file >> _z;
		X.push_back(_x);
		Y.push_back(_y);
		Z.push_back(_z);
	}

	double x_min    = *std::min_element(X.begin(), X.end());
	double x_max    = *std::max_element(X.begin(), X.end());
	double y_min    = *std::min_element(Y.begin(), Y.end());
	double y_max    = *std::max_element(Y.begin(), Y.end());

    // Расширение области
	double left     = x_min - 0.2 * (x_max - x_min);
	double bottom   = y_min - 0.2 * (y_max - y_min);
	double right    = x_max + 0.2 * (x_max - x_min);
	double top      = y_max + 0.2 * (y_max - y_min);
	
    // Расчетная сетка
	auto x = vector<double>(nx * ny);
	auto y = vector<double>(nx * ny);
    // Выходные значения в узлах сетки
	auto z = vector<double>(nx * ny);

	double dx = (right - left) / (nx - 1);
	double dy = (top - bottom) / (ny - 1);

	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j) {
			x[ny * i + j] = left + (right - left) * i / (nx - 1);
			y[ny * i + j] = bottom + (top - bottom) * j / (ny - 1);
		}

    auto start = chrono::high_resolution_clock::now();
	inverse_distances<double>(
        X, Y, Z,
        nx, ny,
        x, y, z,
        beta,
        10, 
        1., 
        true,
        cpu_threads
    );	
	auto end = chrono::high_resolution_clock::now();

	cout << "Total time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;

	ofstream fout("output.dat");
	fout << setprecision(16);
	for (int j = 0; j < ny; ++j)
		for (int i = 0; i < nx; ++i) {
			fout << x[ny * i + j] << '\t' << y[ny * i + j] << '\t' << z[ny * i + j] << endl;
	}
	fout.close();

	//system("pause");
}
