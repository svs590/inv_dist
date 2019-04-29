#pragma once

#include <cstdint>
#include <vector>
#include <cmath>
#include <exception>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <omp.h>

#include "utils.h"
#include "smoothing.h"

template <typename T>
std::vector<T> inverse_distances(
    std::vector<T> data_X,
    std::vector<T> data_Y,
	std::vector<T> data_Z, 
    int Nx, 
    int Ny, 
    std::vector<T> out_x,
	std::vector<T> out_y,
    std::vector<T> &out_z,
	T beta,  
    int neighbors = -1,
	T slope = 0, 
    bool smooth = false, 
    int cpu_threads = CPU_MAX_THREADS
) {

	if (!check_inputs(data_X, data_Y, data_Z, out_x, out_y, out_z))
		return out_z;

	int N = (int)data_X.size();
	beta = abs(beta);
	slope = (N > 1) ? abs(slope) : 0;
	T eps = 100 * std::numeric_limits<T>().epsilon();

	T x_len = *max_element(data_X.begin(), data_X.end()) - *min_element(data_X.begin(), data_X.end());
	T y_len = *max_element(data_Y.begin(), data_Y.end()) - *min_element(data_Y.begin(), data_Y.end());

	T r = 1. / 4 * sqrt(x_len*x_len + y_len*y_len);
	T r2 = r * r;

	T dx = x_len / (Nx - 1);
	T dy = y_len / (Ny - 1);
	T r_not_smooth = 0.5 * sqrt(dx * dx + dy * dy);
	std::vector<bool> to_smooth(out_x.size(), true);

	if (cpu_threads == CPU_MAX_THREADS)
		cpu_threads = omp_get_max_threads();

	std::vector<int> nn_indexes;

	if (neighbors < 0 || neighbors > N) {
		neighbors = N;

		for (int i = 0; i < N; ++i)
			nn_indexes.push_back(i);
	}
	else {
		neighbors = std::max(10, neighbors);
	}

	// Частные производные в окрестности точки данных
    std::vector<T> A;
    std::vector<T> B;
	// Коэффициент наклона
	T nu = 0;

	if (slope > 0) {
		A.resize(N);
		B.resize(N);

		// Расчет частных производных
		for (int p = 0; p < N; ++p) {
			A[p] = 0;
			B[p] = 0;
			double denominator = 0;

			if (neighbors < N)
				nn_indexes = precise_knn(data_X, data_Y, data_X[p], data_Y[p], neighbors);

			for (int j = 0; j < neighbors; ++j) {
				int q = nn_indexes[j];

				if (p != q) {
					double dist = sqare_dist(data_X[q], data_Y[q], data_X[p], data_Y[p]);
					double weight = 1. / pow(dist, beta / 2);

					denominator += weight;

					A[p] += weight * (data_Z[q] - data_Z[p]) *
						(data_X[q] - data_X[p]) / dist;

					B[p] += weight * (data_Z[q] - data_Z[p]) *
						(data_Y[q] - data_Y[p]) / dist;
				}
			}

			A[p] /= denominator;
			B[p] /= denominator;
		}

		// Поиск максимального значения A[i]^2 + B[i]^2
		T max_A2_B2 = 0, tmp;
		for (int p = 0; p < N; ++p) {
			tmp = A[p] * A[p] + B[p] * B[p];
				if (tmp > max_A2_B2)
					max_A2_B2 = tmp;
		}

		nu = slope * (*max_element(data_Z.begin(), data_Z.end()) -
			*min_element(data_Z.begin(), data_Z.end())
			) / sqrt(max_A2_B2);

	}

#pragma omp parallel for num_threads(cpu_threads) private(nn_indexes)
		for (int i = 0; i < out_x.size(); ++i) {
			T numerator = 0;
			T denominator = 0;
			T grid_coord_x = out_x[i];
			T grid_coord_y = out_y[i];
			T dz = 0;

			int c = 0;
			T __r = 0;

			if (neighbors > 0 && neighbors < N) {

				nn_indexes = precise_knn(data_X, data_Y, grid_coord_x, grid_coord_y, neighbors);

				for (int k = 0; k < neighbors; ++k) {

					T neighbor_x = data_X[nn_indexes[k]];
					T neighbor_y = data_Y[nn_indexes[k]];

					T d = sqare_dist(grid_coord_x, grid_coord_y,
						neighbor_x, neighbor_y);

					if (d <= r2)
						c++;
				}
				if (c < 4)
					__r = dist(
						grid_coord_x,
						grid_coord_y,
						data_X[nn_indexes[4]],
						data_Y[nn_indexes[4]]
					);
				else if (c >= 4 && c < 10)
					__r = dist(
						grid_coord_x,
						grid_coord_y,
						data_X[nn_indexes[c]],
						data_Y[nn_indexes[c]]
					);
				else
					__r = r;
			}

			for (int k = 0; k < neighbors; ++k) {
				T neighbor_x = data_X[nn_indexes[k]];
				T neighbor_y = data_Y[nn_indexes[k]];
				T neighbor_z = data_Z[nn_indexes[k]];
				T s;

				T d = dist(grid_coord_x, grid_coord_y,
					neighbor_x, neighbor_y);

				// С учетом соседей другая развесовка
				if (neighbors > 0 && neighbors < N)
					if (d <= __r / 3)
						s = 1. / d;
					else if (d > __r / 3 && d <= __r)
						s = 27. * sqr(d / __r - 1) / 4. / __r;
					else
						s = 0;
				else
					s = 1. / d;

				if (smooth && d <= r_not_smooth)
					to_smooth[i] = false;

				if (d < eps) {
					numerator = neighbor_z;
					denominator = 1;
					break;
				}

				if (slope > 0) {
					dz = slope * (A[nn_indexes[k]] * (grid_coord_x - neighbor_x) +
						B[nn_indexes[k]] * (grid_coord_y - neighbor_y)
						) * nu / (nu + d);
				}

				numerator += (neighbor_z + dz) * pow(s, beta);
				denominator += pow(s, beta);
			}

			out_z[i] = numerator / denominator;
		}

		if (smooth) {
			int iter_num = static_cast<int>(floor(log2(std::max(Nx, Ny)))) - 1;
			iter_num = (iter_num > 5) ? pow(2, iter_num) : 32;

			biharmonic_smoothing(out_z, Nx, Ny, out_x, out_y, to_smooth, iter_num);
		}

	return out_z;
}
