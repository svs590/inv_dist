#ifndef SMOOTHING
#define SMOOTHING

#include <vector>
#include <memory>

template<typename T>
void biharmonic_smoothing(std::vector<T> &interp, int Nx, int Ny,
	const std::vector<T> &grid_x, const std::vector<T> &grid_y,
	const std::vector<bool> &to_smooth, int iter_num) {

	std::vector<T> tmp(interp.size());
	std::copy(interp.begin(), interp.end(), tmp.begin());

	//std::shared_ptr<vector<T>> p1(&interp);
	//std::shared_ptr<vector<T>> p2(&tmp);

	for (int iter = 0; iter < iter_num; ++iter) {

		// Узлы, далекие от границ
		for (int i = 2; i < Nx - 2; ++i) {
			for (int j = 2; j < Ny - 2; ++j) {
				if (to_smooth[i * Ny + j] == true)
					interp[i * Ny + j] = (T)-1. / 20 * (
						interp[(i + 2) * Ny + j]
						+ interp[i * Ny + j + 2]
						+ interp[(i - 2) * Ny + j]
						+ interp[i * Ny + j - 2]
						+ 2 * (
							interp[(i + 1) * Ny + j + 1]
							+ interp[(i - 1) * Ny + j + 1]
							+ interp[(i + 1) * Ny + j - 1]
							+ interp[(i - 1) * Ny + j - 1]
							)
						- 8 * (
							interp[(i + 1) * Ny + j]
							+ interp[(i - 1) * Ny + j]
							+ interp[i * Ny + j + 1]
							+ interp[i * Ny + j - 1]
							)
						);
			}
		}

		// Угловые узлы:
		std::vector<int> ud_i = {1, 1, -1, -1};
		std::vector<int> ud_j = { 1, -1, -1, 1 };
		std::vector<int> pos_i = { 0,	0,		Nx,		Nx - 1 };
		std::vector<int> pos_j = { 0,	Ny - 1,	-1,		0 };

		for (int k = 0; k < pos_i.size(); ++k) {
			int i = pos_i[k];
			int j = pos_j[k];

			if (to_smooth[i * Ny + j] == true)
				interp[i * Ny + j] = (T)0.5 * (
					  2 * interp[i * Ny + j + ud_j[k]]
					+ 2 * interp[(i + ud_i[k]) * Ny + j]
					- interp[i * Ny + j + 2 * ud_j[k]]
					- interp[(i + 2 * ud_i[k]) * Ny + j]
					);
		}

		// Узлы возле угловых, лежащие на диагонали (например, (1,1)):
		ud_i = { 1, 1, -1, -1 };
		ud_j = { 1, -1, -1, 1 };
		pos_i = { 1,	1,		Nx - 1,	Nx - 2 };
		pos_j = { 1,	Ny - 2,	-2,		1 };

		for (int k = 0; k < pos_i.size(); ++k) {
			int i = pos_i[k];
			int j = pos_j[k];

			if (to_smooth[i * Ny + j] == true)
				interp[i * Ny + j] = (T)-1. / 18 * (
					  interp[i * Ny + j + 2*ud_j[k]]
					+ interp[(i + 2*ud_i[k]) * Ny + j]
					+ interp[(i - ud_i[k]) * Ny + j + ud_j[k]]
					+ interp[(i + ud_i[k]) * Ny + j - ud_j[k]]
					+ 2 * interp[(i + ud_i[k]) * Ny + j + ud_j[k]]
					- 8 * (interp[i * Ny + j + ud_j[k]] + interp[(i + ud_i[k]) * Ny + j])
					- 4 * (interp[i * Ny + j - ud_j[k]] + interp[(i - ud_i[k]) * Ny + j])
					);
		}

		// Узел, рядом с угловым, лежащий на границе:
		ud_i = { 1, 1, -1, -1 };
		ud_j = { 1, -1, -1, 1 };
		// На боковых границах:
		pos_i = { 0,	0,			Nx - 1,	Nx - 1	};
		pos_j = { 1,	Ny - 2,		Ny - 2,	1	 };

		for (int k = 0; k < pos_i.size(); ++k) {
			int i = pos_i[k];
			int j = pos_j[k];

			if (to_smooth[i * Ny + j] == true)
				interp[i * Ny + j] = (T)-1. / 6 * (
					  interp[i * Ny + j + 2 * ud_j[k]]
					+ interp[(i + ud_i[k]) * Ny + j + ud_j[k]]
					+ interp[(i + ud_i[k]) * Ny + j - ud_j[k]]
					+ interp[(i + 2 * ud_i[k]) * Ny + j]
					- 2 * interp[i * Ny + j - ud_j[k]]
					- 4 * (interp[(i + ud_i[k]) * Ny + j] + interp[i * Ny + j + ud_j[k]])
					);
		}

		// На верхней и нижней границах:
		pos_i = { 1,	1,			Nx - 2,		Nx - 2 };
		pos_j = { 0,	Ny - 1,		Ny - 1,		0 };

		for (int k = 0; k < pos_i.size(); ++k) {
			int i = pos_i[k];
			int j = pos_j[k];

			if (to_smooth[i * Ny + j] == true)
				interp[i * Ny + j] = (T)-1. / 6 * (
					  interp[i * Ny + j + 2*ud_j[k]]
					+ interp[(i + ud_i[k]) * Ny + j + ud_j[k]]
					+ interp[(i - ud_i[k]) * Ny + j + ud_j[k]]
					+ interp[(i + 2 * ud_i[k]) * Ny + j]
					- 2 * interp[(i - ud_i[k]) * Ny + j]
					- 4 * (interp[(i + ud_i[k]) * Ny + j] + interp[i * Ny + j + ud_j[k]])
					);
		}

		// Граница и вторая граница. На верхней и нижней меняются 
		// знаки инкремента по j, на боковых инкременты по i j меняются местами
		ud_i = {1, 1, 1, -1};
		ud_j = {1, -1, 1, 1};

		std::vector<std::pair<int, int>> edge = {
			{2, Nx - 3},
			{2, Nx - 3},
			{2, Ny - 3},
			{2, Ny - 3}
		};

		for (int k = 0; k < edge.size(); ++k) {
			int i_coeff = (k < 2) ? 1 : 0;
			int j_coeff = 1 - i_coeff;
			int i0 = 0, j0 = 0;
			int i1 = 1, j1 = 1;
			int swap = (k < 2) ? 0 : 1;

			if (k == 1) {
				j0 = Ny - 1;
				j1 = Ny - 2;
			}
			else if (k == 3) {
				i0 = Nx - 1;
				i1 = Nx - 2;
			}

			for (int p = edge[k].first; p < edge[k].second; ++p) {

				int i = i1 + i_coeff * p;
				int j = j1 + j_coeff * p;

				if (to_smooth[i * Ny + j] == true)
					interp[i * Ny + j] = (T)-1. / 19 * (
						interp[(i - 2 * (1 - swap)*ud_i[k]) * Ny + j - 2 * swap*ud_j[k]]
						+ interp[(i + 2 * (1 - swap)*ud_i[k]) * Ny + j + 2 * swap*ud_j[k]]
						+ interp[(i + 2 * swap*ud_i[k]) * Ny + j + 2 * (1 - swap)*ud_j[k]]
						+ 2 * (
							interp[(i - (1 - swap)*ud_i[k] + swap * ud_i[k]) * Ny + j + (1 - swap)*ud_j[k] - swap * ud_j[k]]
							+ interp[(i + (1 - swap)*ud_i[k] + swap * ud_i[k]) * Ny + j + (1 - swap)*ud_j[k] + swap * ud_j[k]]
							)
						+ interp[(i - (1 - swap)*ud_i[k] - swap * ud_i[k]) * Ny + j - (1 - swap)*ud_j[k] - swap * ud_j[k]]
						+ interp[(i + (1 - swap)*ud_i[k] - swap * ud_i[k]) * Ny + j - (1 - swap)*ud_j[k] + swap * ud_j[k]]
						- 8 * (
							interp[(i - (1 - swap)*ud_i[k]) * Ny + j - swap * ud_j[k]]
							+ interp[(i + swap * ud_i[k]) * Ny + j + (1 - swap)*ud_j[k]]
							+ interp[(i + (1 - swap)*ud_i[k]) * Ny + j + swap * ud_j[k]]
							)
						- 4 * interp[(i - swap * ud_i[k]) * Ny + j - (1 - swap)*ud_j[k]]
						);


				i = i0 + i_coeff * p;
				j = j0 + j_coeff * p;

				if (to_smooth[i * Ny + j] == true)
					interp[i * Ny + j] = (T)-1. / 7 * (
						  interp[(i - 2*(1 - swap)*ud_i[k]) * Ny + j - 2*swap*ud_j[k]]
						+ interp[(i + 2*(1 - swap)*ud_i[k]) * Ny + j + 2*swap*ud_j[k]]
						+ interp[(i + 2*swap*ud_i[k]) * Ny + j + 2*(1 - swap)*ud_j[k]]
						+ interp[(i - (1 - swap)*ud_i[k] + swap*ud_i[k]) * Ny + j + (1 - swap)*ud_j[k] - swap*ud_j[k]]
						+ interp[(i + ud_i[k]) * Ny + j + ud_j[k]]
						- 4 * (
								  interp[(i - (1 - swap)*ud_i[k]) * Ny + j - swap*ud_j[k]]
								+ interp[(i + swap*ud_i[k]) * Ny + j + (1 - swap) * ud_j[k]]
								+ interp[(i + (1 - swap)*ud_i[k]) * Ny + j + swap * ud_j[k]]
							)
						);

				
			}
		}

		//interp = tmp;
	}

}


#endif