#pragma once

#include <cmath>
#include <functional>
#include <vector>
#include <numeric>

#define CPU_MAX_THREADS -1
#define NOT_CALCULATED numeric_limits<float>::lowest()

#ifdef _MSC_VER
#define forceinline __forceinline
#elif defined(__GNUC__)
#define forceinline __attribute__((always_inline)) inline
#else
#define forceinline inline
#endif

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#undef small
#elif __linux__
// TODO
#endif

template <typename T>
forceinline T dist(const T &x1, const T &y1, const T &x2, const T &y2) {
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

template <typename T>
forceinline T sqare_dist(T &x1, T &y1, T &x2, T &y2) {
	return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

template <typename T>
forceinline T cube(T x) {
	return x * x * x;
}

template <typename T>
forceinline T sqr(T x) {
	return x * x;
}

template<typename T>
bool check_inputs(const std::vector<T> &data_X, const std::vector<T> &data_Y,
	const std::vector<T> &data_Z, const std::vector<T> &out_x,
	const std::vector<T> &out_y, std::vector<T> &out_z) {

	size_t N = data_X.size();
	if ((N != data_Y.size()) || (N != data_Z.size()))
		throw std::invalid_argument("invalid data points size");
	if (out_x.size() != out_y.size())
		throw std::invalid_argument("invalid grid size");
	if (out_x.size() < 1)
		throw std::invalid_argument("invalid grid size");

	out_z.resize(out_x.size());

	if (N == 1) {
		out_z.resize(out_x.size());
		std::fill(out_z.begin(), out_z.end(), data_Z[0]);
		return false;
	}

	return true;
}

class progress {
protected:
	float percent = 0;
	float min_increase = 5;
	std::function<void(float)> callback = [](float p) {};

	virtual void change_progress(int n) = 0;

public:
	void set_callback(std::function<void(float)> callback) {
		this->callback = callback;
	}

	void set_min_increase(float inc) {
		min_increase = inc;
	}

	virtual void reset() = 0;
};

template <typename T>
std::vector<int> precise_knn(const std::vector<T> &data_X, 
			const std::vector<T> &data_Y, T x, T y, int K) {

	std::vector<T> D(data_X.size());
	std::vector<int> indexes(K);
	T dist2;

	for (size_t i = 0; i < data_X.size(); ++i) {
		T X = data_X[i];
		T Y = data_Y[i];

		dist2 = sqare_dist(x, y, X, Y);

		D[i] = dist2;
	}

	T min, last_min = 0;
	size_t min_index = 0;

	for (int k = 0; k < K; ++k) {
		min = std::numeric_limits<T>::max();
		for (size_t i = 0; i < D.size(); ++i) {
			dist2 = D[i];
			if (dist2 < min && dist2 > last_min) {
				min = dist2;
				min_index = i;
			}
		}

		last_min = min;
		indexes[k] = (int)min_index;
	}

	return indexes;
}
