// algorithms.h
// 


#ifndef ALGORITHMS_H
#define ALGORITHMS_H 

#include "Matrix.h"
#include "fourier.h"

#include <cmath>


template <typename T>
T norm(const T* values,
		int N) {
	if (1 > N) return 0;

	T sum = 0.;
	for (int i = 0; i < N; ++i) {
		sum += values[i] * values[i];
	}

	return sqrt (sum);
}

template <typename T>
void scale(const T* buff, T* out, int len, T factor) {
	for (int i = 0; i < len; ++i) {
		out[i] = buff[i] * factor;
	}
}

template <typename T>
T inner_prod (const T* a, const T* b, int size) {
	// NB: it assumes vector have the same size
	T sum = 0;

	for (int i = 0; i < size; ++i) {
		sum += a[i] * b[i];
	}
	if (std::isnan(sum) || std::isinf(sum)) return 0;
	else return sum;
}


template <typename T>
void make_dictionary (DynamicMatrix<T>& mat, int J, T fdef, T SR) {
	int N = pow (2., J);
	std::vector<T> buff (N);
	T comma = pow (2., 1. / fdef);
	int j = 2;
	while (j <= J) {
		int n = pow (2., j);
		int u = 0;
		while (u <= (N - n)) {
			T fn = SR / (T) n;
			while (fn < (SR / 2)) {
				std::cout << "n = " << n << "; u = " << u << "; fn = " << fn << std::endl;
				memset (&buff[0], 0, sizeof (T) * buff.size ());
				make_window (&buff[u], n, .5, .5, 0.);				
				// gauss_window (&buff[u], n, 8.);
				for (unsigned t = 0; t < n; ++t) {
					buff[t + u] *= cos (2. * M_PI * (T) t / (T) SR  * (T) fn);
				}				
				T nn = norm (&buff[0], N);
				scale (&buff[0], &buff[0], N, 1. / nn);
				mat.push_back (buff);							
				fn *= comma;
			}
			u += n;
		}
		++j;
	}
}

template <typename T>
void pursuit_decomposition (const DynamicMatrix<T>& dictionary,
	int iterations, const std::vector<T>& target,
	T sr, std::vector<T>& decomposition) {

	decomposition.resize (dictionary.size ());

	double no = norm<double> (&target[0], 
		target.size ());
	no *= no;

	std::vector<double> residual (target.size ());
	for (unsigned i = 0; i < target.size (); ++i) {
		residual[i] = target[i];
	}

	T max_modulus = 0;
	T max_prod = 0;
	int max_index = 0;
	std::vector<T> est_copy (target.size ());

	for (unsigned i = 0; i < iterations; ++i) {		
		for (unsigned k = 0; k < dictionary.size (); ++k) {
			std::vector<T> p = dictionary[k];

			double d = (inner_prod(&residual[0], &p[0], residual.size ()));
			double mod = fabs (d);

			if (mod > max_modulus) {
				max_modulus = mod;
				max_prod = d;
				max_index = k;
			}
		}
		for (unsigned k = 0; k < residual.size (); ++k) {
			residual[k] -= (dictionary[max_index][k] * max_prod);
		}
		decomposition[max_index] = max_prod;
		// std::cout << max_index << " " << max_prod << std::endl;
		max_modulus = 0;
	}
}	

template <typename T>
void reconstruct_from_projections (const DynamicMatrix<T>& dictionary, const std::vector<T>& decomposition, 
	std::vector<T>& output) {

	output.resize (dictionary[0].size ());
	memset (&output[0], 0, sizeof (T) * output.size ());
	for (unsigned i = 0; i < decomposition.size (); ++i) {
		if (decomposition[i] != 0) {
			for (unsigned t = 0; t < dictionary[i].size (); ++t) {
				output[t] += decomposition[i] * dictionary[i][t];
			}
		}
	}
}

#endif	// ALGORITHMS_H 

// EOF
