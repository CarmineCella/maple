// algorithms.h
// 


#ifndef ALGORITHMS_H
#define ALGORITHMS_H 

#include "Matrix.h"
#include "fourier.h"
#include "WavFile.h"

#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>

template <typename T>
struct Atom {
	Atom () { position = 0; weight = 1;}	
	Atom (int p, T w) { position = p; weight = w;}
	int position;
	T weight;
};

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
inline T minimum(
		const T* values,
		int N, int& minPos) {
	if (1 > N) return 0;

	T min = values[0];
	for (int i = 0; i < N; ++i) {
		if (min > values[i]) {
			min = values[i];
			minPos = i;
		}
	}

	return min;
}

template <typename T>
inline T maximum(
		const T* values,
		int N, int& maxPos) {
	if (1 > N) return 0;

	T max = values[0];
	for (int i = 0; i < N; ++i) {
		if (max < values[i]) {
			max = values[i];
			maxPos = i;
		}
	}

	return max;
}

template <typename T>
T sum(
		const T* values,
		int N) {
	if (1 > N) return 0;

	T sum = 0.;
	for (int i = 0; i < N; ++i) {
		sum += values[i];
	}

	return sum;
}

template <typename T>
T mean(
		const T* values,
		int N) {
	if (1 > N) return 0;

	return sum(values, N) / N;
}

template <typename T>
T frand(T min, T max) {
	return ((max - min) * ((T)rand() / RAND_MAX) + min);
}
template <typename T>
int wchoice(T* dist, int n) {
	T R = frand<T>(0., 1.);
	for (int i = 0; i < n; ++i) {
		if (dist[i] >= R) {
			return i;
		}
	}
	return (int) frand<double>(0, n);
}

#define FAST_DOT
#ifdef FAST_DOT
	#include <pmmintrin.h>
	float dot_prod (const float *a, const float *b, int len) {
		float total;
		int i;
		__m128 num1, num2, num3, num4;
		num4= _mm_setzero_ps();  //sets sum to zero
		for(i=0; i<len; i+=4) {
			num1 = _mm_loadu_ps(a+i);   //loads unaligned array a into num1  num1= a[3]  a[2]  a[1]  a[0]
			num2 = _mm_loadu_ps(b+i);   //loads unaligned array b into num2  num2= b[3]   b[2]   b[1]  b[0]
			num3 = _mm_mul_ps(num1, num2); //performs multiplication   num3 = a[3]*b[3]  a[2]*b[2]  a[1]*b[1]  a[0]*b[0]
			num3 = _mm_hadd_ps(num3, num3); //performs horizontal addition
		    			                  //num3=  a[3]*b[3]+ a[2]*b[2]  a[1]*b[1]+a[0]*b[0]  a[3]*b[3]+ a[2]*b[2]  a[1]*b[1]+a[0]*b[0]
			num4 = _mm_add_ps(num4, num3);  //performs vertical addition
		}
		num4= _mm_hadd_ps(num4, num4);
		_mm_store_ss(&total,num4);
		return total;
	}
#else
	template <typename T>
	T dot_prod (const T* a, const T* b, int size) {
		// NB: it assumes vector have the same size
		T sum = 0;

		for (int i = 0; i < size; ++i) {
			sum += a[i] * b[i];
		}
		return sum;
		// if (std::isnan(sum) || std::isinf(sum)) return 0;
		// else return sum;
	}
#endif
template <typename T>
void make_fourier_dictionary (DynamicMatrix<T>& mat, int J, T flimit, T SR) {
	int N = pow (2., J);
	std::vector<T> buff (N);

	T f0 = SR / (T) N;
	T fn = f0;
	while (fn < flimit) {
		// std::cout << "fn " << fn << std::endl;
		for (unsigned t = 0; t < N; ++t) {
			buff[t] = cos (2. * M_PI * (T) t / (T) SR  * (T) fn);
		}				
		T nn = norm (&buff[0], N);
		scale<T> (&buff[0], &buff[0], N, 1. / nn);
		mat.push_back (buff);							
		fn += f0;
	}
}

template <typename T>
void make_gabor_dictionary (DynamicMatrix<T>& mat, int J, int minj, T fdef, T flimit, T SR, bool normalize) {
	int N = pow (2., J);
	std::vector<T> buff (N);
	T comma = pow (2., 1. / fdef);
	int j = minj;
	while (j <= J) {
		int n = pow (2., j);
		int u = 0;
		while (u <= (N - n)) {
			T fn = SR / (T) n;
			while (fn < flimit) {
				// std::cout << "n = " << n << "; u = " << u << "; fn = " << fn <<  std::endl;
				memset (&buff[0], 0, sizeof (T) * buff.size ());
				make_window<T> (&buff[u], n, .5, .5, 0.);	
				// gauss_window<T> (&buff[u], n, 4.);
				for (unsigned t = 0; t < n; ++t) {
					buff[t + u] *=  cos ((2. * M_PI * (T) t / (T) SR  * (T) fn));
				}				
				T nn = (T) (2 * J);
				if (normalize) nn = norm (&buff[0], N);
				scale<T> (&buff[0], &buff[0], N, 1. / nn);
				mat.push_back (buff);							
				fn *= comma;
			}
			u += n;
		}
		++j;
	}
}

template <typename T>
void interleave (T* stereo, const T* l, const T* r, int n) {
	for (int i = 0; i < n; ++i) {
		stereo[2 * i] = l[i];
		stereo[2 * i + 1] = r[i];
	}
}

template <typename T>
void pursuit_decomposition (const DynamicMatrix<T>& dictionary,
	int iterations, const std::vector<T>& target,
	T sr, std::vector<Atom<T> >& decomposition) {
	static int fnum = 0;
	std::vector<T> residual (target.size ());
	for (unsigned i = 0; i < target.size (); ++i) {
		residual[i] = target[i];
	}
	int sz = residual.size ();

	T max_modulus = 0;
	T max_prod = 0;
	int max_index = 0;

	for (unsigned i = 0; i < iterations; ++i) {		
		for (unsigned k = 0; k < dictionary.size (); ++k) {
			T d = (dot_prod (&residual[0], &dictionary[k][0], sz));
			T mod = fabs (d);

			if (mod > max_modulus) {
				max_modulus = mod;
				max_prod = d;
				max_index = k;
			}
		}
		for (unsigned k = 0; k < sz; ++k) {
			residual[k] -= (dictionary[max_index][k] * max_prod);
		}
		
		decomposition.push_back (Atom<T> (max_index, max_prod));
// #define DEBUG_DECOMPOSITION
#ifdef DEBUG_DECOMPOSITION
		std::vector<T> outv (2 * sz);
		interleave (&outv[0], &residual[0], &dictionary[max_index][0], sz);

		std::stringstream n;
		n << "iteration_" << fnum << "_" << std::setw (3) << std::setfill ('0') << i << ".wav";
		WavOutFile ooo (n.str ().c_str (), 44100, 16, 2);
		ooo.write (&outv[0], sz * 2);
#endif
		// std::cout << max_index << " " << max_prod << std::endl;
		max_modulus = 0;
	}
	++fnum;
}	

template <typename T>
void reconstruct_from_decomposition (const DynamicMatrix<T>& dictionary, 
	const std::vector<Atom<T> >& decomposition, 
	std::vector<T>& output) {
	int sz = dictionary[0].size ();
	output.resize (sz);
	memset (&output[0], 0, sizeof (T) * output.size ());
	for (unsigned i = 0; i < decomposition.size (); ++i) {
		int p = decomposition[i].position;
		T w = decomposition[i].weight;
		for (unsigned t = 0; t < sz; ++t) {
			output[t] +=  w * dictionary[p][t];
		}
	}
	// std::cout << std::endl;
}

template <typename T>
void export_decomposition_channel (T SR, const DynamicMatrix<T>& dictionary, 
	const DynamicMatrix<Atom<float> >& decompositions, 
	int channel) {
	std::stringstream name;
	name << "decomp_ch_" << std::setw (4) << std::setfill ('0') << channel << ".wav";
	WavOutFile chout (name.str ().c_str (), SR, 16, 1);
	for (unsigned i = 0; i < decompositions.size (); ++i) {
		int p = decompositions[i][channel].position;
		std::vector<T> t (dictionary[p]);
		//t[0] = -.5;
		chout.write(&t[0], t.size ());
	}
}

template <typename T>
void compute_transition_model (
	const DynamicMatrix<T>& dictionary, 
	const DynamicMatrix<Atom<float> >& decompositions, 
	Matrix<T>& transitions) {
	transitions.resize(dictionary.size(), dictionary.size ());
	int channels = decompositions[0].size ();
	for (unsigned ch = 0; ch < channels; ++ch) {
		for (unsigned i = 0; i < decompositions.size () - 1; ++i) {
			int m = decompositions[i][ch].position;
			int n = decompositions[i + 1][ch].position;
			transitions[m][n] += 1;
		}
	}

	// compute probabilities
	for (unsigned i = 0; i < transitions.rows(); ++i) {
		T s = sum<T> (transitions[i], transitions.cols());
		if (s) scale<T> (transitions[i],transitions[i], transitions.cols (), 1. / s);
	}
}

template <typename T>
void generate_probabilistic_decomposition (
	const DynamicMatrix<T>& dictionary, 
	const Matrix<T>& transitions, int frames, 
	int channels, DynamicMatrix<Atom <T> >& decompositions) {

	T** markov = transitions.data();
	decompositions.resize (frames);
	for (unsigned i = 0; i < frames; ++i) {
		decompositions[i].resize (channels);
	}
	for (int ch = 0; ch < channels; ++ch ) {
		int state = rand () % dictionary.size ();
		decompositions[0][ch].position = state;
		for (long i = 1; i < frames; ++i) {
			state =  wchoice (markov[state], dictionary.size ());
			decompositions[i][ch].position = state;
		}
	}
}

#endif	// ALGORITHMS_H 

// EOF
