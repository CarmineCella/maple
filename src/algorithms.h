// algorithms.h
// 


#ifndef ALGORITHMS_H
#define ALGORITHMS_H 

#include "Matrix.h"
#include "fourier.h"
#include "utils.h"
#include "WavFile.h"

#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

#define DOT_PROD_SPEED 1
// #define DEBUG_DECOMPOSITION

// types
template <typename T>
struct Dictionary {
	DynamicMatrix<T> atoms;
	DynamicMatrix<T> parameters;
};
template <typename T>
using Decomposition = std::vector<DynamicMatrix<T> >;

// vector maths
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

#if DOT_PROD_SPEED == 0 // naive dot product
	template <typename T>
	T dot_prod (const T* a, const T* b, int size) {
		// NB: it assumes vector have the same size
		T sum = 0;

		for (int i = 0; i < size; ++i) {
			sum += a[i] * b[i];
		}
		return sum;
		if (std::isinf(sum)) return 0;
		else return sum;
	}
#elif DOT_PROD_SPEED == 1 // SSE3 dot product
	#include <pmmintrin.h>
	float dot_prod (const float *a, const float *b, int len) {
		float total;
		int i;
		__m128 num1, num2, num3, num4;
		num4= _mm_setzero_ps();  			//sets sum to zero
		for(i=0; i<len; i+=4) {
			num1 = _mm_loadu_ps(a+i);   	//loads unaligned array a into num1  
											// num1= a[3]  a[2]  a[1]  a[0]
			num2 = _mm_loadu_ps(b+i);   	//loads unaligned array b into num2 
											// num2= b[3]   b[2]   b[1]  b[0]
			num3 = _mm_mul_ps(num1, num2); 	//performs multiplication   
											// num3 = a[3]*b[3]  a[2]*b[2]  a[1]*b[1]  a[0]*b[0]
			num3 = _mm_hadd_ps(num3, num3); //performs horizontal addition
		    			                  	//num3 =  a[3]*b[3]+ a[2]*b[2]  a[1]*b[1]+a[0]*b[0]  a[3]*b[3]+ a[2]*b[2]  a[1]*b[1]+a[0]*b[0]
			num4 = _mm_add_ps(num4, num3);  //performs vertical addition
		}
		num4= _mm_hadd_ps(num4, num4);
		_mm_store_ss(&total,num4);
		return total;
	}
#else
	#error DOT_PROD_SPEED_INVALID
#endif

template <typename T>
void gram_schmidt (const Matrix<T>& matrix, Matrix<T>& base) {
	int dim = matrix.cols (); // must be squared
   	Matrix<T> r (matrix);
    Matrix<T> v (dim, dim);

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            v[i][j] = matrix[i][j];
        }
    }

    for (int i = 0; i < dim; ++i) {
        r[i][i] = norm (&v[i][0], v.size ());
        for (int j = 0; j < dim; ++j) {
            base[i][j] = v[i][j] / r[i][i];
        }
        for (int k = i + 1;  k < dim; ++k) {
            r[i][k] = dot_product (&base[i][0], &v[k][0], base[i].size ());
            for (int j=0; j < dim; ++j) {
                v[k][j] = v[k][j] - r[i][k] * base[i][j];
            }
        }
    }
}

// others
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

// matching pursuit
template <typename T>
void store_vector (std::vector<T>& atom, std::vector<T> params,
	Dictionary<T>& dict) {
		int N = atom.size ();
		T nn = norm (&atom[0], N);
		scale<T> (&atom[0], &atom[0], (T) N, 1. / nn);
		dict.atoms.push_back (atom);
		dict.parameters.push_back (params);
}
template <typename T>
void make_dictionary (const Parameters<T>& p, Dictionary<T>& dict) {
	int N = pow (2., p.J);
	std::vector<T> buff (N);
	std::vector<T> cbuff (2 * N, 0);
	if (p.dictionary_type == "cosine") {
		T f0 = p.SR / (T) N;
		T fn = f0;
		while (fn < p.freq_limit) {
			for (unsigned t = 0; t < N; ++t) {
				buff[t] = cos (2. * M_PI * (T) t / (T) p.SR  * (T) fn);
			}				
			
			store_vector(buff, std::vector<T> {fn, (T) N, 0.}, dict);
			fn += f0;
		}
	} else if (p.dictionary_type == "gabor" || p.dictionary_type == "gammatone") {
		T comma = pow (2., 1. / p.oct_div);

		int j = p.minj;
		while (j <= p.J) {
			int n = pow (2., j);
			int u = 0;
			while (u <= (N - n)) {
				T fn = p.SR / (T) n;
				while (fn < p.freq_limit) {
					T phi_incr = 2. * M_PI / p.phi_slices;
					T phi = 0.;
					 while (phi < 2. * M_PI) {
						memset (&buff[0], 0, sizeof (T) * buff.size ());
						if (p.dictionary_type == "gabor") {
							gauss_window<T> (&buff[u], n, 4.);
						} else if (p.dictionary_type == "gammatone") {				
							gamma_window<T>(&buff[u], n, 1.5);
						}
						for (unsigned t = 0; t < n; ++t) {
							buff[t + u] *=  cos ((2. * M_PI * (T) t / (T) p.SR  * (T) fn) + phi);
						}	
						phi += phi_incr;
						store_vector(buff, std::vector<T>  {fn, (T) n, (T) u}, dict);
					} 
					fn *= comma;
				}
				u += n;
			}
			++j;
		}
	} else {
		throw std::runtime_error ("invalid dictionary type");
	}
}
template <typename T>
void decompose_frame (T sr, int components, const Dictionary<T>& dictionary,
	const std::vector<T>& target, DynamicMatrix<T>& frame) {
	std::vector<T> residual (target.size (), 0);
	for (unsigned i = 0; i < target.size (); ++i) {
		residual[i] = target[i];
	}
	T* input = &residual[0];

	T max_modulus = 0;
	T max_prod = 0;
	int max_index = 0;
	static int fnum = 0;
	for (unsigned i = 0; i < components; ++i) {		
		for (unsigned k = 0; k < dictionary.atoms.size (); ++k) {
			int n = dictionary.parameters[k][1];
			int u = dictionary.parameters[k][2];
			const T* proj = &dictionary.atoms[k][u];
			T d = (dot_prod (input + u, proj, n));
			T mod = fabs (d);

			if (mod > max_modulus) {
				max_modulus = mod;
				max_prod = d;
				max_index = k;
			}
		}
		int n = dictionary.parameters[max_index][1];
		int u = dictionary.parameters[max_index][2];
		for (unsigned k = 0; k < n; ++k) {
			residual[k + u] -= (dictionary.atoms[max_index][k + u] * max_prod);
		}
		frame.push_back (std::vector<T> {(T) max_index, max_prod});

#ifdef DEBUG_DECOMPOSITION
		int sz = residual.size ();
		std::vector<T> outv (2 * sz);
		interleave (&outv[0], &residual[0], &dictionary.atoms[max_index][0], sz);

		std::stringstream name;
		name << "iteration_" << fnum << "_" << std::setw (3) << std::setfill ('0') << i << ".wav";
		WavOutFile ooo (name.str ().c_str (), 44100, 16, 2);
		ooo.write (&outv[0], sz * 2);
#endif
		// std::cout << max_index << " " << max_prod << std::endl;
		max_modulus = 0;
	}
	++fnum;
}	

template <typename T>
void pursuit_decomposition (const Parameters<T>& p, const Dictionary<T>& dictionary, 
	const std::vector<T>& target, Decomposition<T>& decomposition, std::ostream& out) {
	int ptr = 0;
	int N = pow (2., p.J);
	int hop = (int) ((T) (N / p.overlap) * p.stretch);	
	std::vector<T> buffer (N);
	int r = target.size ();
	while (ptr < r) {
		int perc = (int) ((T) ptr / r * 100.); 
		if (perc % 10 == 0) {
			out << "analysing........... " << perc << "%\r"; out.flush ();
		}		
		for (int i  = 0; i < N; ++i) {
			if (i + ptr >= r) buffer[i] = 0;
			else buffer[i] = target[i + ptr];
		}
		DynamicMatrix<T> frame;
		decompose_frame<T>(p.SR, p.comp, dictionary, buffer, frame);
		ptr += hop;
		decomposition.push_back(frame);
	}
}

template <typename T>
void reconstruct_frame (T ratio, const Dictionary<T>& dictionary, 
	const DynamicMatrix<T>& decomposition, std::vector<T>& output) {
	int sz = dictionary.atoms[0].size ();
	output.resize (sz);
	memset (&output[0], 0, sizeof (T) * output.size ());
	std::vector<T> buff;
	std::vector<T> window (sz);
	for (unsigned i = 0; i < decomposition.size (); ++i) {
		int p = decomposition[i][0]; // position
		T w = decomposition[i][1];  // weight
		T* ptr = (T*) &dictionary.atoms[p][0];
		if (ratio != 1.) {
			int nsamp = (int) ((T) sz * ratio);
			buff.resize (nsamp + 1, 0);
			T phi = 0;
			for (unsigned t = 0; t < nsamp; ++t) {
				int index = (int) phi;
				T frac = phi - index;
				int next = index == sz - 1 ? 0 : index + 1;
				buff[t] = dictionary.atoms[p][index] * (1. - frac) + dictionary.atoms[p][next] * frac;
				phi += ratio;
				if (phi >= sz) phi -= sz;
			}
			ptr = &buff[0];
		} 

		for (unsigned t = 0; t < sz; ++t) {
			output[t] +=  w * ptr[t];
		}
	}
	// std::cout << std::endl;
}

template <typename T> 
void pursuit_reconstruction (const Parameters<T>& p, const Dictionary<T>& dictionary, 
	const Decomposition<T>& decomposition, 	std::vector<T>& output) {
	int N = pow (2., p.J);
	int hop = N / p.overlap;	
	std::vector<T> buffer (N);
	std::vector<T> hann (N, 1.);
	if (p.overlap > 1) {
		cosine_window<T>(&hann[0], N, .5, .5, 0.);
	}	
	int samples = (int) ((T) decomposition.size () * hop);
	output.resize (samples + N, 0);
	memset (&output[0], 0, sizeof (T) * samples);
	int ptr = 0;
	for (unsigned i = 0; i < decomposition.size (); ++i) {
		reconstruct_frame (p.ratio, dictionary, decomposition[i], buffer);
		for (int i  = 0; i < N; ++i) output[i + ptr] += (buffer[i] * hann[i] / (T) p.overlap);
		ptr += hop;
	}
}

template <typename T>
void export_decomposition_channel (T SR, int channel, const Dictionary<T>& dictionary, 
	const Decomposition<T>& decompositions) {
	std::stringstream name;
	name << "decomp_ch_" << std::setw (4) << std::setfill ('0') << channel << ".wav";
	WavOutFile chout (name.str ().c_str (), SR, 16, 1);
	for (unsigned i = 0; i < decompositions.size (); ++i) {
		int p = decompositions[i][channel][0]; // position
		std::vector<T> t (dictionary.atoms[p]);
		chout.write(&t[0], t.size ());
	}
}


#endif	// ALGORITHMS_H 

// EOF
