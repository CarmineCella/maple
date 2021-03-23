// fourier.h
// 

#ifndef FOURIER_H
#define FOURIER_H

#include <algorithm>
#include <iostream>
#include <cmath>

#ifndef TWOPI_H
#define TWOPI_H
static const double TWOPI = 8. * atan ((double) 1.);
#endif

template <typename T>
void cosine_window (T* out, int N, T a0, T a1, T a2) {
    // .5, .5, 0     --> hanning
    // .54, .46, 0   --> hamming
    // .42, .5, 0.08 --> blackman
    for (int i = 0; i < N; ++i) {
        out[i] = a0 - a1 * cos ((TWOPI * (T) i) / (N - 1)) + a2 *
                 cos ((2 * TWOPI * (T) i) / (N - 1)); // hamming, hann or blackman
    }
}

template <typename T>
void gauss_window (T* win, int N, T alpha) {
	int N2 = (N  / 2);
	for (int n = -N2; n < N2; ++n) {
		T r = alpha * (T) n / (T)(N / 2.);
		win[n + N2] = exp (-0.5 * r * r);
	}
}

template <typename T>
void gamma_window (T* win, int N, T order) {
	for (unsigned t = 0; t < N; ++t) {
		win[t] = pow ((T)t / N, order - 1) * exp (-2. * M_PI * (T) t / (T) N);
	}
	// int SMOOTH = 512 > N ? N : 512;
	// T env = 1;
	// T decr = 1. / (T) SMOOTH;
	// for (unsigned  i = SMOOTH; i > 0; --i) {
	// 	win[N - i] *= env;
	// 	env -= decr;
	// }
}

//! Metatemplate sine/cosine generation
template<unsigned M, unsigned N, unsigned B, unsigned A>
	struct SinCosSeries {
	static double value () {
		return 1 - (A * M_PI / B) * (A * M_PI / B) / M / (M + 1)
			*SinCosSeries<M + 2, N, B, A>::value ();
	}
};

template<unsigned N, unsigned B, unsigned A>
	struct SinCosSeries<N, N, B, A> {
	static double value () { return 1.; }
};

template<unsigned B, unsigned A, typename T = double>
	struct Sin;

template<unsigned B, unsigned A>
	struct Sin<B,A,float> {
	static float value () {
		return (A * M_PI / B) * SinCosSeries<2, 24, B, A>::value ();
	}
};

template<unsigned B, unsigned A>
	struct Sin<B, A, double> {
	static double value () {
		return (A * M_PI / B) * SinCosSeries<2, 34, B, A>::value ();
	}
};

template<unsigned B, unsigned A, typename T=double>
	struct Cos;

template<unsigned B, unsigned A>
	struct Cos<B, A, float> {
	static float value () {
		return SinCosSeries<1, 23, B, A>::value ();
	}
};
template<unsigned B, unsigned A>
	struct Cos<B, A, double> {
	static double value () {
		return SinCosSeries<1, 33, B, A>::value ();
	}
};

//! Metatemplate version of the Danielson Lanczos Algorithm
template<unsigned N, typename T = double>
	class DanielsonLanczos {
	DanielsonLanczos<N/2, T> next;
	public:
	void apply (T* data) {
		next.apply (data);
		next.apply (data + N);

		T wtemp, tempr, tempi, wr, wi, wpr, wpi;
		//      wtemp = sin (M_PI / N);
		wtemp = -Sin<N, 1, T>::value ();
		wpr = -2.0 * wtemp * wtemp;
		//      wpi = -sin (2 * M_PI / N);
		wpi = -Sin<N, 2, T>::value ();
		wr = 1.0;
		wi = 0.0;
		for (unsigned i = 0; i < N; i += 2) {
			tempr = data[i + N] * wr - data[i + N + 1] * wi;
			tempi = data[i + N] * wi + data[i + N + 1] * wr;
			data[i + N] = data[i] - tempr;
			data[i + N + 1] = data[i + 1] - tempi;
			data[i] += tempr;
			data[i + 1] += tempi;

			wtemp = wr;
			wr += wr * wpr - wi * wpi;
			wi += wi * wpr + wtemp * wpi;
		}
	}
};

template<typename T>
class DanielsonLanczos<4, T> {
public:
   void apply (T* data) {
	  T tr = data[2];
	  T ti = data[3];
	  data[2] = data[0] - tr;
	  data[3] = data[1] - ti;
	  data[0] += tr;
	  data[1] += ti;
	  tr = data[6];
	  ti = data[7];
	  data[6] = data[5] - ti;
	  data[7] = tr - data[4];
	  data[4] += tr;
	  data[5] += ti;

	  tr = data[4];
	  ti = data[5];
	  data[4] = data[0] - tr;
	  data[5] = data[1] - ti;
	  data[0] += tr;
	  data[1] += ti;
	  tr = data[6];
	  ti = data[7];
	  data[6] = data[2] - tr;
	  data[7] = data[3] - ti;
	  data[2] += tr;
	  data[3] += ti;
   }
};

template<typename T>
class DanielsonLanczos<2, T> {
public:
   void apply (T* data) {
	  T tr = data[2];
	  T ti = data[3];
	  data[2] = data[0] - tr;
	  data[3] = data[1] - ti;
	  data[0] += tr;
	  data[1] += ti;
   }
};

//! Base class for the all FFTs
template<typename T>
class AbstractFFT {
public:
	virtual ~AbstractFFT () {}
	virtual void forward (T*) = 0;
	virtual void inverse (T*) = 0;
};

//! Fast-Fourier transform based on metatemplates
template<unsigned P, typename T = double>
class FFT : public AbstractFFT<T> {
	enum { N = 1 << P };
	DanielsonLanczos<N, T> recursion;
public:
	FFT<P, T> () {
#ifdef USE_VDSP
		m_fftSize = N;
		m_fftSizeOver2 = m_fftSize / 2;
		m_log2n = P; //log2f (m_fftSize); // bins
		m_log2nOver2 = m_log2n / 2;

		m_splitData.realp = new T[m_fftSize];
		m_splitData.imagp = new T[m_fftSize];

		m_scale = 1.0f / (T) (4.0f * m_fftSize);

		// allocate the fft object once
		m_fftSetup = vDSP_create_fftsetup (m_log2n, FFT_RADIX2);
		if (m_fftSetup == NULL || /*m_inReal == NULL || m_outReal == NULL || */
				m_splitData.realp == NULL || m_splitData.imagp == NULL /*|| m_window == NULL*/) {
			throw std::runtime_error ("FFT_Setup failed to allocate enough memory");
		}
#endif
	}

	virtual ~FFT () {
#ifdef USE_VDSP
		delete [] (m_splitData.realp);
		delete [] (m_splitData.imagp);
#endif
	}
	void conjugate (T* data) {
		for (int i = 0; i < N; ++i) {
			data[2 * i + 1] *= -1;
		}
	}
	void scramble (T* data) {
		int i, m, j = 1;
		for (i = 1; i < 2 * N; i += 2) {
			if (j > i) {
				std::swap (data[j - 1], data[i - 1]);
				std::swap (data[j], data[i]);
			}
			m = N;
			while (m >= 2 && j > m) {
				j -= m;
				m >>= 1;
			}
			j += m;
		}
	}
	public:
	enum { id = P };
	void forward (T* data) {
#ifdef USE_VDSP
		//convert to split complex format with evens in real and odds in imag
		vDSP_ctoz ((COMPLEX *) data, 2, &m_splitData, 1, m_fftSize);

		//calc fft
		vDSP_fft_zip (m_fftSetup, &m_splitData, 1, m_log2n, FFT_FORWARD);

		//m_splitData.imagp[0] = 0.0;

		vDSP_ztoc (&m_splitData, 1, (COMPLEX*) data, 2, m_fftSize);
#else
		scramble (data);
		recursion.apply (data);
#endif
	}
	void inverse (T* data) {
#ifdef USE_VDSP
		//convert to split complex format with evens in real and odds in imag
		vDSP_ctoz ((COMPLEX *) data, 2, &m_splitData, 1, m_fftSize);

		vDSP_fft_zip (m_fftSetup, &m_splitData, 1, m_log2n, FFT_INVERSE);
		vDSP_ztoc (&m_splitData, 1, (COMPLEX*) data, 2, m_fftSize);

		//vDSP_vsmul (buffer, 1, &m_scale, buffer, 1, m_fftSize);
#else
		conjugate (data);
		scramble (data);
		recursion.apply (data);

#endif
	}
private:
#ifdef USE_VDSP
	int	m_fftSize;
	int m_fftSizeOver2;
	int m_log2n;
	int m_log2nOver2;
	int m_windowSize;

	T* m_inReal;
	T* m_outReal;
	T* m_window;

	T m_scale;

	FFTSetup m_fftSetup;
	COMPLEX_SPLIT m_splitData;
#endif
};

template <typename T>
	AbstractFFT<T>* createFFT (int N) {
		switch (N) {
		case 4:
			return new FFT<2, T> ();
			break;
		case 8:
			return new FFT<3, T> ();
			break;
		case 16:
			return new FFT<4, T> ();
			break;
		case 32:
			return new FFT<5, T> ();
			break;
		case 64:
			return new FFT<6, T> ();
			break;
		case 128:
			return new FFT<7, T> ();
			break;
		case 256:
			return new FFT<8, T> ();
			break;
		case 512:
			return new FFT<9, T> ();
			break;
		case 1024:
			return new FFT<10, T> ();
			break;
		case 2048:
			return new FFT<11, T> ();
			break;
		case 4096:
			return new FFT<12, T> ();
			break;
		case 8192:
			return new FFT<13, T> ();
			break;
		case 16384:
			return new FFT<14, T> ();
			break;
		case 32768:
			return new FFT<15, T> ();
			break;
		case 65536:
			return new FFT<16, T> ();
			break;
		case 131072:
			return new FFT<17, T> ();
			break;
		case 262144:
			return new FFT<18, T> ();
			break;
		case 524288:
			return new FFT<19, T> ();
			break;
		case 1048576:
			return new FFT<20, T> ();
			break;
		case 2097152:
			return new FFT<21, T> ();
			break;
		case 4194304:
			return new FFT<22, T> ();
			break;
		case 8388608:
			return new FFT<23, T> ();
			break;
		case 16777216:
			return new FFT<24, T> ();
			break;
		default:
			throw std::runtime_error ("invalid size requested for fft");
		}
}

// ---------------------------------------------------------------------- //

template <typename T>
int locmax (const T* amp, int N, std::vector<int>& max) {
	T maxPeak = amp[1];
	if (N < 2) return 0;
	if (amp[0] > amp[1]) {
		max.push_back(0);
	}

	for (int i = 1; i < N - 1; ++i) {
		T magCurr = amp[i];
		T magPrev = amp[(i - 1)];
		T magNext = amp[(i + 1)];

		if (magCurr > magPrev && magCurr > magNext) {
			max.push_back(i);
			if (magCurr > maxPeak) maxPeak = magCurr;
		}
	}
	return max.size ();
}

template <typename T>
inline T specflux (T* amplitudes, T* oldAmplitudes, int N) {
	T sf = 0; // spectral flux
	T a = 0;
	for (int i = 0; i < N; ++i) {
		a = (amplitudes[i] - oldAmplitudes[i]);
		oldAmplitudes[i] = amplitudes[i];
		sf += a < 0 ? 0 : a; // rectification
		//sf += a;
		//sf += fabs (a);
	}

	return sf;
}

#endif	// FOURIER_H 

// EOF
