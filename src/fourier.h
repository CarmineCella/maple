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
struct Peak {
	Peak() {
		amp = 0;
		freq = 0;
	}
	T amp;
	T freq;
};

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

template <typename T>
void fft(T* data, const int n, const int isign) {
	int nn,mmax,m,j,istep,i;
	T wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	//if (n<2 || n&(n-1)) throw("n must be power of 2 in four1");
	nn = n << 1;
	j = 1;
	for (i=1;i<nn;i+=2) {
		if (j > i) {
			std::swap(data[j-1],data[i-1]);
			std::swap(data[j],data[i]);
		}
		m=n;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

template <typename T>
void amp_freq (const T* cbuffer, T* amp, T* freq, int N,  double R) {
    T freqPerBin = (R ) / (T) N;
	T min = 0;
    for (int i = 0; i < N; ++i) {
        amp[i] = sqrt (cbuffer[2 * i] * cbuffer[2 * i] + cbuffer[2 * i + 1] * cbuffer[2 * i + 1]);
    }
	
	freq[0] = freqPerBin;
	for (int i = 1; i < N - 1; ++i) {
		freq[i] = cubicInterp (freqPerBin * (i - 1), freqPerBin * (i), freqPerBin * (i + 1),
    		amp[i - 1], amp[i], amp[i + 1], &min);
	}
}

template <typename T>
T locmax_amp_freq (const T* amp, const T* freq, int size, std::vector<int>& max, T FUNDAMENTAL) {
	T maxPeak = amp[2];
	for (int i = 2; i < size - 2; ++i) {
		T magCurr = amp[i]; //20 * log10 (amp[i] + .00000001);
		T magPrev = amp[i - 1]; //20 * log10 (amp[i - 1] + .00000001);
        T magPrevPrev = amp[i - 2]; //20 * log10 (amp[i - 1] + .00000001);
		T magNext = amp[i + 1]; // 20 * log10 (amp[i + 1] + .00000001);
		T magNextNext = amp[i + 2]; // 20 * log10 (amp[i + 1] + .00000001);
        
        if (magCurr > magPrev && magCurr > magNext &&
			magCurr > magPrevPrev && magCurr > magNextNext) {
            max.push_back (i);
            if (magCurr > maxPeak) maxPeak = magCurr;
        }        
	}

	return maxPeak;
}

template<typename T>
void sort (Peak<T>* peaks, int k) {
	Peak<T> p;
	for (int i = 0; i < k; ++i) {
		for (int j = i + 1; j < k; ++j) {
			//if (peaks[i].freq > peaks[j].freq) {
            if (peaks[i].amp < peaks[j].amp) {
				p = peaks[i];
				peaks[i] = peaks[j];
				peaks[j] = p;
			}
		}	
	}	
}

#endif	// FOURIER_H 

// EOF
