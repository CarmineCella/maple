// fourier.h
// 

#ifndef FOURIER_H
#define FOURIER_H

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
void makeWindow (T* out, int N, T a0, T a1, T a2) {
    // .5, .5, 0     --> hanning
    // .54, .46, 0   --> hamming
    // .42, .5, 0.08 --> blackman
    for (int i = 0; i < N; ++i) {
        out[i] = a0 - a1 * cos ((TWOPI * (T) i) / (N - 1)) + a2 *
                 cos ((2 * TWOPI * (T) i) / (N - 1)); // hamming, hann or blackman
    }
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
T cubicInterp (T x1, T x2, T x3, T y1, T y2, T y3, T *min) {
    T a, b, c;
    T pos;
    a= ((y1-y2)/(x1-x2)-(y2-y3)/(x2-x3))/(x1-x3);
    b= (y1-y2)/(x1-x2) - a*(x1+x2);
    c= y1-a*x1*x1-b*x1;

    *min= c;

    // dy/dx = 2a*x + b = 0

    pos= -b/2.0/a;

    return pos;

}

template <typename T>
void ampFreq (const T* cbuffer, T* amp, T* freq, int N,  double R) {
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
T locmaxAmpFreq (const T* amp, const T* freq, int size, std::vector<int>& max, T FUNDAMENTAL) {
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
