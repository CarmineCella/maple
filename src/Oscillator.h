// Oscillator.h
// 

#ifndef OSCILLATOR_H
#define OSCILLATOR_H 

#include <cmath>

template <typename T>
//! Digital oscillator with linear interpolation
class Oscillator {
public:
	Oscillator (T* tab, long len, T sr = 44100) :
		m_phi (0), m_amp (1), m_freq (440) {
		m_tab = tab;
		m_len = len;
		if (sr > 2)
			m_sr = sr;
		else
			m_sr = 44100;
		m_sicvt = (T) m_len / m_sr;
	}
	void amplitude (T amp) {
		if (amp >= 0 && amp <= 1) m_amp = amp;
	}
	void frequency (T freq) {
		m_freq = freq;
	}
	T frequency () const {
		return m_freq;
	}
	void phase (T phi) {
		if (phi >= 0 && phi <= 1) m_phi = phi * m_len;
	}
	void process (T* output, int len) {
		T incr = m_sicvt * m_freq;

		for (int i = 0; i < len; ++i) {
			m_phi += incr;
			while (m_phi >= m_len) {
				m_phi -= m_len;
			}
	
			while (m_phi < 0) {
				m_phi += m_len;
			}
	
			long index = (long) m_phi;
			long next = index + 1;
			next %= m_len;
			
			T frac = m_phi - index;
			T lo = m_tab[index];
			T hi = m_tab[next];
	
			output[i] = (lo + frac * (hi - lo)) * m_amp;
		}
	}
	static void gen (T* tab, long len) {
		for (long l = 0; l < len; l++) {
			tab[l] = sin (2 * M_PI * (T) l / len);
		}
	}

private:
	T m_sr;
	T m_sicvt;
	T m_phi;
	T m_amp;
	T m_freq;
	T* m_tab;
	long m_len;
};

#endif	// OSCILLATOR_H 

// EOF
