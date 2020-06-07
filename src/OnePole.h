// OnePole.h
// 

#ifndef ONEPOLE_H
#define ONEPOLE_H 

#include <cmath>

//! A simple one-pole filter
template <typename T>
class OnePole {
private:
	OnePole& operator= (OnePole&);
	OnePole (const OnePole&);
public: 
	OnePole (T sr, T cutoff, bool lp = true) : m_y1 (0) {
		m_sr = sr;
		m_y1 = 0;
		setup (cutoff, lp);
	}
	virtual ~OnePole () {
	}
	void reset () {
		m_y1 = 0;	
	}   
	void setup (T cutoff, bool lp = true) {
		T omega = 2. * M_PI * cutoff;	
		T b = 2. - cos (omega / m_sr);
		
		if (lp) {
			m_a1 = cos (omega / m_sr) - 2. + 
			sqrt ((b * b)  - 1.);
			m_b0 = 1. + (m_a1);
		} else {
			m_a1 = b - sqrt ((b * b)  - 1.);
			m_b0 = 1. - (m_a1);
		}
	}
	T* process (const T* input, T* output, int size) {
		T y = 0;
		T* tmp = (T*) input;
		for (int i = 0; i < size; ++i) {
			y = (m_b0 * tmp[i]) - (m_a1 * m_y1);
			output[i] = m_y1 = y;
		}		
		return output;
	}
private:
	T m_sr;
	T m_y1;
	T m_a1;
	T m_b0;
};

#endif	// ONEPOLE_H 

// EOF
