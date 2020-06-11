// TwoPoles.h
//

#ifndef TWOPOLES_H
#define TWOPOLES_H

#include <cmath>

template <typename T>
//! Two-poles resonant filter: for each buffer returns a buffer of same size
class TwoPoles {
public:
    TwoPoles (T sr) {
        m_sr = sr;
		reset (110, 1.5);
    }
    TwoPoles (T sr, T frequency, T tau) {
        m_sr = sr;
        reset (frequency, tau);
    }
    virtual ~TwoPoles () {}
    void reset (T frequency, T tau) {
        T om = 2 * M_PI * (frequency / m_sr);
        T B = 1. / (tau / 2);
        T t = 1. / m_sr;
        T radius = exp (-1. * M_PI * B * t);
        m_a1 = -2 * radius * cos (om);
        m_a2 = radius * radius;

        m_gain = radius * sin (om);

        m_x1 = m_y1 = m_y2 = 0;
    }
    void process (const T* in, T* out, int len) {
        T* ptr = (T*) in;
        T* optr = out;
        T v = 0;
        while (len--) {
            v = m_gain * m_x1 - (m_a1 * m_y1) - (m_a2 * m_y2);
            m_x1 = *ptr++;
            m_y2 = m_y1;
            m_y1 = v;

            *optr++ = v;
        }
    }
private:
	T m_sr;
    T m_gain;
    T m_a1;
    T m_a2;
    T m_y1;
    T m_y2;
    T m_x1;
};

#endif	// TWOPOLES_H

// EOF
