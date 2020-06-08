// Matrix.h
// 

#ifndef MATRIX_H
#define MATRIX_H 

#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>

#include <iostream>

template <class T> class Matrix {
public:
	typedef T value_type;

    Matrix () : m_rows (0), m_cols (0), m_data (0) {};
    Matrix (unsigned int rows, unsigned int cols) : m_rows (rows), m_cols (cols) {
		if (rows <= 0 || cols <= 0) throw std::runtime_error ("invalid dimensions for matrix allocation");
		alloc ();
    };

    Matrix (const Matrix& M) : m_rows (M.m_rows), m_cols (M.m_cols) {
		alloc ();
		for (unsigned int i = 0; i < m_rows; ++i) { 
			for (unsigned int j = 0; j < m_cols; ++j) {
				m_data[i][j] = M.m_data[i][j];
			}
		}
    };
    Matrix (int rows, int cols) : m_rows (rows), m_cols (cols) {
		alloc ();
		for (unsigned int i = 0; i < m_rows; ++i) { 
			for (unsigned int j = 0; j < m_cols; ++j) {
				m_data[i][j] = 0;
			}
		}
    };
    Matrix (const T** data, int rows, int cols) : m_rows (rows), m_cols (cols) {
		alloc ();
		for (unsigned int i = 0; i < m_rows; ++i) { 
			for (unsigned int j = 0; j < m_cols; ++j) {
				m_data[i][j] = data[i][j];
			}
		}
    };

    virtual ~Matrix () {
        free ();
    };

    void null () {
		for (unsigned int i = 0; i < m_rows; ++i) { 
			for (unsigned int j = 0; j < m_cols; ++j) {
				m_data[i][j] = 0;
			}
		}
    }
    void rand () {
		for (unsigned int i = 0; i < m_rows; ++i) { 
			for (unsigned int j = 0; j < m_cols; ++j) {
				m_data[i][j] = (T) std::rand () / RAND_MAX;
			}
		}
    }
	void unit () {
		int row = std::min (m_rows, m_cols);
		m_rows = m_cols = row;
		
		for (unsigned int i = 0; i < m_rows; ++i) {
			for (unsigned int j = 0; j < m_cols; ++j) {
				m_data[i][j] = i == j ? T(1) : T(0);
			}
		}
	}
    unsigned cols () {
        return m_cols;
    }

    unsigned rows () {
        return m_rows;
    }
	int size () const {
		return m_rows * m_cols;
	}
	void clear () {
		delete [] m_data;
		m_rows = 0;
		m_cols = 0;
		m_data = 0;
	}
	void resize (int i, int j) {
		Matrix<T> temp = (*this);

		free ();
		m_rows = i;
		m_cols = j;
		alloc ();
			
		for (int ii = 0; ii < std::min<int> (temp.rows (), i); ++ii) {
			for (int jj = 0; jj < std::min<int> (temp.cols (), j); ++jj) {
				m_data[ii][jj] = temp[ii][jj];
			}
		}
	}	
    T* operator[] (unsigned i) {
        return m_data[i];
    }
	
	T operator() (int i, int j) const {
		return m_data[i][j];
	}
	T& operator() (int i, int j) {
		return m_data[i][j];
	}
	T** data () const {
		return m_data;
	}
	// algebra
	Matrix<T> inv () {
		T a1, a2, *rowptr;
		
		if (m_rows != m_cols) {
			throw std::runtime_error ("cannot invert a non-square matrix");
		}
		Matrix<T> temp (m_rows, m_cols);
		temp.unit ();
		for (unsigned int k = 0; k < m_rows; ++k) {
			int indx = pivot (k);
			if (indx == -1) {
				throw std::runtime_error ("cannot invert a singular matrix");
			}
		
			if (indx != 0) {
				rowptr = temp.m_data[k];
				temp.m_data[k] = temp.m_data[indx];
				temp.m_data[indx] = rowptr;
			}
			a1 = m_data[k][k];
			for (unsigned int j = 0; j < m_rows; ++j) {
				m_data[k][j] /= a1;
				temp.m_data[k][j] /= a1;
			}
			for (unsigned int i = 0; i < m_rows; ++i) {
				if (i != k) {
					a2 = m_data[i][k];
					for (unsigned int j = 0; j < m_rows; ++j) {
						m_data[i][j] -= a2 * m_data[k][j];
						temp.m_data[i][j] -= a2 * temp.m_data[k][j];
					}
				}
			}
		}
		return temp;
	}		
	// calculate the determinant of a matrix
	T det () const {
		T piv, detVal = T(1);
		
		if (m_rows != m_cols) {
			throw std::runtime_error ("cannot compute determinant a non-square matrix");
		}
		
		Matrix<T> temp (*this);
		for (int k = 0; k < m_rows; ++k) {
			int indx = temp.pivot (k);
			if (indx == -1) return 0;
			if (indx != 0) detVal = - detVal;
			detVal = detVal * temp.m_data[k][k];
			for (int i = k + 1; i < m_rows; ++i) {
				piv = temp.m_data[i][k] / temp.m_data[k][k];
				for (int j = k + 1; j < m_rows; ++j)
					temp.m_data[i][j] -= piv * temp.m_data[k][j];
			}
		}
		return detVal;
	}
	// calculate the norm of a matrix (Frobenius)
	inline T norm () {
		T retVal = T(0);
		
		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j) {
				retVal += m_data[i][j] * m_data[i][j];
			}
		}
		retVal = sqrt (retVal);
		return retVal;
	}	
private:
    unsigned int m_rows;
	unsigned int m_cols;
    T** m_data;

	void alloc () {
		m_data = new T*[m_rows];

		for (unsigned int i = 0; i < m_rows; ++i) { 
			m_data[i] = new T[m_cols];
			for (unsigned int j = 0; j < m_cols; ++j) {
				m_data[i][j] = 0;
			}
		}
	}
	void free () {
		for (unsigned int i = 0; i < m_rows; ++i) {
			if (m_data[i]) delete [] m_data[i];
		}
	
		if (m_data) delete [] m_data;
	}
	// private partial pivoting method
	int pivot (int row) {
		int k = row;
		T amax, temp;
		
		amax = -1;
		for (unsigned int i = row; i < m_rows; i++) {
			if ((temp = fabs ((T) m_data[i][row])) > amax && temp != 0.0) {
				amax = temp;
				k = i;
			}
		}
		if (m_data[k][row] == T(0)) return -1;
		if (k != row) {
			T* rowptr = m_data[k];
			m_data[k] = m_data[row];
			m_data[row] = rowptr;
			return k;
		}
		return 0;
	}
};

// overloaded operators
// unary transpose operator
template <typename T>
	inline Matrix<T> operator~ (Matrix<T>& m) {
	Matrix<T> temp (m.cols (), m.rows ());
	
	for (unsigned int i = 0; i < m.rows (); ++i)
		for (unsigned int j = 0; j < m.cols (); ++j) {
			T x = m (i, j);
			temp (j, i) = x;
		}
	return temp;
}
// unary inversion operator
template <typename T>
	inline Matrix<T> operator! (const Matrix<T>& m) {
	Matrix<T> temp = m;
	return temp.inv ();
}	
template <typename T>	
// logical equal-to operator
	inline bool operator== (const Matrix<T>& m1, const Matrix<T>& m2) {
	if (m1.rows () != m2.rows () || m1.cols () != m2.cols ()) {
		return false;
	}
	
	for (int i = 0; i < m1.rows (); ++i) {
		for (int j = 0; j < m1.cols (); ++j) {
			if (m1 (i,j) != m2 (i,j)) {
				return false;
			}
		}
	}
	
	return true;
}
// logical no-equal-to operator
template <typename T>
	inline bool	operator != (const Matrix<T>& m1, const Matrix<T>& m2) {
	return (m1 == m2) ? false : true;
}

template <typename T>
	std::ostream& operator<< (std::ostream& s, Matrix<T>& m_a) {
	s << m_a.rows () << "x" << m_a.cols () <<
		std::endl << "[" << std::endl;
	for (unsigned int i = 0; i < m_a.rows (); ++i) {
		s << "[";
		for (unsigned int j = 0; j < m_a.cols (); ++j) {
			s << m_a (i, j);
			if (j != m_a.cols () - 1) s << " ";
		}
	
		s << "]" << std::endl;
	}
	s << "]" << std::endl;
	return s;
}

template <typename T>
using DynamicMatrix = std::vector<std::vector <T> >;

#endif	// MATRIX_H 

// EOF
