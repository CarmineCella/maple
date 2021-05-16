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



template <typename T>
void complex_mul (
    const T* src1, const T* src2, T* dest, int num) {
    while (num--) {
        T r1 = *src1++;
        T r2 = *src2++;
        T i1 = *src1++;
        T i2 = *src2++;

        *dest++ = r1*r2 - i1*i2;
        *dest++ = r1*i2 + r2*i1;
    }
}

					for (unsigned i = 0; i < N; ++i) {
						cbuff[2 * i] = buff[i];
						cbuff[2 * i + 1] = 0;
					}
					fft(&cbuff[0], N, -1); // forward		
					dict.spectra.push_back (cbuff);		

#elif DOT_PROD_SPEED == 2 // frequency domain dot product
	template <typename T>
	T dot_prod (const T* complex_a, const T* complex_b, int len) {
		std::cout << "*** " << std::endl;
		T* out = (T*) complex_a;
		complex_mul(complex_a, complex_b, out, len);
		//fft(out, len, 1); // inverse
		return out[0] / len;
	}
#else


				#if DOT_PROD_SPEED == 2
				for (unsigned i = 0; i < sz; ++i) {
					cbuff[2 * i] = residual[i];
					cbuff[2 * i + 1] = 0;
				}
				fft(&cbuff[0], sz, -1); // forward
				input = &cbuff[0];
				proj = &dictionary.spectra[k][0];
			#endif

				
				
%FUNCTION algo_omp solves y = Ax, takes the input parameters y, A, k where
%y is the output field, A is the dataset field and k is the sparsity. It
%return the solution of x.

function  x = algo_omp(k, A, y)
    xbeg = zeros(size(A,2),1);
    support=[];
    temp=y;
    count = 1;
    while count < k+1
        ST = abs(A' * temp);
        [a, b] = max(ST);
        support = [support b];
        xfinal = A(:, support)\y;
        temp = y-A(:,support) * xfinal;
        count = count + 1;
    end
    x = xbeg;
    t = support';
    x(t) = xfinal;
end









		for (map<Prefix, vector<int> >::iterator i = tab.begin (); i != tab.end (); ++i) {
			cout << "[";
			for (unsigned l = 0; l < i->first.size (); ++l) {
				cout << i->first.at (l) << " ";
			}
			cout << "] ";
			for (unsigned j = 0; j < i->second.size (); ++j) {
				cout << i->second.at (j) << " ";
			}
			cout << endl;
		}

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
