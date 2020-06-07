// mptypes.cpp
// 

#include <stdexcept>
#include <iostream>

#include "algorithms.h"
#include "WavFile.h"

using namespace std;

int main (int argc, char* argv[]) {
	try {
		float SR = 44100;
		float J = (12);
		int N = pow (2., J);
		int NCOMP = 100;
		float fdef = 6;
		float flimit = 19000;

		DynamicMatrix<float> dictionary;

		cout << "make dictionary..."; cout.flush();
		make_gabor_dictionary<float>(dictionary, J, fdef, flimit, SR);
		// make_fourier_dictionary<float>(dictionary, J, flimit, SR);
		std::cout << "done (" << dictionary.size () << " atoms)" << std::endl;
		cout << "save dictionary..."; cout.flush();
		WavOutFile dict ("dictionary.wav", SR, 16, 1);
		for (unsigned i = 0; i < dictionary.size (); ++i) {
			dict.write (&dictionary[i][0], dictionary[i].size ());
		}
		std::cout << "done" << endl << endl;

		WavInFile in ("../../samples/Jarrett_Vienna_cut.wav");
		WavOutFile out ("reconstruction.wav", SR, 16, 1);

		std::cout << "analysing ["; cout.flush ();
		vector<float> buffer (N);
		vector<float> hann (N);
		make_window<float>(&hann[0], N, .5, .5, 0.);
		vector<float> target (in.getNumSamples());
		vector<float> rebuild (in.getNumSamples() + N);
		memset (&rebuild[0], 0, sizeof (float) * in.getNumSamples());

		int r = in.read (&target[0], in.getNumSamples());
		int p = 0;
		int hop = N;
		vector<float> decomposition (dictionary.size ());
		while (p < r) {
			for (int i  = 0; i < N; ++i) buffer[i] = target[i + p];
			memset(&decomposition[0], 0, sizeof (float) * decomposition.size ());
			pursuit_decomposition<float>(dictionary, NCOMP, buffer, SR, decomposition);
			reconstruct_from_decomposition (dictionary, decomposition, buffer);
			for (int i  = 0; i < N; ++i) rebuild[i + p] += buffer[i]; // * hann[i];
			cout << "*"; cout.flush ();
			p += hop;
		}
		cout << "]" << endl;
		out.write (&rebuild[0], r);
	}
	catch (exception& e) {
		cout << "Error: " << e.what () << endl;
	}
	catch (...) {
		cout << "Fatal error: unknown exception" << endl;
	}
	return 0;
}

// EOF

