// mptypes.cpp
// 

#include <stdexcept>
#include <iostream>

#include "algorithms.h"
#include "utils.h"
#include "WavFile.h"

using namespace std;

int main (int argc, char* argv[]) {
	cout << "[mptypes]" << endl << endl;
	cout << "matching pursuit decomposition" << endl;
	cout << "(C) 2020 www.carminecella.com" << endl << endl;
	try {
		if (argc != 3) {
			throw runtime_error ("syntax is 'mptypes input.wav params.txt");
		}
		Parameters<float> p;
		p.read(argv[1]);

		int N = pow (2., p.J);

		DynamicMatrix<float> dictionary;

		cout << "make dictionary..."; cout.flush();
		if (p.dictionary == "gabor") {
			make_gabor_dictionary<float>(dictionary, p.J, p.minj, p.oct_div, p.freq_limit, p.SR, p.normalize);
		} else if (p.dictionary == "fourier") {
			make_fourier_dictionary<float>(dictionary, p.J, p.freq_limit, p.SR);			
		} else throw runtime_error ("invalid dictionary type");
		std::cout << "done (" << dictionary.size () << " atoms)" << std::endl;
		cout << "save dictionary..."; cout.flush();
		WavOutFile dict ("dictionary.wav", p.SR, 16, 1);
		for (unsigned i = 0; i < dictionary.size (); ++i) {
			dict.write (&dictionary[i][0], dictionary[i].size ());
		}
		std::cout << "done" << endl << endl;

		WavInFile in (argv[2]);
		if (in.getSampleRate() != p.SR) throw runtime_error ("incompatible sampling rate");

		WavOutFile out ("reconstruction.wav", p.SR, 16, 1);
		std::cout << "analysing ["; cout.flush ();
		vector<float> buffer (N);
		vector<float> hann (N, 1.);
		if (p.windowing) {
			make_window<float>(&hann[0], N, .5, .5, 0.);
		}
		vector<float> target (in.getNumSamples());
		vector<float> rebuild (in.getNumSamples() + N);
		memset (&rebuild[0], 0, sizeof (float) * in.getNumSamples());

		int r = in.read (&target[0], in.getNumSamples());
		int ptr = 0;
		int olap = p.windowing ? 4 : 1;
		int hop = N / olap;
		vector<float> decomposition (dictionary.size ());
		while (ptr < r) {
			for (int i  = 0; i < N; ++i) buffer[i] = target[i + ptr];
			memset(&decomposition[0], 0, sizeof (float) * decomposition.size ());
			pursuit_decomposition<float>(dictionary, p.comp, buffer, p.SR, decomposition);
			reconstruct_from_decomposition (dictionary, decomposition, buffer);
			for (int i  = 0; i < N; ++i) rebuild[i + ptr] += (buffer[i] * hann[i] / (float) olap);
			cout << "*"; cout.flush ();
			ptr += hop;
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

