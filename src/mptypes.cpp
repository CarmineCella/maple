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
			throw runtime_error ("syntax is 'mptypes params.txt input.wav");
		}
		Parameters<float> p;
		p.read(argv[1]);

		int N = pow (2., p.J);

		DynamicMatrix<float> dictionary;

		cout << "making dictionary..."; cout.flush();
		if (p.dictionary == "gabor") {
			make_gabor_dictionary<float>(dictionary, p.J, p.minj, p.oct_div, p.freq_limit, p.SR, p.normalize);
		} else if (p.dictionary == "fourier") {
			make_fourier_dictionary<float>(dictionary, p.J, p.freq_limit, p.SR);			
		} else throw runtime_error ("invalid dictionary type");
		std::cout << "done (" << dictionary.size () << " atoms)" << std::endl;
		cout << "saving dictionary..."; cout.flush();

		WavOutFile dict ("dictionary.wav", p.SR, 16, 1);
		for (unsigned i = 0; i < dictionary.size (); ++i) {
			dict.write (&dictionary[i][0], dictionary[i].size ());
		}
		std::cout << "done" << endl;

		WavInFile in (argv[2]);
		if (in.getSampleRate() != p.SR) throw runtime_error ("incompatible sampling rate");
		vector<float> buffer (N);
		vector<float> hann (N, 1.);
		if (p.windowing) {
			make_window<float>(&hann[0], N, .5, .5, 0.);
		}
		vector<float> target (in.getNumSamples());

		int r = in.read (&target[0], in.getNumSamples());
		int ptr = 0;
		int olap = p.windowing ? 4 : 1;
		int hop = N / olap;
		DynamicMatrix<Atom<float> > decompositions;
		
		clock_t tic = clock ();
		while (ptr < r) {
			for (int i  = 0; i < N; ++i) buffer[i] = target[i + ptr];
			vector<Atom <float> > decomposition;
			pursuit_decomposition<float>(dictionary, p.comp, buffer, p.SR, decomposition);
			int perc = (int) ((float) ptr / r * 100.); 
			if (perc % 10 == 0) {
				cout << "analysing..........." << perc << "%\r"; cout.flush ();
			}
			ptr += hop;
			decompositions.push_back(decomposition);
		}
		clock_t toc = clock ();
		cout << "analysing...........done (" << (float) (toc - tic) / CLOCKS_PER_SEC << " sec.)" << endl; cout.flush ();
		
		vector<float> rebuild (in.getNumSamples() + N);
		memset (&rebuild[0], 0, sizeof (float) * in.getNumSamples());
		ptr = 0;
		cout << "rebuilding.........."; cout.flush();
		for (unsigned i = 0; i < decompositions.size (); ++i) {
			reconstruct_from_decomposition (dictionary, decompositions[i], buffer);
			for (int i  = 0; i < N; ++i) rebuild[i + ptr] += (buffer[i] * hann[i] / (float) olap);
			ptr += hop;
		}
		WavOutFile reconstruction ("reconstruction.wav", p.SR, 16, 1);
		reconstruction.write (&rebuild[0], r);
		cout << "done"<< endl;

		cout << "modelling..........."; cout.flush();
		Matrix<float> transitions;
		compute_transition_model(dictionary, decompositions, transitions);
		ofstream outmat ("transitions.txt");
		for (unsigned i = 0; i < transitions.rows(); ++i) {
			for (unsigned j = 0; j < transitions.cols(); ++j) {
				outmat << transitions[i][j] << " ";
			}
			outmat << endl;
		}
		cout << "done"<< endl;
		cout << "generating.........."; cout.flush();
		int gframes = in.getNumSamples() / hop;
		DynamicMatrix<Atom<float> > gendec;
		generate_probabilistic_decomposition<float>(dictionary, transitions, 
			gframes, p.comp, gendec);
		rebuild.resize(gframes * hop + N);
		memset (&rebuild[0], 0, sizeof (float) * in.getNumSamples());
		ptr = 0;
		for (unsigned i = 0; i < gendec.size (); ++i) {
			reconstruct_from_decomposition (dictionary, gendec[i], buffer);
			for (int i  = 0; i < N; ++i) rebuild[i + ptr] += (buffer[i] * hann[i] / (float) olap);
			ptr += hop;
		}
		WavOutFile generation ("generation.wav", p.SR, 16, 1);
		generation.write (&rebuild[0], rebuild.size ());
		cout << "done"<< endl << endl;

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

