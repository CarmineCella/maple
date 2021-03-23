// maple.cpp
// 

#include <stdexcept>
#include <iostream>

#include "algorithms.h"
#include "utils.h"
#include "WavFile.h"

using namespace std;

#define VERSION 0.2

int main (int argc, char* argv[]) {
	cout << "[ma.p.l.e, ver. " << VERSION << "]" << endl << endl;
	cout << "matching pursuit linaer expansion" << endl;
	cout << "(C) 2021 www.carminecella.com" << endl << endl;
	try {
		if (argc != 4) {
			throw runtime_error ("syntax is 'maple params.txt input.wav output.wav");
		}
		Parameters<float> p (argv[1]);
		p.print(cout);

		cout << "reading input......."; cout.flush ();
		WavInFile in (argv[2]);
		if (in.getSampleRate() != p.SR || in.getNumChannels () != 1) {
			throw runtime_error ("incompatible sampling rate or number of channels");
		}
		vector<float> target (in.getNumSamples());
		in.read (&target[0], in.getNumSamples());
		cout << " done" << endl;

		// dictionary
		Dictionary<float> dictionary;
		cout << "making dictionary..."; cout.flush();
		make_dictionary<float>(p, dictionary);
		std::cout << " done (" << dictionary.atoms.size () << " atoms)" << std::endl;
		cout << "saving dictionary..."; cout.flush();
		WavOutFile dict ("dictionary_atoms.wav", p.SR, 16, 1);
		for (unsigned i = 0; i < dictionary.atoms.size (); ++i) {
			dict.write (&dictionary.atoms[i][0], dictionary.atoms[i].size ());
		}
		ofstream dict_params ("dictionary_params.txt");
		for (unsigned i = 0; i < dictionary.parameters.size (); ++i) {
			for (unsigned j = 0; j < dictionary.parameters[i].size (); ++j) {
				dict_params << dictionary.parameters[i][j] << " ";
			}
			dict_params << endl;
		}
		dict_params.close ();
		cout << " done" << endl;

		// analysis
		Decomposition<float> decomposition;
		clock_t tic = clock ();
		pursuit_decomposition(p, dictionary, target, decomposition, cout);
		clock_t toc = clock ();
		cout << "analysing........... done (" << 
			(float) (toc - tic) / CLOCKS_PER_SEC << " sec.)" << endl; cout.flush ();
		
		cout << "saving analysis..... "; cout.flush ();
		ofstream decomp_out ("decomposition.txt");
		for (unsigned i = 0; i < p.comp; ++i) {
			for (unsigned j = 0; j < decomposition.size (); ++j) {
				decomp_out << decomposition[j][i][0] << " " << decomposition[j][i][1] << " " << decomposition[j][i][2] << endl;
			}
			decomp_out << endl;
		}
		cout << "done" << endl;

		// synthesis
		cout << "reconstructing......"; cout.flush ();
		vector<float> rebuild;
		tic = clock ();
		pursuit_reconstruction (p, dictionary, decomposition, rebuild);
		toc = clock ();
		cout << " done (" << (float) (toc - tic) / CLOCKS_PER_SEC << " sec.)" 
			<< endl; cout.flush ();

		cout << "saving output......."; cout.flush ();
		WavOutFile reconstruction (argv[3], p.SR, 16, 1);
		reconstruction.write (&rebuild[0], rebuild.size ());
		cout << " done"<< endl << endl;
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

