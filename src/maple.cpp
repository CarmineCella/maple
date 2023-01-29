// maple.cpp
// 

#include <stdexcept>
#include <iostream>

#include "algorithms.h"
#include "utils.h"
#include "WavFile.h"

using namespace std;

#define VERSION 0.3

// TODO:
// - file-based decomposition
// - filters for file-based dictionary

int main (int argc, char* argv[]) {
	srand (time (NULL));
	cout << "[ma.p.l.e, ver. " << VERSION << "]" << endl << endl;
	cout << "matching pursuit linear expansion" << endl;
	cout << "(C) 2021-2023 www.carminecella.com" << endl << endl;
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

#ifdef DEBUG_DICTIONARY				
		cout << "saving dictionary..."; cout.flush();
		WavOutFile dict ("dictionary_atoms.wav", p.SR, 16, 1);
		for (unsigned i = 0; i < dictionary.atoms.size (); ++i) {
			dict.write (&dictionary.atoms[i][0], dictionary.atoms[i].size ());
		}
		ofstream dict_params ("dictionary_params.txt");
		if (!dict_params.good ()) throw runtime_error ("cannot create file for dictionary parameters");
		for (unsigned i = 0; i < dictionary.parameters.size (); ++i) {
			for (unsigned j = 0; j < dictionary.parameters[i].size (); ++j) {
				dict_params << dictionary.parameters[i][j] << " ";
			}
			dict_params << endl;
		}
		dict_params.close ();
		cout << " done" << endl;
#endif

		// analysis
		Decomposition<float> decomposition;
		clock_t tic = clock ();
		pursuit_decomposition(p, dictionary, target, decomposition, cout);
		clock_t toc = clock ();
		cout << "analysing........... done (" << 
			(float) (toc - tic) / CLOCKS_PER_SEC << " sec.)" << endl; cout.flush ();

		if (p.npref != 0) {
			vector<StateTab> tabs;
			Prefix prefix; // current input prefix
			cout << "building transitions"; cout.flush();
			build_transitions (p, decomposition, p.npref, prefix, tabs);
			cout << " done (order " << p.npref << ")"  << endl << endl;		
			cout << "saving transitions.. "; cout.flush();
			ofstream trans_out ("transitions.txt");
			if (!trans_out.good ()) throw runtime_error ("cannot create file for transitions");		
			export_transitions (trans_out, tabs);
			cout << " done"<< endl;
			cout << "generating.......... "; cout.flush();
			probabilistic_generation (p, tabs, dictionary, in.getNumSamples (), decomposition);
			cout << " done"<< endl;
		}		
		cout << "saving analysis..... "; cout.flush ();
		ofstream decomp_out ("decomposition.txt");
		if (!decomp_out.good ()) throw runtime_error ("cannot create file for decomposition");
		save_decomposition (p, decomp_out, decomposition);
		cout << "done" << endl;

#ifdef DEBUG_DECOMPOSITION
		for (int i = 0; i < p.comp; ++i) {
			export_decomposition_channel(p.SR, i, dictionary, decomposition);
		}
#endif

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

 