// mptypes.cpp
// 

#include <stdexcept>
#include <iostream>

#include "algorithms.h"
#include "WavFile.h"

using namespace std;

int main (int argc, char* argv[]) {
	try {
		double SR = 44100;
		double J = (13);
		int N = pow (2., J);
		int NCOMP = 30;
		double fdef = 2;
		DynamicMatrix<double> dictionary;

		cout << "make dictionary..."; cout.flush();
		make_dictionary (dictionary, J, fdef, SR);
		std::cout << "done (" << dictionary.size () << " atoms)" << std::endl;
		cout << "save dictionary..."; cout.flush();
		WavOutFile dict ("dictionary.wav", SR, 16, 1);
		for (unsigned i = 0; i < dictionary.size (); ++i) {
			dict.write (&dictionary[i][0], dictionary[i].size ());
		}
		std::cout << "done" << endl;

		WavInFile in ("Vox.wav");

		WavOutFile out ("reconstruction.wav", SR, 16, 1);

		std::cout << "analysing "; cout.flush ();
		vector<double> buffer (N);
		while (!in.eof ()) {
			int r = in.read (&buffer[0], N);
			vector<double> decomposition;
			pursuit_decomposition(dictionary, NCOMP, buffer, SR, decomposition);
			reconstruct_from_projections (dictionary, decomposition, buffer);
			out.write (&buffer[0], r);
			cout << "*"; cout.flush ();
		}
		cout << endl;
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

