// utils.h
// 


#ifndef UTILS_H
#define UTILS_H 

#include <string>
#include <fstream>
#include <deque>
#include <sstream>
#include <stdexcept>
#include <cmath>

std::string trim (std::string const& source,
	char const* delims = " \t\r\n") {
	std::string result (source);
	std::string::size_type index = result.find_last_not_of (delims);
	if (index != std::string::npos)
		result.erase (++index);

	index = result.find_first_not_of (delims);
	if (index != std::string::npos)
		result.erase (0, index);
	else
		result.erase ();
	return result;
}

template <typename T>
struct Parameters {
	Parameters () {
		setup ();
	}
	Parameters (const char* config_file) {
		setup ();
		read (config_file);
	}
	void setup () {
		window_type = "gabor"; SR = 44100;  J = 12; minj = 8; comp = 100; 
		harm_coeff = 1;
		geom_coeff = 0;
		phi_slices = 4;	
		overlap = 4; freq_limit = 17000;

		ratio = 1.;
		stretch = 1.;
		threshold = 0.;
	}
	void print (std::ostream& out) {
		out << "dictionary type..... " << window_type << std::endl;
		out << "sampling rate....... " << SR << " Hz" << std::endl;
 		out << "lowest frequency.... " << SR / pow (2., J) << " Hz" << std::endl;
		out << "smallest time....... " << (T) pow (2, minj) / SR * 1000. << " ms" << std::endl;
		out << "harmonic coeff...... " << harm_coeff << std::endl;
		out << "geometric coeff..... " << geom_coeff << std::endl;
		out << "phase step.......... " << 2. * M_PI / phi_slices << "/" << 2. * M_PI << std::endl;
		out << "highest frequency... " << freq_limit << " Hz" << std::endl << std::endl;
		out << "components.......... " << comp << std::endl;
		out << "overlap............. " << overlap << std::endl;
		out << "pitch shift......... " << ratio << std::endl;
		out << "time stretch........ " << stretch << std::endl << std::endl;
 	}
	void read (const char* config_file) {
		std::ifstream config (config_file);
		if (!config.good ()) {
			throw std::runtime_error ("cannot open configuration file");
		}

	    int line = 0;
	    while (!config.eof ()) {
	        std::string inp;
	        std::string opcode;

	        ++line;
	        std::getline (config, inp, '\n');

	        inp = trim (inp);
	        if (inp.size () == 0) continue;

	        std::istringstream istr (inp, std::ios_base::out);

	        std::deque <std::string> tokens;
	        while (!istr.eof ()) {
	            istr >> opcode;
	            tokens.push_back (opcode);
	        }

	        if (tokens[0][0] == ';') continue;
	        if (tokens.size () < 2) {
	            std::stringstream err;
	            err << "invalid syntax at line " << line;
	            throw std::runtime_error (err.str ());
	        }

	        set_parameter (tokens);
	    }
	}

	void set_parameter (std::deque<std::string>& tokens) {
   		if (tokens[0] == "SR") {
        	SR = atof (tokens[1].c_str ());
        } else if (tokens[0] == "J") {
        	J = atol (tokens[1].c_str ());
        } else if (tokens[0] == "minj") {
        	minj = atol (tokens[1].c_str ());
        } else if (tokens[0] == "components") {
        	comp = atol (tokens[1].c_str ());
        } else if (tokens[0] == "harm_coeff") {
        	harm_coeff = atof (tokens[1].c_str ());
        } else if (tokens[0] == "geom_coeff") {
        	geom_coeff = atof (tokens[1].c_str ());
        } else if (tokens[0] == "phi_slices") {
        	phi_slices = atol (tokens[1].c_str ());
        } 
        else if (tokens[0] == "overlap") {
        	overlap = atol (tokens[1].c_str ());
        } else if (tokens[0] == "freq_limit") {
        	freq_limit = atof (tokens[1].c_str ());
        } else if (tokens[0] == "window") {
			window_type = tokens[1];
        } else if (tokens[0] == "ratio") {
        	ratio = atof (tokens[1].c_str ());
        } else if (tokens[0] == "stretch") {
        	stretch = atof (tokens[1].c_str ());
        } else if (tokens[0] == "threshold") {
        	threshold = atof (tokens[1].c_str ());
        } else {
            std::stringstream err;
            err << "invalid parameter " << tokens[0];
            throw std::runtime_error (err.str ());
        }
	}

	T SR;
	int J;
	int minj;
	int comp;
	T harm_coeff;
	T geom_coeff;
	int phi_slices;
	int overlap;
	T freq_limit; 
	T threshold;
	std::string window_type;

	T ratio;
	T stretch;
};

// -----------------------------------------------------------------------------
template <typename T>
void interleave (T* stereo, const T* l, const T* r, int n) {
	for (int i = 0; i < n; ++i) {
		stereo[2 * i] = l[i];
		stereo[2 * i + 1] = r[i];
	}
}

template <typename T>
void deinterleave (const T* stereo, T* l, T* r, int n) {
	for (int i = 0; i < n; ++i) {
		l[i] = stereo[2 * i];
		r[i] = stereo[2 * i + 1];
	}
}


#endif	// UTILS_H 

// EOF
