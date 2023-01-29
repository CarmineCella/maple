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
#include <cstring>

#include <sys/types.h>
#include <dirent.h>

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
		npref = 0;
		dictionary_type = "gabor"; SR = 44100;  J = 12; minj = 8; comp = 100; oct_div = 12;
		phi_slices = 4;	
		overlap = 4; freq_limit = 17000;

		ratio = 1.;
		stretch = 1.;
	}
	void print (std::ostream& out) {
		out << "dictionary type..... " << dictionary_type;
		if (dictionary_type == "frames" || dictionary_type == "onsets" || dictionary_type == "files") {
			out << " (" << dictionary_path << ") ";
			if (dictionary_type == "onsets") out << db_onset_threshold << " " << db_onset_timegate; 
		}
		out << std::endl;
		if (dictionary_type == "onsets") {
		out << "segmentation........ " << target_onset_threshold << " " << target_onset_timegate << std::endl;
		}
		out << "sampling rate....... " << SR << " Hz" << std::endl;
 		out << "lowest frequency.... " << SR / pow (2., J) << " Hz" << std::endl;
		if (dictionary_type == "gabor" || dictionary_type == "gammatone") {
			out << "smallest time....... " << (T) pow (2, minj) / SR * 1000. << " ms" << std::endl;
			out << "frequency factor.... " << (T) pow (2., 1. / (T) oct_div) 
				<< " (1/" << oct_div <<  " oct)" << std::endl;
			out << "phase factor........ " << 2. * M_PI / phi_slices << std::endl;
		}
		if (dictionary_type == "gabor" || dictionary_type == "gammatone") {
			out << "highest frequency... " << freq_limit << " Hz" << std::endl << std::endl;
		}
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
   		if (tokens[0] == "npref") {
        	npref = atol (tokens[1].c_str ());
			if (SR < 8000) throw std::runtime_error ("invaild SR");
        }else if (tokens[0] == "SR") {
        	SR = atof (tokens[1].c_str ());
			if (SR < 8000) throw std::runtime_error ("invaild SR");
        } else if (tokens[0] == "J") {
        	J = atol (tokens[1].c_str ());
			if (J < 5) throw std::runtime_error ("invaild J");
        } else if (tokens[0] == "minj") {
        	minj = atol (tokens[1].c_str ());
			if (minj < 1 || minj > J) throw std::runtime_error ("invaild minj");
        } else if (tokens[0] == "components") {
        	comp = atol (tokens[1].c_str ());
			if (comp < 1) throw std::runtime_error ("invaild number of components");
        } else if (tokens[0] == "oct_divisions") {
        	oct_div = atol (tokens[1].c_str ());
			if (oct_div < 1) throw std::runtime_error ("invaild octave divisions");
        } 
        else if (tokens[0] == "phi_slices") {
        	phi_slices = atol (tokens[1].c_str ());
			if (phi_slices < 1) throw std::runtime_error ("invaild phase slices");
        } 
        else if (tokens[0] == "overlap") {
        	overlap = atol (tokens[1].c_str ());
			if (overlap < 1) throw std::runtime_error ("invaild overlap");
        } else if (tokens[0] == "freq_limit") {
        	freq_limit = atof (tokens[1].c_str ());
			if (freq_limit < 1 || freq_limit > SR / 2) throw std::runtime_error ("invaild frequency limit");
        } else if (tokens[0] == "dictionary") {
			dictionary_type = tokens[1];
			if (tokens.size () < 3) std::runtime_error ("invalid number of parameters in config file (frames)");		
			if (dictionary_type == "frames") {
				dictionary_path = tokens[2];
			}
			if (dictionary_type == "onsets") {
				if (tokens.size () < 4) std::runtime_error ("invalid number of parameters in config file (onsets)");
				dictionary_path = tokens[2];
				db_onset_threshold = atof (tokens[3].c_str ());
				db_onset_timegate = atof (tokens[4].c_str ());
				if (db_onset_threshold < 0 || db_onset_timegate < 0) throw std::runtime_error ("invaild db onset parameters");
			}	
			if (dictionary_type == "files") {
				dictionary_path = tokens[2];
			}		
        } else if (tokens[0] == "segmentation") {
			if (tokens.size () < 3) std::runtime_error ("invalid number of parameters in config file (segmentation)");		
        	target_onset_threshold = atof (tokens[1].c_str ());
			target_onset_timegate = atof (tokens[2].c_str ());
			if (target_onset_threshold < 0 || target_onset_timegate < 0) throw std::runtime_error ("invaild target onset parameters");
        } else if (tokens[0] == "ratio") {
        	ratio = atof (tokens[1].c_str ());
			if (ratio <= 0) throw std::runtime_error ("invaild ratio");
        } else if (tokens[0] == "stretch") {
        	stretch = atof (tokens[1].c_str ());
			if (stretch <= 0) throw std::runtime_error ("invaild stretch");
        } else {
            std::stringstream err;
            err << "invalid parameter " << tokens[0];
            throw std::runtime_error (err.str ());
        }
	}

	std::string dictionary_type;
	std::string dictionary_path;

	int npref;

	T SR;
	int J;
	T freq_limit; 

	int minj;	
	int oct_div;
	int phi_slices;
	
	T db_onset_threshold;
	T db_onset_timegate;
	T target_onset_threshold;
	T target_onset_timegate;

	int comp;
	int overlap;

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

// ------------------------------------------------------------------------------

void list_files (const char *path, std::vector<std::string>& files) {
    struct dirent *dp;
    DIR *dir = opendir (path);
    if (!dir) return; 

    while ((dp = readdir (dir)) != NULL) {
		if (strcmp(dp->d_name, ".") != 0 && strcmp(dp->d_name, "..") != 0) { // skip . and ..
        	std::string n (dp->d_name);
			files.push_back (n);
		}
    }

    closedir (dir);
}


#endif	// UTILS_H 

// EOF
