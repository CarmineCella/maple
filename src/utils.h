// utils.h
// 


#ifndef UTILS_H
#define UTILS_H 

#include <string>
#include <fstream>
#include <deque>
#include <sstream>
#include <stdexcept>

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
struct Params {

};

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
		SR = 44100; J = 12; minj = 5; comp = 100; oct_div = 1; normalize = 0; 
		windowing = false; freq_limit = 15000;
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
        } else if (tokens[0] == "oct_divisions") {
        	oct_div = atol (tokens[1].c_str ());
        } else if (tokens[0] == "normalize") {
        	normalize = (bool) atol (tokens[1].c_str ());
        } else if (tokens[0] == "windowing") {
        	windowing = (bool) atol (tokens[1].c_str ());
        } else if (tokens[0] == "freq_limit") {
        	freq_limit = atof (tokens[1].c_str ());
        } else if (tokens[0] == "dictionary") {
			dictionary = tokens[1];
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
	int oct_div;
	bool normalize; 
	bool windowing;
	T freq_limit; 
	std::string dictionary;
};

#endif	// UTILS_H 

// EOF
