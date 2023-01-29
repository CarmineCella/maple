// algorithms.h
// 


#ifndef ALGORITHMS_H
#define ALGORITHMS_H 

#include "Matrix.h"
#include "fourier.h"
#include "utils.h"
#include "WavFile.h"

#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <map>
#include <deque>
#include <cmath>

// TODO: 
// - save dictionary and analysis to reuse
// - add probabilistic generation
// - implement shift/stretch in frequency domain
// - improve segmentation 

//#define DEBUG_DECOMPOSITION
//#define USE_FAST_DOTPROD

#define ONSET_FADE_MS 10

// types
template <typename T>
struct Dictionary {
	DynamicMatrix<T> atoms;
	DynamicMatrix<T> parameters;
};
template <typename T>
using Decomposition = std::vector<DynamicMatrix<T> >;

// vector maths
template <typename T>
T norm(const T* values,
		int N) {
	if (1 > N) return 0;

	T sum = 0.;
	for (int i = 0; i < N; ++i) {
		sum += values[i] * values[i];
	}

	return sqrt (sum);
}
template <typename T>
void scale(const T* buff, T* out, int len, T factor) {
	for (int i = 0; i < len; ++i) {
		out[i] = buff[i] * factor;
	}
}
template <typename T>
inline T minimum(
		const T* values,
		int N, int& minPos) {
	if (1 > N) return 0;

	T min = values[0];
	for (int i = 0; i < N; ++i) {
		if (min > values[i]) {
			min = values[i];
			minPos = i;
		}
	}

	return min;
}
template <typename T>
inline T maximum(
		const T* values,
		int N, int& maxPos) {
	if (1 > N) return 0;

	T max = values[0];
	for (int i = 0; i < N; ++i) {
		if (max < values[i]) {
			max = values[i];
			maxPos = i;
		}
	}

	return max;
}
template <typename T>
T sum(
		const T* values,
		int N) {
	if (1 > N) return 0;

	T sum = 0.;
	for (int i = 0; i < N; ++i) {
		sum += values[i];
	}

	return sum;
}
template <typename T>
T mean(
		const T* values,
		int N) {
	if (1 > N) return 0;

	return sum(values, N) / N;
}

typedef float (*dot_function)(const float*, const float*, int);

template <typename T>
T dot_prod (const T* a, const T* b, int size) {
	// NB: it assumes vector have the same size
	T sum = 0;

	for (int i = 0; i < size; ++i) {
		sum += a[i] * b[i];
	}
	return sum;
	if (std::isinf(sum)) return 0;
	else return sum;
}

#ifdef USE_FAST_DOTPROD
#include <pmmintrin.h>
float dot_prod_sse (const float *a, const float *b, int len) {
	float total;
	int i;
	__m128 num1, num2, num3, num4;
	num4= _mm_setzero_ps();  			//sets sum to zero
	for(i=0; i<len; i+=4) {
		num1 = _mm_loadu_ps(a+i);   	//loads unaligned array a into num1  
										// num1= a[3]  a[2]  a[1]  a[0]
		num2 = _mm_loadu_ps(b+i);   	//loads unaligned array b into num2 
										// num2= b[3]   b[2]   b[1]  b[0]
		num3 = _mm_mul_ps(num1, num2); 	//performs multiplication   
										// num3 = a[3]*b[3]  a[2]*b[2]  a[1]*b[1]  a[0]*b[0]
		num3 = _mm_hadd_ps(num3, num3); //performs horizontal addition
										//num3 =  a[3]*b[3]+ a[2]*b[2]  a[1]*b[1]+a[0]*b[0]  a[3]*b[3]+ a[2]*b[2]  a[1]*b[1]+a[0]*b[0]
		num4 = _mm_add_ps(num4, num3);  //performs vertical addition
	}
	num4= _mm_hadd_ps(num4, num4);
	_mm_store_ss(&total,num4);
	return total;
}
#else

#define dot_prod_sse dot_prod
#endif

// others
template <typename T>
T frand (T min, T max) {
	return ((max - min) * ((T)rand() / RAND_MAX) + min);
}
template <typename T>
int wchoice (T* dist, int n) {
	T R = frand<T>(0., 1.);
	for (int i = 0; i < n; ++i) {
		if (dist[i] >= R) {
			return i;
		}
	}
	return (int) frand<double>(0, n);
}

// segmentation
template <typename T>
void get_onsets (const T* buffer, int samples,
	int bsize, int hop, T sr, T threshold, T timegate, std::vector<T>& onsets) {
	T* cdata = new T[bsize * 2];
	T* spectrum = new T[bsize];
	T* old_spectrum = new T[bsize];
	memset (old_spectrum, 0, sizeof (T) * bsize);
	AbstractFFT<T>* fft = createFFT<T>(bsize);

	T* win = new T[bsize];
	cosine_window<T>(win, bsize, .5, .5, 0.); // hanning

	std::vector<T> flux;
	int frames = 0;
	for (unsigned i = 0; i < samples; i += hop) {
		memset(cdata, 0, sizeof(T) * bsize * 2);

		int rsize = i + bsize > samples ? samples - i : bsize;
		for (unsigned j = 0; j < rsize; ++j) {
			cdata[2 * j] = buffer[i + j] * win[j]; // windowing
		}

		fft->forward (cdata);
		for (unsigned j = 0; j < bsize; ++j) {
			spectrum[j] = sqrt (cdata[j * 2] * cdata[j * 2] +
				cdata[j * 2 + 1] * cdata[j * 2 + 1]);

		}

		T v = specflux (spectrum, old_spectrum, bsize);
		flux.push_back (v);
		++frames;
	}

	int mpos = 0;
	T ma = maximum (&flux[0], flux.size(), mpos);
	scale<T>(&flux[0], &flux[0], flux.size (), 1. / ma);
	// save_vector("flux.txt", flux);

	std::vector<int> fluxpeaks;
	locmax(&flux[0], flux.size (), fluxpeaks);

	T prev_onset = 0;
	for (unsigned i = 0; i < fluxpeaks.size (); ++i) {
		if (flux[fluxpeaks[i]] > threshold) {
			T pos = fluxpeaks[i] * hop / sr;
			T dist = fabs (pos - prev_onset);
			if (dist > timegate || i == 0) {
				onsets.push_back(pos);
				prev_onset = pos;
			}
		}
	}

	if (onsets.size () == 0) onsets.push_back (0);

	delete [] cdata;
	delete [] spectrum;
	delete [] old_spectrum;
	delete [] win;
	delete fft;
}

// matching pursuit
template <typename T>
void store_vector (std::vector<T>& atom, std::vector<T> params,
	Dictionary<T>& dict) {
		int N = atom.size ();
		T nn = norm (&atom[0], N);
		scale<T> (&atom[0], &atom[0], (T) N, 1. / nn);
		dict.atoms.push_back (atom);
		dict.parameters.push_back (params);
}
template <typename T>
bool check_file (const std::string& name, WavInFile& in, const Parameters<T>& p) {
	if (in.getNumChannels () != 1) {
		std::stringstream msg;
		std::cerr << "cannot use " << name.c_str ()  << " (only mono files supported)" << std::endl;
		return false;
	}
	if (in.getSampleRate () != p.SR) {
		std::stringstream msg;
		std::cerr << "invalid SR for " << name << std::endl;
		return false;
	}
	return true;
}

template <typename T>
void make_dictionary (const Parameters<T>& p, Dictionary<T>& dict) {
	int N = pow (2., p.J);
	std::vector<T> buff (N);
	std::vector<T> cbuff (2 * N, 0);
	std::vector<T> hann (N, 1.);
	cosine_window<T>(&hann[0], N, .5, .5, 0.);

	if (p.dictionary_type == "frames") {
		std::vector<std::string> files;
		list_files (p.dictionary_path.c_str (), files);

		for (unsigned i = 0; i < files.size (); ++i) {
			if (files.at (i).find (".wav") != files.at (i).size () - 4) continue; // check WAV extension
			std::string name = p.dictionary_path + "/" + files.at (i);
			WavInFile in (name.c_str ());
			if (!check_file (name, in, p)) continue;
			while (!in.eof ()) {
				memset (&buff[0], 0, sizeof (T) * N);
				int r = in.read (&buff[0], N);
				if (r == N) {
					for (unsigned j = 0; j < r; ++j) buff[j] *= hann[j];
					store_vector(buff, std::vector<T> {0., (T) N, 0.}, dict);
				}
			}
		}
	} else if (p.dictionary_type == "onsets") {
		std::vector<std::string> files;
		list_files (p.dictionary_path.c_str (), files);

		for (unsigned i = 0; i < files.size (); ++i) {
			if (files.at (i).find (".wav") != files.at (i).size () - 4) continue; // check WAV extension
			std::string name = p.dictionary_path + "/" + files.at (i);
			WavInFile in (name.c_str ());
			if (!check_file (name, in, p)) continue;
			long nsamp = in.getNumSamples ();
			std::vector<T> data (nsamp);
			in.read (&data[0], nsamp);
			std::vector<T> onsets;
			get_onsets (&data[0], nsamp, N, N / p.overlap, p.SR, p.db_onset_threshold, p.db_onset_timegate, onsets);
			for (unsigned i = 0; i < onsets.size (); ++i) {
				int start =  (int) (onsets[i] * p.SR);
				int len = (int) (i == onsets.size () - 1 ? nsamp - start
					: (int) ((onsets[i + 1] - onsets[i]) * p.SR));
				std::vector<T> o (len);
				for (unsigned i = 0; i < len; ++i) o[i] = data[i + start];
				int fade_samps = (int) ((T) ONSET_FADE_MS / 1000.  * p.SR);
				if (len < 2 * fade_samps) { // onset too short
					 continue;
				}
				T env = 0.;

				T incr = 1. / fade_samps;
				for (int h = 0; h < fade_samps; ++h) { // fade in
					o[h] *= env;
					env += incr;
				}
				for (int h = len - fade_samps; h < len; ++h) { // fade out
					o[h] *= env;
					env -= incr;
				}

				store_vector(o, std::vector<T> {0., (T) len, 0.}, dict);
			}
		}
	}  else if (p.dictionary_type == "files") {
		std::vector<std::string> files;
		list_files (p.dictionary_path.c_str (), files);

		for (unsigned i = 0; i < files.size (); ++i) {
			if (files.at (i).find (".wav") != files.at (i).size () - 4) continue; // check WAV extension
			std::string name = p.dictionary_path + "/" + files.at (i);
			WavInFile in (name.c_str ());
			if (!check_file (name, in, p)) continue;
			long nsamp = in.getNumSamples ();
			std::vector<T> data (nsamp);
			in.read (&data[0], nsamp);
			store_vector(data, std::vector<T> {0., (T) nsamp, 0.}, dict);
		}
	} else if (p.dictionary_type == "gabor" || p.dictionary_type == "gammatone") {
		T comma = pow (2., 1. / p.oct_div);

		int j = p.minj;
		while (j <= p.J) {
			int n = pow (2., j);
			int u = 0;
			while (u <= (N - n)) {
				T fn = p.SR / (T) n;
				while (fn < p.freq_limit) {
					T phi_incr = 2. * M_PI / p.phi_slices; // PI?
					T phi = 0.;
					 while (phi < 2. * M_PI) {
						memset (&buff[0], 0, sizeof (T) * buff.size ());
						if (p.dictionary_type == "gabor") {
							gauss_window<T> (&buff[u], n, 4.);
						} else if (p.dictionary_type == "gammatone") {				
							gamma_window<T>(&buff[u], n, 1.5);
						}
						for (unsigned t = 0; t < n; ++t) {
							buff[t + u] *=  cos ((2. * M_PI * (T) t / (T) p.SR  * (T) fn) + phi);
						}	
						phi += phi_incr;
						store_vector(buff, std::vector<T>  {fn, (T) n, (T) u}, dict);
					} 
					fn *= comma;
				}
				u += n; // IS IT ENOUGH?
			}
			++j;
		}
	} else {
		throw std::runtime_error ("invalid dictionary type");
	}
}
template <typename T>
void decompose_segment (T sr, int components, const Dictionary<T>& dictionary,
	const std::vector<T>& target, int time_pos, dot_function dot,  DynamicMatrix<T>& frame) {
	std::vector<T> residual (target.size (), 0);
	for (unsigned i = 0; i < target.size (); ++i) {
		residual[i] = target[i];
	}
	T* input = &residual[0];

	T max_modulus = 0;
	T max_prod = 0;
	int max_index = 0;
	static int fnum = 0;
	for (unsigned i = 0; i < components; ++i) {		
		for (unsigned k = 0; k < dictionary.atoms.size (); ++k) {
			int n = dictionary.parameters[k][1];
			n = n > target.size () ? target.size () : n;
			int u = dictionary.parameters[k][2];
			const T* proj = &dictionary.atoms[k][u];
			T d = (dot (input + u, proj, n));
			T mod = fabs (d);

			if (mod > max_modulus) {
				max_modulus = mod;
				max_prod = d;
				max_index = k;
			}
		}
		int n = dictionary.parameters[max_index][1];
		n = n > target.size () ? target.size () : n;
		int u = dictionary.parameters[max_index][2];
		for (unsigned k = 0; k < n; ++k) {
			residual[k + u] -= (dictionary.atoms[max_index][k + u] * max_prod);
		}
		frame.push_back (std::vector<T> {(T) max_index, max_prod, (T) time_pos});
#ifdef DEBUG_DECOMPOSITION
		int sz = residual.size ();
		std::vector<T> outv (2 * sz);
		interleave (&outv[0], &residual[0], &dictionary.atoms[max_index][0], sz);

		std::stringstream name;
		name << "iteration_" << fnum << "_" << std::setw (3) << std::setfill ('0') << i << ".wav";
		WavOutFile ooo (name.str ().c_str (), 44100, 16, 2);
		ooo.write (&outv[0], sz * 2);
#endif
		max_modulus = 0;
	}
	++fnum;
}	

template <typename T>
void pursuit_decomposition (const Parameters<T>& p, const Dictionary<T>& dictionary, 
	const std::vector<T>& target, Decomposition<T>& decomposition, std::ostream& out) {
	if (dictionary.atoms.size () == 0) {
		throw std::runtime_error ("no atoms found in dictionary");
	}
	int ptr = 0;
	int N = pow (2., p.J);
	int hop = (int) ((T) (N / p.overlap));	
	std::vector<T> buffer (N);
	int r = target.size ();
	std::vector<T> onsets;	
	if (p.dictionary_type == "onsets" || p.dictionary_type == "files") {
		get_onsets (&target[0], r, N, N / p.overlap, p.SR, p.target_onset_threshold, p.target_onset_timegate, onsets);
		if (onsets.size () == 0) {
			throw std::runtime_error ("no onsets found in target sound");
		}
	}
	int ocntr = 0;
	while (ptr < r) {
		DynamicMatrix<T> frame;
		int perc = (int) ((T) ptr / r * 100.); 
		if (perc % 10 == 0) {
			out << "analysing........... " << perc << "%\r"; out.flush ();
		}		
		if (p.dictionary_type == "onsets" || p.dictionary_type == "files") {
			int start =  (int) (onsets[ocntr] * p.SR);
			int len = (int) (ocntr == onsets.size () - 1 ? r - start
				: (int) ((onsets[ocntr + 1] - onsets[ocntr]) * p.SR));

			buffer.resize (len, 0);
			for (unsigned i = 0; i < len; ++i) {
				if (i + len >= r) buffer[i] = 0;
				buffer[i] = target[i + start];
			}

			decompose_segment<T>(p.SR, p.comp, dictionary, buffer, ptr, &dot_prod, frame);			
			ptr += len;
			++ocntr;
			decomposition.push_back(frame);
			if (ocntr == onsets.size ()) break; 
		} else {
			for (int i  = 0; i < N; ++i) {
				if (i + ptr >= r) buffer[i] = 0;
				else buffer[i] = target[i + ptr];
			}
			decompose_segment<T>(p.SR, p.comp, dictionary, buffer, ptr, &dot_prod_sse, frame);			
			ptr += hop;
			decomposition.push_back(frame);
		}
	}

	if (decomposition.size () == 0) {
		throw std::runtime_error ("no elements found in the decomposition");
	}
}

template <typename T>
int reconstruct_segment (T ratio, const Dictionary<T>& dictionary, 
	const DynamicMatrix<T>& decomposition, std::vector<T>& output) {
	
	int max_sz = 0;
	for (unsigned i = 0; i < decomposition.size (); ++i) {
		int p = decomposition[i][0]; // position in dictionary
		if (max_sz < dictionary.atoms[p].size ()) max_sz = dictionary.atoms[p].size ();
	}
	max_sz /= ratio;
	output.resize (max_sz);
	memset (&output[0], 0, sizeof (T) * output.size ());	

	for (unsigned i = 0; i < decomposition.size (); ++i) {
		int p = decomposition[i][0]; // position in dictionary
		T w = decomposition[i][1];  // weight
		T* ptr = (T*) &dictionary.atoms[p][0];
		int sz = dictionary.atoms[p].size ();

		std::vector<T> buff;
		if (ratio != 1.) {
			int nsamp = (int) ((T) sz / ratio);
			buff.resize (nsamp);
			memset (&buff[0], 0, sizeof (T) * buff.size ());				
			T phi = 0;
			for (unsigned t = 0; t < nsamp; ++t) {
				int index = (int) phi;
				T frac = phi - index;
				int next = index + 1;
				if (next >= sz) break;
				buff[t] = dictionary.atoms[p][index] * (1. - frac) + dictionary.atoms[p][next] * frac;	
				phi += ratio;
				if (phi >= sz) break;
			}
			ptr = &buff[0];
			sz = nsamp;
		} 
		for (unsigned t = 0; t < sz; ++t) {
			output[t] +=  w * ptr[t];
		}
	}
	return max_sz;
}

template <typename T> 
void pursuit_reconstruction (const Parameters<T>& p, const Dictionary<T>& dictionary, 
	const Decomposition<T>& decomposition, 	std::vector<T>& output) {
	int samples = (int) ((T) p.stretch * decomposition.at (decomposition.size () - 1)[0][2]);
	int last_segment_len = (int) ((T) decomposition.at (decomposition.size () - 1)[0][0]);
	output.resize (samples + last_segment_len, 0);
	memset (&output[0], 0, sizeof (T) * samples);
	for (unsigned i = 0; i < decomposition.size (); ++i) {
		
		std::vector<T> buffer;
		int sz = reconstruct_segment (p.ratio, dictionary, decomposition[i], buffer);
		int ptr = (int) ((T) decomposition[i][0][2] * p.stretch); // time position (all components for each segment start together)
		for (int i  = 0; i < sz; ++i) {
			if (i + ptr >= output.size ()) break;
			output[i + ptr] += (buffer[i] / p.overlap);
		}
	}
}

template <typename T>
void export_decomposition_channel (T SR, int channel, const Dictionary<T>& dictionary, 
	const Decomposition<T>& decompositions) {
	std::stringstream name;
	name << "decomp_ch_" << std::setw (4) << std::setfill ('0') << channel << ".wav";
	WavOutFile chout (name.str ().c_str (), SR, 16, 1);
	for (unsigned i = 0; i < decompositions.size (); ++i) {
		int p = decompositions[i][channel][0]; // position
		std::vector<T> t (dictionary.atoms[p]);
		chout.write(&t[0], t.size ());
	}
}

// ------------------------------------------------------------------
typedef std::deque<int> Prefix;
typedef std::map<Prefix, std::vector<int> > StateTab;

// add: add word to suffix deque, update prefix
void add_state (Prefix& prefix, int s, int npref, StateTab& tab) {
	if (prefix.size() == npref) { 
		tab[prefix].push_back(s);
		prefix.pop_front();
	}
	prefix.push_back(s);
}

// build: read input words, build state table 
template <typename T>
void build_channel_transitions (const Decomposition<T>& decomposition, int channel, int npref, 
	Prefix& prefix, StateTab& tab) {
	for (unsigned i = 0; i < decomposition.size (); ++i) {
		const DynamicMatrix<float>& m = decomposition.at (i);
		if (channel < 0 || channel > m.size ()) {
			throw std::runtime_error ("invalid channel requested in generation");
		}
		int w = (int) m.at (channel).at (0); // 0-th is the atom index
		add_state (prefix, w, npref, tab);
	}	
}

template <typename T>
void build_transitions (const Parameters<T>& p, const Decomposition<T>& decomposition, int npref, 
	Prefix& prefix, std::vector<StateTab>& tabs) {
	for (unsigned i = 0; i < p.comp; ++i) {
		StateTab tab;
		build_channel_transitions (decomposition, i, npref, prefix, tab);
		tabs.push_back (tab);
	}
}

// generate: produce output
template <typename T>
int generate_channel (const Parameters<T>& p,  StateTab& tab, const Dictionary<T>& dict,
	Prefix& current,int channel, int& time_pos, DynamicMatrix<T>& frame) {
	std::vector<int>& suf = tab[current]; 
	while (suf.size () == 0) {
		StateTab::iterator item = tab.begin();
		int random_index = rand() % tab.size();
		std::advance(item, random_index);			
		suf = item->second;
		current = item->first;
	}
	assert (suf.size () > 0);
	int rp = rand() % suf.size();
	int w = suf[rp]; 
	int len = dict.parameters.at (w).at (1); // 1st is length
	T max_prod = 1;
	frame.push_back (std::vector<T> {(T) w, max_prod, (T) time_pos});
	current.pop_front();
	current.push_back(w);
	if (p.dictionary_type != "onsets") {
		int N = pow (2., p.J);
		int hop = (int) ((T) N / p.overlap);	
		return hop;
	} else return len;	
}
template <typename T>
void probabilistic_generation (const Parameters<T>& p, std::vector<StateTab>& tabs, 
	const Dictionary<T>& dict, int target_len,
	Decomposition<T>& generation) {
	generation.clear ();
	std::vector<Prefix> currents (p.comp);
	int time_pos = 0;
	
	while (time_pos < target_len) {
		DynamicMatrix<T> frame; // FIXME
		int max_len = 0;
		for (unsigned i = 0; i < p.comp; ++i) {	
			int l = generate_channel (p, tabs.at (i), dict, currents.at (i), i, time_pos, frame);
			if (l > max_len) max_len = l;
		}
		time_pos += max_len;
		if (frame.size () == p.comp) generation.push_back (frame);	
	}
}

// ----------------------------------------------------------------------------------
void export_transitions (std::ostream& out, std::vector<StateTab>& tabs) {
	for (unsigned h = 0; h < tabs.size (); ++h) {
		out << "CHANNEL " << h << std::endl;
		StateTab& tab = tabs.at (h);
		for (std::map<Prefix, std::vector<int> >::iterator i = tab.begin (); i != tab.end (); ++i) {
			out << "[";
			for (unsigned l = 0; l < i->first.size (); ++l) {
				out << i->first.at (l) << " ";
			}
			out << "] ";
			for (unsigned j = 0; j < i->second.size (); ++j) {
				out << i->second.at (j) << " ";
			}
			out << std::endl;
		}
		out << std::endl;
	}	
}

template <typename T>
void save_decomposition (const Parameters<T>& p, std::ostream& out, Decomposition<T>& decomposition) {
	for (unsigned j = 0; j < decomposition.size (); ++j) {
		for (unsigned i = 0; i < p.comp; ++i) {
			out << decomposition[j][i][0] << " " << decomposition[j][i][1] << " " << decomposition[j][i][2] << " ";
		}
		out << std::endl;
	}
	out << std::endl;
}

#endif	// ALGORITHMS_H 

// EOF
