// WavFile.h
//


#ifndef WAVFILE_H
#define WAVFILE_H

#include <stdexcept>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <limits.h>

#ifndef uint
	typedef unsigned int uint;
#endif

const static char riffStr[] = "RIFF";
const static char waveStr[] = "WAVE";
const static char fmtStr[]  = "fmt ";
const static char dataStr[] = "data";

typedef struct {
	char riff_char[4];
	int  package_len;
	char wave[4];
} WavRiff;

typedef struct {
	char  fmt[4];
	int   format_len;
	short fixed;
	short channel_number;
	int   sample_rate;
	int   byte_rate;
	short byte_per_sample;
	short bits_per_sample;
} WavFormat;

typedef struct {
	char  data_field[4];
	uint  data_len;
} WavData;


typedef struct {
	WavRiff   riff;
	WavFormat format;
	WavData   data;
} WavHeader;


#ifdef BYTE_ORDER
	#if BYTE_ORDER == BIG_ENDIAN
		#define _BIG_ENDIAN_
	#endif
#endif

#ifdef _BIG_ENDIAN_
	static inline void _swap32 (unsigned int &dwData) {
		dwData = ((dwData >> 24) & 0x000000FF) |
				 ((dwData >> 8)  & 0x0000FF00) |
				 ((dwData << 8)  & 0x00FF0000) |
				 ((dwData << 24) & 0xFF000000);
	}
	
	static inline void _swap16 (unsigned short &wData) {
		wData = ((wData >> 8) & 0x00FF) |
				((wData << 8) & 0xFF00);
	}
	
	static inline void _swap16Buffer (unsigned short *pData, unsigned int dwNumWords) {
		unsigned long i;
	
		for (i = 0; i < dwNumWords; i ++) {
			_swap16 (pData[i]);
		}
	}

#else   // BIG_ENDIAN
	static inline void _swap32 (unsigned int &dwData) {
		// do nothing
	}
	
	static inline void _swap16 (unsigned short &wData) {
		// do nothing
	}
	
	
	static inline void _swap16Buffer (unsigned short *pData, unsigned int dwNumBytes) {
		// do nothing
	}

#endif  // BIG_ENDIAN

// test if character code is between a white space ' ' and little 'z'
static int isAlpha (char c) {
	return (c >= ' ' && c <= 'z') ? 1 : 0;
}


// test if all characters are between a white space ' ' and little 'z'
static int isAlphaStr (char *str) {
	int c;

	c = str[0];
	while (c) {
		if (isAlpha(c) == 0) return 0;
		str ++;
		c = str[0];
	}

	return 1;
}

class WavInFile {
public:	
	WavInFile (const char *fileName) {
		int hdrsOk;
	
		// Try to open the file for reading
		fptr = fopen(fileName, "rb");
		if (fptr == NULL) {
			// didn't succeed
			std::string msg = "unable to open file \"";
			msg += fileName;
			msg += "\" for reading.";
			throw std::runtime_error(msg);
		}
	
		// Read the file headers
		hdrsOk = readWavHeaders();
		if (hdrsOk != 0) {
			// Something didn't match in the wav file headers
			std::string msg = "File \"";
			msg += fileName;
			msg += "\" is corrupt or not a WAV file";
			throw std::runtime_error(msg);
		}
	
		if (header.format.fixed != 1) {
			std::string msg = "File \"";
			msg += fileName;
			msg += "\" uses unsupported encoding.";
			throw std::runtime_error(msg);
		}
	
		dataRead = 0;
	}
	
	
	
	~WavInFile() {
		close();
	}
	
	
	
	void rewind() {
		int hdrsOk;
	
		fseek(fptr, 0, SEEK_SET);
		hdrsOk = readWavHeaders();
		assert(hdrsOk == 0);
		dataRead = 0;
	}
	
	
	int checkCharTags() {
		// header.format.fmt should equal to 'fmt '
		if (memcmp(fmtStr, header.format.fmt, 4) != 0) return -1;
		// header.data.data_field should equal to 'data'
		if (memcmp(dataStr, header.data.data_field, 4) != 0) return -1;
	
		return 0;
	}
	
	
	int read(char *buffer, int maxElems) {
		int numBytes;
		uint afterDataRead;
	
		// ensure it's 8 bit format
		if (header.format.bits_per_sample != 8) {
			throw std::runtime_error("read(char*, int) works only with 8bit samples.");
		}
		assert(sizeof(char) == 1);
	
		numBytes = maxElems;
		afterDataRead = dataRead + numBytes;
		if (afterDataRead > header.data.data_len) {
			// Don't read more samples than are marked available in header
			numBytes = header.data.data_len - dataRead;
			assert(numBytes >= 0);
		}
	
		numBytes = fread(buffer, 1, numBytes, fptr);
		dataRead += numBytes;
	
		return numBytes;
	}
	
	
	int read(short *buffer, int maxElems) {
		unsigned int afterDataRead;
		int numBytes;
		int numElems;
	
		if (header.format.bits_per_sample == 8) {
			// 8 bit format
			char *temp = new char[maxElems];
			int i;
	
			numElems = read(temp, maxElems);
			// convert from 8 to 16 bit
			for (i = 0; i < numElems; i ++) {
				buffer[i] = temp[i] << 8;
			}
			delete[] temp;
		} else {
			// 16 bit format
			assert(header.format.bits_per_sample == 16);
			assert(sizeof(short) == 2);
	
			numBytes = maxElems * 2;
			afterDataRead = dataRead + numBytes;
			if (afterDataRead > header.data.data_len) {
				// Don't read more samples than are marked available in header
				numBytes = header.data.data_len - dataRead;
				assert(numBytes >= 0);
			}
	
			numBytes = fread(buffer, 1, numBytes, fptr);
			dataRead += numBytes;
			numElems = numBytes / 2;
	
			// 16bit samples, swap byte order if necessary
			_swap16Buffer((unsigned short *)buffer, numElems);
		}
	
		return numElems;
	}
	
	
	
	int read(float *buffer, int maxElems) {
		short *temp = new short[maxElems];
		int num;
		int i;
		double fscale;
	
		num = read(temp, maxElems);
	
		fscale = 1.0 / 32768.0;
		// convert to floats, scale to range [-1..+1[
		for (i = 0; i < num; i ++) {
			buffer[i] = (float)(fscale * (double)temp[i]);
		}
	
		delete[] temp;
	
		return num;
	}
	
	int read(double *buffer, int maxElems) {
		short *temp = new short[maxElems];
		int num;
		int i;
		double fscale;
	
		num = read(temp, maxElems);
	
		fscale = 1.0 / 32768.0;
		// convert to doubles, scale to range [-1..+1[
		for (i = 0; i < num; i ++) {
			buffer[i] = (double)(fscale * (double)temp[i]);
		}
	
		delete[] temp;
	
		return num;
	}
	
	int eof() const {
		// return true if all data has been read or file eof has reached
		return (dataRead == header.data.data_len || feof(fptr));
	}
	
	
	void close() {
		fclose(fptr);
		fptr = NULL;
	}
	
	
	int readRIFFBlock() {
		fread(&(header.riff), sizeof(WavRiff), 1, fptr);
	
		// swap 32bit data byte order if necessary
		_swap32((unsigned int &)header.riff.package_len);
	
		// header.riff.riff_char should equal to 'RIFF');
		if (memcmp(riffStr, header.riff.riff_char, 4) != 0) return -1;
		// header.riff.wave should equal to 'WAVE'
		if (memcmp(waveStr, header.riff.wave, 4) != 0) return -1;
	
		return 0;
	}
	
	int readHeaderBlock() {
		char label[5];
		std::string sLabel;
	
		// lead label std::string
		fread(label, 1, 4, fptr);
		label[4] = 0;
	
		if (isAlphaStr(label) == 0) return -1;    // not a valid label
	
		// Decode blocks according to their label
		if (strcmp(label, fmtStr) == 0) {
			int nLen, nDump;
	
			// 'fmt ' block
			memcpy(header.format.fmt, fmtStr, 4);
	
			// read length of the format field
			fread(&nLen, sizeof(int), 1, fptr);
			// swap byte order if necessary
			_swap32((unsigned int &)nLen); // int format_len;
			header.format.format_len = nLen;
	
			// calculate how much length differs from expected
			nDump = nLen - (sizeof(header.format) - 8);
	
			// if format_len is larger than expected, read only as much data as we've space for
			if (nDump > 0) {
				nLen = sizeof(header.format) - 8;
			}
	
			// read data
			fread(&(header.format.fixed), nLen, 1, fptr);
	
			// swap byte order if necessary
			_swap16((unsigned short &)header.format.fixed);            // short int fixed;
			_swap16((unsigned short &)header.format.channel_number);   // short int channel_number;
			_swap32((unsigned int   &)header.format.sample_rate);      // int sample_rate;
			_swap32((unsigned int   &)header.format.byte_rate);        // int byte_rate;
			_swap16((unsigned short &)header.format.byte_per_sample);  // short int byte_per_sample;
			_swap16((unsigned short &)header.format.bits_per_sample);  // short int bits_per_sample;
	
			// if format_len is larger than expected, skip the extra data
			if (nDump > 0) {
				fseek(fptr, nDump, SEEK_CUR);
			}
	
			return 0;
		} else if (strcmp(label, dataStr) == 0) {
			// 'data' block
			memcpy(header.data.data_field, dataStr, 4);
			fread(&(header.data.data_len), sizeof(uint), 1, fptr);
	
			// swap byte order if necessary
			_swap32((unsigned int &)header.data.data_len);
	
			return 1;
		} else {
			uint len, i;
			uint temp;
			// unknown block
	
			// read length
			fread(&len, sizeof(len), 1, fptr);
			// scan through the block
			for (i = 0; i < len; i ++) {
				fread(&temp, 1, 1, fptr);
				if (feof(fptr)) return -1;   // unexpected eof
			}
		}
		return 0;
	}


	int readWavHeaders() {
		int res;
	
		memset(&header, 0, sizeof(header));
	
		res = readRIFFBlock();
		if (res) return 1;
		// read header blocks until data block is found
		do {
			// read header blocks
			res = readHeaderBlock();
			if (res < 0) return 1;  // error in file structure
		} while (res == 0);
		// check that all required tags are legal
		return checkCharTags();
	}
	
	
	uint getNumChannels() const {
		return header.format.channel_number;
	}
	
	
	uint getNumBits() const {
		return header.format.bits_per_sample;
	}
	
	
	uint getBytesPerSample() const {
		return getNumChannels() * getNumBits() / 8;
	}
	
	
	uint getSampleRate() const {
		return header.format.sample_rate;
	}
	
	
	
	uint getDataSizeInBytes() const {
		return header.data.data_len;
	}
	
	
	uint getNumSamples() const {
		return header.data.data_len / header.format.byte_per_sample;
	}
	
	
	uint getLengthMS() const {
		uint numSamples;
		uint sampleRate;
	
		numSamples = getNumSamples();
		sampleRate = getSampleRate();
	
		assert(numSamples < UINT_MAX / 1000);
		return (1000 * numSamples / sampleRate);
	}
private:
	FILE *fptr;
	uint dataRead;
	WavHeader header;
};


/// Class for writing WAV audio files.
class WavOutFile {
private:
	/// Pointer to the WAV file
	FILE *fptr;

	/// WAV file header data.
	WavHeader header;

	/// Counter of how many bytes have been written to the file so far.
	int bytesWritten;

public:
	//////////////////////////////////////////////////////////////////////////////
	//
	// Class WavOutFile
	//

	WavOutFile(const char *fileName, int sampleRate, int bits, int channels) {
		bytesWritten = 0;
		fptr = fopen(fileName, "wb");
		if (fptr == NULL) {
			std::string msg = "unable to open file \"";
			msg += fileName;
			msg += "\" for writing.";
			//pmsg = msg.c_str;
			throw std::runtime_error(msg);
		}
	
		fillInHeader(sampleRate, bits, channels);
		writeHeader();
	}
	
	
	
	~WavOutFile() {
		close();
	}
	
	
	
	void fillInHeader(uint sampleRate, uint bits, uint channels) {
		// fill in the 'riff' part..
	
		// copy std::string 'RIFF' to riff_char
		memcpy(&(header.riff.riff_char), riffStr, 4);
		// package_len unknown so far
		header.riff.package_len = 0;
		// copy std::string 'WAVE' to wave
		memcpy(&(header.riff.wave), waveStr, 4);
	
	
		// fill in the 'format' part..
	
		// copy std::string 'fmt ' to fmt
		memcpy(&(header.format.fmt), fmtStr, 4);
	
		header.format.format_len = 0x10;
		header.format.fixed = 1;
		header.format.channel_number = (short)channels;
		header.format.sample_rate = sampleRate;
		header.format.bits_per_sample = (short)bits;
		header.format.byte_per_sample = (short)(bits * channels / 8);
		header.format.byte_rate = header.format.byte_per_sample * sampleRate;
		header.format.sample_rate = sampleRate;
	
		// fill in the 'data' part..
	
		// copy std::string 'data' to data_field
		memcpy(&(header.data.data_field), dataStr, 4);
		// data_len unknown so far
		header.data.data_len = 0;
	}
	
	
	void finishHeader() {
		// supplement the file length into the header structure
		header.riff.package_len = bytesWritten + 36;
		header.data.data_len = bytesWritten;
	
		writeHeader();
	}
	
	
	
	void writeHeader() {
		WavHeader hdrTemp;
	
		// swap byte order if necessary
		hdrTemp = header;
		_swap32((unsigned int   &)hdrTemp.riff.package_len);
		_swap32((unsigned int   &)hdrTemp.format.format_len);
		_swap16((unsigned short &)hdrTemp.format.fixed);
		_swap16((unsigned short &)hdrTemp.format.channel_number);
		_swap32((unsigned int   &)hdrTemp.format.sample_rate);
		_swap32((unsigned int   &)hdrTemp.format.byte_rate);
		_swap16((unsigned short &)hdrTemp.format.byte_per_sample);
		_swap16((unsigned short &)hdrTemp.format.bits_per_sample);
		_swap32((unsigned int   &)hdrTemp.data.data_len);
	
		// write the supplemented header in the beginning of the file
		fseek(fptr, 0, SEEK_SET);
		fwrite(&hdrTemp, sizeof(hdrTemp), 1, fptr);
		// jump back to the end of the file
		fseek(fptr, 0, SEEK_END);
	}
	
	
	
	void close() {
		finishHeader();
		fclose(fptr);
		fptr = NULL;
	}
	
	
	void write(const char *buffer, int numElems) {
		int res;
	
		if (header.format.bits_per_sample != 8) {
			throw std::runtime_error("write(const char*, int) accepts only 8bit samples.");
		}
		assert(sizeof(char) == 1);
	
		res = fwrite(buffer, 1, numElems, fptr);
		if (res != numElems) {
			throw std::runtime_error("problem while writing to a wav file.");
		}
	
		bytesWritten += numElems;
	}
	
	
	void write(const short *buffer, int numElems) {
		int res;
	
		// 16 bit samples
		if (numElems < 1) return;   // nothing to do
	
		if (header.format.bits_per_sample == 8) {
			int i;
			char *temp = new char[numElems];
			// convert from 16bit format to 8bit format
			for (i = 0; i < numElems; i ++) {
				temp[i] = buffer[i] >> 8;
			}
			// write in 8bit format
			write(temp, numElems);
			delete[] temp;
		} else {
			// 16bit format
			unsigned short *pTemp = new unsigned short[numElems];
	
			assert(header.format.bits_per_sample == 16);
	
			// allocate temp buffer to swap byte order if necessary
			memcpy(pTemp, buffer, numElems * 2);
			_swap16Buffer(pTemp, numElems);
	
			res = fwrite(pTemp, 2, numElems, fptr);
	
			delete[] pTemp;
	
			if (res != numElems) {
				throw std::runtime_error("problem while writing to a wav file.");
			}
			bytesWritten += 2 * numElems;
		}
	}
	
	
	void write(const float *buffer, int numElems) {
		int i;
		short *temp = new short[numElems];
		int iTemp;
	
		// convert to 16 bit integer
		for (i = 0; i < numElems; i ++) {
			// convert to integer
			iTemp = (int)(32768.0f * buffer[i]);
	
			// saturate
			if (iTemp < -32768) iTemp = -32768;
			if (iTemp > 32767)  iTemp = 32767;
			temp[i] = (short)iTemp;
		}
	
		write(temp, numElems);
	
		delete[] temp;
	}
	
	void write(const double *buffer, int numElems) {
		int i;
		short *temp = new short[numElems];
		int iTemp;
	
		// convert to 16 bit integer
		for (i = 0; i < numElems; i ++) {
			// convert to integer
			iTemp = (int)(32768.0f * buffer[i]);
	
			// saturate
			if (iTemp < -32768) iTemp = -32768;
			if (iTemp > 32767)  iTemp = 32767;
			temp[i] = (short)iTemp;
		}
	
		write(temp, numElems);
	
		delete[] temp;
	}

};

#endif
