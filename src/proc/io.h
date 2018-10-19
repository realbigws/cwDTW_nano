#ifndef IO_H__
#define IO_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

namespace g{
namespace io{

bool ReadATCG(const char* name, std::vector<char>& genomes);

bool WriteATCG(const char* name, const std::vector<char>& genomes);

bool ReadSignalSequence(const char* name, std::vector<double>& signals);

bool ReadSignalSequence_int(const char* name, std::vector<int>& signals);

bool WriteSignalSequence(const char* name, const std::vector<double>& signals);

bool WriteSignalSequence_int(const char* name, const std::vector<int>& signals);

//---- write signal with name ----//
bool WriteSignalSequence_withName(const char* name,
        const std::vector<double>& signals,
        const std::vector<std::string>& kmer_rec);
bool WriteSignalSequence_int_withName(const char* name, 
	const std::vector<int>& signals,
	const std::vector<std::string>& kmer_rec);

}

}

#endif
