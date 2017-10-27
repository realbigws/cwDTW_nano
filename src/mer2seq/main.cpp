#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <proc/io.h>
#include <map>
#include "util/exception.h"
#include "5mer/5mer_index.h"

using namespace std;

bool Genomes2SignalSequence(const std::vector<char>& genomes, std::vector<double>& signals, int scale)
{
	size_t bound = genomes.size()-5;//genomes.size()%5;
	for(size_t i = 0; i < bound; i++){
		int idx = g::Mer2Signal::FiveMer2Index(genomes[i], genomes[i+1], genomes[i+2], genomes[i+3], genomes[i+4]);
		double sigval = g::Mer2Signal::AvgSignalAt(idx);
		
		for(int c = scale; c--;){
			signals.push_back(sigval);
		}
	}
	
	for(size_t i = bound; i < genomes.size(); i++){
		for(int c = scale; c--;){
			signals.push_back(100);
		}
	}
}


//----------------------- 5mer index mapping --------------------//
int Char_To_Int(char c)
{
        switch(c)
        {
                case 'A':return 0;
                case 'C':return 1;
                case 'G':return 2;
                case 'T':return 3;
        }
}
char Int_To_Char(int i)
{
        switch(i)
        {
                case 0:return 'A';
                case 1:return 'C';
                case 2:return 'G';
                case 3:return 'T';
        }
}


int FiveMer_Mapping_Generate(map<string, int > &nam_mapping,vector <string> &anti_mapping)
{
        map<string, int>::iterator iter;
        nam_mapping.clear();
        anti_mapping.clear();
        int count=0;

        int i;
        int size=1024;
        for(i=0;i<1024;i++)
        {
                int a1 = i/256;
                int a2 = (i - a1*256)/64;
                int a3 = (i - a1*256 - a2*64)/16;
                int a4 = (i - a1*256 - a2*64 - a3*16)/4;
                int a5 = (i - a1*256 - a2*64 - a3*16 - a4*4);
                string fivemer="";
                fivemer=fivemer+Int_To_Char(a1)+Int_To_Char(a2)+Int_To_Char(a3)+Int_To_Char(a4)+Int_To_Char(a5);
		int idx = g::Mer2Signal::FiveMer2Index(fivemer);
                iter = nam_mapping.find(fivemer);
                if(iter != nam_mapping.end())
                {
                        fprintf(stderr,"duplicagted mapping %s \n",fivemer.c_str());
                        exit(-1);
                }
                count++;
                nam_mapping.insert(map < string, int >::value_type(fivemer, count));
                anti_mapping.push_back(fivemer);
		//---- output ----//
		printf("%s %d %d \n",fivemer.c_str(),idx,count);
	}
}


//------------------- main ----------------//
int main(int argc, char **argv)
{
	std::string input="";
	std::string output="";



//------ generate 5mer_index_mapping --------//
/*
	map<string, int > nam_mapping;
	vector <string> anti_mapping;
	FiveMer_Mapping_Generate(nam_mapping,anti_mapping);
	exit(-1);
*/

	struct options opts;
	if(GetOpts(argc, argv, &opts) < 0){
		return 0;
	}

	input=opts.input;
	output=opts.output;
	if(input=="" || output=="")
	{
		fprintf(stderr,"input or output is NULL \n");
		exit(-1);
	}

	
	EX_TRACE("Transform genomes to signal sequence...\n");
	
	std::vector<char> genomes;
	std::vector<double> signals;
	
	g::io::ReadATCG(opts.input, genomes);
	EX_TRACE("%ld genomes are readed.\n", genomes.size());
	
	Genomes2SignalSequence(genomes, signals, 1);
	
	g::io::WriteSignalSequence(opts.output, signals);
	
// 	genome::Mer2Signal::FiveMer2Index("ATAAA");
	
	return 0;
}

