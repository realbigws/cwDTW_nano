#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/cwdtw.h"
#include "util/exception.h"
#include "kmer/kmer_index.h"
#include <malloc.h>
#include <cmath>
#include <iomanip>

using namespace std;


//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}
//---------------------- utility ------//over



//--------- pore_models: from genome to expected signal ----------//
bool Genomes2SignalSequence(const std::vector<char>& genomes, 
	std::vector<int>& index, std::vector<double>& signals, int scale, int FIVE_or_SIX)
{
	long bound;
	if(FIVE_or_SIX==0) //-> 5mer model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-4;//genomes.size()%5;
		for(long i = 0; i < bound; i++){
			double sigval = g::Mer2Signal::AvgSignalAt_5mer(index[i]);
			sigval = 5.7*sigval+14;
			for(int c = scale; c--;){
				signals.push_back(sigval);
			}
		}
	}
	else
	{
		g::Mer2Signal::Genome2Index_6mer(genomes, index);
		bound = genomes.size()-5;//genomes.size()%5;
		for(long i = 0; i < bound; i++){
			double sigval = g::Mer2Signal::AvgSignalAt_6mer(index[i]);
			sigval = 5.7*sigval+14;
			for(int c = scale; c--;){
				signals.push_back(sigval);
			}
		}
	}

	//---- tail five_mer ------//	
	for(long i = bound; i < genomes.size(); i++){
		for(int c = scale; c--;){
			signals.push_back(5.7*100+14);
		}
	}

}


//--------- pore_models: from genome to expected signal (RNA) ----------//
bool Genomes2SignalSequence_RNA(const std::vector<char>& genomes, 
	std::vector<int>& index, std::vector<double>& signals, int scale, int mv200_or_mv180)
{
	long bound;
	if(mv200_or_mv180==1) //-> 200mv model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-4;//genomes.size()%5;
		for(long i = 0; i < bound; i++){
			double sigval = g::Mer2Signal::RnaSignalAt_5mer_200mv(index[i]);
			sigval = 5.7*sigval+14;
			for(int c = scale; c--;){
				signals.push_back(sigval);
			}
		}
	}
	else                  //-> 180mv model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-4;//genomes.size()%5;
		for(long i = 0; i < bound; i++){
			double sigval = g::Mer2Signal::RnaSignalAt_5mer_180mv(index[i]);
			sigval = 5.7*sigval+14;
			for(int c = scale; c--;){
				signals.push_back(sigval);
			}
		}
	}

	//---- tail five_mer ------//	
	for(long i = bound; i < genomes.size(); i++){
		for(int c = scale; c--;){
			signals.push_back(5.7*100+14);
		}
	}

}

//----------- write alignment to file (nano) -------------//__2017.10.15__//(Sheng modified)
void WriteSequenceAlignment_nano(const char* output,int KMER, 
	const std::vector<double>& reference, const std::vector<double>& peer,
	const std::vector<int>& refer_orig, const std::vector<int>& peer_orig,
	vector<pair<long,long> >& alignment, int swap, std::string &refer_str)
{
	vector <std::string> tmp_rec;
	double diff;
	for(long i = 0; i < alignment.size(); i++)
	{
		//----- output to string ----//
		std::ostringstream o;
		std::string sub_str;
		diff = std::fabs(reference[alignment[i].first]-peer[alignment[i].second]);
		if(swap==0)
		{
			o<<setw(5)<<peer_orig[alignment[i].second]<<" "<<setw(5)<<refer_orig[alignment[i].first]<<" | ";
			o<<setw(10)<<alignment[i].second+1<<" "<<setw(10)<<alignment[i].first+1<<" | ";
			o<<setw(15)<<peer[alignment[i].second]<<" "<<setw(15)<<reference[alignment[i].first];
			//-- judge --//
			if(alignment[i].first>refer_str.size()-KMER)break;
			sub_str=refer_str.substr(alignment[i].first,KMER);
		}
		else
		{
			o<<setw(5)<<peer_orig[alignment[i].first]<<" "<<setw(5)<<refer_orig[alignment[i].second]<<" | ";
			o<<setw(10)<<alignment[i].first+1<<" "<<setw(10)<<alignment[i].second+1<<" | ";
			o<<setw(15)<<reference[alignment[i].first]<<" "<<setw(15)<<peer[alignment[i].second];
			//-- judge --//
			if(alignment[i].second>refer_str.size()-KMER)break;
			sub_str=refer_str.substr(alignment[i].second,KMER);
		}
		o<<"          diff:"<<setw(15)<<diff;
		o<<"   "<<sub_str;
		//----- record string -----//
		std::string s=o.str();
		tmp_rec.push_back(s);
	}
	//----- output to file ------//
	FILE *fp=fopen(output,"wb");
	for(long i=0;i<(long)tmp_rec.size();i++)fprintf(fp,"%s\n",tmp_rec[i].c_str());
	fclose(fp);
}


//----------------- main -------------------//
int main(int argc, char **argv)
{
	struct options opts;
	opts.radius  = 50;
	opts.level   = 3;
	opts.scale0  = sqrt(2);
	opts.verbose = 0;       //-> [0] no verbose; 1 verbose
	opts.test    = 0;       //-> [0] not use test mode; 1 equal_ave, 2 peak_ave, 3 Fast_DTW
	opts.mode    = 0;       //-> [0] block bound; 1 diagonol bound
	opts.kmer    = 1;       //-> [1] to use 6mer; 0 to use 6mer
	opts.rna     = 0;       //-> [0] to use DNA; 1 to use 200mv RNA; -1 to use 180mv RNA
	

	//----- parse arguments -----//
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return -1;
	}

	std::string input1=opts.input;
	std::string input2=opts.peer;
	std::string output=opts.output;
	if(input1=="" || input2=="")
	{
		fprintf(stderr,"input1 or input2 is NULL \n");
		return -1;
	}

	//----- determine Kmer -------//
	int KMER;
	if(opts.kmer==0)KMER=5;
	else if(opts.kmer==1)KMER=6;
	else
	{
		fprintf(stderr,"kmer should be 5 or 6 \n");
		return -1;
	}

	//------ if RNA model is use, then Kmer should be fixed to 5mer -----//
	if(opts.rna!=0)
	{
		KMER=5;
	}

	//======================= START Procedure ===================================//

	//=========================================//
	//----- 1. read genome translated signal -----//
	std::vector<char> genomes;   //genome sequence
	if(!g::io::ReadATCG(opts.input, genomes)){
		EX_TRACE("Cannot open %s.\n", opts.input);
		return -1;
	}
	std::string genome_str(genomes.begin(), genomes.end() );
	//----- 1.1 pore_model: transform genome sequence to expected signal -------//
	std::vector<int> refer_orig;
	std::vector<double> reference;  //reference: genome signal
	if(opts.rna==0)Genomes2SignalSequence(genomes, refer_orig, reference, 1, opts.kmer);
	else Genomes2SignalSequence_RNA(genomes, refer_orig, reference, 1, opts.rna);
	//---- get nanopore reference name ----//start
	std::string genom_name_orig=opts.input;
	std::string genom_name;
	getBaseName(genom_name_orig,genom_name,'/','.');
	//---- get nanopore reference name ----//end


	//========================================//
	//----- 2. read nanopore raw signal -----//
	std::vector<int> peer_orig;             //peer: nanopore signal
	if(!g::io::ReadSignalSequence_int(opts.peer, peer_orig)){
		EX_TRACE("Cannot open %s.\n", opts.peer);
		return -1;
	}
	std::vector<double> peer(peer_orig.begin(), peer_orig.end());
	//---- get nanopore raw signal name ----//start
	std::string signal_name_orig=opts.peer;
	std::string signal_name;
	getBaseName(signal_name_orig,signal_name,'/','.');
	//---- get nanopore raw signal name ----//end


	//----- length check ------//
	int swap=0;
	if(reference.size()>peer.size())
	{
		fprintf(stderr," reference DNA sequence length %d larger than raw signal length %d \n",
			reference.size(),peer.size());
		return -1;
	}


	//==================================================//
	//------3. process initial input signals ----------//

	//----- 3.1 Zscore normaliza on both signals -----//
	g::proc::ZScoreNormalize(reference);
	g::proc::ZScoreNormalize(peer);

	//----- 3.2 calculate length ratio between input signals -----//	
	double alpha = (double)peer.size()/reference.size();


	//====================================================//
	//----- 4. continous wavelet transform --------------//
	std::vector<std::vector<double> > rcwt, pcwt;

	if(opts.verbose ==1){
		EX_TRACE("CWT Analysis...\n");
	}
	
	long npyr = opts.level;          // default: 3
	double scale0 = opts.scale0;	// default: sqrt(2)
	double dscale = 1;              // default: 1

	g::cwdtw::CWTAnalysis(reference, rcwt, scale0, dscale, npyr);	
	g::cwdtw::CWTAnalysis(peer, pcwt, scale0*alpha, dscale, npyr);

	//------ 4.1 Zscore normaliza on both CWT signals -----//	
	//if multiscale is used, pyr logical should be added.
	for(long i = 0; i < npyr; i++){
		g::proc::ZScoreNormalize(rcwt[i]);
		g::proc::ZScoreNormalize(pcwt[i]);
	}


	//============================================//
	//------ 5. multi-level WaveletDTW ----------//
	std::vector<std::pair<long,long> > cosali;
	double tdiff;

	if(opts.verbose ==1){
		EX_TRACE("Coarse Alignment...\n");
	}
	g::cwdtw::MultiLevel_WaveletDTW(reference, peer, rcwt, pcwt, cosali, opts.radius, opts.test, opts.mode, &tdiff);
	if(opts.verbose ==1){
		EX_TRACE("Average Deviation (%.1lf/%ld=%.3lf)\n", tdiff, cosali.size(), tdiff/cosali.size());
	}

	//------ 5.1 generate final boundary -------//
	std::vector<std::pair<long,long> > bound;
	g::cwdtw::BoundGeneration(cosali, opts.radius, bound, opts.mode);

	//------ 5.2 generate final alignment via cDTW ------//
	std::vector<std::pair<long,long> > alignment;
	tdiff = g::proc::BoundDynamicTimeWarping(reference, peer, bound, alignment);  //-> NOT use restrict version !!!
	fprintf(stderr,"%s %s %lf %d %lf\n",signal_name.c_str(),genom_name.c_str(),tdiff,alignment.size(),tdiff/alignment.size());


	//=================================================//
	//------ 6. output final alignment to file -------//
	if(output!="")
		WriteSequenceAlignment_nano(opts.output, KMER, reference, peer, refer_orig, peer_orig, alignment, swap, genome_str);

	//----- exit -----//	
	return 0;
}

