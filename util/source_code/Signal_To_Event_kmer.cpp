#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
//#include "FiveMer_Data.h"
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

//----- Z-score ----//
void Z_score(vector <double> &in,double &mean,double &vari)
{
	//init
	mean=0;
	vari=1;
	//proc
	int i;
	int size=(int)in.size();
	if(size==0)return;
	//-> calculate mean
	for(i=0;i<size;i++)mean+=in[i];
	mean/=size;
	//-> calculate vari
	vari=0;
	for(i=0;i<size;i++)vari+=(in[i]-mean)*(in[i]-mean);
	vari=1.0*sqrt(vari/size);
	//-> calculate Z-score
	for(i=0;i<size;i++)in[i]=(in[i]-mean)/vari;
}


//==================== load KMer data =======================//
//-> format (5mer or 6mer)
/*
0 AAAAA 85.083612
1 AAAAC 76.635809
2 AAAAG 83.489101
3 AAAAT 75.716163
4 AAACA 87.429392
5 AAACC 84.244828
....
*/

int Load_KMer_Data(string &input,vector <double> &kmer_data)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",input.c_str());
		exit(-1);
	}
	//load
	int count=0;
	kmer_data.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		double value;
		www>>temp>>temp>>value;
		kmer_data.push_back(value);
		count++;
	}
	//return
	return count;
}


//==================== kmer mapping ==================//

//-> generate kmer mapping
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
	//mapping
	map<string, int>::iterator iter;
	nam_mapping.clear();
	anti_mapping.clear();
	int count=0;
	//process
	int i;
	int size=1024;
	for(i=0;i<1024;i++)
	{
		//decode
		int a1 = i/256;
		int a2 = (i - a1*256)/64;
		int a3 = (i - a1*256 - a2*64)/16;
		int a4 = (i - a1*256 - a2*64 - a3*16)/4;
		int a5 = (i - a1*256 - a2*64 - a3*16 - a4*4);
		//get string
		string fivemer="";
		fivemer=fivemer+Int_To_Char(a1)+Int_To_Char(a2)+Int_To_Char(a3)+Int_To_Char(a4)+Int_To_Char(a5);
		//mapping
		iter = nam_mapping.find(fivemer);
		if(iter != nam_mapping.end())
		{
			fprintf(stderr,"duplicagted mapping %s \n",fivemer.c_str());
			exit(-1);
		}
		//count++
		count++;
		nam_mapping.insert(map < string, int >::value_type(fivemer, count));
		anti_mapping.push_back(fivemer);
	}
	//return
	return count;
}
int SixMer_Mapping_Generate(map<string, int > &nam_mapping,vector <string> &anti_mapping)
{
	//mapping
	map<string, int>::iterator iter;
	nam_mapping.clear();
	anti_mapping.clear();
	int count=0;
	//process
	int i;
	int size=4096;
	for(i=0;i<4096;i++)
	{
		//decode
		int a1 = i/1024;
		int a2 = (i - a1*1024)/256;
		int a3 = (i - a1*1024 - a2*256)/64;
		int a4 = (i - a1*1024 - a2*256 - a3*64)/16;
		int a5 = (i - a1*1024 - a2*256 - a3*64 - a4*16)/4;
		int a6 = (i - a1*1024 - a2*256 - a3*64 - a4*16 - a5*4);
		//get string
		string sixmer="";
		sixmer=sixmer+Int_To_Char(a1)+Int_To_Char(a2)+Int_To_Char(a3)+Int_To_Char(a4)+Int_To_Char(a5)+Int_To_Char(a6);
		//mapping
		iter = nam_mapping.find(sixmer);
		if(iter != nam_mapping.end())
		{
			fprintf(stderr,"duplicagted mapping %s \n",sixmer.c_str());
			exit(-1);
		}
		//count++
		count++;
		nam_mapping.insert(map < string, int >::value_type(sixmer, count));
		anti_mapping.push_back(sixmer);
	}
	//return
	return count;
}

//-> get mapping
int Get_Mapping(map<string, int > &nam_mapping,string &kmer)
{
	map<string, int>::iterator iter;
	iter = nam_mapping.find(kmer);
	if(iter == nam_mapping.end())
	{
		fprintf(stderr,"failed mapping %s \n",kmer.c_str());
		return -1;
	}
	int retv=nam_mapping[kmer]-1;
	return retv;
}
int Get_Mapping_NoWarn(map<string, int > &nam_mapping,string &kmer)
{
	map<string, int>::iterator iter;
	iter = nam_mapping.find(kmer);
	if(iter == nam_mapping.end())return -1;
	int retv=nam_mapping[kmer]-1;
	return retv;
}


//===================== process signal_label data ================//
//-------- read signal_label ----------//
//-> format
/*
478 703
571 703
531 703
542 767
536 767
541 767
534 767
549 767
545 767
548 767
534 767
528 767
521 767
546 767
.....
*/
int Read_Gold_Standard(string &input, vector <double> &signal, vector <int> &label)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",input.c_str());
		exit(-1);
	}
	//load
	signal.clear();
	label.clear();
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		double raw_signal;
		int code;
		www>>raw_signal>>code;
		signal.push_back(raw_signal);
		label.push_back(code);
		//count++
		count++;
	}
	//return
	return count;
}
int Remove_Blank_HeadTail(vector <double> &signal, vector <int> &label,
	vector <double> &signal_, vector <int> &label_)
{
	int i;
	int size=(int)signal.size();
	//check head
	int head=0;
	for(i=0;i<size;i++)
	{
		if(label[i]!=-1)
		{
			head=i;
			break;
		}
	}
	//check tail
	int tail=0;
	for(i=size-1;i>=0;i--)
	{
		if(label[i]!=-1)
		{
			tail=i;
			break;
		}
	}
	//assign
	signal_.resize(tail-head+1);
	label_.resize(tail-head+1);
	int cur=0;
	for(i=head;i<=tail;i++)
	{
		signal_[cur]=signal[i];
		label_[cur]=label[i];
		cur++;
	}
	return cur;
}

//------------------ Signal_To_Event ------------------//
void Signal_To_Event(vector <double> &signal,vector <int> &label, 
	vector <int> &MoveOrNot, vector <int> &label_,
	vector <double> &signal_, vector <double> &signal_var, 
	vector <double> &head_mean, vector <double> &head_vari,
	vector <int> &durate_ )
{
	int i;
	int size=(int)label.size();
	MoveOrNot.resize(size);
	label_.clear();
	label_.push_back(label[0]);
	int count=0;
	//-> signal mean and vari
	signal_.clear();
	signal_var.clear();
	//-> head mean and vari
	head_mean.clear();
	head_vari.clear();
	//-> duration
	durate_.clear();
	double signal_cur=signal[0];
	int signal_num=1;
	for(i=0;i<size-1;i++)
	{
		if(label[i]==label[i+1])
		{
			MoveOrNot[i]=0;
			signal_cur+=signal[i+1];
			signal_num++;
		}
		else
		{
			MoveOrNot[i]=1;
			label_.push_back(label[i+1]);
			//-> calculate mean
			signal_cur/=signal_num;
			signal_.push_back(signal_cur);
			//-> calculate vari
			double vari=0;
			for(int k=0;k<signal_num;k++)
			{
				vari+=(signal[i-k]-signal_cur)*(signal[i-k]-signal_cur);
			}
			vari=1.0*sqrt(vari/signal_num);
			signal_var.push_back(vari);
			//-> calculate head mean
			if(count>0)
			{
				int prev_len=durate_[count-1];
				int head_len=1;
				//-> head mean
				double head_m=signal[i-signal_num];
				for(int k=1;k<=prev_len/4;k++)
				{
					head_m+=signal[i-signal_num-k];
					head_len++;
				}
				for(int k=1;k<=signal_num/4;k++)
				{
					head_m+=signal[i-signal_num+k];
					head_len++;
				}
				head_m/=head_len;
				head_mean.push_back(head_m);
				//-> head vari
				double head_v=(signal[i-signal_num]-head_m)*(signal[i-signal_num]-head_m);
				for(int k=1;k<=prev_len/4;k++)
				{
					head_v+=(signal[i-signal_num-k]-head_m)*(signal[i-signal_num-k]-head_m);
				}
				for(int k=1;k<=signal_num/4;k++)
				{
					head_v+=(signal[i-signal_num+k]-head_m)*(signal[i-signal_num+k]-head_m);
				}
				head_v=1.0*sqrt(head_v/head_len);
				head_vari.push_back(head_v);
			}
			else
			{
				head_mean.push_back(0);
				head_vari.push_back(1);
			}

			//-> calculate length
			durate_.push_back(signal_num);
			//--- next position ---//
			signal_cur=signal[i+1];
			signal_num=1;
			count++;
		}
	}
	MoveOrNot[size-1]=0;
	//-> calculate mean
	signal_cur/=signal_num;
	signal_.push_back(signal_cur);
	//-> calculate vari
	i=size-1;
	double vari=0;
	for(int k=0;k<signal_num;k++)
	{
		vari+=(signal[i-k]-signal_cur)*(signal[i-k]-signal_cur);
	}
	vari=1.0*sqrt(vari/signal_num);
	signal_var.push_back(vari);
	//-> calculate head mean
	if(count>0)
	{
		int prev_len=durate_[count-1];
		int head_len=1;
		//-> head mean
		double head_m=signal[i-signal_num];
		for(int k=1;k<=prev_len/4;k++)
		{
			head_m+=signal[i-signal_num-k];
			head_len++;
		}
		for(int k=1;k<=signal_num/4;k++)
		{
			head_m+=signal[i-signal_num+k];
			head_len++;
		}
		head_m/=head_len;
		head_mean.push_back(head_m);
		//-> head vari
		double head_v=(signal[i-signal_num]-head_m)*(signal[i-signal_num]-head_m);
		for(int k=1;k<=prev_len/4;k++)
		{
			head_v+=(signal[i-signal_num-k]-head_m)*(signal[i-signal_num-k]-head_m);
		}
		for(int k=1;k<=signal_num/4;k++)
		{
			head_v+=(signal[i-signal_num+k]-head_m)*(signal[i-signal_num+k]-head_m);
		}
		head_v=1.0*sqrt(head_v/head_len);
		head_vari.push_back(head_v);
	}
	else
	{
		head_mean.push_back(0);
		head_vari.push_back(1);
	}
	//-> calculate length
	durate_.push_back(signal_num);
}

//============================== From Label to Signal ====================//
//-> prepare genome data (transfer to Z-score)
void Genome_Data_Generate(string &genome_seq, 
	vector <double> &k_mer_zsco, vector <double> &k_mer_orig, int KMER,
	vector <double> &genome_orig, vector <double> &genome_data,int Genome_Scale)
{
	//init
	map<string, int > nam_mapping;
	vector <string> anti_mapping;
	int kmer_size;
	if(KMER==5)
	{
		kmer_size=FiveMer_Mapping_Generate(nam_mapping,anti_mapping);
	}
	else if(KMER==6)
	{
		kmer_size=SixMer_Mapping_Generate(nam_mapping,anti_mapping);
	}
	else
	{
		fprintf(stderr,"only support KMER=5 or 6 \n");
		exit(-1);
	}
	//proc
	int i,k;
	int len=(int)genome_seq.length();
	genome_orig.resize(len*Genome_Scale);
	genome_data.resize(len*Genome_Scale);
	for(i=0;i<len-KMER;i++)
	{
		string curr_str=genome_seq.substr(i,KMER);
		int curr_pos=Get_Mapping(nam_mapping,curr_str);
		for(k=0;k<Genome_Scale;k++)
		{
			genome_orig[i*Genome_Scale+k]=k_mer_zsco[curr_pos];
			genome_data[i*Genome_Scale+k]=k_mer_orig[curr_pos];
		}
	}
	for(i=len-KMER;i<len;i++)
	{
		for(k=0;k<Genome_Scale;k++)genome_orig[i*Genome_Scale+k]=0;
		for(k=0;k<Genome_Scale;k++)genome_data[i*Genome_Scale+k]=100;
	}
	//Z-score
	double mean,vari;
	Z_score(genome_data,mean,vari);
}
void Genome_Data_Generate(vector <int> &genome_label, 	
  vector <double> &k_mer_zsco, vector <double> &k_mer_orig, int KMER,
	vector <double> &genome_orig, vector <double> &genome_data,int Genome_Scale)
{
	//proc
	int i,k;
	int len=(int)genome_label.size();
	genome_orig.resize(len*Genome_Scale);
	genome_data.resize(len*Genome_Scale);
	for(i=0;i<len-KMER;i++)
	{
		int curr_pos=genome_label[i];
		for(k=0;k<Genome_Scale;k++)
		{
			genome_orig[i*Genome_Scale+k]=k_mer_zsco[curr_pos];
			genome_data[i*Genome_Scale+k]=k_mer_orig[curr_pos];
		}
	}
	for(i=len-KMER;i<len;i++)
	{
		for(k=0;k<Genome_Scale;k++)genome_orig[i*Genome_Scale+k]=0;
		for(k=0;k<Genome_Scale;k++)genome_data[i*Genome_Scale+k]=100;
	}
	//Z-score
	double mean,vari;
	Z_score(genome_data,mean,vari);
}

//---------- final deconding -----------//
void Final_Decoding(vector <int> &decode_label,vector <string> &anti_mapping,string &decode_str, int KMER)
{
	//determine KMER
	int kmer_size;
	if(KMER==5)
	{
		kmer_size=1024;
	}
	else if(KMER==6)
	{
		kmer_size=4096;
	}
	else
	{
		fprintf(stderr,"only support KMER=5 or 6 \n");
		exit(-1);
	}
	//proc
	int i;
	int size=(int)decode_label.size();
	decode_str="";
	//-> first position
	{
		i=0;
		int pos=decode_label[i];
		if(pos<0 || pos>=kmer_size)
		{
			fprintf(stderr,"over_range error in FiveMer_To_String !! %d \n",pos);
			exit(-1);
		}
		char code=anti_mapping[pos][0];
		decode_str.push_back(code);
	}
	//-> later position
	for(i=1;i<size;i++)
	{
		if(decode_label[i]!=decode_label[i-1])
		{
			int pos=decode_label[i];
			if(pos<0 || pos>=kmer_size)
			{
				fprintf(stderr,"over_range error in FiveMer_To_String !! %d \n",pos);
				exit(-1);
			}
			char code=anti_mapping[pos][0];
			decode_str.push_back(code);
		}
	}
}


//---------- main ----------//
int main(int argc,char **argv)
{
	//---- Signal_To_Event ----//
	{
		if(argc<3)
		{
			fprintf(stderr,"./Signal_To_Event_kmer <signal_label_input> <KMer_poremodel> \n");
//			fprintf(stderr,"./Signal_To_Event_kmer <signal_label_input> \n");
			exit(-1);
		}
		string signal_label_input=argv[1];
		string kmer_poremodel=argv[2];
		//load
		//-> load original signal_label
		vector <double> signal_;
		vector <int> label_;
		Read_Gold_Standard(signal_label_input,signal_,label_);
		//-> load pore model
/*
		int kmer_count=1024;
		int KMER=5;
		vector <double> kmer_data(kmer_count,0);
		for(int i=0;i<kmer_count;i++)kmer_data[i]=FiveMer_MeanVariLen[i][0];
*/

		vector <double> kmer_data;
		int kmer_count=Load_KMer_Data(kmer_poremodel,kmer_data);
		int KMER;
		if(kmer_count==1024)KMER=5;
		else if(kmer_count==4096)KMER=6;
		else
		{
			fprintf(stderr,"KMER should be 5 or 6 \n");
			exit(-1);
		}

		//-> remove head & tail blank region
		vector <double> signal_r;
		vector <int> label_r;
		Remove_Blank_HeadTail(signal_,label_,signal_r,label_r);

		//-> transform to Zsco
		double mean;
		double vari;
		Z_score(signal_r,mean,vari);
		//-> transform 1024 5mer to Zsco
		vector <double> k_mer_value(kmer_count,0);
		for(int i=0;i<kmer_count;i++)k_mer_value[i]=kmer_data[i];
		Z_score(k_mer_value,mean,vari);

		//-> merge data from L1 to L2
		vector <int> MoveOrNot;
		vector <int> label;
		vector <double> signal;
		vector <double> signal_var;
		vector <double> head_mean;
		vector <double> head_vari;
		vector <int> durate;		
		Signal_To_Event(signal_r,label_r,MoveOrNot,label,signal,signal_var,
			head_mean,head_vari,durate);

		//proc
		vector <double> genome_orig;
		vector <double> genome_data;
		Genome_Data_Generate(label, k_mer_value, kmer_data, KMER, genome_orig, genome_data,1);
		//name
		string name;
		getBaseName(signal_label_input,name,'/','.');
		//output
		{
			//-> get fivemer mapping
			map<string, int > nam_mapping;
			vector <string> anti_mapping;
			int kmer_size;
			if(KMER==5)
			{
				kmer_size=FiveMer_Mapping_Generate(nam_mapping,anti_mapping);
			}
			else if(KMER==6)
			{
				kmer_size=SixMer_Mapping_Generate(nam_mapping,anti_mapping);
			}
			else
			{
				fprintf(stderr,"only support KMER=5 or 6 \n");
				exit(-1);
			}
			//-> start output
			string native_str;
			Final_Decoding(label, anti_mapping, native_str, KMER);
			//-> output
			printf(">%s\n",name.c_str());
			printf("%s\n",native_str.c_str());
			for(int i=0;i<(int)label.size()-KMER;i++)
			{
				string str=native_str.substr(i,KMER);
//				printf("%s %5d | %10.7f %10.7f | %10.7f %5d\n",
//					str.c_str(),label[i],genome_orig[i],genome_data[i],signal[i], durate[i]);
				printf("%s %5d | %10.7f %10.7f | %10.7f %10.7f | %10.7f %10.7f | %5d \n",
					str.c_str(),label[i],genome_orig[i], genome_data[i],
					head_mean[i], head_vari[i], signal[i], signal_var[i], durate[i]);
			}
		}
		//exit
		exit(0);
	}
}



