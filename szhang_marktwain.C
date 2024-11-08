//Author: Shuaixiang Zhang; Mar 21, 2022---
//iminate the method used in yw_marktwain.C by Prof. Liu---

#include <algorithm> //next_permutation() is in this headfile--
#include <TH1.h>
using namespace std;

int szhang_marktwain(){
	const int ndata = 18;
	const int n1 = 8;
	const int n2 = 10;	
    const double data[18] = {
    0.225, 0.262, 0.217, 0.240, 0.230, 0.229, 0.235, 0.217,
    0.209,0.205,0.196,0.210, 0.202,0.207,0.224, 0.223, 0.220, 0.201
    };
    TH1F *hDiff = new TH1F("hDiff", "diffs", 200, -0.1, 0.1);
    int label[18];
    
    for(int i=0; i<8; ++i){//initialization of the label---
    	label[i] = 1;
	}
	for(int i=8; i<18; ++i){
		label[i] = 0;
	}
	for(int i=0; i<18; ++i){//to test the label[]---
		cout<<label[i]<<"  ";
	}
	cout<<"\n"<<endl;
	
	//permutation: C^8_{18}---
	double sum1 = 0.0;
	double sum2 = 0.0;
	double diff = 0.0;//to show (sum1/8.-sum2/10.)---
	double diff0 = 0.0;//to record the initial value---
	int irun = 0;//to record the number of permutations---
	int np = 0;//to record the case larger than initial array---
	
	do{
		for(int i=0; i<ndata; ++i){
			if(label[i]==1)
			    sum1 += data[i];
			else
			    sum2 += data[i];
		}
		diff = sum1/float(n1) - sum2/float(n2);
		
		if(irun==0){//Show initial value without permutation---
		    diff0 = diff;
			cout<<"------------------------------"<<endl;
			cout<<"The initial value is: "<<diff0<<endl;
			cout<<"sum1/8.0 = "<<sum1/float(n1);
			cout<<"; sum2/10.0 = "<<sum2/float(n2)<<endl;
			cout<<"------------------------------"<<endl;
		}
		
		if(fabs(diff) > fabs(diff0)){
			++np;
			cout<<"diff: "<<diff<<endl;
			for(int i=0; i<ndata; ++i){
				if(label[i]==1)//to show the position of '1'---
				    cout<<i<<"    ";
			}
			cout<<"\n"<<endl;
			for(int i=0; i<ndata; ++i){//new selected first eight values---
				if(label[i]==1)
				    cout<<data[i]<<" ";
			}
			cout<<"\n"<<endl;
		}
		
		hDiff->Fill(diff);//fill the histo---
		
		irun++;//to record the num of permutations---
		sum1 = 0;
		sum2 = 0;
		diff = 0.0;
	}while(std::prev_permutation(label, label + ndata));
	
	//Final result---
	cout<<"============Final Result============"<<endl;
	cout<<"np/ntrials: "<<double(np)<<"/"<<double(irun);
	cout<<" = "<<double(np)/double(irun)<<endl;
	
	return 0;
}
