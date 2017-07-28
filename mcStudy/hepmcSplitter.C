
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

int hepmcSplitter(int num=0)
{
	string filename = "rapgapBatches/hepmc_HI_jetSel_b15.dat"; //input
	ifstream infile(filename.c_str());
	char outfilename[100];
	sprintf(outfilename,"hepmc/PbPb2015_jet_hepmc_b15/rapgap_%d.dat",num+1);
	ofstream output(outfilename);
	if (! infile.is_open()) { cout << "\t ERROR: I can not open \"" << filename << "\"" << endl; return 0; }
	string temp_string, temp;
	istringstream curstring;

	Int_t M = 10000; // N_events	
	Int_t K = num*M+1;  // first event 
	Int_t evt_n;
	Int_t eventsInNewFile = 0;

//	getline(infile,temp_string); // The very first line is useless
	output <<  "HepMC::IO_Ascii-START_EVENT_LISTING" << endl;

	TDatabasePDG *massFinder = new TDatabasePDG();

	while (getline(infile,temp_string)) {

		curstring.clear(); 
		curstring.str(temp_string);
		if(strstr(temp_string.c_str(),"E")){
			curstring >> temp >> evt_n;
			if(evt_n >=K && evt_n < K+M){
				eventsInNewFile++;
				output << temp_string << endl;
			}

		} else if(strstr(temp_string.c_str(),"V"))	{
			if(evt_n >=K && evt_n < K+M){
				output << temp_string << endl;
			}
		} else if(strstr(temp_string.c_str(),"P"))	{
			if(evt_n >=K && evt_n < K+M){
				output << temp_string << endl;
			}
		}
	}
	output << "HepMC::IO_Ascii-END_EVENT_LISTING" << endl;
	infile.close();
	output.close();
	cout << eventsInNewFile << " events written in " << outfilename << endl;
	return num;
}

void makefiles(int number_of_files=100) {
	int particle_number = hepmcSplitter(0);
	for(int i=0; i<number_of_files; i++){
		cout<<i<<endl;
		particle_number = hepmcSplitter(i);
	}
}

/* Explaination of the format :
 * +++ Event +++
 * E 1 -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 20 0 1 0 0
 *   1 : event number  			<-------
 *   -1 : event scale
 *   -1 : alpha_QCD
 *   -1 : alpha_QED
 *   20 : signal process ID
 *   0 : signal process vertex barcode
 *   1 : number of vertices 		<-------
 *   0 : list of vertices
 *   0 : ?
 *
 * +++ Vertex +++
 * V -1 0 0 0 0 0 0 4 0
 *   -1 : vertex barcode (unique)       <-------
 *    0 : vertex id
 *    0 0 0 0 : vertex x,y,z,t
 *    0 : number of orphans
 *    4 : number of out particles       <-------
 *    0 : weights
 *
 * +++ Particle +++
 * P 5 2212 -2.0 1.2 8.1 5.2 1 0 0 0 0   
 *    5 : barcode			<-------
 *    0 : pdg_id			<-------
 *   -2.0 : px				<-------
 *    1.2 : py				<-------
 *    8.1 : pz				<-------
 *    5.2 : e				<-------
 *    1 : status			<-------
 *    0 0  : polarization eta , phi
 *    0 0  : vertex and ?
 */
