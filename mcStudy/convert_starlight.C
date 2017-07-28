
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

int makeEventsFile(int num=0,int particle_number = 0)
{
	string filename = "slight.out"; //input
	ifstream infile(filename.c_str());
	char outfilename[100];
	sprintf(outfilename,"hepmc/dpmjet/starlight%d.dat",num+1);
	ofstream output(outfilename);
	if (! infile.is_open()) { cout << "\t ERROR: I can not open \"" << filename << "\"" << endl; return 0; }

	string temp_string, temp;
	istringstream curstring;
	Int_t M = 100; // N_events	
	Int_t K = num*M+1;  // first event 
	const double MU = 0.105658369; // muon mass [GeV]
	Int_t evt_n=0; // event_number, read from the input-file
	Int_t nn=0; // event_counter, in the output file

	getline(infile,temp_string); // The very first line is useless
	output <<  "HepMC::IO_Ascii-START_EVENT_LISTING" << endl;

	Int_t N = 0;
	Int_t N_piece = 0;

	TDatabasePDG *massFinder = new TDatabasePDG();


	while (getline(infile,temp_string)) {

		curstring.clear(); 
		curstring.str(temp_string);
		if(strstr(temp_string.c_str(),"EVENT"))	{
			curstring >> temp >> evt_n;

			if(evt_n >=K && evt_n < K+M)  {
				output << "E " << evt_n << " -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 20 " << evt_n << " 1 0 0" << endl;
				nn++;
			}

		} else if(strstr(temp_string.c_str(),"VERTEX"))	{
			float x,y,z,t,a,b,c;
			curstring >> temp >> x >> y >> z >> t >> a >> b >> c >> N;
			if(evt_n >=K && evt_n < K+M) output << "V " << evt_n << " 0 " << x << " " << y << " " << z << " " << t << " 0 " << N <<  " 0" << endl;
		} else if(strstr(temp_string.c_str(),"TRACK"))	{
			int useless, part_n, pdg_id;
			float px, py, pz;

			curstring >> temp >> useless >> px >> py >> pz >> part_n >> useless >> useless >> pdg_id;
			N_piece++;
			particle_number++;

			TParticlePDG *particleProspect= massFinder->GetParticle(pdg_id);
			double mass = particleProspect->Mass();
			if(evt_n >=K && evt_n < K+M) output << "P " << N_piece << " " << pdg_id << " " << px << " " << py << " " << pz << " " << sqrt(mass*mass + px*px + py*py + pz*pz) << " 1 0 0 0 0\n";
		}
	}
	output << "HepMC::IO_Ascii-END_EVENT_LISTING" << endl;
	infile.close();
	output.close();
	cout << nn << " events written in " << outfilename << endl;
	return particle_number++;
}
void makefiles(int number_of_files=100) {
	int particle_number = makeEventsFile(0,0);
	for(int i=0; i<number_of_files; i++)
		particle_number = makeEventsFile(i,particle_number);
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
