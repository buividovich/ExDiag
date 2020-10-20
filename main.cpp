#include <iostream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>
#include <thread>

#include "spin_chain.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
	//XXZ model parameters
	uint   L  = 8;    //Spin chain length
	double h  = 0.0;  //Strength of random magnetic field
	uint   nh = 1;    //Number of random magnetic field replicas
	
	string datadir = "G:/LAT/ExDiag/data/";
	
	//Option parsing
	po::options_description desc("Simulation of a spin chain which features a chaos-order transition");
	desc.add_options()
	     ("help,h", "produce this help message")
	     ("L",   		po::value<uint>(          &(L))->default_value( 8), "Spin chain length")
	     ("datadir", 	po::value<string>( &datadir ),     "Directory for data output" )
	     ("h",   po::value<double>(        &(h))->default_value( 0.0), "Random magnetic field length")
		 ("nh",  po::value<uint>(         &(nh))->default_value(   1), "Number of random magnetic field replicas");
	     
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	// Check if help is needed
	if(vm.count( "help" )){cout<<desc<<endl; return 1;};
	//if(vm.count("min-h")) quantum = true; 
	
	printf("\nHere we perform exact diagonalization of Eq. (30) in 1803.08050 in sector with total Sz = 0\n");
	printf("The MBL phase transition happens at h=3.75, as discussed in 1411.0660\n");
	printf("At this transition min(dE_i,dE_{i+1})/max(dE_i,dE_{i+1}) changes from 0.53 to 0.39 (at the centrum of the spectrum at E=0.5)\n\n");

	cout << "Datadir:\t" << datadir << endl << endl;
	
	SpinChain* SC0 = new SpinChain(L, h, false, false);
	int NS0 = SC0->NS0;
	double* ES = new double [NS0*nh];
	delete SC0;
	
	printf("Calculating spectra for %u replicas using %i threads...\n", nh, omp_get_max_threads());
	
	#pragma omp parallel for
	for(uint ih=0; ih<nh; ih++)
	{
		int ithread = omp_get_thread_num();
		SpinChain* SC = new SpinChain(L, h);
		
		bool ordered_check = true;
		for(uint ie=0; ie<NS0; ie++)
		{
			if(ie<NS0-1) ordered_check = ordered_check && (SC->E[ie+1] > SC->E[ie]);
			ES[ih*NS0 + ie] = SC->E[ie];
		};
		if(!ordered_check){ fprintf(stderr, "\x1b[1;31m Eigenvalues NOT in order \x1b[0m\n"); fflush(stderr);};
		
		delete SC;
		fprintf(stdout, "."); fflush(stdout);
	};
	printf("\n");
	
	char Ename[512];
	sprintf(Ename, "%s/E_L%u_h%2.4lf.dat", datadir.c_str(), L, h);
	printf("\nOutput file for energies:\t%s\n", Ename);
	FILE* f = fopen(Ename, "a");
	if(f==NULL){cerr << "The file " << Ename << " could not be opened" << endl; };
	
	printf("Calculating the characteristic ratio...\n", nh);	
	//Calculating the characteristic ratio for all energy levels between -0.5 and 0.5
	double aR = 0.0, dR = 0.0; uint ncR = 0;	
	
	for(uint ih=0; ih<nh; ih++)
	{
		double R = 0.0; uint nR = 0;
		double* E = &(ES[ih*NS0]);
		for(uint ie=0; ie<NS0; ie++)
		{
			if(ie>0 && ie<NS0-1 && E[ie]>-0.5 && ES[ie]<0.5)
			{
				double dE1 = E[ie  ] - E[ie-1];
				double dE2 = E[ie+1] - E[ie  ];
				R += min(dE1, dE2)/max(dE1, dE2);
				nR ++;
			};
			if(f!=NULL) fprintf(f, "%+2.8E\n", E[ie]);
		};
			
		if(nR > 0){	R = R/(double)nR; aR += R; dR += R*R; ncR++;};
	};
	if(f!=NULL) fclose(f);
	
	aR /= (double)ncR;
	dR /= (double)ncR;
	dR = sqrt((dR - aR*aR)/(double)(ncR-1));
	
	printf("\n\t >> \x1b[1;36m Characteristic ratio:\t \x1b[1;35m %2.2lf +/- %2.2lf \x1b[1;37m (over %u points out of %u)\x1b[0m\n", aR, dR,ncR, nh);
	
	delete [] ES;
	
	return EXIT_SUCCESS;
}

