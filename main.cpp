#include <iostream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>
#include <thread>

#include <OpenBLAS/cblas.h>
#include <OpenBLAS/lapacke.h>

#include "spin_chain.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
	//XXZ model parameters
	uint   L  = 8;    //Spin chain length
	double h  = 0.0;  //Strength of random magnetic field
	uint   nh = 1;    //Number of random magnetic field replicas
	
	//Option parsing
	po::options_description desc("Simulation of a spin chain which features a chaos-order transition");
	desc.add_options()
	     ("help,h", "produce this help message")
	     ("L",   po::value<uint>(          &(L))->default_value(   8), "Spin chain length")
	     ("h",   po::value<double>(        &(h))->default_value( 0.0), "Random magnetic field length")
		 ("nh",  po::value<uint>(         &(nh))->default_value(   1), "Number of random magnetic field replicas");
	     
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	// Check if help is needed
	if(vm.count( "help" )){cout<<desc<<endl; return 1;};
	//if(vm.count("min-h")) quantum = true; 
	
	char suffix[512];
	sprintf(suffix, "L%u_h%2.4lf", L, h);
	
	for(uint ih=0; ih<nh; ih++)
	{
		printf("\tIteration %u out of %u\n", ih, nh);
		SpinChain* SC = new SpinChain(L, h);
		
		FILE* f = fopen("G:\\LAT\\ExDiag\\hi.dat", "w");
		for(uint i=0; i<SC->L; i++)
			fprintf(f, "%+2.8E\n", SC->hi[i]);
		fclose(f);
		
		//Saving the whole matrix to output
		double* HM = new double[SC->N2];
		SC->get_H_matrix(HM);
		
		f = fopen("G:\\LAT\\ExDiag\\H.dat", "w");
		for(uint i=0; i<SC->N; i++)
		{
			for(uint j=0; j<SC->N; j++)
				fprintf(f, "%+2.8E ", HM[i*SC->N + j]);
			fprintf(f, "\n");
		};
		fclose(f);
		
		delete [] HM;
		
		printf("Eigenspectrum:\n");
		f = fopen("G:\\LAT\\ExDiag\\E.dat", "w");
		for(uint i=0; i<SC->NS0; i++)
		{
			fprintf(f, "%+2.8E\n", SC->E[i]);
			printf("%+2.4lf\n", SC->E[i]);
		};
		fclose(f);
		fflush(stdout);
		
		delete SC;
	};
	
	return EXIT_SUCCESS;
}

