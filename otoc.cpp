#include <iostream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>

#include "spin_chain.hpp"
#include "timing.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
	cout << endl;
	
	cout << sizeof(float) << endl;
	cout << sizeof(double) << endl;
	cout << sizeof(long double) << endl;
	
	double x = 22.0/7.0;
	double y = sqrt(x);
	double r = 7.0*y*y/22.0 - 1.0;
    printf("%2.4E", r); 
	
	return EXIT_SUCCESS;
	
	uint nthreads = 0; //If 0, the number of threads will be set automatically
	double a_time = 0.0; //Timing variable
	//XXZ model parameters
	uint   L    = 8;    //Spin chain length
	double h    = 0.0;  //Strength of random magnetic field
	uint   nh   = 0;    //Number of random magnetic field replicas, if 0 = num of OpenMP threads
	//Trotter decomposition parameters
	double dt   = 0.01; //Trotter discretization step
	double tmax = 1.0;  //Total evolution time
	double beta = 1.0;  //Inverse temperature
	
	string datadir = "./data/";
	
	//Option parsing
	po::options_description desc("Simulation of a spin chain which features a chaos-order transition");
	desc.add_options()
	     ("help,h", "produce this help message")
	     ("nthreads",	po::value<uint>( &(nthreads))->default_value(  0), "Number of OpenMP threads to use, 0 = automatic choice")
	     ("L",   		po::value<uint>(        &(L))->default_value(  8), "Spin chain length")
	     ("datadir", 	po::value<string>(                    &datadir  ), "Directory for data output" )
	     ("h",          po::value<double>(    &(h))->default_value(  0.0), "Random magnetic field length")
		 ("nh",         po::value<uint>(     &(nh))->default_value(    0), "Number of random magnetic field replicas")
		 ("tmax",       po::value<double>( &(tmax))->default_value(  1.0), "Number of Trotter evolution steps")
		 ("beta",       po::value<double>( &(beta))->default_value(  1.0), "Inverse temperature")
		 ("dt",         po::value<double>(   &(dt))->default_value( 0.01), "Trotter decomposition step");
	     
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	// Check if help is needed
	if(vm.count( "help" )){cout<<desc<<endl; return 1;};
	//if(vm.count("min-h")) quantum = true; 
	
	uint   nt   = (int)ceil( tmax/dt);  //Number of steps
	
	if(nthreads==0)
	{
		nthreads = omp_get_max_threads();
	}
	else
	{
		omp_set_num_threads(nthreads);
		openblas_set_num_threads(nthreads);
	};
	cout << ansi::white << "Using " << ansi::magenta << nthreads << ansi::white << " OpenMP threads, openblas_num_threads = " << ansi::magenta << openblas_get_num_threads() << ansi::reset << endl;
	
	nh = (nh==0? nthreads : nh);
	
	printf("\nHere we perform exact diagonalization of Eq. (30) in 1803.08050 in sector with total Sz = 0\n");
	printf("The MBL phase transition happens at h=3.75, as discussed in 1411.0660\n");
	printf("Here we calculate the real-time evolution operators with Trotter steps of %2.2E up to time %2.2E (%u steps)\n", dt, (dt*nt), nt);
	printf("And calculate the out-of-time order correlators for the operator sz_0\n\n");
	
	cout << ansi::green << "L  =   " << ansi::magenta << L    << ansi::reset << endl;
	cout << ansi::green << "h  =   " << ansi::magenta << h    << ansi::reset << endl;
	cout << ansi::green << "nh =   " << ansi::magenta << nh   << ansi::reset << endl;
	cout << ansi::green << "dt =   " << ansi::magenta << dt   << ansi::reset << endl;
	cout << ansi::green << "tmax = " << ansi::magenta << tmax << ansi::reset << endl;
	cout << ansi::green << "beta = " << ansi::magenta << beta << ansi::reset << endl;
	cout << endl;
	
	//Diagonalizing the spin chain without any magnetic field ...
	SpinChain* SC0 = new SpinChain(L, 0.0, false, false);
	uint NS0 = SC0->NS0; uint NS02 = NS0*NS0;
	cout << ansi::cyan << "Space for storing the eigensystem in spin-0 sector: " << ansi::red << SC0->storage_size << " Gb" << ansi::reset << endl << endl;
		
	std::cout << "Diagonalizing h=0 Hamiltonian ..." << endl << endl;
	double* E0     = new double[NS0];
	double* psi0   = new double[NS02];
	double* HS0_h0 = new double[NS02];
	SC0->get_HS0_matrix(HS0_h0, false);
	
	TIMING_START;
	SC0->diagonalize_HS0_h0(E0, psi0);
	TIMING_FINISH;
	std::cout << "... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl << endl;
	
	std::cout << "Sorting and checking the eigensystem..." << endl;
	TIMING_START;
	 sort_eigensystem(        E0, psi0, NS0);
	check_eigensystem(HS0_h0, E0, psi0, NS0);
	TIMING_END;
	std::cout << "... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl << endl;
	delete [] HS0_h0;
	
	std::cout << "Calculating Trotter evolution for " << nh << " replicas using " << nthreads << " threads, nt = " << nt << "..." << endl;
	
	double* otocs = new double[nh*nt];
	
	TIMING_START;
	#pragma omp parallel for
	for(uint ih=0; ih<nh; ih++)
	{
		int ithread = omp_get_thread_num();
		int rng_seed = (int)omp_get_wtime() + ithread;
		SpinChain* SC = new SpinChain(L, h, false, false, false, rng_seed);
		SC->init_magnetic_hamiltonian();	 //This is a diagonal matrix in our spin eigenbasis, thus we only store diagonal elements
		
		//First, we initialize the density matrix
		
		t_complex* rho  = SC->trotter_density_matrix(    E0, psi0, beta, dt);
		t_complex* Udt  = SC->trotter_evolution_operator(E0, psi0,       dt);
		t_complex* Udtc = hermitian_conjugate(Udt, NS0);
		
		//Now the OTOC operators A and B
		t_complex* A = SC->s3(0); t_complex* B = SC->s3(0);
		
		double t = 0.0;
		t_complex* tmp  = new t_complex[NS02]; t_complex* tmp1 = new t_complex[NS02];
		//Now the real-time evolution
		for(uint it=0; it<nt; it++)
		{
			//Real-time evolution of B
			A_eq_B_mult_C(tmp,    B,  Udt, NS0);
			A_eq_B_mult_C(  B, Udtc,  tmp, NS0);

			//Forming the OTOC
			commutator(tmp, A, B, NS0);
			A_eq_B_mult_C(tmp1, tmp, tmp, NS0);
			A_eq_B_mult_C(tmp, rho, tmp1, NS0);
			otocs[ih*nt + it] = -1.0*real(tr(tmp, NS0));
			//
			t += dt;			   
			if(ithread==0){cout << "."; fflush(stdout);};
		};
		
		delete [] rho;
		delete [] tmp;
		delete [] Udt;
		delete [] Udtc;
		delete [] tmp1;
		
		delete SC;
	};//End of loop over magnetic field replicas
	cout << endl;
	TIMING_FINISH;
	std::cout << ansi::green << "... Time per Trotter step:\t" << ansi::magenta << (a_time/(double)nt) << ansi::green << " sec.\n" << ansi::reset << endl << endl;
	
	char fname[512];
	sprintf(fname, "%s/OTOC/s3s3_L%i_beta%2.2E_dt%2.2E_h%2.4lf.dat", datadir.c_str(), L, beta, dt, h);
	FILE* f = fopen(fname, "w");
	if(f==NULL) fprintf(stderr, "Cannot open the file %s for writing!\n", fname);
	
	for(uint it=0; it<nt; it++)
	{
		double aotoc = 0.0, dotoc = 0.0;
		for(uint ih=0; ih<nh; ih++)
		{
			aotoc += otocs[ih*nt + it];
			dotoc += otocs[ih*nt + it]*otocs[ih*nt + it];
		};
		aotoc = aotoc/(double)nh; 
		if(nh>1) dotoc = sqrt((dotoc/(double)nh - aotoc*aotoc)/(double)(nh-1));
		
		cout << ansi::cyan << ((it+1)*dt) << "\t" << ansi::green;
		printf("%2.4E\t%2.4E", aotoc, dotoc);
		cout << ansi::reset << endl;
		
		if(f!=NULL) fprintf(f, "%4.4lf\t%2.4E\t%2.4E\n", ((it+1)*dt), aotoc, dotoc);   
	};
	
	if(f!=NULL)
	{
		fclose(f);
		fprintf(stdout, "Data saved to the file %s", fname);
	};
	
	delete [] otocs;
	
	cout << endl;
	return EXIT_SUCCESS;
}

