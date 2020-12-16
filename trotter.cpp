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
	
	uint nthreads = 0; //If 0, the number of threads will be set automatically
	double a_time = 0.0; //Timing variable
	//XXZ model parameters
	uint   L    = 8;    //Spin chain length
	double h    = 0.0;  //Strength of random magnetic field
	uint   nh   = 0;    //Number of random magnetic field replicas, if 0 = num of OpenMP threads
	//Trotter decomposition parameters
	double dt   = 0.01; //Trotter discretization step
	double tmax = 1.0;  //Total evolution time
	
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
		 ("dt",         po::value<double>(   &(dt))->default_value( 0.01), "Trotter decomposition step");
	     
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	// Check if help is needed
	if(vm.count( "help" )){cout<<desc<<endl; return 1;};
	//if(vm.count("min-h")) quantum = true; 
	
	uint   nt   = (int)ceil(tmax/dt);  //Number of steps
	
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
	printf("Here we calculate the real-time evolution operators with Trotter steps of %2.2E and %2.2E up to time %2.2E (%u steps)\n", dt, 0.5*dt, (dt*nt), nt);
	printf("And compare the Frobenius norm of the difference\n\n");
	
	cout << ansi::green << "L  =   " << ansi::magenta << L    << ansi::reset << endl;
	cout << ansi::green << "h  =   " << ansi::magenta << h    << ansi::reset << endl;
	cout << ansi::green << "nh =   " << ansi::magenta << nh   << ansi::reset << endl;
	cout << ansi::green << "dt =   " << ansi::magenta << dt   << ansi::reset << endl;
	cout << ansi::green << "tmax = " << ansi::magenta << tmax << ansi::reset << endl;
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
	
	double* trotter_errs = new double[nh*nt];
	double* norm_errs    = new double[nh*nt];
	
	TIMING_START;
	#pragma omp parallel for
	for(uint ih=0; ih<nh; ih++)
	{
		int ithread = omp_get_thread_num();
		int rng_seed = (int)omp_get_wtime() + ithread;
		SpinChain* SC = new SpinChain(L, h, false, false, false, rng_seed);
		
		SC->init_magnetic_hamiltonian();	 //This is a diagonal matrix in our spin eigenbasis, thus we only store diagonal elements
		double t = 0.0; 
		
		/*t_complex* psi_in = double2complex(&(psi0[502*NS0]), NS0); //SC->rand_vec(NS0, true);
		t_complex* in1  = new t_complex[NS0];
		t_complex* out1 = new t_complex[NS0];
		t_complex* in2  = new t_complex[NS0];
		t_complex* out2 = new t_complex[NS0];
		
		std::copy(psi_in, psi_in + NS0, in1);
		std::copy(psi_in, psi_in + NS0, in2);*/
		
		t_complex* U1   = identity_matrix(NS0);
		t_complex* U2   = identity_matrix(NS0);
		t_complex* Udt1 = SC->trotter_evolution_operator(E0, psi0, 1.0*dt);
		t_complex* Udt2 = SC->trotter_evolution_operator(E0, psi0, 0.5*dt);
		t_complex* tmp  = new t_complex[NS02];
				
		for(uint it=0; it<nt; it++)
		{
			//Trotter evolution with step dt for in1
			A_eq_B_mult_C(tmp, U1, Udt1, NS0);
			std::copy(tmp, tmp + NS02, U1);
			
			/*SC->trotter_evolve(in1, out1, E0, psi0, dt);
			std::copy(out1, out1 + NS0, in1);*/
			
			//Trotter evolution with step dt/2 for in2
			A_eq_B_mult_C(tmp, U2, Udt2, NS0);
			std::copy(tmp, tmp + NS02, U2);
			A_eq_B_mult_C(tmp, U2, Udt2, NS0);
			std::copy(tmp, tmp + NS02, U2);
			
			/*SC->trotter_evolve(in2, out2, E0, psi0, 0.5*dt);
			std::copy(out2, out2 + NS0, in2);
			SC->trotter_evolve(in2, out2, E0, psi0, 0.5*dt);
			std::copy(out2, out2 + NS0, in2);*/
			
			t += dt;
			
			/*trotter_errs[nt*ih + it] = norm_diff(in1, in2, NS0);
			   norm_errs[nt*ih + it] = std::max(abs(norm(in1, NS0) - 1.0), abs(norm(in2, NS0) - 1.0));*/
			trotter_errs[nt*ih + it] = norm_diff(U1, U2, NS02)/(double)NS0;
			   norm_errs[nt*ih + it] = std::max(unitarity_err(U1, NS0), unitarity_err(U2, NS0));  
			   
			if(ithread==0){cout << "."; fflush(stdout);};
		};
		
		/*delete [] psi_in; 
		delete [] in1;	delete [] out1;	
		delete [] in2;	delete [] out2;	*/
		delete [] U1; delete [] U2;
		delete [] Udt1;	delete [] Udt2;
		delete [] tmp;
		
		delete SC;
	};//End of loop over magnetic field replicas
	cout << endl;
	TIMING_FINISH;
	std::cout << ansi::green << "... Time per Trotter step:\t" << ansi::magenta << (a_time/(double)nt) << ansi::green << " sec.\n" << ansi::reset << endl << endl;
	
	//Statistical averaging and output
	char fname[512];
	sprintf(fname, "%s/trotter/trotter_errs_L%i_dt%2.2E_h%2.4lf_frobenius.dat", datadir.c_str(), L, dt, h);
	FILE* f = fopen(fname, "w");
	if(f==NULL) fprintf(stderr, "Cannot open the file %s for writing!\n", fname);
	
	for(uint it=0; it<nt; it++)
	{
		double t = dt*(double)it;
		double a_trotter_err = 0.0;	double d_trotter_err = 0.0;
		double a_norm_err = 0.0;	double d_norm_err = 0.0;
		for(uint ih=0; ih<nh; ih++)
		{
			a_trotter_err += trotter_errs[nt*ih + it];
			d_trotter_err += trotter_errs[nt*ih + it]*trotter_errs[nt*ih + it];
			a_norm_err    +=    norm_errs[nt*ih + it];
			d_norm_err    +=    norm_errs[nt*ih + it]*norm_errs[nt*ih + it];
		};
		a_trotter_err /= (double)nh; d_trotter_err /= (double)nh;
		   a_norm_err /= (double)nh;    d_norm_err /= (double)nh;
		d_trotter_err = (nh>1? sqrt(abs(d_trotter_err - a_trotter_err*a_trotter_err)/(double)(nh-1)) : 0.0);
	   	d_norm_err =    (nh>1? sqrt(abs(d_norm_err -    a_norm_err*a_norm_err   )/(double)(nh-1))    : 0.0);
		
		if(f!=NULL) fprintf(f, "%4.4lf\t%2.4E\t%2.4E\t%2.4E\t%2.4E\n", t, a_trotter_err, d_trotter_err, a_norm_err, d_norm_err);   
	};
	if(f!=NULL)
	{
		fclose(f);
		fprintf(stdout, "Data saved to the file %s", fname);
	};
	
	cout << endl;
	return EXIT_SUCCESS;
}

