#ifndef _SPIN_CHAIN_HPP_
#define _SPIN_CHAIN_HPP_

#include<complex>
#include<chrono>
#include<random>
#include<limits>

#include<OpenBLAS/cblas.h>
#include<OpenBLAS/lapacke.h>

using namespace std; 

typedef complex<double> t_complex;
typedef unsigned int    uint;

#define SQR(_x) ((_x)*(_x))

uint count_bits(uint a, uint L); //Counts the number of binary 1's in the binary notation of a, up to L 

class SpinChain
{
	private:
		std::ranlux48 rng_engine;
		std::uniform_real_distribution<double> rng_uniform_dist;
		std::normal_distribution<double> rng_normal_dist{0.0, 1.0};
	public:
		SpinChain(uint L, double h, bool diagonalize=true, bool check=true, bool noisy=false);
		~SpinChain(){};
		uint L     = 0;										//Spin chain length
		uint N     = 1;										//Size of the Hilbert space
		uint N2	   = 1;										//N*N
		uint NS0   = 1;										//Number of spin zero states
		uint NS02  = 1;										//NS0*NS0
		uint*   S0_basis  = NULL;							//uints which encode basis vectors with total spin 0
		uint*   S0_lookup = NULL;							//S0_lookup[ips] produces the index of the corresponding state vector in S0_basis if it's total spin is zero
		double* hi = NULL;									//A particular instance of a random magnetic field
		/* Hamiltonian definition */
		void    H(double* in, double* out);	 				//Hamiltonian, in and out are real since the Hamiltonian is real
		void    HS0(double* in, double* out);				//Hamiltonian in spin0 sector
		void    sz(double* in, double* out);       			//global SigmaZ operator
		void    sz(double* in, double* out, uint i);       	//local  SigmaZ operator
		void    diagonalize_HS0(bool check=true, bool noisy=false);
		double* E  = NULL;									//Eigenenergies
		double* U  = NULL;									//The real-value matrix of eigenvectors, H = O.E.O^T, E is treated as diagonal matrix
		double  max_orthogonality_err	=	0.0;
		double  max_eigensystem_err		=	0.0;
		//Spectral properties
		void    get_HS0_matrix(double* HM);
		//Some useful stuff
		void    rand_vec(double* in, uint n, bool normalize=true);	//Random real vector, filled with Gaussian numbers with unit dispersion
		double  scalar_prod(double* psi1, double* psi2, uint n);
		void    aA_plus_bB(double* out, double a, double* A, double b, double* B, uint n);
		double  norm_diff(double* psi1, double* psi2, uint n);
};

#endif
