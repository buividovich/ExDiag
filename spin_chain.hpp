#ifndef _SPIN_CHAIN_HPP_
#define _SPIN_CHAIN_HPP_

#include<random>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<algorithm>
#include<vector>
#include<climits>

#include "linalg.hpp"
#include "ansi_io.hpp"
#include "timing.hpp"

uint  count_bits(uint a, uint L); //Counts the number of binary 1's in the binary notation of a, up to L 
uint rotate_bits(uint a, uint L, uint s=1); //Cyclic shift of L bits by s bits

class SpinChain
{
	private:
		std::ranlux48 rng_engine;
		std::uniform_real_distribution<double> rng_uniform_dist;
		std::normal_distribution<double> rng_normal_dist{0.0, 1.0};
	public:
		SpinChain(uint L, double h, bool diagonalize=true, bool check=true, bool noisy=false, int rng_engine_seed=0, string eigensystem_file_name="");
		~SpinChain(){};
		uint        L     = 0;												//Spin chain length
		uint        N     = 1;												//Size of the Hilbert space
		uint        N2	   = 1;												//N*N
		uint        NS0   = 1;												//Number of spin zero states
		uint        NS02  = 1;										//NS0*NS0
		uint*       S0_basis  = NULL;							//uints which encode basis vectors with total spin 0
		uint*       S0_lookup = NULL;							//S0_lookup[ips] produces the index of the corresponding state vector in S0_basis if it's total spin is zero
		double*     hi = NULL;									//A particular instance of a random magnetic field
		/* Hamiltonian definition */
		void        H(double* in, double* out);	 				//Hamiltonian, in and out are real since the Hamiltonian is real
		//void        HS0(double*    in, double*    out);				//Hamiltonian in spin0 sector
		//void        HS0(t_complex* in, t_complex* out);				//Hamiltonian in spin0 sector
		template<typename T> 
		void        HS0(T* in, T* out, bool magnetic_field_on=true);
		void        sz(double* in, double* out);       			//global SigmaZ operator
		void        sz(double* in, double* out, uint i);       	//local  SigmaZ operator
		void        diagonalize_HS0_h0(double* E0, double* psi0, bool check=true, bool noisy=false);
		void        diagonalize_HS0(bool check=true, bool noisy=false);
		t_complex*  evolution_operator(t_complex dt);
		void        init_magnetic_hamiltonian();  //returns HI is the term in the Hamiltonian with external magnetic field. The data is an array of NS0 diagonal matrix elements in the basis of S=0 states
		void        trotter_evolve(t_complex* in, t_complex* out, double* E0, double* psi0, t_complex dt);
		t_complex*  trotter_evolution_operator(double* E0, double* psi0, t_complex dt);
		void        exact_evolve(t_complex* in, t_complex* out, double t);
		double      storage_size = 0.0;								//size in Gb required to store U
		double*     E    = NULL;									//Eigenenergies
		double*     U    = NULL;									//The real-value matrix of eigenvectors, H = O.E.O^T, E is treated as diagonal matrix
		double*     HI	 = NULL;									//Diagonal elements of the magnetic hamiltonian
		double      max_orthogonality_err	=	0.0;
		double      max_eigensystem_err		=	0.0;
		//Spectral properties
		void        get_HS0_matrix(double* HM, bool magnetic_field_on=true);
		//Saving/reading eigensystems
		void        read_eigensystem(string filename, bool noisy=true);
		void       write_eigensystem(string filename, bool noisy=true);
		//Misc
		void        rand_vec(double*    out, uint n, bool normalize=true);	//Random real vector, filled with Gaussian numbers with unit dispersion
		void        rand_vec(t_complex* out, uint n, bool normalize=true);
		t_complex*  rand_vec(uint n, bool normalize=true);
};

#endif
