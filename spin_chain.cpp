#include "spin_chain.hpp"

uint count_bits(uint a, uint L) //Counts the number of binary 1's in the binary notation of a, up to L 
{
	uint r = 0;
	for(uint ib=0; ib<L; ib++)
		r += (a & (1<<ib))>>ib;
	return r;
}

uint rotate_bits(uint a, uint L, uint s) //Cyclic shift of L bits by s bits
{
	uint bm = (1<<L)-1;   //Bit mask with 1 for all bits at positions 0 ... L-1
	uint r = (a<<(s%L));  //Shifting by (s%L)
	uint q = (r>>L);      //What goes out of range?
	return ((r&bm)|q);    //(r&bm) sets to zero except for 0 ... L-1, and | q puts the bits which went out of range during ordinary shift back to the sequence 
}

SpinChain::SpinChain(uint L, double h, bool diagonalize, bool check, bool noisy, int rng_engine_seed, string eigensystem_file_name)
{
	if(L%2!=0)
	{
		fprintf(stderr, "We work in spin 0 sector and hence the code supports only even L!!!\n");
		fflush(stderr);
		exit(EXIT_FAILURE);
	};
	
	this->L   = L;
	this->N   = 1 << L;
	this->N2  = this->N*this->N; 

	NS0 = 0;
	for(uint ips=0; ips<N; ips++)
		NS0 += (count_bits(ips, L)==L/2? 1 : 0);
	NS02 = NS0*NS0;
	
	storage_size = (double)(NS02*sizeof(double))/(double)(1024*1024*1024);

	S0_basis = new uint[NS0];

	S0_lookup = new uint[N];
	for(uint ips=0; ips<N; ips++)
		S0_lookup[ips] = UINT_MAX;

	uint is = 0;
	for(uint ips=0; ips<N; ips++)
		if(count_bits(ips, L)==L/2)
		{
			S0_basis[is] = ips;
			S0_lookup[ips] = is;
			is++;
		};

	if(noisy){printf("\n >>>>>> L = %u, N=%u, h = %2.4lf, NS0=%u >>>>>>>>>>>>> \n\n", this->L, this->N, h, NS0); fflush(stdout);};
	
	hi = new double[L];
	//Initialize the random number generator
	if(rng_engine_seed==0)
		rng_engine.seed((int)round(omp_get_wtime()));
	else
		rng_engine.seed(rng_engine_seed);
		
	if(noisy) printf("\n hi = [");
	for(uint i=0; i<L; i++)
	{
		hi[i] = -0.5*h + h*rng_uniform_dist(rng_engine);
		if(noisy) printf("%+2.4lf, ", hi[i]);
	};
	if(noisy){ printf("]\n\n"); fflush(stdout); };
	
	if(diagonalize)
	{
		read_eigensystem(eigensystem_file_name);
		if(E==NULL || U==NULL)
		{
			diagonalize_HS0(check, false);
			write_eigensystem(eigensystem_file_name);
		};
	};
}

void        SpinChain::read_eigensystem(string filename, bool noisy)
{
	if(filename.length()>0)
	{
		std::ifstream afile(filename.c_str(), std::ios::binary);
		if(afile.is_open())
		{
			E = new double[NS0]; U = new double[NS02];
			afile.seekg(0, std::ios::beg);
			afile.read((char*)E, NS0*sizeof(double));
			if(afile)
				afile.read((char*)U, NS02*sizeof(double));
			if(!afile)
			{
				cerr << "Cannot read eigenenergies from " << filename << ", only " << afile.gcount() << "bytes out of " << ((NS0 + NS02)*sizeof(double)) << endl;
				delete [] E; delete [] U;
				E = NULL; U = NULL;
			};
			afile.close();
			if(noisy) cout << endl << "\t Data successfully loaded from the file " << filename << endl << endl;
		}
		else
			cerr << "Cannot open the file " << filename << " for binary reading!" << endl;
	};
}

void       SpinChain::write_eigensystem(string filename, bool noisy)
{
	if(filename.length()>0)
	{
		std::ofstream afile(filename.c_str(), std::ios::binary | std::ios::trunc);
		if(afile.is_open())
		{
			afile.seekp(0, std::ios::beg);
			afile.write((char*)E, NS0*sizeof(double));
			if(afile)
				afile.write((char*)U, NS02*sizeof(double));
			if(!afile) cerr << "Cannot write eigenenergies to" << filename << endl;
			afile.close();
			if(noisy) cout << endl << "\t Data successfully written to the file " << filename << endl << endl;
		}
		else
			cerr << "Cannot open the file " << filename << " for binary writing!" << endl;
	};
}

void SpinChain::H(double* in, double* out)
{
	for(uint ips=0; ips<N; ips++)
		out[ips] = 0.0;
		
	//This implements Eq. (30) in arXiv:1803.08050	
	for(uint i1=0; i1<L; i1++)
	{
		uint i2 = (i1+1)%L;
		uint bm1 = 1 << i1; //Bit at position i1 set to one
		uint bm2 = 1 << i2; //Bit at position i2 set to one
		for(uint ips1=0; ips1<N; ips1++)
		{
			bool C1 = (bool)(ips1 & bm1); //In the state ips1 bit at i1 is one
			bool C2 = (bool)(ips1 & bm2); //In the state ips1 bit at i2 is one
			bool C  = (C1 != C2);
			out[ips1] += ((C? -0.25 : +0.25 ) + (C1? -hi[i1] : hi[i1]))*in[ips1]; 
			if(C) out[ips1] += 0.5*in[(bm2^(bm1^ips1))];
		};
	};
}

template<typename T>
void SpinChain::HS0(T* in, T* out, bool magnetic_field_on) //Hamiltonian in spin0 sector
{
	for(uint is=0; is<NS0; is++)
		out[is] = 0.0;
		
	//This implements Eq. (30) in arXiv:1803.08050	
	for(uint i1=0; i1<L; i1++)
	{
		uint i2 = (i1+1)%L;
		uint bm1 = 1 << i1;
		uint bm2 = 1 << i2;
		for(uint is1=0; is1<NS0; is1++)
		{
			uint ips1 = S0_basis[is1];
			bool C1 = (bool)(ips1 & bm1); 
			bool C2 = (bool)(ips1 & bm2);
			bool C  = (C1 != C2);
			out[is1] += (C? -0.25 : +0.25 )*in[is1];
			
			if(magnetic_field_on)
				out[is1] += (C1? -hi[i1] : hi[i1])*in[is1];
			
			uint ips2 = (bm2^(bm1^ips1));
			uint is2 = S0_lookup[ips2];
			 
			if(C && (is2<NS0)) out[is1] += 0.5*in[is2];
		};
	};
}

void SpinChain::sz(double* in, double* out, uint i)       //local SigmaZ operator
{
	uint bm = 1 << i;
	for(uint ips=0; ips<N; ips++)
		out[ips] =  (bm&ips? -in[ips] : in[ips]);
}

void SpinChain::sz(double* in, double* out)		 //global SigmaZ operator
{
	for(uint ips=0; ips<N; ips++)
	{
		//Count the number of 1's in the bit string of ips - this will give the total spin
		uint n1=0;
		for(uint i=0; i<L; i++)	n1 += ((ips >> i) & 1);
		//Multiply by the total spin
		out[ips] =  (2.0*(double)n1 - (double)L)*in[ips];
	};	
}

void SpinChain::rand_vec(double* out, uint n, bool normalize)
{
	for(int ips=0; ips<n; ips++)
		out[ips] = rng_normal_dist(rng_engine);
	if(normalize)
	{
		double nn = norm(out, n);
		rescale(out, 1.0/nn, n);
	};
}

void SpinChain::rand_vec(t_complex* out, uint n, bool normalize)
{
	for(int ips=0; ips<n; ips++)
		out[ips] = rng_normal_dist(rng_engine) + 1.0i*rng_normal_dist(rng_engine);
	if(normalize)
	{
		double nn = norm(out, n);
		rescale(out, 1.0/nn + 0.0i, n);
	};
}

t_complex* SpinChain::rand_vec(uint n, bool normalize)
{
	t_complex* res = new t_complex[n];
	rand_vec(res, n, normalize);
	return res;
}

void    SpinChain::get_HS0_matrix(double* HM, bool magnetic_field_on)
{
	double* in  = new double[NS0];
	double* out = new double[NS0];
	
	for(uint ib=0; ib<NS0; ib++)
	{
		for(uint is=0; is<NS0; is++) in[is] = 0.0;
		in[ib] = 1.0;
		this->HS0(in, out);
		for(uint is=0; is<NS0; is++) HM[NS0*is + ib] = out[is];
	};
	
	delete [] in;
	delete [] out;
}

void    SpinChain::diagonalize_HS0(bool check, bool noisy)					//find the eigenspectrum of the Hamiltonian
{
	//Initialize storage for evecs
	E	= new double[NS0];
	U   = new double[NS02];
	
	get_HS0_matrix(U);
	
	double* HM  = NULL; 
	if(check)
	{
		HM = new double[NS02];
		std::copy(U, U + NS02, HM);
	};
	
	int res = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', NS0, U, NS0, E);
	if(res!=0){ fprintf(stderr, "\nSomething went wrong, LAPACKE_dsyev returned %i !!!\n", res); fflush(stderr);};
	
	//Sorting the eigensystem
	for(int istep=0; istep<NS0-1; istep++)
	{
		//Find maximal eigenvalue and its index
		int minie = -1; double minE = 1.0E+16;
		for(int ie=istep; ie<NS0; ie++)
			if(E[ie] < minE){minE = E[ie]; minie=ie;};
		//Place it on position istep
		double tmp = E[istep]; E[istep] = E[minie]; E[minie] = tmp;
		//Exchange also the eigenvectors 
		for(int i=0; i<NS0; i++)
		{
			tmp = U[i*NS0 + istep]; U[i*NS0 + istep] = U[i*NS0 + minie]; U[i*NS0 + minie] = tmp;
		};
	};
	
	if(check)
	{	
		double* C = new double[NS0];
		
		//Checking the eigensystem	
		max_eigensystem_err = 0.0;
		for(uint ie=0; ie<NS0; ie++)
		{
			for(uint i=0; i<NS0; i++)
			{
				C[i] = 0.0;
				for(uint j=0; j<NS0; j++)
					C[i] += HM[i*NS0 + j]*U[j*NS0 + ie];
			};
			double err = 0.0;
			for(uint i=0; i<NS0; i++)
				err += SQR(C[i] - E[ie]*U[i*NS0 + ie]);
			max_eigensystem_err = max(max_eigensystem_err, sqrt(err));
		};
		if(noisy){ printf("\t >> Max. eigenspectrum err:\t %2.4E\n", max_eigensystem_err); };
	
		//Orthogonality check
		max_orthogonality_err = 0.0;
		for(uint ie1=0; ie1<NS0; ie1++)
			for(uint ie2=0; ie2<NS0; ie2++)
			{
				double sp = 0.0;
				for(uint i=0; i<NS0; i++)
					sp += U[i*NS0 + ie1]*U[i*NS0 + ie2];
				sp = (ie1==ie2? sp - 1.0 : sp);
				max_orthogonality_err = max(max_orthogonality_err, fabs(sp));
			};
		if(noisy){ printf("\t >> Max. orthogonality err:\t %2.4E\n", max_orthogonality_err); };
		
		delete [] HM;
		delete [] C;
	};
}

t_complex* SpinChain::evolution_operator(t_complex dt)
{
	using namespace std::complex_literals;
	t_complex* res = NULL;
	if((E!=NULL) && (U!=NULL))
	{
		res = new t_complex[NS02];
		for(int i=0; i<NS0; i++)
			for(int j=0; j<NS0; j++)
			{
				res[i*NS0 + j] = 0.0 + 0.0i;
				for(int ie=0; ie<NS0; ie++)
					res[i*NS0 + j] += U[i*NS0 + ie]*exp(1i*dt*E[ie])*U[j*NS0 + ie];
			};
	};
	return res;
}

void SpinChain::init_magnetic_hamiltonian()
{
	HI  = new double[NS0];
	for(uint is=0; is<NS0; is++)
	{
		uint ips = S0_basis[is]; //Get the full bit mask for state number is
		HI[is] = 0.0;
		//First we save in HI the diagonal elements of the magnetic Hamiltonian
		for(uint x=0; x<L; x++)
		{
			bool C = (bool)(ips & (1 << x)); //If spin up at position x
			HI[is] += (C? -hi[x] : hi[x]);
		};
	};
}

void SpinChain::diagonalize_HS0_h0(double* E0, double* psi0, bool check, bool noisy)
{
	double a_time;
	
	bool* state_used  = new bool[NS0];
	uint  state_count = 0;
	double max_eigenstate_err    = 0.0;
	double max_orthogonality_err = 0.0;
	
	t_complex* psi = new t_complex[NS0]; //This will be our work vectors
	t_complex* chi = new t_complex[NS0];
	
	for(uint p=0; p<=L/2; p++)
	{
		std::vector<t_complex*> momentum_eigenstates;
		//Diagonalizing h=0 Hamiltonian in the basis of states with lattice momentum 2*pi*p/L
		std::fill_n(state_used, NS0, false);
		uint is=0;									 //We start with the first state in a list, loop over all states obtained by shifts, and mark these states as used
		while(is<NS0)
		{
			uint ips = S0_basis[is]; 				 //Bit representation of our unused basis state
			t_complex* mstate = new t_complex[NS0]; 	 //psi will contain momentum eigenstate constructed from |ips>
			std::fill_n(mstate, NS0, 0.0 + 0.0i);
			//Build up all momentum states with this state, some might be zero
			for(uint x=0; x<L; x++)
			{
				uint is1 = S0_lookup[ips];
				state_used[is1] = true;
				mstate[is1] += exp(2.0i*M_PI*(double)(p*x)/(double)L);
				ips = rotate_bits(ips, L);
			};
			//Check whether the new state has nonzero norm?
			double mstate_norm = real(scalar_prod(mstate, mstate, NS0));
			if(mstate_norm>1.0E-10)
			{
			    rescale(mstate, 1.0/sqrt(mstate_norm) + 0.0i, NS0);
			    momentum_eigenstates.push_back(mstate);	
			};
			while(is<NS0 && state_used[is]) is++; //Move is to the next state which was not used yet...
		};//End of loop over spin-0 states (iterator is)
		//Now the vector momentum_eigenstates contains the basis of momentum - p states
		uint ns = momentum_eigenstates.size();
		
		//Form the matrix of HS0 in the basis of p-eigenstates
		t_complex* HS0p = new t_complex[ns*ns];
		double*     ES  = new double[ns];

		for(uint j=0; j<ns; j++)
		{
			HS0(momentum_eigenstates[j], chi, false);
			for(uint i=0; i<ns; i++)
				HS0p[i*ns + j] = scalar_prod(momentum_eigenstates[i], chi, NS0);
		};
		
		//Now use LAPACK to diagonalize HS0p
		cout << ansi::green << "Diagonalizing (" << ansi::magenta << ns << " x " << ns << ansi::green << ") Hamiltonian for momentum " << ansi::magenta << p << ansi::green << " ... ";
		TIMING_START;
		int res = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', ns, HS0p, ns, ES);
		TIMING_END;
		cout << ansi::green << "Done in " << ansi::magenta << a_time << ansi::green << " sec." << ansi::reset << endl;
		if(res!=0) cerr << endl << ansi::white << "Something went wrong, LAPACKE_dsyev returned" << ansi::red << res << ansi::white << "!!!" << endl << endl;

		TIMING_START;
		if(noisy) cout << "Eigenstates:\t" << endl;
		for(uint ie=0; ie<ns; ie++)
		{
			std::fill_n(psi, NS0, 0.0 + 0.0i);
			
			for(uint i=0; i<ns; i++)
				A_pluseq_bB(psi, HS0p[i*ns + ie], momentum_eigenstates[i], NS0);
			//psi is now the eigenstate
			
			if(check)
			{
				//Norm of H|psi> - E|psi>
				HS0(psi, chi, false);
				A_pluseq_bB(chi, -ES[ie] + 0.0i, psi, NS0);
				max_eigenstate_err = std::max(max_eigenstate_err, norm(chi, NS0));
			};	
			if(noisy) cout << "\t E = " << std::showpos << std::setw(6) << ansi::cyan << ES[ie] << std::noshowpos << ansi::reset << ", |psi| = " << norm(psi, NS0) << ", ||H |psi> - E |psi>|| = " << norm(chi, NS0) << endl;
			
			//Add real part of psi to U0 and ES to E0
			for(uint i=0; i<NS0; i++)
				psi0[NS0*state_count + i] = real(psi[i]);
			double npsi0 = norm(&(psi0[NS0*state_count]), NS0);
			rescale(&(psi0[NS0*state_count]), 1.0/npsi0, NS0);
			E0[state_count] = ES[ie];
			state_count ++;
			
			if(p>0 && p<L/2)
			{
				//Also add imaginary part of the state vector, as it forms a linearly independent state for p->L-p
				for(uint i=0; i<NS0; i++)
					psi0[NS0*state_count + i] = imag(psi[i]);
				double npsi0 = norm(&(psi0[NS0*state_count]), NS0);
				rescale(&(psi0[NS0*state_count]), 1.0/npsi0, NS0);
				E0[state_count] = ES[ie];
				state_count ++;
			};
		}; //End of loop over eigenstates within the fixed-momentum block
		TIMING_END;
		cout << ansi::green << "Eigenstate prep: " << ansi::magenta << a_time << ansi::green << " sec." << ansi::reset << endl;

		delete [] HS0p;	delete [] ES;
		//Free all eigenstates from memory and reset the vector
		for(const auto& state: momentum_eigenstates) delete [] state;
		momentum_eigenstates.clear();
	}; //End of loop over all momenta
	
	if(state_count!=NS0) cerr << ansi::red << "Number of states " << state_count << "is not equal to the expected number" << NS0 << endl;
	
	if(check) cout << ansi::yellow << "Max. norm of |H|psi> - E|psi>| among " << ansi::magenta << state_count << ansi::yellow << " out of " << ansi::magenta << NS0 << ansi::yellow << " states is " << ansi::magenta << max_eigenstate_err << ansi::reset << endl;

	delete [] state_used;
	delete [] psi;
	delete [] chi;
}

void	SpinChain::trotter_evolve(t_complex* in, t_complex* out, double* E0, double* psi0, t_complex dt)
{
	std::fill(out, out + NS0, 0.0 + 0.0i);
	for(uint ie=0; ie<NS0; ie++)
	{
		t_complex psi0in = 0.0 + 0.0i;
		for(uint is=0; is<NS0; is++)
			psi0in += psi0[ie*NS0 + is]*exp(0.5i*HI[is]*dt)*in[is];
		psi0in *= exp(1.0i*dt*E0[ie]);
		for(uint is=0; is<NS0; is++)
			out[is] += exp(0.5i*HI[is]*dt)*psi0[ie*NS0 + is]*psi0in;
	};
}

t_complex*  SpinChain::trotter_evolution_operator(double* E0, double* psi0, t_complex dt)
{
	t_complex* U = new t_complex[NS02];
	for(uint is=0; is<NS0; is++)
		for(uint js=0; js<NS0; js++)
		{
			U[is*NS0 + js] = 0.0 + 0.0i;
			for(uint ie=0; ie<NS0; ie++)
			{
				t_complex fl = exp(0.5i*HI[is]*dt)*psi0[ie*NS0 + is];
				t_complex fr = exp(0.5i*HI[js]*dt)*psi0[ie*NS0 + js];
				U[is*NS0 + js] += fl*exp(1.0i*E0[ie]*dt)*fr;
			};
		};
	return U;
}

void SpinChain::exact_evolve(t_complex* in, t_complex* out, double t)
{
	if((E!=NULL) && (U!=NULL))
	{
		std::fill(out, out + NS0, 0.0 + 0.0i);
		for(uint ie=0; ie<NS0; ie++)
		{
			t_complex psi0in = 0.0 + 0.0i;
			for(uint is=0; is<NS0; is++)
				psi0in += U[is*NS0 + ie]*in[is];
			psi0in *= exp(1.0i*t*E[ie]);
			for(uint is=0; is<NS0; is++)
				out[is] += U[is*NS0 + ie]*psi0in;
		};
	};
}
