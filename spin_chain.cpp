#include "spin_chain.hpp"

uint count_bits(uint a, uint L) //Counts the number of binary 1's in the binary notation of a, up to L 
{
	uint r = 0;
	for(uint ib=0; ib<L; ib++)
		r += (a & (1<<ib))>>ib;
	return r;
}

SpinChain::SpinChain(uint L, double h, bool diagonalize, bool check)
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

	printf("\n >>>>>> L = %u, N=%u, h = %2.4lf, NS0=%u >>>>>>>>>>>>> \n\n", this->L, this->N, h, NS0); fflush(stdout);
	
	hi = new double[L];
	//Initialize the random number generator
	rng_engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
	printf("\n hi = [");
	for(uint i=0; i<L; i++)
	{
		hi[i] = -0.5*h + h*rng_uniform_dist(rng_engine);
		printf("%+2.4lf, ", hi[i]);
	};
	printf("]\n\n"); fflush(stdout);
	
	//Test hermiticity and commutativity with sz
	if(check)
	{
		double* psi1 = new double[N];
		double* psi2 = new double[N];
		double* psi3 = new double[N];
		double* tmp  = new double[N];
		
		rand_vec(psi1, N); rand_vec(psi2, N);
		
		//Hermiticity test
		H(psi1, tmp);
		double r1 = scalar_prod(psi2, tmp, N);
		H(psi2, tmp);
		double r2 = scalar_prod(psi1, tmp, N);
		printf("\t >> Hermiticity error (full):        \t %2.4E\n", fabs(r1 - r2));
		
		//Hermiticity test for spin-0 Hamiltonian
		rand_vec(psi1, NS0); rand_vec(psi2, NS0);
		HS0(psi1, tmp);
		r1 = scalar_prod(psi2, tmp, NS0);
		HS0(psi2, tmp);
		r2 = scalar_prod(psi1, tmp, NS0);
		printf("\t >> Hermiticity error (spin zero):    \t %2.4E\n", fabs(r1 - r2));
		
		//Commuting with total spin
		sz(psi1, tmp); H( tmp, psi2);
		H( psi1, tmp); sz(tmp, psi3);
		printf("\t >> SZ commutativity error:\t %2.4E\n", norm_diff(psi2, psi3));
		
		delete [] psi1;
		delete [] psi2;
		delete [] psi3;
		delete [] tmp;
	};
	
	if(diagonalize)
		diagonalize_HS0(check);
}

void SpinChain::H(double* in, double* out)
{
	for(uint ips=0; ips<N; ips++)
		out[ips] = 0.0;
		
	//This implements Eq. (30) in arXiv:1803.08050	
	for(uint i1=0; i1<L; i1++)
	{
		uint i2 = (i1+1)%L;
		uint bm1 = 1 << i1;
		uint bm2 = 1 << i2;
		for(uint ips1=0; ips1<N; ips1++)
		{
			bool C1 = (bool)(ips1 & bm1); 
			bool C2 = (bool)(ips1 & bm2);
			bool C  = (C1 != C2);
			out[ips1] += ((C? -0.25 : +0.25 ) + (C1? -hi[i1] : hi[i1]))*in[ips1]; 
			if(C) out[ips1] += 0.5*in[(bm2^(bm1^ips1))];
		};
	};
}

void SpinChain::HS0(double* in, double* out) //Hamiltonian in spin0 sector
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
			out[is1] += ((C? -0.25 : +0.25 ) + (C1? -hi[i1] : hi[i1]))*in[is1];
			
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

void SpinChain::rand_vec(double* in, uint n, bool normalize)
{
	for(int ips=0; ips<n; ips++)
		in[ips] = rng_normal_dist(rng_engine);
	if(normalize)
	{
		double nn = 0.0;
		for(int ips=0; ips<n; ips++)
			nn += in[ips]*in[ips];
		nn = 1.0/sqrt(nn);
		for(int ips=0; ips<n; ips++)
			in[ips] *= nn;
	};
}

double  SpinChain::scalar_prod(double* psi1, double* psi2, uint n)
{
	return cblas_ddot(n, psi1, 1, psi2, 1);
}

void    SpinChain::aA_plus_bB(double* out, double a, double* A, double b, double* B, uint n)
{
	for(int ips=0; ips<n; ips++)
		out[ips] = a*A[ips] + b*B[ips];
}

double    SpinChain::norm_diff(double* psi1, double* psi2, uint n)
{
	double nn = 0.0;
	for(int ips=0; ips<n; ips++)
		nn += (psi1[ips] - psi2[ips])*(psi1[ips] - psi2[ips]);
	return sqrt(nn);
}

void    SpinChain::get_H_matrix(double* HM)
{
	double* in  = new double[N];
	double* out = new double[N];
	
	for(uint ib=0; ib<N; ib++)
	{
		for(uint ips=0; ips<N; ips++) in[ips] = 0.0;
		in[ib] = 1.0;
		this->H(in, out);
		for(uint ips=0; ips<N; ips++) HM[N*ips + ib] = out[ips];
	};
	
	delete [] in;
	delete [] out;
}

void    SpinChain::get_HS0_matrix(double* HM)
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

void    SpinChain::diagonalize_HS0(bool check)					//find the eigenspectrum of the Hamiltonian
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
	if(res!=0) printf("\nSomething went wrong, LAPACKE_dsyev returned %i !!!\n", res);
	
	if(check)
	{
		double max_err = 0.0;	
		double* C = new double[NS0];
		
		//Checking the eigensystem	
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
			max_err = max(max_err, sqrt(err));
		};
		printf("\t >> Max. eigenspectrum err:\t %2.4E\n", max_err);
	
		//Orthogonality check
		max_err = 0.0;
		for(uint ie1=0; ie1<NS0; ie1++)
			for(uint ie2=0; ie2<NS0; ie2++)
			{
				double sp = 0.0;
				for(uint i=0; i<NS0; i++)
					sp += U[i*NS0 + ie1]*U[i*NS0 + ie2];
				sp = (ie1==ie2? sp - 1.0 : sp);
				max_err = max(max_err, fabs(sp));
			};
		printf("\t >> Max. orthogonality err:\t %2.4E\n", max_err);
		
		delete [] HM;
		delete [] C;
	};
}
