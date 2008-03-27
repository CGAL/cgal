//=============================================================================
//
//  CLASS Taucs_solver:
//  direct solver for symmetric positive definite sparse systems
//
//=============================================================================


#ifndef CGAL_TAUCS_SOLVER_H
#define CGAL_TAUCS_SOLVER_H

#include <vector>
#include <set>

#include <CGAL/Taucs_fix.h>

// Uncomment the next line to see libraries selected by auto-link
//#define CGAL_LIB_DIAGNOSTIC
#include <CGAL/auto_link/TAUCS.h>

#include <CGAL/surface_reconstruction_assertions.h>

CGAL_BEGIN_NAMESPACE


class Taucs_solver
{
public:
	Taucs_solver()
		: PAP(0),
		  L(0),
			SL(0),
			n_rows(0),
			perm(0),
			invperm(0),
			supernodal_(true)
	{
		m_io_handle = NULL;

#ifdef DEBUG_TRACE
    // Turn on TAUCS trace
    std::cerr.flush();
    taucs_logfile((char*)"stderr");
#endif
	}

	~Taucs_solver()
	{
		delete_matrices();
	}

	void begin_row()
	{
		if (colptr.empty() || colptr.back() != (int)values.size())
		{
			colptr.push_back((int)values.size());
			n_rows = (int)colptr.size()-1;
		}
	}
	void add_value(int _i,
		             double _val,
								 const bool symmetric = true)
	{
		if(symmetric)
		{
			if(_i <= n_rows)
			{
				values.push_back(_val);
				rowind.push_back(_i);
			}
		}
		else
		{
			values.push_back(_val);
			rowind.push_back(_i);
		}
	}


	void end_row()
	{
		if (colptr.empty() || colptr.back() != (int)values.size())
		{
			colptr.push_back((int)values.size());
			n_rows = (int)(colptr.size()-1);
		}
	}

	void finalize_matrix()
	{
		// setup ccs matrix
		A.n        = (int)(colptr.size()-1);
		A.m        = (int)(colptr.size()-1);
		A.flags    = (TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
		A.colptr   = &colptr[0];
		A.rowind   = &rowind[0];
		A.values.d = &values[0];
	}

	bool factorize(bool _use_supernodal = true)
	{
		supernodal_ = _use_supernodal;

		// delete old matrices
		delete_matrices();

		finalize_matrix();

		// bandlimitation
		taucs_ccs_order(&A, &perm, &invperm, (char*)"metis");
		PAP = taucs_ccs_permute_symmetrically(&A, perm, invperm);
		if (!PAP)
		{
			std::cerr << "Taucs_solver: permutation failed\n";
			return false;
		}

		// Cholesky factorization
		if (supernodal_)  SL = taucs_ccs_factor_llt_mf (PAP);
		else               L = taucs_ccs_factor_llt    (PAP, 0, 0);
		if (!(L || SL))
		{
			std::cerr << "Taucs_solver: factorization failed\n";
			return false;
		}

		return true;
	}

	bool factorize_ooc()
	{
		// delete old matrices
		delete_matrices();

		finalize_matrix();

		// bandlimitation
		taucs_ccs_order(&A, &perm, &invperm, (char*)"metis");
		PAP = taucs_ccs_permute_symmetrically(&A, perm, invperm);
		if (!PAP)
		{
			std::cerr << "Taucs_solver: permutation failed\n";
			return false;
		}

		// out-of-core Cholesky factorization.
        // LS 03/2008: ooc file opening will fail if 2 instances of the application
        //             run at the same time. Better use tempnam().
		m_io_handle = taucs_io_create_multifile((char*)"taucs-ooc");
		if(m_io_handle == NULL)
		{
			//CGAL_surface_reconstruction_assertion(false);
			std::cerr << "ooc file opening failed\n";
			return false;
		}

		double available_memory = taucs_available_memory_size();
		int result = taucs_ooc_factor_llt(PAP,m_io_handle,available_memory);
		if(result != TAUCS_SUCCESS)
		{
			//CGAL_surface_reconstruction_assertion(false);
			std::cerr << "ooc factorization failed\n";
			return false;
		}
		return true;
	}



	bool solve(std::vector<double>& _b, std::vector<double>& _x)
	{
		const unsigned int N = A.n;

		if (N != _b.size() || N != _x.size())
		{
			std::cerr << "Taucs_solver: matrix size doesn't match vector size\n";
			return false;
		}

		std::vector<double>  PB(N), PX(N);


		// permute rhs
		for (unsigned int i=0; i<N; ++i)
			PB[i] = _b[perm[i]];


		// solve by back-substitution
		if ((supernodal_ ?
				taucs_supernodal_solve_llt(SL, &PX[0], &PB[0]) :
				taucs_ccs_solve_llt(L, &PX[0], &PB[0]))
				!= TAUCS_SUCCESS)
		{
			std::cerr << "Taucs_solver: back-substitution failed\n";
			return false;
		}


		// re-permute x
		for (unsigned int i=0; i<N; ++i)
			_x[i] = PX[invperm[i]];

		return true;
	}

	bool solve_ooc(std::vector<double>& _b,
		             std::vector<double>& _x)
	{
		const unsigned int N = A.n;

		if (N != _b.size() || N != _x.size())
		{
			std::cerr << "Taucs_solver: matrix size doesn't match vector size\n";
			return false;
		}

		std::vector<double>  PB(N), PX(N);


		// permute rhs
		for (unsigned int i=0; i<N; ++i)
			PB[i] = _b[perm[i]];


		if(taucs_ooc_solve_llt(m_io_handle,&PX[0],&PB[0]) != TAUCS_SUCCESS)
		{
			std::cerr << "Taucs_solver: back-substitution failed\n";
			return false;
		}

		// re-permute x
		for (unsigned int i=0; i<N; ++i)
			_x[i] = PX[invperm[i]];

		return true;
	}


	bool solve_minres(std::vector<double>& _b,
		                std::vector<double>& _x)
	{
		double* x = new double[_x.size()];
		double* b = new double[_b.size()];
		for(unsigned int i=0;i<_b.size();i++)
		{
			x[i] = _x[i];
			b[i] = _b[i];
		}
		int result = taucs_minres(&A,NULL,NULL,x,b,1000,0.01);
		delete [] x;
		delete [] b;
		if(result == TAUCS_SUCCESS)
			return true;
		return false;
	}

	bool solve_linear(std::vector<double>& _b,
		                std::vector<double>& _x)
	{
		double* x = new double[_x.size()];
		double* b = new double[_b.size()];
		for(unsigned int i=0;i<_b.size();i++)
		{
			x[i] = _x[i];
			b[i] = _b[i];
		}
		// char* options[] = { "taucs.factor.LU=true", NULL };
		char* options[] = {NULL};
		int result = taucs_linsolve(&A,NULL,1,x,b,options,NULL);
		delete [] x;
		delete [] b;
		if(result == TAUCS_SUCCESS)
			return true;
		return false;
	}

	bool solve_conjugate_gradient(std::vector<double>& b,
	                              std::vector<double>& x,
																const int itermax,
																const double convergetol)
	{
		finalize_matrix();
		int result = taucs_conjugate_gradients(&A,NULL,NULL,&x[0],&b[0],itermax,convergetol);
		return (result == TAUCS_SUCCESS) ? true : false;
	}


	bool solve(std::vector<double>& b,
             std::vector<double>& x,
             const unsigned int dim)
	{
		const unsigned int N = A.n;
		CGAL_surface_reconstruction_assertion(dim >= 1);

		if (N != b.size()/dim || N != x.size()/dim)
		{
			std::cerr << "Taucs_solver: matrix size doesn't match vector size\n";
			return false;
		}

		std::vector<double>	PB(N);
		std::vector<double>	PX(N);

		// solver component-wise
		for (unsigned int c=0;c<dim;c++)
		{
			// permute rhs
			for(unsigned int i=0; i<N; ++i)
				PB[i] = b[c*N+perm[i]];

			// solve by back-substitution
			if((supernodal_ ?
		      taucs_supernodal_solve_llt(SL, &PX[0], &PB[0]) :
		      taucs_ccs_solve_llt(L, &PX[0], &PB[0])) != TAUCS_SUCCESS)
			{
				std::cerr << "Taucs_solver: back-substitution failed\n";
				return false;
			}

			// re-permute x
			for (unsigned int i=0; i<N; ++i)
				x[c*N + i] = PX[invperm[i]];
		}
		return true;
	}

private:

	void delete_matrices()
	{
		if (PAP)      { taucs_ccs_free(PAP);               PAP     = 0; }
		if (L)        { taucs_ccs_free(L);                 L       = 0; }
		if (SL)       { taucs_supernodal_factor_free(SL);  SL      = 0; }
		if (perm)     { free(perm);                        perm    = 0; }
		if (invperm)  { free(invperm);                     invperm = 0; }

		if(m_io_handle != NULL)
		{
			taucs_io_delete(m_io_handle);
			// taucs_io_close(m_io_handle);
			m_io_handle = NULL;
		}
	}


private:

  taucs_ccs_matrix           A, *PAP, *L;
  void                       *SL;
  std::vector<double>        values;
  std::vector<int>           colptr;
  std::vector<int>           rowind;
  int                        n_rows;
  int                        *perm, *invperm;
  bool                       supernodal_;
	taucs_io_handle *m_io_handle;
};


CGAL_END_NAMESPACE

#endif // CGAL_TAUCS_SOLVER_H
