/*
 *  OpenNL<T>: Generic Numerical Library
 *  Copyright (C) 2004 Bruno Levy
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  In addition, as a special exception, the INRIA gives permission to link the 
 *  code of this program with the CGAL library, and distribute linked combinations 
 *  including the two. You must obey the GNU General Public License in all respects 
 *  for all of the code used other than CGAL. 
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ISA-ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 *
 *  Laurent Saboret 01/2005: Change for CGAL:
 *		- Added OpenNL namespace
 *		- solve() returns true on success
 *		- test divisions by zero
 */
 
#ifndef __CONJUGATE_GRADIENT__
#define __CONJUGATE_GRADIENT__

#include "blas.h"

#include <cassert>

namespace OpenNL {


/**
 *  The Conjugate Gradient algorithm without preconditioner:
 *  Ashby, Manteuffel, Saylor
 *     A taxononmy for conjugate gradient methods
 *     SIAM J Numer Anal 27, 1542-1568 (1990)
 *
 * This implementation is inspired by the lsolver library,
 * by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de/IAM/Research/projectskr/lin_solver/
 *
 * @param A generic matrix, a function
 *   mult(const MATRIX& M, const double* x, double* y)
 *   needs to be defined.
 * @param b right hand side of the system.
 * @param x initial value.
 * @param eps threshold for the residual.
 * @param max_iter maximum number of iterations.
 */

template< class MATRIX, class VECTOR > class Solver_CG {
public:
    typedef MATRIX Matrix ;
    typedef VECTOR Vector ;
    typedef typename Vector::CoeffType CoeffType ;

    Solver_CG() {
        epsilon_ = 1e-6 ;
        max_iter_ = 0 ;
    }

    void set_epsilon(CoeffType eps) { epsilon_ = eps ; }
    void set_max_iter(unsigned int max_iter) { max_iter_ = max_iter ; }

 	// Solve the sparse linear system "A*x = b" for A SYMMETRIC POSITIVE. Return true on success.
	// 
	// Preconditions:
	// * A.dimension() == b.dimension()
	// * A.dimension() == x.dimension()
	bool solve(const MATRIX &A, const VECTOR& b, VECTOR& x)
	{
        assert(A.dimension() == b.dimension()) ;
        assert(A.dimension() == x.dimension()) ;
        assert (A.dimension() > 0);

        unsigned int n = A.dimension() ;						// (Square) matrix dimension
        unsigned int max_iter = max_iter_ ;						// Max number of iterations
        if(max_iter == 0) {
            max_iter = 5 * n ;
        }
        Vector g(n) ;
        Vector r(n) ;
        Vector p(n) ;
        CoeffType gg;
        unsigned int its=0;										// Loop counter
        CoeffType t, tau, sig, rho, gam;
        CoeffType bnorm2 = BLAS<Vector>::dot(b,b) ; 
        CoeffType err=epsilon_*epsilon_*bnorm2 ;				// Error to reach
        mult(A,x,g);
        BLAS<Vector>::axpy(-1,b,g);
        BLAS<Vector>::scal(-1,g);
        BLAS<Vector>::copy(g,r);
        gg=BLAS<Vector>::dot(g,g);								// Current error (g|g)

        while ( gg>err && its < max_iter) 
		{
			// Debug trace
			//if (its % 25 == 0)
			//	std::cerr << "Solver_CG<>::solve: gg(=" << gg << ") > err(=" << err << ")" << std::endl;

            mult(A,r,p);
            rho=BLAS<Vector>::dot(p,p);
            sig=BLAS<Vector>::dot(r,p);
            tau=BLAS<Vector>::dot(g,r);
			assert( fabs(sig)>1e-40 );
            t=tau/sig;
            BLAS<Vector>::axpy(t,r,x);
            BLAS<Vector>::axpy(-t,p,g);
			assert( fabs(tau)>1e-40 );
            gam=(t*t*rho-tau)/tau;
            BLAS<Vector>::scal(gam,r);
            BLAS<Vector>::axpy(1,g,r);
            gg=BLAS<Vector>::dot(g,g);
            its++;
        }

		bool success = (its < max_iter);
		if ( ! success )
			std::cerr << "Solver_CG<>::solve: failure: gg(=" << gg << ") > err(=" << err << ")" << std::endl;
		return success;
    }

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


}; // namespace OpenNL

#endif

