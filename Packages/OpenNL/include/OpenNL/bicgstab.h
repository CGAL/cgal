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
 *		- added comments and traces
 */

#ifndef __OPENNL_BICGSTAB__
#define __OPENNL_BICGSTAB__

#include <cmath>
#include <cfloat>
#include <climits>
#include <cassert>

namespace OpenNL {


// Utility macro to display a variable's value
// Usage: x=3.7; cerr << OPENNL_STREAM_TRACE(x) << endl;
//        prints
//        x=3.7
#define OPENNL_STREAM_TRACE(var) #var << "=" << var << " "


/**
 *  The BICGSTAB algorithm without preconditioner:
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

template < class MATRIX, class VECTOR> class Solver_BICGSTAB {
public:
    typedef MATRIX Matrix ;
    typedef VECTOR Vector ;
    typedef typename Vector::CoeffType CoeffType ;

    Solver_BICGSTAB() {
		// LS 03/2005: change epsilon from 1e-6 to 1e-4 to parameterize venus-loop.off w/ authalic/square method
		// LS 04/2005: change epsilon back to 1e-6 to increase precision
        epsilon_ = 1e-6 ;
        max_iter_ = 0 ;
    }

    void set_epsilon(CoeffType eps) { epsilon_ = eps ; }
    void set_max_iter(unsigned int max_iter) { max_iter_ = max_iter ; }

	// Solve the sparse linear system "A*x = b". Return true on success.
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

//#ifndef NDEBUG 
//		// Debug trace
//		fprintf(stderr, "\n");
//		if (n < 20)	// if small matrix, print it entirely
//		{
//			fprintf(stderr, "******************  A:  ******************\n");
//			for (int i=0; i<n; i++)  {
//				for (int j=0; j<n; j++)
//					fprintf(stderr, "%lf\t", (double)A.get_coef(i,j));
//				fprintf(stderr, "\n");
//			}
//			fprintf(stderr, "******************  B:  ******************\n");
//			for (int j=0; j<n; j++)
//				fprintf(stderr, "%lf\t", (double)b[j]);
//			fprintf(stderr, "\n");
//			fprintf(stderr, "******************************************\n");
//		}
//		else		// if large matrix, print only not null elements
//		{
//			fprintf(stderr, "******************  A*x=b  ******************\n");
//			for (int i=0; i<n; i++)  {
//				for (int j=0; j<n; j++)
//					if ( ! IsZero(A.get_coef(i,j)) )
//						fprintf(stderr, "A[%d][%d] = %lf\t", i, j, (double)A.get_coef(i,j));
//				fprintf(stderr, "\n");
//			}
//			for (int j=0; j<n; j++)
//				if ( ! IsZero(b[j]) )
//					fprintf(stderr, "b[%d] = %lf\t", j, (double)b[j]);
//			fprintf(stderr, "\n");
//			fprintf(stderr, "******************************************\n");
//		}
//#endif

        unsigned int max_iter = max_iter_ ;						// Max number of iterations
        if(max_iter == 0) {
			// LS 03/2005: change max_iter from 5*n to 10*n to parameterize venus-loop.off w/ authalic/square method
            max_iter = 10*n ;	
        }
        Vector rT(n) ;											// Initial residue rT=Ax-b
        Vector d(n) ;
        Vector h(n) ;
        Vector Ad(n) ;
        Vector t(n) ;
        CoeffType rTh, rTAd=0, rTr, alpha=0, beta=0, omega=0, ht=0, tt=0;
        unsigned int its=0;										// Loop counter
        CoeffType err=epsilon_*epsilon_*BLAS<Vector>::dot(b,b);	// Error to reach

        Vector r(n) ;											// Current residue r=A*x-b
        mult(A,x,r);
        BLAS<Vector>::axpy(-1,b,r);

	    // Initially, d=h=rT=r=A*x-b
        BLAS<Vector>::copy(r,d);								// d = r
        BLAS<Vector>::copy(d,h);								// h = d
        BLAS<Vector>::copy(h,rT);								// rT = h
        assert( ! IsZero( BLAS<Vector>::dot(rT,rT) ) );

		rTh=BLAS<Vector>::dot(rT,h);							// rTh = (rT|h)
		assert( ! IsZero(rTh) );
        rTr=BLAS<Vector>::dot(r,r);								// Current error rTr = (r|r)

		while ( rTr>err && its < max_iter) 
		{
			mult(A,d,Ad);										// Ad = A*d
            rTAd=BLAS<Vector>::dot(rT,Ad);						// rTAd = (rT|Ad)
			assert( ! IsZero(rTAd) );
			alpha=rTh/rTAd;										// alpha = rTh/rTAd
            BLAS<Vector>::axpy(-alpha,Ad,r);					// r = r - alpha*Ad
            BLAS<Vector>::axpy(-alpha,Ad,h);					// h = h - alpha*Ad
            mult(A,h,t);										// t = A*h
            ht=BLAS<Vector>::dot(h,t);							// ht = (h|t)
            tt=BLAS<Vector>::dot(t,t);							// tt = (t|t)
			if ( IsZero(ht) || IsZero(tt) )
                omega = 0 ;
		    else
	      		omega = ht/tt;									// omega = ht/tt
            BLAS<Vector>::axpy(-alpha,d,x);						// x = x - alpha*d
            BLAS<Vector>::axpy(-omega,h,x);						// x = x - omega*h
            BLAS<Vector>::axpy(-omega,t,r);						// r = r - omega*t
	        rTr=BLAS<Vector>::dot(r,r);							// Current error rTr = (r|r)
            BLAS<Vector>::axpy(-omega,t,h);						// h = h - omega*t
//#ifndef NDEBUG 
//			// Debug trace
//			std::cerr << "Solver_BICGSTAB<>::solve: " << OPENNL_STREAM_TRACE(its) << OPENNL_STREAM_TRACE(rTr) 
//				      << OPENNL_STREAM_TRACE(alpha) << OPENNL_STREAM_TRACE(beta) << OPENNL_STREAM_TRACE(omega) 
//					  << OPENNL_STREAM_TRACE(rTh) 
//					  << OPENNL_STREAM_TRACE(rTAd) << OPENNL_STREAM_TRACE(ht) << OPENNL_STREAM_TRACE(tt) 
//					  << std::endl;
//#endif
			if (IsZero(omega)) {								// LS 03/2005: break to avoid division by zero (see Laspack implementation)
				std::cerr << "Solver_BICGSTAB<>::solve: warning: omega = 0" << std::endl;
				break;		
			}
			if (IsZero(rTh)) {									// LS 04/2005: don't know what do do if division by zero 
				std::cerr << "Solver_BICGSTAB<>::solve: error: rTh = 0" << std::endl;
				break;	
			}
            beta=(alpha/omega)/rTh; 
			rTh=BLAS<Vector>::dot(rT,h); 						// rTh = (rT|h)
			beta*=rTh;											// beta = (rTh / previous rTh) * (alpha / omega)
            BLAS<Vector>::scal(beta,d);							// d = beta*d
            BLAS<Vector>::axpy(1,h,d);							// d = d + h
            BLAS<Vector>::axpy(-beta*omega,Ad,d);				// d = d - beta*omega*Ad
            its++ ;
        }

		bool success = (rTr <= err);
#ifndef NDEBUG 
		if ( ! success )
			std::cerr << "Solver_BICGSTAB<>::solve failure: "
			          << OPENNL_STREAM_TRACE(its) << OPENNL_STREAM_TRACE(max_iter) 
				      << OPENNL_STREAM_TRACE(rTr) << OPENNL_STREAM_TRACE(err) << OPENNL_STREAM_TRACE(alpha) << OPENNL_STREAM_TRACE(beta) << OPENNL_STREAM_TRACE(omega) 
					  << OPENNL_STREAM_TRACE(rTh) << OPENNL_STREAM_TRACE(rTAd) << OPENNL_STREAM_TRACE(ht) << OPENNL_STREAM_TRACE(tt) 
				      << std::endl;
#endif
		return success;
    }

private:
	// Test if a floating point number is (close to) 0.0
	static inline bool IsZero(CoeffType a) 
	{
		// LS 03/2005: replace Lsolver test by Laspack test to parameterize venus-loop.off w/ authalic/square method
		//#define IsZero(a) (fabs(a) < 1e-40)
		return (fabs(a) < 10.0 * std::numeric_limits<CoeffType>::min());
	}

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


}; // namespace OpenNL

#endif

