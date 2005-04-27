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
 *  Laurent Saboret 01/2005-04/2005: Changes for CGAL:
 *      - Added OpenNL namespace
 *      - solve() returns true on success
 *      - test divisions by zero with IsZero() method
 *      - added comments and traces
 */

#ifndef __OPENNL_BICGSTAB__
#define __OPENNL_BICGSTAB__

#include <OpenNL/blas.h>

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
    bool solve(const MATRIX &A, const VECTOR& b, VECTOR& x) {
        assert(A.dimension() == b.dimension()) ;
        assert(A.dimension() == x.dimension()) ;
        assert (A.dimension() > 0);

        unsigned int n = A.dimension() ;                        // (Square) matrix dimension
        unsigned int max_iter = max_iter_ ;                     // Max number of iterations
        if(max_iter == 0) {
            max_iter = 5 * n ;
        }
        Vector rT(n) ;                                          // Initial residue rT=Ax-b
        Vector d(n) ;
        Vector h(n) ;
        Vector u(n) ;
        Vector Ad(n) ;
        Vector t(n) ;
        Vector& s = h ;
        CoeffType rTh, rTAd, rTr, alpha, beta, omega, st, tt;
        unsigned int its=0;                                     // Loop counter
        CoeffType err=epsilon_*epsilon_*BLAS<Vector>::dot(b,b); // Error to reach
        Vector r(n) ;                                           // Current residue r=A*x-b
        mult(A,x,r);
        BLAS<Vector>::axpy(-1,b,r);
        BLAS<Vector>::copy(r,d);
        BLAS<Vector>::copy(d,h);
        BLAS<Vector>::copy(h,rT);
        assert( ! IsZero( BLAS<Vector>::dot(rT,rT) ) );
        rTh=BLAS<Vector>::dot(rT,h);                            // rTh = (rT|h)
        assert( ! IsZero(rTh) );
        rTr=BLAS<Vector>::dot(r,r);                             // Current error rTr = (r|r)

        while ( rTr>err && its < max_iter) {
            mult(A,d,Ad);
            rTAd=BLAS<Vector>::dot(rT,Ad);
            assert( ! IsZero(rTAd) );
            alpha=rTh/rTAd;
            BLAS<Vector>::axpy(-alpha,Ad,r);
            BLAS<Vector>::copy(h,s);
            BLAS<Vector>::axpy(-alpha,Ad,s);
            mult(A,s,t);
            BLAS<Vector>::axpy(1,t,u);
            BLAS<Vector>::scal(alpha,u);
            st=BLAS<Vector>::dot(s,t);
            tt=BLAS<Vector>::dot(t,t);
            if ( IsZero(st) || IsZero(tt) )
                omega = 0 ;
            else
                omega = st/tt;
            BLAS<Vector>::axpy(-omega,t,r);
            BLAS<Vector>::axpy(-alpha,d,x);
            BLAS<Vector>::axpy(-omega,s,x);
            BLAS<Vector>::copy(s,h);
            BLAS<Vector>::axpy(-omega,t,h);
            assert( ! IsZero(omega) );
            assert( ! IsZero(rTh) );
            beta=(alpha/omega)/rTh; rTh=BLAS<Vector>::dot(rT,h); beta*=rTh;
            BLAS<Vector>::scal(beta,d);
            BLAS<Vector>::axpy(1,h,d);
            BLAS<Vector>::axpy(-beta*omega,Ad,d);
            rTr=BLAS<Vector>::dot(r,r);
            its++ ;
        }

        bool success = (rTr <= err);
#ifndef NDEBUG
        // Trace on error
        if ( ! success )
            std::cerr << "Solver_BICGSTAB<>::solve failure: "
                      << "(" << OPENNL_STREAM_TRACE(its) << OPENNL_STREAM_TRACE(max_iter)
                             << OPENNL_STREAM_TRACE(rTr) << OPENNL_STREAM_TRACE(err)
                      << ")" << std::endl;
#endif
        return success;
    }

private:
    // Test if a floating point number is (close to) 0.0
    static inline bool IsZero(CoeffType a)
    {
        return (fabs(a) < 10.0 * std::numeric_limits<CoeffType>::min());
    }

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


}; // namespace OpenNL

#endif // __OPENNL_BICGSTAB__

