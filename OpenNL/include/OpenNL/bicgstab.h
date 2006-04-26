/*
 * author:  Bruno Levy, INRIA, project ALICE
 * website: http://www.loria.fr/~levy/software
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License version 2.1 as published by the Free Software Foundation
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Scientific work that use this software can reference the website and
 * the following publication:
 *
 * @INPROCEEDINGS {levy:NMDGP:05,
 *    AUTHOR = Bruno Levy,
 *    TITLE  = Numerical Methods for Digital Geometry Processing,
 *    BOOKTITLE =Israel Korea Bi-National Conference,
 *    YEAR=November 2005,
 *    URL=http://www.loria.fr/~levy/php/article.php?pub=../publications/papers/2005/Numerics
 * }
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
            if (IsZero(rTAd))
                break;                                          // stop if bad conditioning
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
            if (IsZero(omega))
                break;                                          // stop if bad conditioning
            if (IsZero(rTh))
                break;                                          // stop if bad conditioning
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
        return (std::fabs(a) < 10.0 * std::numeric_limits<CoeffType>::min());
    }

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


}; // namespace OpenNL

#endif // __OPENNL_BICGSTAB__

