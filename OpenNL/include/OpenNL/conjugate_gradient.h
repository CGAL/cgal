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

#ifndef __OPENNL_CONJUGATE_GRADIENT__
#define __OPENNL_CONJUGATE_GRADIENT__

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

    // Solve the sparse linear system "A*x = b" for A SYMMETRIC POSITIVE
    // Return true on success
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
        Vector g(n) ;
        Vector r(n) ;
        Vector p(n) ;
        unsigned int its=0;                                     // Loop counter
        CoeffType t, tau, sig, rho, gam;
        CoeffType bnorm2 = BLAS<Vector>::dot(b,b) ;
        CoeffType err=epsilon_*epsilon_*bnorm2 ;                // Error to reach
        // Current residue g=b-A*x
        mult(A,x,g);
        BLAS<Vector>::axpy(-1,b,g);
        BLAS<Vector>::scal(-1,g);
        // Initially, r=g=b-A*x
        BLAS<Vector>::copy(g,r);                                // r = g
        CoeffType gg=BLAS<Vector>::dot(g,g);                    // Current error gg = (g|g)

        while ( gg>err && its < max_iter)
        {
            mult(A,r,p);
            rho=BLAS<Vector>::dot(p,p);
            sig=BLAS<Vector>::dot(r,p);
            tau=BLAS<Vector>::dot(g,r);
            if (IsZero(sig))
                break;                                          // stop if bad conditioning
            t=tau/sig;
            BLAS<Vector>::axpy(t,r,x);
            BLAS<Vector>::axpy(-t,p,g);
            if (IsZero(tau))
                break;                                          // stop if bad conditioning
            gam=(t*t*rho-tau)/tau;
            BLAS<Vector>::scal(gam,r);
            BLAS<Vector>::axpy(1,g,r);
            gg=BLAS<Vector>::dot(g,g);                          // Current error gg = (g|g)
            its++;
        }

        bool success = (gg <= err);
#ifndef NDEBUG
        // Trace on error
        if ( ! success )
            std::cerr << "Solver_CG<>::solve: failure: "
                      << "(" << OPENNL_STREAM_TRACE(its) << OPENNL_STREAM_TRACE(max_iter)
                             << OPENNL_STREAM_TRACE(gg) << OPENNL_STREAM_TRACE(err)
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

    // Test if 2 floating point numbers are very close
    static inline bool AreEqual(CoeffType a, CoeffType b)
    {
        if (IsZero(a))
            return IsZero(b);
        else
            return IsZero(b/a - 1.0);
    }

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


}; // namespace OpenNL

#endif // __OPENNL_CONJUGATE_GRADIENT__

