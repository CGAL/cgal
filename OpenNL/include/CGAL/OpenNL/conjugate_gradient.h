// Copyright (c) 2005-2008  Inria Loria (France).
/*
 * author:  Bruno Levy, INRIA, project ALICE
 * website: http://www.loria.fr/~levy/software
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
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
 *  Laurent Saboret 2005-2006: Changes for CGAL:
 *      - Added OpenNL namespace
 *      - solve() returns true on success
 *      - check divisions by zero
 *      - added comments
 *      - copied Conjugate Gradient algorithm WITH preconditioner from Graphite 1.9 code
 */

#ifndef __OPENNL_CONJUGATE_GRADIENT__
#define __OPENNL_CONJUGATE_GRADIENT__

#include <CGAL/OpenNL/blas.h>
#include <CGAL/assertions.h>

#include <cmath>
#include <cfloat>
#include <climits>

namespace OpenNL {


/**
 *  The Conjugate Gradient algorithm WITHOUT preconditioner:
 *  Ashby, Manteuffel, Saylor
 *     A taxononmy for conjugate gradient methods
 *     SIAM J Numer Anal 27, 1542-1568 (1990)
 *
 * This implementation is inspired by the lsolver library,
 * by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de/IAM/Research/projectskr/lin_solver/
 *
 * @param A generic square matrix; a function
 *   mult(const MATRIX& M, const double* x, double* y)
 * and a member function 
 *   int dimension() const
 * must to be defined.
 * @param b right hand side of the system.
 * @param x initial value.
 * @param eps threshold for the residual.
 * @param max_iter maximum number of iterations.
 */

template<class MATRIX, class VECTOR> class Solver_CG {
public:
    typedef MATRIX Matrix ;
    typedef VECTOR Vector ;
    typedef typename Vector::CoeffType CoeffType ;

public:
    Solver_CG() {
        epsilon_ = 1e-6 ;
        max_iter_ = 0 ;
    }

    // Default copy constructor, operator =() and destructor are fine

    void set_epsilon(CoeffType eps) { epsilon_ = eps ; }
    void set_max_iter(unsigned int max_iter) { max_iter_ = max_iter ; }

    // Solve the sparse linear system "A*x = b" for A symmetric positive definite
    // Return true on success
    bool solve(const MATRIX &A, const VECTOR& b, VECTOR& x) 
    {
#ifdef DEBUG_TRACE
        std::cerr << "  Call Conjugate Gradient" << std::endl;
#endif
        CGAL_assertion(A.dimension() > 0);
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
        // residue g=b-A*x
        mult(A,x,g);
        BLAS<Vector>::axpy(-1,b,g);
        BLAS<Vector>::scal(-1,g);
        // Initially, r=g=b-A*x
        BLAS<Vector>::copy(g,r);                                // r = g

        CoeffType gg=BLAS<Vector>::dot(g,g);                    // error gg = (g|g)
        while ( gg>err && its < max_iter) {
            mult(A,r,p);
            rho=BLAS<Vector>::dot(p,p);
            sig=BLAS<Vector>::dot(r,p);
            tau=BLAS<Vector>::dot(g,r);
            CGAL_assertion(sig != 0.0);
            t=tau/sig;
            BLAS<Vector>::axpy(t,r,x);
            BLAS<Vector>::axpy(-t,p,g);
            CGAL_assertion(tau != 0.0);
            gam=(t*t*rho-tau)/tau;
            BLAS<Vector>::scal(gam,r);
            BLAS<Vector>::axpy(1,g,r);
            gg=BLAS<Vector>::dot(g,g);                          // Update error gg = (g|g)
            ++its;
        }

        bool success = (gg <= err);
        return success;
    }

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


/**
 *  The Conjugate Gradient algorithm WITH preconditioner:
 *  Ashby, Manteuffel, Saylor
 *     A taxononmy for conjugate gradient methods
 *     SIAM J Numer Anal 27, 1542-1568 (1990)
 *
 * This implementation is inspired by the lsolver library,
 * by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de/IAM/Research/projectskr/lin_solver/
 *
 * @param A generic square matrix; a function
 *   mult(const MATRIX& M, const double* x, double* y)
 * and a member function 
 *   int dimension() const
 * must to be defined.
 * @param C preconditioner; a function
 *   mult(const PC_MATRIX& C, const double* x, double* y)
 * needs to be defined.
 * @param b right hand side of the system.
 * @param x initial value.
 * @param eps threshold for the residual.
 * @param max_iter maximum number of iterations.
 */

template< class MATRIX, class PC_MATRIX, class VECTOR > 
class Solver_preconditioned_CG 
{
public:
    typedef MATRIX Matrix ;
    typedef PC_MATRIX Preconditioner ;
    typedef VECTOR Vector ;
    typedef typename Vector::CoeffType CoeffType ;

public:
    Solver_preconditioned_CG() {
        epsilon_ = 1e-6 ;
        max_iter_ = 0 ;
    }

    // Default copy constructor, operator =() and destructor are fine

    void set_epsilon(CoeffType eps) { epsilon_ = eps ; }
    void set_max_iter(unsigned int max_iter) { max_iter_ = max_iter ; }

    // Solve the sparse linear system "A*x = b" for A symmetric positive definite
    // Return true on success
    bool solve(const MATRIX &A, const PC_MATRIX &C, const VECTOR& b, VECTOR& x) 
    {
#ifdef DEBUG_TRACE
        std::cerr << "  Call Conjugate Gradient with preconditioner" << std::endl;
#endif
        CGAL_assertion(A.dimension() > 0);
        unsigned int n = A.dimension() ;                        // (Square) matrix dimension

        unsigned int max_iter = max_iter_ ;                     // Max number of iterations
        if(max_iter == 0) {
            max_iter = 5 * n ;
        }

        Vector r(n) ;                                           // residue r=Ax-b
        Vector d(n) ;
        Vector h(n) ;
        Vector Ad = h ;
        unsigned int its=0;                                     // Loop counter
        CoeffType rh, alpha, beta;
        CoeffType bnorm2 = BLAS<Vector>::dot(b,b) ;
        CoeffType err=epsilon_*epsilon_*bnorm2 ;                // Error to reach
        mult(A,x,r);
        BLAS<Vector>::axpy(-1,b,r);
        mult(C,r,d);
        BLAS<Vector>::copy(d,h);
        rh=BLAS<Vector>::dot(r,h);

        CoeffType rr=BLAS<Vector>::dot(r,r);                    // error rr = (r|r)
        while ( rr>err && its < max_iter) {
            mult(A,d,Ad);
            CGAL_assertion(BLAS<Vector>::dot(d,Ad) != 0.0);
            alpha=rh/BLAS<Vector>::dot(d,Ad);
            BLAS<Vector>::axpy(-alpha,d,x);
            BLAS<Vector>::axpy(-alpha,Ad,r);
            mult(C,r,h);
            CGAL_assertion(rh != 0.0);
            beta=1/rh; rh=BLAS<Vector>::dot(r,h); beta*=rh;
            BLAS<Vector>::scal(beta,d);
            BLAS<Vector>::axpy(1,h,d);
            rr=BLAS<Vector>::dot(r,r);                          // Update error rr = (r|r)
            ++its;
        }

        bool success = (rr <= err);
        return success;
    }

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


} // namespace OpenNL

#endif // __OPENNL_CONJUGATE_GRADIENT__
