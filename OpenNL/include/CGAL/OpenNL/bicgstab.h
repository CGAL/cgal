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
 *  Laurent Saboret 2005-2008: Changes for CGAL:
 *      - Added OpenNL namespace
 *      - solve() returns true on success
 *      - check divisions by zero
 *      - added comments and traces
 *      - copied BICGSTAB algorithm WITH preconditioner from Graphite 1.9 code
 */

#ifndef __OPENNL_BICGSTAB__
#define __OPENNL_BICGSTAB__

#include <CGAL/OpenNL/blas.h>
#include <CGAL/assertions.h>

#include <cmath>
#include <cfloat>
#include <climits>

namespace OpenNL {


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

template <class MATRIX, class VECTOR> class Solver_BICGSTAB {
public:
    typedef MATRIX Matrix ;
    typedef VECTOR Vector ;
    typedef typename Vector::CoeffType CoeffType ;

public:
    Solver_BICGSTAB() {
        epsilon_ = 1e-6 ;
        max_iter_ = 0 ;
    }

    // Default copy constructor, operator =() and destructor are fine

    void set_epsilon(CoeffType eps) { epsilon_ = eps ; }
    void set_max_iter(unsigned int max_iter) { max_iter_ = max_iter ; }

    // Solve the sparse linear system "A*x = b". Return true on success.
    bool solve(const MATRIX &A, const VECTOR& b, VECTOR& x) 
    {
#ifdef DEBUG_TRACE
        std::cerr << "  Call BICGSTAB" << std::endl;
#endif
        CGAL_assertion(A.dimension() > 0);
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
        Vector r(n) ;                                           // residue r=A*x-b
        mult(A,x,r);
        BLAS<Vector>::axpy(-1,b,r);
        BLAS<Vector>::copy(r,d);
        BLAS<Vector>::copy(d,h);
        BLAS<Vector>::copy(h,rT);
        //CGAL_assertion(BLAS<Vector>::dot(rT,rT) == 0.0);      // may happen for small systems
        rTh=BLAS<Vector>::dot(rT,h);                            // rTh = (rT|h)
        rTr=BLAS<Vector>::dot(r,r);                             // error rTr = (r|r)

        while ( rTr>err && its < max_iter) {
            mult(A,d,Ad);
            rTAd=BLAS<Vector>::dot(rT,Ad);
            CGAL_assertion(rTAd != 0.0);
            alpha=rTh/rTAd;
            BLAS<Vector>::axpy(-alpha,Ad,r);
            BLAS<Vector>::copy(h,s);
            BLAS<Vector>::axpy(-alpha,Ad,s);
            mult(A,s,t);
            BLAS<Vector>::axpy(1,t,u);
            BLAS<Vector>::scal(alpha,u);
            st=BLAS<Vector>::dot(s,t);
            tt=BLAS<Vector>::dot(t,t);
            if (st == 0.0 || tt == 0.0)
                omega = 0 ;
            else
                omega = st/tt;
            BLAS<Vector>::axpy(-omega,t,r);
            BLAS<Vector>::axpy(-alpha,d,x);
            BLAS<Vector>::axpy(-omega,s,x);
            rTr=BLAS<Vector>::dot(r,r);
            BLAS<Vector>::copy(s,h);
            BLAS<Vector>::axpy(-omega,t,h);
            if (omega == 0.0)                                   // stop if omega==0 (success)
            {
#ifdef DEBUG_TRACE
                std::cerr << "  BICGSTAB: omega=0" << std::endl;
#endif
                break;                                          
            }
            if (rTh == 0.0)                                     // stop if rTh==0 (failure)
            {
#ifdef DEBUG_TRACE
                std::cerr << "  BICGSTAB: rTh=0" << std::endl;
#endif
                break;                                          
            }
            beta=(alpha/omega)/rTh;                             // beta = (rTh/"old rTh") * (alpha/omega)
            rTh=BLAS<Vector>::dot(rT,h); 
            beta*=rTh;
            BLAS<Vector>::scal(beta,d);
            BLAS<Vector>::axpy(1,h,d);
            BLAS<Vector>::axpy(-beta*omega,Ad,d);
            ++its;
        }

        bool success = (rTr <= err);
        return success;
    }

private:
    // Test if a floating point number is (close to) 0.0
    static bool IsZero(CoeffType a)
    {
        return (std::fabs(a) < 10.0 * (std::numeric_limits<CoeffType>::min)());
    }

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


/**
 *  The BICGSTAB algorithm WITH preconditioner:
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
class Solver_preconditioned_BICGSTAB 
{
public:
    typedef MATRIX Matrix ;
    typedef PC_MATRIX Preconditioner ;
    typedef VECTOR Vector ;
    typedef typename Vector::CoeffType CoeffType ;

public:
    Solver_preconditioned_BICGSTAB() {
        epsilon_ = 1e-6 ;
        max_iter_ = 0 ;
    }

    // Default copy constructor, operator =() and destructor are fine

    void set_epsilon(CoeffType eps) { epsilon_ = eps ; }
    void set_max_iter(unsigned int max_iter) { max_iter_ = max_iter ; }

    // Solve the sparse linear system "A*x = b". Return true on success.
    bool solve(const MATRIX &A, const PC_MATRIX &C, const VECTOR& b, VECTOR& x) 
    {
#ifdef DEBUG_TRACE
        std::cerr << "  Call BICGSTAB with preconditioner" << std::endl;
#endif
        CGAL_assertion(A.dimension() > 0);
        unsigned int n = A.dimension() ;                        // (Square) matrix dimension

        unsigned int max_iter = max_iter_ ;                     // Max number of iterations
        if(max_iter == 0) {
            max_iter = 5 * n ;
        }

        Vector rT(n) ;                                          // Initial residue rT=Ax-b
        Vector d(n) ;
        Vector h(n) ;
        Vector u(n) ;
        Vector Sd(n) ;
        Vector t(n) ;
        Vector aux(n) ;
        Vector& s = h ;
        CoeffType rTh, rTSd, rTr, alpha, beta, omega, st, tt;
        unsigned int its=0;                                     // Loop counter
        CoeffType err=epsilon_*epsilon_*BLAS<Vector>::dot(b,b); // Error to reach
        Vector r(n) ;                                           // residue r=A*x-b
        mult(A,x,r);
        BLAS<Vector>::axpy(-1,b,r);
        mult(C,r,d);
        BLAS<Vector>::copy(d,h);
        BLAS<Vector>::copy(h,rT);
        //CGAL_assertion(BLAS<Vector>::dot(rT,rT) == 0.0);      // may happen for small systems
        rTh=BLAS<Vector>::dot(rT,h);                            // rTh = (rT|h)
        rTr=BLAS<Vector>::dot(r,r);                             // error rTr = (r|r)

        while (rTr>err && its < max_iter) {
            mult(A,d,aux);
            mult(C,aux,Sd);
            rTSd=BLAS<Vector>::dot(rT,Sd);
            if (rTSd == 0.0)                                     // stop if rTSd==0 (failure)
            {
#ifdef DEBUG_TRACE
                std::cerr << "  BICGSTAB with preconditioner: rTSd=0" << std::endl;
#endif
                break;                                          
            }
            alpha=rTh/rTSd;
            BLAS<Vector>::axpy(-alpha,aux,r);
            BLAS<Vector>::copy(h,s);
            BLAS<Vector>::axpy(-alpha,Sd,s);
            mult(A,s,aux);
            mult(C,aux,t);
            BLAS<Vector>::axpy(1,t,u);
            BLAS<Vector>::scal(alpha,u);
            st=BLAS<Vector>::dot(s,t);
            tt=BLAS<Vector>::dot(t,t);
            if (st == 0.0 || tt == 0.0)
                omega = 0 ;
            else
                omega = st/tt;
            BLAS<Vector>::axpy(-omega,aux,r);
            BLAS<Vector>::axpy(-alpha,d,x);
            BLAS<Vector>::axpy(-omega,s,x);
            rTr=BLAS<Vector>::dot(r,r);
            BLAS<Vector>::copy(s,h);
            BLAS<Vector>::axpy(-omega,t,h);
            if (omega == 0.0)                                   // stop if omega==0 (success)
            {
#ifdef DEBUG_TRACE
                std::cerr << "  BICGSTAB with preconditioner: omega=0" << std::endl;
#endif
                break;                                          
            }
            if (rTh == 0.0)                                     // stop if rTh==0 (failure)
            {
#ifdef DEBUG_TRACE
                std::cerr << "  BICGSTAB with preconditioner: rTh=0" << std::endl;
#endif
                break;                                          
            }
            beta=(alpha/omega)/rTh;                             // beta = (rTh/"old rTh") * (alpha/omega)
            rTh=BLAS<Vector>::dot(rT,h); 
            beta*=rTh;
            BLAS<Vector>::scal(beta,d);                         // d = h + beta * (d - omega * Sd);
            BLAS<Vector>::axpy(1,h,d);
            BLAS<Vector>::axpy(-beta*omega,Sd,d);
            ++its;
        }

        bool success = (rTr <= err);
        return success;
    }

private:
    // Test if a floating point number is (close to) 0.0
    static bool IsZero(CoeffType a)
    {
        return (std::fabs(a) < 10.0 * (std::numeric_limits<CoeffType>::min)());
    }

private:
    CoeffType epsilon_ ;
    unsigned int max_iter_ ;
} ;


} // namespace OpenNL

#endif // __OPENNL_BICGSTAB__
