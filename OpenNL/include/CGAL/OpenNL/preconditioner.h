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
 *  Laurent Saboret 2006: Changes for CGAL:
 *      - copied Jacobi preconditioner from Graphite 1.9 code
 *      - Added OpenNL namespace
 */

#ifndef __OPENNL_PRECONDITIONER__
#define __OPENNL_PRECONDITIONER__

#include <CGAL/OpenNL/blas.h>
#include <CGAL/OpenNL/sparse_matrix.h>
#include <CGAL/OpenNL/full_vector.h>

#include <CGAL/assertions.h>
#include <CGAL/tss.h>

namespace OpenNL {


/**
 * Base class for some preconditioners.
 */

template <class T> 
class Preconditioner {
public:
    typedef T CoeffType ;

public:
    /**
     * The matrix A should be square.
     */
    Preconditioner(
        const SparseMatrix<T>& A, CoeffType omega = 1.0
    ) ;
    
    const SparseMatrix<T>& A() const { return A_ ; }
    CoeffType omega() const { return omega_ ; }

    /**
     * To use this function, the matrix A should be symmetric.
     */
    void mult_upper_inverse(const FullVector<T>& x, FullVector<T>& y) const ;

    /**
     * To use this function, the matrix A should be symmetric.
     */
    void mult_lower_inverse(const FullVector<T>& x, FullVector<T>& y) const ;

    void mult_diagonal(FullVector<T>& xy) const ;
    void mult_diagonal_inverse(FullVector<T>& xy) const ;
    
private:
    const SparseMatrix<T>& A_ ;
    CoeffType omega_ ;
} ;

template <class T> 
Preconditioner<T>::Preconditioner(
    const SparseMatrix<T>& A, CoeffType omega
) : A_(A), omega_(omega) {
    //CGAL_assertion(A.is_square()) ;
}

template <class T> 
void Preconditioner<T>::mult_lower_inverse(
    const FullVector<T>& x, FullVector<T>& y
) const {
    //CGAL_assertion(A_.has_symmetric_storage()) ;
    //CGAL_assertion(A_.rows_are_stored()) ;
    int n = A_.dimension() ;
    for(int i=0; i<n; i++) {
        double S = 0 ;
        const typename SparseMatrix<T>::Row& Ri = A_.row(i) ;
        for(int ij=0; ij < Ri.size(); ij++) {
            const typename SparseMatrix<T>::Coeff& c = Ri[ij] ;
            if (c.index < i) // traverse only lower half matrix
                S += c.a * y[c.index] ;
        }
        y[i] = (x[i] - S) * omega_ / A_.get_coef(i,i) ;
    }
}
template <class T> 
void Preconditioner<T>::mult_upper_inverse(
    const FullVector<T>& x, FullVector<T>& y
) const {
    //CGAL_assertion(A_.has_symmetric_storage()) ;
    //CGAL_assertion(A_.columns_are_stored()) ;
    int n = A_.dimension() ;
    for(int i=n-1; i>=0; i--) {
        double S = 0 ;
        const typename SparseMatrix<T>::Row& Ci = A_.row(i) ; // column i == row i
        for(int ij=0; ij < Ci.size(); ij++) {
            const typename SparseMatrix<T>::Coeff& c = Ci[ij] ;
            if(c.index > i) // traverse only upper half matrix
                S += c.a * y[c.index] ;
        }
        y[i] = (x[i] - S) * omega_ / A_.get_coef(i,i) ;
    }
}

template <class T> 
void Preconditioner<T>::mult_diagonal(FullVector<T>& xy) const {
    int n = A_.dimension() ;
    for(int i=0; i<n; i++) {
        xy[i] *= ( A_.get_coef(i,i) / omega_ ) ;
    }
}

template <class T> 
void Preconditioner<T>::mult_diagonal_inverse(FullVector<T>& xy) const {
    int n = A_.dimension() ;
    for(int i=0; i<n; i++) {
        xy[i] *= ( omega_ / A_.get_coef(i,i) ) ;
    }
}
    
/**
 * Jacobi preconditioner
 */

template <class T> 
class Jacobi_Preconditioner : public Preconditioner<T> {
public:
    typedef T CoeffType ;

public:
    Jacobi_Preconditioner(
        const SparseMatrix<T>& A, CoeffType omega = 1.0
    ) ;
} ;

template <class T> 
Jacobi_Preconditioner<T>::Jacobi_Preconditioner(
    const SparseMatrix<T>& A, CoeffType omega
) : Preconditioner<T>(A, omega) {
}

template <class T> 
void mult(const Jacobi_Preconditioner<T>& M, const FullVector<T>& x, FullVector<T>& y) {
    BLAS< FullVector<T> >::copy(x, y) ;
    M.mult_diagonal_inverse(y) ;
}


/**
 * The SSOR preconditioner, sharing storage with the matrix.
 */
 
template <class T> 
class SSOR_Preconditioner : public Preconditioner<T> {
public:
    typedef T CoeffType ;

public:
    /**
     * The matrix A should be symmetric.
     */
    SSOR_Preconditioner(
        const SparseMatrix<T>& A, CoeffType omega = 1.0
    ) ;
} ;

template <class T> 
SSOR_Preconditioner<T>::SSOR_Preconditioner(
    const SparseMatrix<T>& A, CoeffType omega
) : Preconditioner<T>(A, omega) {
}

/** y <- M*x */
template <class T> 
void mult(const SSOR_Preconditioner<T>& M, const FullVector<T>& x, FullVector<T>& y) {

    CGAL_STATIC_THREAD_LOCAL_VARIABLE(FullVector<T>, work,0) ;

    const SparseMatrix<T>& A = M.A() ;
    int n = A.dimension() ;

    if(work.dimension() != n) {
        work = FullVector<T>(n) ;
    }

    M.mult_lower_inverse(x, work) ;
    M.mult_diagonal(work) ;
    M.mult_upper_inverse(work, y) ;

    BLAS< FullVector<T> >::scal(2 - M.omega(), y) ;
}


} // namespace OpenNL

#endif // __OPENNL_PRECONDITIONER__
