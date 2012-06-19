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
 *  Laurent Saboret 01/2005: Change for CGAL:
 *      - Added OpenNL namespace
 *      - FullVector is now a model of the SparseLinearAlgebraTraits_d::Vector concept
 *      - Coefficients are initialized with zeros
 */


#ifndef __OPENNL_FULL_VECTOR__
#define __OPENNL_FULL_VECTOR__

#include <CGAL/OpenNL/blas.h>
#include <CGAL/assertions.h>

#include <cstdlib>

namespace OpenNL {


// Class FullVector
// Model of the SparseLinearAlgebraTraits_d::Vector concept
template <class T> class FullVector
{
// Public types
public:
    typedef T CoeffType ;

    // for SparseLinearAlgebraTraits_d::Vector concept
    typedef T NT;

    // Create a vector initialized with zeros
    FullVector(unsigned int dim) {
        dimension_ = dim ; coeff_ = new T[dim] ;

        // init with zeros (for SparseLinearAlgebraTraits_d::Vector concept)
        for (unsigned int i=0; i < dimension_; i++)
            coeff_[i] = 0;
    }
    // Copy constructor
    FullVector(const FullVector& toCopy) {
        dimension_ = toCopy.dimension_ ;
        coeff_ = new T[dimension_] ;
        for (unsigned int i=0; i < dimension_; i++)
            coeff_[i] = toCopy.coeff_[i];
    }
    // operator =()
    FullVector& operator =(const FullVector& toCopy) {
        delete[] coeff_ ;

        dimension_ = toCopy.dimension_ ;
        coeff_ = new T[dimension_] ;
        for (unsigned int i=0; i < dimension_; i++)
            coeff_[i] = toCopy.coeff_[i];

        return *this;
    }

    ~FullVector() {
        delete[] coeff_ ;
        coeff_ = NULL ;
    }

    // Return the vector's number of coefficients
    unsigned int dimension() const {
        return dimension_ ;
    }
    // Read/write access to 1 vector coefficient.
    //
    // Preconditions:
    // * 0 <= i < dimension()
    const T& operator[](unsigned int i) const {
        CGAL_assertion(i < dimension_) ;
        return coeff_[i] ;
    }
    T& operator[](unsigned int i) {
        CGAL_assertion(i < dimension_) ;
        return coeff_[i] ;
    }

private:
    unsigned int dimension_ ;
    T* coeff_ ;
} ;

template <class T> class BLAS< FullVector<T> > {
public:
    typedef FullVector<T> VectorType ;
    typedef T CoeffType ;

    /** y <- y + a*x  */
    static void axpy(CoeffType a, const VectorType& x, VectorType& y) {
        CGAL_assertion(x.dimension() == y.dimension()) ;
        for(unsigned int i=0; i<x.dimension(); i++) {
            y[i] += a * x[i] ;
        }
    }

    /** x <- a*x */
    static void scal(CoeffType a, VectorType& x) {
        for(unsigned int i=0; i<x.dimension(); i++) {
            x[i] *= a ;
        }
    }

    /** y <- x */
    static void copy(const VectorType& x, VectorType& y) {
        CGAL_assertion(x.dimension() == y.dimension()) ;
        for(unsigned int i=0; i<x.dimension(); i++) {
            y[i] = x[i] ;
        }
    }

    /** returns x^t * y */
    static CoeffType dot(const VectorType& x, const VectorType& y) {
        CoeffType result = 0 ;
        CGAL_assertion(x.dimension() == y.dimension()) ;
        for(unsigned int i=0; i<x.dimension(); i++) {
            result += y[i] * x[i] ;
        }
        return result ;
    }
} ;


} // namespace OpenNL

#endif // __OPENNL_FULL_VECTOR__
