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
 *		- FullVector is now a model of the SparseLinearAlgebraTraits_d::Vector concept
 *		- Coefficients are initialized with zeros
 */


#ifndef __FULL_VECTOR__
#define __FULL_VECTOR__

#include "blas.h"
#include <cstdlib>
#include <cassert>

namespace OpenNL {


// Class FullVector
// Model of the SparseLinearAlgebraTraits_d::Vector concept
template <class T> class FullVector {
public:
    typedef T CoeffType ;
	typedef T NT;								// for SparseLinearAlgebraTraits_d::Vector concept

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
	// * 0 <= row < dimension()
    const T& operator[](unsigned int i) const {
        assert(i < dimension_) ;
        return coeff_[i] ;
    }
    T& operator[](unsigned int i) {
        assert(i < dimension_) ;
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
        assert(x.dimension() == y.dimension()) ;
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
        assert(x.dimension() == y.dimension()) ;
        for(unsigned int i=0; i<x.dimension(); i++) {
            y[i] = x[i] ;
        }
    }

    /** returns x^t * y */
    static CoeffType dot(const VectorType& x, const VectorType& y) {
        CoeffType result = 0 ;
        assert(x.dimension() == y.dimension()) ;
        for(unsigned int i=0; i<x.dimension(); i++) {
            result += y[i] * x[i] ;
        }
        return result ;
    }

} ;


}; // namespace OpenNL

#endif
