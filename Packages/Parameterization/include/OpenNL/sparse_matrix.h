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
 *		- SparseMatrix is now a model of the SparseLinearAlgebraTraits_d::Matrix concept
*/

#ifndef ___SPARSE_MATRIX__
#define ___SPARSE_MATRIX__

#include "full_vector.h"
#include <vector>
#include <cstdlib>
#include <cassert>

namespace OpenNL {


//________________________________________________________________
// Class SparseMatrix
// Model of the SparseLinearAlgebraTraits_d::Matrix concept
template <class T> class SparseMatrix {
public:

    typedef T CoeffType ;
	typedef T NT;			// for SparseLinearAlgebraTraits_d::Matrix concept

    struct Coeff {
        Coeff() { }
        Coeff(unsigned int i, T val) : index(i), a(val) { }
        unsigned int index ;
        T a ;
    } ;
    
//__________________________________________________

   /**
    * A row or a column of a SparseMatrix. The row/column is
    * compressed, and stored in the form of a list of
    * (value,index) couples.
    */
    class Row : public std::vector<Coeff> {
        typedef typename std::vector<Coeff> superclass ;
    public:
        /** a_{index} <- a_{index} + val */
        void add(unsigned int index, T val) 
		{
			// search for coefficient in superclass vector
            for(typename superclass::iterator it = superclass::begin() ; 
                it != superclass::end() ; 
				it++)
            {
                if(it->index == index) {
                    it->a += val ;						// +=
                    return ;
				}
            }
            superclass::push_back(Coeff(index, val)) ;	// coefficient doesn't exist yet if we reach this point
        }

        // a_{index} <- val
		// (added for SparseLinearAlgebraTraits_d::Matrix concept)
        void set_coef(unsigned int index, T val) 
		{
			// search for coefficient in superclass vector
            for(typename superclass::iterator it = superclass::begin() ; 
                it != superclass::end() ; 
				it++)
            {
                if(it->index == index) {
                    it->a = val ;						// =
                    return ;
				}
            }
            superclass::push_back(Coeff(index, val)) ;	// coefficient doesn't exist yet if we reach this point
        }

        // return a_{index} (0 by default)
		// (added for SparseLinearAlgebraTraits_d::Matrix concept)
        T get_coef(unsigned int index) const
		{
			// search for coefficient in superclass vector
            for(typename superclass::const_iterator it = superclass::begin() ; 
                it != superclass::end() ; 
				it++)
            {
                if(it->index == index)
                    return it->a ;						// return value
            }
            return 0 ;									// coefficient doesn't exist if we reach this point
        }
    } ;

    //__________ constructors / destructor _____

	// Create a square matrix initialized with zeros
    SparseMatrix(unsigned int dim) {
		assert(dim > 0);
		dimension_ = dim ;
        row_ = new Row[dimension_] ;
    }
	// Create a matrix initialized with zeros (for SparseLinearAlgebraTraits_d::Matrix concept)
	// WARNING: this class supports square matrices only
	SparseMatrix (unsigned int rows, unsigned int columns ) {
		assert(rows == columns);
		assert(columns > 0);
		dimension_ = columns ;
		row_ = new Row[dimension_] ;
    }
    
    ~SparseMatrix() {
        delete[] row_ ;
        row_ = NULL ;
    }
    
    //___________ access ________________________
    
	// Return the matrix dimension
    unsigned int dimension() const {  return dimension_ ;  }
	// For SparseLinearAlgebraTraits_d::Matrix concept:
	unsigned int row_dimension() const    { return dimension(); }
	unsigned int column_dimension() const { return dimension(); }

    Row& row(unsigned int i) {
        assert(i < dimension_) ;
        return row_[i] ;
    }

    const Row& row(unsigned int i) const {
        assert(i < dimension_) ;
        return row_[i] ;
    }

	// Write access to 1 matrix coefficient: aij <- aij + val
    void add(unsigned int i, unsigned int j, T val) {
        assert(i < dimension_) ;
        assert(j < dimension_) ;
        row(i).add(j, val) ;
    }

	// Write access to 1 matrix coefficient: aij <- val
	//(for SparseLinearAlgebraTraits_d::Matrix concept)
    void set_coef(unsigned int i, unsigned int j, NT  val) {
        assert(i < dimension_) ;
        assert(j < dimension_) ;
        row(i).set_coef(j, val) ;
    }

	// Read access to 1 matrix coefficient 
	// (for SparseLinearAlgebraTraits_d::Matrix concept)
	NT  get_coef (unsigned int i, unsigned int j) const {
        assert(i < dimension_) ;
        assert(j < dimension_) ;
        return row(i).get_coef(j) ;
	}
            
    /**
     * removes all the coefficients and frees the allocated
     * space.
     */
    void clear() {
        for(unsigned int i=0; i<dimension_; i++) {
            row(i).clear() ;
        }
    }

private:
    unsigned int dimension_ ;
    Row* row_ ;
    
    // SparseMatrix cannot be copied (for the moment, could be implemented if needed).
    SparseMatrix(const SparseMatrix& rhs) ;
    SparseMatrix& operator=(const SparseMatrix& rhs) ;
} ;

template <class T> void mult(const SparseMatrix<T>& M, const FullVector<T>& x, FullVector<T>& y) {
    unsigned int N = M.dimension() ;
    assert(x.dimension() == N) ;
    assert(y.dimension() == N) ;
    for(unsigned int i=0; i<N; i++) {
        y[i] = 0 ;
        const typename SparseMatrix<T>::Row& R = M.row(i) ;
        for(unsigned int jj=0; jj<R.size(); jj++) {
            unsigned int j = R[jj].index ;
            y[i] += R[jj].a * x[j] ;
        }
    }
}


}; // namespace OpenNL

#endif
