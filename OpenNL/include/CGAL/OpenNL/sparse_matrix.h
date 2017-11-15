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
 *      - SparseMatrix is now a model of the SparseLinearAlgebraTraits_d::Matrix concept
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0+
*/

#ifndef __OPENNL_SPARSE_MATRIX__
#define __OPENNL_SPARSE_MATRIX__

#include <CGAL/OpenNL/full_vector.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

#include <vector>
#include <cstdlib>

namespace OpenNL {


//________________________________________________________________
// Class SparseMatrix
// Model of the SparseLinearAlgebraTraits_d::Matrix concept
template <class T> class SparseMatrix
{
// Public types
public:

    typedef T CoeffType ;

    // added for SparseLinearAlgebraTraits_d::Matrix concept
    typedef T NT;

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
        void add_coef(unsigned int index, T val)
        {
            // search for coefficient in superclass vector
            for(typename superclass::iterator it = superclass::begin() ;
                it != superclass::end() ;
                it++)
            {
                if(it->index == index) {
                    it->a += val ;                      // +=
                    return ;
                }
            }
            // coefficient doesn't exist yet if we reach this point
            superclass::push_back(Coeff(index, val)) ;
        }

        // a_{index} <- val
        // (added for SparseLinearAlgebraTraits_d::Matrix concept)
        //
        // Optimization:
        // - Caller can optimize this call by setting 'new_coef' to true
        //   if the coefficient does not already exists in the matrix. 
        void set_coef(unsigned int index, T val, bool new_coef)
        {
            if (!new_coef)
            {
              // search for coefficient in superclass vector
              for(typename superclass::iterator it = superclass::begin() ;
                  it != superclass::end() ;
                  it++)
              {
                  if(it->index == index) {
                      it->a = val ;                       // =
                      return ;
                  }
              }
            }
            // coefficient doesn't exist yet if we reach this point
            superclass::push_back(Coeff(index, val)) ;
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
                    return it->a ;                      // return value
            }
            // coefficient doesn't exist if we reach this point
            return 0 ;
        }
    } ;

// Public operations
public:

    //__________ constructors / destructor _____

    // Create a square matrix initialized with zeros
    SparseMatrix(unsigned int dim) {
        CGAL_assertion(dim > 0);
        dimension_ = dim ;
        row_ = new Row[dimension_] ;
    }
    // Create a rectangular matrix initialized with zeros
    // (added for SparseLinearAlgebraTraits_d::Matrix concept)
    // WARNING: this class supports square matrices only
    SparseMatrix (unsigned int rows, unsigned int columns ) {
        CGAL_USE(rows);
        CGAL_assertion(rows == columns);
        CGAL_assertion(columns > 0);
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

    // Added for SparseLinearAlgebraTraits_d::Matrix concept:
    // Return the matrix number of rows
    unsigned int row_dimension() const    { return dimension(); }
    // Return the matrix number of columns
    unsigned int column_dimension() const { return dimension(); }

    Row& row(unsigned int i) {
        CGAL_assertion(i < dimension_) ;
        return row_[i] ;
    }

    const Row& row(unsigned int i) const {
        CGAL_assertion(i < dimension_) ;
        return row_[i] ;
    }

    // Read access to 1 matrix coefficient
    // (added for SparseLinearAlgebraTraits_d::Matrix concept)
    //
    // Preconditions:
    // * 0 <= i < row_dimension()
    // * 0 <= j < column_dimension()
    NT  get_coef (unsigned int i, unsigned int j) const {
        CGAL_assertion(i < dimension_) ;
        CGAL_assertion(j < dimension_) ;
        return row(i).get_coef(j) ;
    }

    // Write access to 1 matrix coefficient: a_ij <- a_ij + val
    //
    // Preconditions:
    // * 0 <= i < row_dimension()
    // * 0 <= j < column_dimension()
    void add_coef(unsigned int i, unsigned int j, T val) {
        CGAL_assertion(i < dimension_) ;
        CGAL_assertion(j < dimension_) ;
        row(i).add_coef(j, val) ;
    }

    // Write access to 1 matrix coefficient: a_ij <- val
    //(added for SparseLinearAlgebraTraits_d::Matrix concept)
    //
    // Optimization:
    // - Caller can optimize this call by setting 'new_coef' to true
    //   if the coefficient does not already exists in the matrix. 
    //
    // Preconditions:
    // - 0 <= i < row_dimension().
    // - 0 <= j < column_dimension().
    void set_coef(unsigned int i, unsigned int j, NT  val, bool new_coef = false) {
        CGAL_assertion(i < dimension_) ;
        CGAL_assertion(j < dimension_) ;
        row(i).set_coef(j, val, new_coef) ;
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

    // SparseMatrix cannot be copied
    // (for the moment, could be implemented if needed).
    SparseMatrix(const SparseMatrix& rhs) ;
    SparseMatrix& operator=(const SparseMatrix& rhs) ;
} ;

/** y <- M*x */
template <class T> 
void mult(const SparseMatrix<T>& M, const FullVector<T>& x, FullVector<T>& y) {
    unsigned int N = M.dimension() ;
    CGAL_assertion(x.dimension() == N) ;
    CGAL_assertion(y.dimension() == N) ;
    for(unsigned int i=0; i<N; i++) {
        y[i] = 0 ;
        const typename SparseMatrix<T>::Row& R = M.row(i) ;
        for(unsigned int jj=0; jj<R.size(); jj++) {
            unsigned int j = R[jj].index ;
            y[i] += R[jj].a * x[j] ;
        }
    }
}


} // namespace OpenNL

#endif // __OPENNL_SPARSE_MATRIX__
