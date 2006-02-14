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
 */

#ifndef __OPENNL_BLAS__
#define __OPENNL_BLAS__

#include <cassert>

namespace OpenNL {


/** Basic Linear Algebra Subroutines */
template <class VECTOR> class BLAS {
public:
    typedef VECTOR VectorType ;
    typedef typename VECTOR::CoeffType CoeffType ;
    
    /** y <- y + a*x  */
    static void axpy(CoeffType a, const VectorType& x, VectorType& y) {
        assert(false) ;
    }

    /** x <- a*x */
    static void scal(CoeffType a, VectorType& x) {
        assert(false) ;
    }

    /** y = x */
    static void copy(const VectorType& x, VectorType& y) {
        assert(false) ;
    }

    /** returns x^t * y */
    static CoeffType dot(const VectorType& x, const VectorType& y) {
        assert(false) ;
    }
} ;


}; // namespace OpenNL

#endif
