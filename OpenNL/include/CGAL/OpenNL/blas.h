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
 *  Andreas Meyer 2007 changes for CGAL:
 *      - replaced assert with CGAL_assertion/CGAL_error etc.
*/

#ifndef __OPENNL_BLAS__
#define __OPENNL_BLAS__

#include <CGAL/assertions.h>

namespace OpenNL {


/** Basic Linear Algebra Subroutines */
template <class VECTOR> class BLAS {
public:
    typedef VECTOR VectorType ;
    typedef typename VECTOR::CoeffType CoeffType ;

    /** y <- y + a*x  */
    static void axpy(CoeffType /*a*/, const VectorType& /*x*/, VectorType& /*y*/) {
        CGAL_error();
    }

    /** x <- a*x */
    static void scal(CoeffType /*a*/, VectorType& /*x*/) {
        CGAL_error();
    }

    /** y <- x */
    static void copy(const VectorType& /*x*/, VectorType& /*y*/) {
        CGAL_error();
    }

    /** returns x^t * y */
    static CoeffType dot(const VectorType& /*x*/, const VectorType& /*y*/) {
        CGAL_error();
    }
} ;


} // namespace OpenNL

#endif
