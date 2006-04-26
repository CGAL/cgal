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
 *  Laurent Saboret 01/2005: Change for CGAL:
 *      - Added OpenNL namespace
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
