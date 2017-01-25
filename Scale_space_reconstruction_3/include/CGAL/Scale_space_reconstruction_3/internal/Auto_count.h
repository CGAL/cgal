//A unary operator to that counts instances.
//Copyright (C) 2013  INRIA - Sophia Antipolis
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s):      Thijs van Lankveld


#ifndef CGAL_INTERNAL_AUTO_COUNT_H
#define CGAL_INTERNAL_AUTO_COUNT_H

#include <CGAL/license/Scale_space_reconstruction_3.h>


#include <functional>
#include <utility>


namespace CGAL {

namespace internal {

//  constructs a pair containing the object and the number of these pairs previously constructed.
/*  \ingroup PkgScaleSpaceReconstruction3Auxiliary
 *  \tparam T is the type of object to count.
 *  \tparam C is the number type of the counter.
 *  \todo Make thread-safe.
 */
template < class T, class C = unsigned int >
class Auto_count
: public std::unary_function< const T&, std::pair< T, C > > {
    mutable C i; // Note, not thread-safe.
public:
/// \name Constructors
/// \{
    /// starts a new count.
    Auto_count(): i(0) {}
/// \}

/// \name Operations
/// \{
    /// constructs a pair with the object and the number of pairs previously constructed.
    /** \param t is the current object.
     *  \return a pair containing the object and the number of pairs previously constructed.
     */
    std::pair< T, C > operator()( const T& t ) const { return std::make_pair( t, i++ ); }
/// \}
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_INTERNAL_AUTO_COUNT_H
