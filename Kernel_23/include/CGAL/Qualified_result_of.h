// Copyright (c) 1999,2003  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion
 

#ifndef CGAL_QUALIFIED_RESULT_OF_H
#define CGAL_QUALIFIED_RESULT_OF_H

#include <boost/type_traits.hpp>

namespace CGAL {

struct Has_qrt {};

template < typename F, bool b >
struct qrt_or_not
{ typedef typename F::result_type type; };

template < typename F >
struct qrt_or_not <F, true>
{ typedef const typename F::result_type & type; };

template <typename Functor, typename P1 = void, typename P2 = void,
                            typename P3 = void, typename P4 = void,
                            typename P5 = void, typename P6 = void,
                            typename P7 = void, typename P8 = void>
struct Qualified_result_of
  : qrt_or_not<Functor, boost::is_base_and_derived<Has_qrt, Functor>::value>
{};

} //namespace CGAL

#endif // CGAL_QUALIFIED_RESULT_OF_H
