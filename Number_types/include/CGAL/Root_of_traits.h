// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France)
// All rights reserved.
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
// Author(s)     : Sylvain Pion, Monique Teillaud, Athanasios Kakargias

#ifndef CGAL_ROOT_OF_TRAITS_H
#define CGAL_ROOT_OF_TRAITS_H

#include <CGAL/number_type_basic.h>
#include <CGAL/Root_of_2_fwd.h>
#include <CGAL/Quotient_fwd.h>

namespace CGAL {

namespace CGALi {

// Dispatcher for the case Has_sqrt==Tag_true or not.
template < typename RT,
           typename Has_sqrt = typename Number_type_traits<RT>::Has_sqrt >
struct Root_of_traits_helper
{
  typedef Quotient< RT >   RootOf_1;
  typedef Root_of_2< RT >  RootOf_2;
};

// Specialization for Has_sqrt==Tag_true.
template < typename RT >
struct Root_of_traits_helper < RT, Tag_true >
{
  typedef RT  RootOf_1;
  typedef RT  RootOf_2;
};

} // namespace CGALi


// Default Traits class for RT types
template < typename RT >
struct Root_of_traits
  : public CGALi::Root_of_traits_helper<RT> {};

} // namespace CGAL

#endif // CGAL_ROOT_OF_TRAITS_H
