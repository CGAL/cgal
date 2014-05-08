// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_MIX_VECTOR_H
#define CGAL_KD_MIX_VECTOR_H
#include <CGAL/Dimension.h>
namespace CGAL {

template <class Static_, class Dynamic_, class NT_ ,class Dim_, class Max_dim_ = Dim_>
struct Mix_vector
: Dynamic_::template Rebind_dimension<Dim_, Max_dim_>::Other
{
  template <class D2, class D3 = D2>
  struct Rebind_dimension {
    typedef Mix_vector<Static_, Dynamic_, NT_, D2, D3> Other;
  };
};

template <class Static_, class Dynamic_, class NT_, int d, class Max_dim_>
struct Mix_vector<Static_, Dynamic_, NT_, Dimension_tag<d>, Max_dim_>
: Static_::template Rebind_dimension<Dimension_tag<d>, Max_dim_>::Other
{
  template <class D2, class D3 = D2>
  struct Rebind_dimension {
    typedef Mix_vector<Static_, Dynamic_, NT_, D2, D3> Other;
  };
};
}
#endif

