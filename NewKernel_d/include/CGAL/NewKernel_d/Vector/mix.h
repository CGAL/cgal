// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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

// FIXME: shouldn't we dispatch based on Max_dim_ instead?
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

