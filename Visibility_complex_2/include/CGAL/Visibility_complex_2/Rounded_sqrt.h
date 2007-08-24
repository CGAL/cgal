// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_ROUNDED_SQRT_H
#define CGAL_VISIBILITY_COMPLEX_2_ROUNDED_SQRT_H

#include <CGAL/basic.h>
#include <CGAL/Algebraic_structure_traits.h>

CGAL_BEGIN_NAMESPACE

namespace Visibility_complex_2_details {




template <class FT,class Sqrt> struct Rounded_sqrt_aux {
  FT operator()(const FT&x) {
    return Sqrt()(x);
  }
};


template <class FT> struct Rounded_sqrt_aux<FT,Null_functor> {
  FT operator() (const FT& x) {
    typename Algebraic_structure_traits<double>::Sqrt sqrt;
    return static_cast<FT>(sqrt(CGAL_NTS to_double(x)));
  }
};

template <class FT> struct Rounded_sqrt {
  typedef
  Rounded_sqrt_aux<FT,
                   typename Algebraic_structure_traits<FT>::Sqrt>
  Sqrt;
};

  typedef int ploum;

}
CGAL_END_NAMESPACE
#endif // CGAL_VISIBILITY_COMPLEX_2_ROUNDED_SQRT_H
