// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_MIXED_COMPLEX_FILTERED_TRAITS_3_H
#define CGAL_MIXED_COMPLEX_FILTERED_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Mixed_complex_traits_3.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Filtered_predicate.h>

CGAL_BEGIN_NAMESPACE 

template <class K>
class Mixed_complex_filtered_traits_3
  : public Mixed_complex_traits_base_3<K>
{
  // Exact traits is based on the exact kernel.
  typedef Mixed_complex_traits_3<typename K::EK>
                                                   Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Mixed_complex_traits_3<typename K::FK>
                                                   Filtering_traits;

  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

  typedef Mixed_complex_traits_base_3<K>           Base;
public:
  typedef Filtered_predicate<
            typename Exact_traits::Side_of_mixed_cell_3,
            typename Filtering_traits::Side_of_mixed_cell_3,
            Weighted_converter_3<C2E>,
            Weighted_converter_3<C2F> >  Side_of_mixed_cell_3;

  Mixed_complex_filtered_traits_3() {}
  Mixed_complex_filtered_traits_3(typename Base::FT s) : Base(s) {}

  // Only make the predicates filtered, not the constructions:
  Side_of_mixed_cell_3 
  side_of_mixed_cell_3_object() const 
  { 
    return Side_of_mixed_cell_3(Base::get_shrink()); 
  }

};

CGAL_END_NAMESPACE
#endif // CGAL_MIXED_COMPLEX_FILTERED_TRAITS_3_H
