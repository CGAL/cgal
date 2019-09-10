// Copyright (c) 2011  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Olivier Devillers

#ifndef CGAL_HILBERT_SORT_2_H
#define CGAL_HILBERT_SORT_2_H

#include <CGAL/Hilbert_policy_tags.h>
#include <CGAL/Hilbert_sort_median_2.h>
#include <CGAL/Hilbert_sort_middle_2.h>

namespace CGAL {

template <class K,  class Hilbert_policy >
  class Hilbert_sort_2;

template <class K>  
  class Hilbert_sort_2<K, Hilbert_sort_median_policy >
  : public Hilbert_sort_median_2<K>
{
 public:
  Hilbert_sort_2 (const K &k=K() , std::ptrdiff_t limit=1 )
   : Hilbert_sort_median_2<K> (k,limit)
    {}
};

template <class K>
  class Hilbert_sort_2<K, Hilbert_sort_middle_policy >
  : public Hilbert_sort_middle_2<K>
{
 public:
 Hilbert_sort_2 (const K &k=K() , std::ptrdiff_t limit=1 )
   : Hilbert_sort_middle_2<K> (k,limit)
    {}
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_2_H
