// Copyright (c) 2012 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: 
// $Id: 
// 
//

#ifndef CGAL_RANDOM_SHUFFLE_H
#define CGAL_RANDOM_SHUFFLE_H

// needed for iter_swap
#include <algorithm>
#include <CGAL/config.h>

namespace CGAL {

// a random_shuffle that does not vary across different stdlib
// implementation and thus reliably reproduces the same result when
// called with the same seed.

// only the 3 argument version is provided.

// this is implementation found in the gcc 4.6.2 release

template<typename RandomAccessIterator, typename RandomNumberGenerator>
void
random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
#if !defined CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
               RandomNumberGenerator&& rand)
#else
               RandomNumberGenerator& rand)
#endif
{
  if (first == last) return;
  for (RandomAccessIterator i = first + 1; i != last; ++i)
    std::iter_swap(i, first + rand((i - first) + 1));
}

} // CGAL


#endif /* CGAL_RANDOM_SHUFFLE_H */
