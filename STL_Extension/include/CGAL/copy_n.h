// Copyright (c) 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>


// This file is obsolete and exists only for backwards-compatibility.
// Include <CGAL/algorithm.h> instead.

#ifndef CGAL_COPY_N_H
#define CGAL_COPY_N_H 1
#include <cstddef>

// copy_n is usually in the STL as well, but not in the official
// standard. We provide our own copy_n.

CGAL_BEGIN_NAMESPACE

template <class InputIterator, class Size, class OutputIterator>
OutputIterator copy_n( InputIterator first,
                       Size n,
                       OutputIterator result) {
  // copies the first `n' items from `first' to `result'. Returns
  // the value of `result' after inserting the `n' items.
  while( n--) {
    *result = *first;
    first++;
    result++;
  }
  return result;
}

CGAL_END_NAMESPACE
#endif // CGAL_COPY_N_H //
