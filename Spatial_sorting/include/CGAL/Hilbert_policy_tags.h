// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.6-branch/Spatial_sorting/include/CGAL/hilbert_sort.h $
// $Id: hilbert_sort.h 53867 2010-01-28 12:18:19Z lrineau $
//
// Author(s)     : Olivier Devillers

#ifndef CGAL_HILBERT_POLICY_H
#define CGAL_HILBERT_POLICY_H


namespace CGAL {

struct Middle {};
struct Median {};


// A policy to select the sorting strategy.

template < typename Tag >
struct Hilbert_policy {};

typedef Hilbert_policy<Middle>      Hilbert_sort_middle_policy;
typedef Hilbert_policy<Median>      Hilbert_sort_median_policy;

} // namespace CGAL

#endif // CGAL_HILBERT_POLICY_H
