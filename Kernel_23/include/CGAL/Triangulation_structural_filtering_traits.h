// Copyright (c) 2010  GeometryFactory Sarl (France).
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
//
// Author(s)     : Laurent Rineau <Laurent.Rineau__CGAL@normalesup.org>
// 

#ifndef CGAL_TRIANGULATION_STRUCTURAL_FILTERING_TRAITS_H
#define CGAL_TRIANGULATION_STRUCTURAL_FILTERING_TRAITS_H

#include <CGAL/tags.h>

namespace CGAL {

template <typename Geom_traits>
struct Triangulation_structural_filtering_traits {
  typedef Tag_false Use_structural_filtering_tag;
};

} // namespace CGAL

#endif // CGAL_TRIANGULATION_STRUCTURAL_FILTERING_TRAITS_H
