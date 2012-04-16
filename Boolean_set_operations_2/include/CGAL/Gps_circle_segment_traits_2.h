// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_CIRCLE_SEGMENT_TRAITS_2_H
#define CGAL_GPS_CIRCLE_SEGMENT_TRAITS_2_H

#include <CGAL/Gps_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>

namespace CGAL {

template <class Kernel_, bool Filer_ = true>
class Gps_circle_segment_traits_2 : 
  public Gps_traits_2<Arr_circle_segment_traits_2<Kernel_, Filer_> >
{
public:
  Gps_circle_segment_traits_2<Kernel_, Filer_>(bool use_cache = false): 
      Gps_traits_2<Arr_circle_segment_traits_2<Kernel_, Filer_> >()
  {
    this->m_use_cache = use_cache;
  }

};

} //namespace CGAL

#endif
