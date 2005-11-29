// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef GPS_CIRCLE_SEGMENT_TRAITS_2_H
#define GPS_CIRCLE_SEGMENT_TRAITS_2_H

#include <CGAL/Gps_traits_adaptor_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Gps_circle_segment_traits_2 : 
  public Gps_traits_adaptor_2<Arr_circle_segment_traits_2<Kernel_> >
{};

CGAL_END_NAMESPACE

#endif
