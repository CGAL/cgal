// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_POLYGON_FUNCTION_OBJECTS_H
#define CGAL_POLYGON_FUNCTION_OBJECTS_H

#include <CGAL/license/Partition_2.h>


#include <CGAL/is_y_monotone_2.h>

namespace CGAL {

template <class Traits>
class Is_vacuously_valid 
{
  public:

     Is_vacuously_valid(Traits ) {}

     template <class ForwardIterator>
     bool operator()(ForwardIterator, ForwardIterator)
     {  return true; }

};


template <class Traits>
class Is_convex_2
{
  public:
     Is_convex_2(Traits t): traits(t) {}
  
     template <class ForwardIterator>
     bool operator()(ForwardIterator first, ForwardIterator last)
     {  return is_convex_2(first, last, traits); }

  private:
     Traits  traits;
};

template <class Traits>
class Is_y_monotone_2
{
  public:
     Is_y_monotone_2(Traits t): traits(t) {}
  
     template <class ForwardIterator>
     bool operator()(ForwardIterator first, ForwardIterator last)
     {  return is_y_monotone_2(first, last, traits); }

  private:
     Traits  traits;
};

}

#endif // CGAL_POLYGON_FUNCTION_OBJECTS_H
