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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_TRIANGULATION_INDIRECT_TRAITS_2_H
#define CGAL_TRIANGULATION_INDIRECT_TRAITS_2_H

#include <CGAL/license/Partition_2.h>


#include <CGAL/Kernel/function_objects.h>

namespace CGAL {

template <class Circulator>
class Indirect_segment 
{
public:
   Indirect_segment() {}
   Indirect_segment(Circulator s, Circulator t) : _source_ref(s), 
                                                  _target_ref(t)
   {}
   Circulator source() {return _source_ref;}
   Circulator target() {return _target_ref;}

private:
   Circulator _source_ref;
   Circulator _target_ref;
};

template <class Circulator>
class Indirect_triangle 
{
public:
   Indirect_triangle() {}
   Indirect_triangle(Circulator p0, Circulator p1, Circulator p2): 
       _p0(p0), _p1(p1), _p2(p2) 
   {}

private:
   Circulator _p0, _p1, _p2;
};

template <class Traits>
class Indirect_compare_x_2
{
   typedef typename Traits::Compare_x_2 Compare_x_2;
   typedef typename Traits::Point_2     Point;
public:
   Indirect_compare_x_2(const Compare_x_2& compare_x_2)
     : _compare_x_2(compare_x_2)
   {}

   template <class Point_2_ptr>
   Comparison_result operator()(Point_2_ptr p1, Point_2_ptr p2)
   {
      return _compare_x_2(Point(*p1), Point(*p2));
   }

private:
   Compare_x_2 _compare_x_2;
};

template <class Traits>
class Indirect_compare_y_2
{
   typedef typename Traits::Compare_y_2 Compare_y_2;
   typedef typename Traits::Point_2     Point;
public:
   Indirect_compare_y_2(const Compare_y_2& compare_y_2)
     : _compare_y_2(compare_y_2)
   {}

   template <class Point_2_ptr>
   Comparison_result operator()(Point_2_ptr p1, Point_2_ptr p2)
   {
      return _compare_y_2(Point(*p1), Point(*p2));
   }

private:
   Compare_y_2 _compare_y_2;
};

template <class Traits>
class Indirect_orientation_2
{
   typedef typename Traits::Orientation_2 Orientation_2;
   typedef typename Traits::Point_2 Point;
public:
   Indirect_orientation_2(const Orientation_2& orientation_2)
     : _orientation_2(orientation_2)
   {}

   template <class Point_2_ptr>
   Orientation operator()(Point_2_ptr p1, Point_2_ptr p2, Point_2_ptr p3)
   {
      return _orientation_2(Point(*p1), Point(*p2), Point(*p3));
   }

private:
   Orientation_2 _orientation_2;
};

template <class Circulator>
class Construct_circulator_2
{
public:
   Circulator operator()(Circulator p1) const { return p1; }
};

template <class Circulator>
class Construct_indirect_segment_2
{
public:
   typedef Indirect_segment<Circulator>   I_segment;

   I_segment operator()(Circulator p1, Circulator p2)
   {
      return I_segment(p1, p2);
   }
};

template <class Circulator, class Traits>
class Triangulation_indirect_traits_2 
{
public:

  typedef Circulator                      Point_2;
  typedef Indirect_segment<Circulator>    Segment_2;
  typedef Indirect_triangle<Circulator>   Triangle_2;

  typedef Indirect_orientation_2<Traits>  Orientation_2;
  typedef Indirect_compare_x_2<Traits>    Compare_x_2;
  typedef Indirect_compare_y_2<Traits>    Compare_y_2;
  typedef Construct_indirect_segment_2<Circulator>      Construct_segment_2;
  typedef Construct_circulator_2<Circulator>            Construct_point_2;

  // constructor
  Triangulation_indirect_traits_2 (const Traits& traits)
    : _traits(traits)
  { }

   Compare_x_2 compare_x_2_object() const
   {
     return Compare_x_2(_traits.compare_x_2_object());
   }  

   Compare_y_2 compare_y_2_object() const
   {
     return Compare_y_2(_traits.compare_y_2_object());
   }  

   Orientation_2 orientation_2_object() const
   {
     return Orientation_2(_traits.orientation_2_object());
   }

   Construct_segment_2
   construct_segment_2_object() const
   { return Construct_segment_2(); }

   Construct_point_2 construct_point_2_object() const
   {
     return Construct_point_2();
   }

private:
   const Traits& _traits;
};

}

#endif // CGAL_TRIANGULATION_INDIRECT_TRAITS_2_H

// For the Emacs editor
// Local Variables:
// c-basic-offset: 3
// End:
