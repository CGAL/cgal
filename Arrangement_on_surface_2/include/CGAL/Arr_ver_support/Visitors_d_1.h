// Copyright (c) 2010  Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
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
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

#ifndef CGAL_VISITORS_D_1_H
#define CGAL_VISITORS_D_1_H

#include <CGAL/tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>

#include "boost/variant.hpp"


namespace CGAL {
namespace Arr_vertical_rational_arc {

template < class Non_vertical_arc,class Vertical_Segment >
class Source_infinite_in_x_visitor
  : public boost::static_visitor <Arr_parameter_space>
{
public:
  Arr_parameter_space operator() (const Non_vertical_arc & arc) const
  {
    return arc.source_infinite_in_x();
  }
  Arr_parameter_space operator() (const Vertical_Segment & ver) const
  {
    return CGAL::ARR_INTERIOR;
  }
};  //Source_infinite_in_x_visitor

template < class Non_vertical_arc,class Vertical_Segment >
class Source_infinite_in_y_visitor
  : public boost::static_visitor <Arr_parameter_space>
{
public:
  Arr_parameter_space operator() (const Non_vertical_arc & arc) const
  {
    return arc.source_infinite_in_y();
  }
  Arr_parameter_space operator() (const Vertical_Segment & ver) const
  {
    return ver.min_parameter_space();
  }
};  //Source_infinite_in_y_visitor
template < class Non_vertical_arc,class Vertical_Segment >
class Target_infinite_in_x_visitor
  : public boost::static_visitor <Arr_parameter_space>
{
public:
  Arr_parameter_space operator() (const Non_vertical_arc & arc) const
  {
    return arc.target_infinite_in_x();
  }
  Arr_parameter_space operator() (const Vertical_Segment & ver) const
  {
    return CGAL::ARR_INTERIOR;
  }
};  //Target_infinite_in_x_visitor

template < class Non_vertical_arc,class Vertical_Segment >
class Target_infinite_in_y_visitor
  : public boost::static_visitor <Arr_parameter_space>
{
public:
  Arr_parameter_space operator() (const Non_vertical_arc & arc) const
  {
    return arc.target_infinite_in_y();
  }
  Arr_parameter_space operator() (const Vertical_Segment & ver) const
  {
    return ver.min_parameter_space();
  }
};  //Target_infinite_in_y_visitor

template < class Non_vertical_arc,class Vertical_Segment, class Point_2 >
class Source_visitor
  : public boost::static_visitor <const Point_2&>
{
public:
  const Point_2& operator() (const Non_vertical_arc & arc) const
  {
    return arc.source();
  }
  const Point_2& operator() (const Vertical_Segment & ver) const
  {
    return ver.min();
  }
};  //Source_visitor

template < class Non_vertical_arc,class Vertical_Segment, class Algebraic >
class Source_x_visitor
  : public boost::static_visitor <Algebraic>
{
public:
  Algebraic operator() (const Non_vertical_arc & arc) const
  {
    return arc.source_x();
  }
  Algebraic operator() (const Vertical_Segment & ver) const
  {
    return ver.x();
  }
};  //Source_x_visitor

template < class Non_vertical_arc,class Vertical_Segment, class Point_2  >
class Target_visitor
  : public boost::static_visitor <const Point_2&>
{
public:
  const Point_2& operator() (const Non_vertical_arc & arc) const
  {
    return arc.target();
  }
  const Point_2& operator() (const Vertical_Segment & ver) const
  {
    return ver.max();
  }
};  //Target_visitor
template < class Non_vertical_arc,class Vertical_Segment, class Algebraic >
class Target_x_visitor
  : public boost::static_visitor <Algebraic>
{
public:
  Algebraic operator() (const Non_vertical_arc & arc) const
  {
    return arc.target_x();
  }
  Algebraic operator() (const Vertical_Segment & ver) const
  {
    return ver.x();
  }
};  //Target_x_visitor
template < class Non_vertical_arc,class Vertical_Segment >
class Left_infinite_in_x_visitor
  : public boost::static_visitor <Arr_parameter_space>
{
public:
  Arr_parameter_space operator() (const Non_vertical_arc & arc) const
  {
    return arc.left_infinite_in_x();
  }
  Arr_parameter_space operator() (const Vertical_Segment & ver) const
  {
    return CGAL::ARR_INTERIOR;
  }
};  //Left_infinite_in_x_visitor
template < class Non_vertical_arc,class Vertical_Segment >
class Left_infinite_in_y_visitor
  : public boost::static_visitor <Arr_parameter_space>
{
public:
  Arr_parameter_space operator() (const Non_vertical_arc & arc) const
  {
    return arc.left_infinite_in_y();
  }
  Arr_parameter_space operator() (const Vertical_Segment & ver) const
  {
    return ver.min_parameter_space();
  }
};  //Left_infinite_in_x_visitor

template < class Non_vertical_arc,class Vertical_Segment >
class Right_infinite_in_x_visitor
  : public boost::static_visitor <Arr_parameter_space>
{
public:
  Arr_parameter_space operator() (const Non_vertical_arc & arc) const
  {
    return arc.right_infinite_in_x();
  }
  Arr_parameter_space operator() (const Vertical_Segment & ver) const
  {
    return CGAL::ARR_INTERIOR;
  }
};  //right_infinite_in_x_visitor

template < class Non_vertical_arc,class Vertical_Segment >
class Right_infinite_in_y_visitor
  : public boost::static_visitor <Arr_parameter_space>
{
public:
  Arr_parameter_space operator() (const Non_vertical_arc & arc) const
  {
    return arc.right_infinite_in_y();
  }
  Arr_parameter_space operator() (const Vertical_Segment & ver) const
  {
    return ver.max_parameter_space();
  }
};  //Right_infinite_in_y_visitor

template < class Non_vertical_arc,class Vertical_Segment, class Point_2 >
class Left_visitor
  : public boost::static_visitor <const Point_2&>
{
public:
  const Point_2& operator() (const Non_vertical_arc & arc) const
  {
    return arc.left();
  }
  const Point_2& operator() (const Vertical_Segment & ver) const
  {
    return ver.min();
  }
};  //Left_visitor

template < class Non_vertical_arc,class Vertical_Segment, class Algebraic >
class Left_x_visitor
  : public boost::static_visitor <Algebraic>
{
public:
  Algebraic operator() (const Non_vertical_arc & arc) const
  {
    return arc.left_x();
  }
  Algebraic operator() (const Vertical_Segment & ver) const
  {
    return ver.x();
  }
};  //Left_x_visitor

template < class Non_vertical_arc,class Vertical_Segment, class Point_2 >
class Right_visitor
  : public boost::static_visitor <const Point_2&>
{
public:
  const Point_2& operator() (const Non_vertical_arc & arc) const
  {
    return arc.right();
  }
  const Point_2& operator() (const Vertical_Segment & ver) const
  {
    return ver.max();
  }
};  //Right_visitor

template < class Non_vertical_arc,class Vertical_Segment,class Algebraic >
class Right_x_visitor
  : public boost::static_visitor <Algebraic>
{
public:
  Algebraic operator() (const Non_vertical_arc & arc) const
  {
    return arc.right_x();
  }
  Algebraic operator() (const Vertical_Segment & ver) const
  {
    return ver.x();
  }
};  //Right_x_visitor

template < class Non_vertical_arc,class Vertical_Segment >
class Is_continuous_visitor
  : public boost::static_visitor <bool>
{
public:
  bool operator() (const Non_vertical_arc & arc) const
  {
    return arc.is_directed_right();
  }
  bool operator() (const Vertical_Segment & ver) const
  {
    return true;
  }
};  //Is_continuous_visitor

template < class Non_vertical_arc,class Vertical_Segment >
class Is_directed_right_visitor
  : public boost::static_visitor <bool>
{
public:
  bool operator() (const Non_vertical_arc & arc) const
  {
    return arc.is_directed_right();
  }
  bool operator() (const Vertical_Segment & ver) const
  {
    return true;
  }
};  //Is_directed_right_visitor


template < class Non_vertical_arc,class Vertical_Segment >
class Print_visitor
  : public boost::static_visitor <std::ostream&>
{
public:
  Print_visitor (std::ostream& os) : _os (&os) {}
  std::ostream& operator() (const Non_vertical_arc & arc) const
  {
    *_os << arc;
    return *_os;
  }
  std::ostream& operator() (const Vertical_Segment & ver) const
  {
    *_os << ver;
    return *_os;
  }
private:
  std::ostream* _os;

};  //Print_visitor

}   //name_space Arr_vertical_rational_arc
}       //namespace CGAL {
#endif //CGAL_VISITORS_D_1_H
