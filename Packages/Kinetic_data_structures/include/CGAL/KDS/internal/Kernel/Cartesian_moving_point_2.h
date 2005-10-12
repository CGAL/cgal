// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_CARTESIAN_MOVING_POINT_2_H_
#define CGAL_KDS_CARTESIAN_MOVING_POINT_2_H_
#include <CGAL/KDS/basic.h>
#include <iostream>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE;


template <class Coordinate_t>
class Cartesian_moving_point_2 {
protected:
  typedef Cartesian_moving_point_2<Coordinate_t> This;
public:
  //! The type for coordinate values
  typedef Coordinate_t Coordinate;
    
  //! What should I do for this
  typedef typename Coordinate::NT NT;

  //! initialize it from polys
  Cartesian_moving_point_2(const Coordinate &x, const Coordinate &y){
    _coords[0]=x;
    _coords[1]=y;
  }

  //! initialize it from a still point
  template <class Static_point>
  Cartesian_moving_point_2(const Static_point &pt){
    _coords[0]=pt.x();
    _coords[1]=pt.y();

  }
   
  //! null
  Cartesian_moving_point_2(){}

  //! homogeneous x
  const Coordinate &hx() const {
    return _coords[0];
  }

  //! homogeneous y
  const Coordinate &hy() const {
    return _coords[1];
  }

  //! homogeneous z
  const Coordinate hz() const {
    return Coordinate(0);
  }

  //! homogeneous w
  const Coordinate hw() const {
    return Coordinate(1);
  }

  //! x
  const Coordinate &x() const {
    return _coords[0];
  }

  //! y
  const Coordinate &y() const {
    return _coords[1];
  }

  //! z
  const Coordinate z() const {
    return Coordinate(0);
  }

  bool operator==(const This &o) const {
    return x()==o.x() && y()==o.y();
  }

  bool operator<(const This &o) const {
    return x()<o.x() || (x()==o.x() && y() < o.y());
  }

 bool operator>(const This &o) const {
  return x()>o.x() || (x()==o.x() && y() > o.y());
  }

#if 0    
  //! Returns the value at time t. 
  /*!
   */
  Static_point operator()(const NT &t) const {
    return Static_point(x()(t), y()(t));
  }

  //! Non-operator version of operator()
  Static_point value_at( NT time){
    return operator()(time);
  }

#endif

  template <class SK>
  struct Static_traits {
    typedef typename SK::Point_2 Static_type;
    
    static Static_type to_static(const This &o, const NT &t, const SK &) {
      return Static_type(o.x()(t), o.y()(t));
    } 
  }; 

  template <class Converter>
  struct Coordinate_converter {
    Coordinate_converter(const Converter &c): c_(c){}
    typedef Cartesian_moving_point_2<typename Converter::argument_type> argument_type;
    typedef Cartesian_moving_point_2<typename Converter::result_type> result_type;
    
    result_type operator()(const argument_type &i) const {
      return result_type(c_(i.x()), c_(i.y()));
    }

    Converter c_;
  };

  //! Reverse the motion, time must be negated also
  template <class NV>
  This transformed_coordinates(NV n) const {
    return This(n(_coords[0]), n(_coords[1]));
  }
    
  void write(std::ostream &out) const {
    out << x() << ", " << y();
  }
protected:
  Coordinate _coords[2];
};

template <class Coordinate>
std::ostream &operator<<(std::ostream &out,
			 const Cartesian_moving_point_2<Coordinate> &point){
  point.write(out);
  return out;
}

template <class Coordinate>
std::istream &operator>>(std::istream &in,
			 Cartesian_moving_point_2<Coordinate> &point){
  Coordinate x, y;
  in >> x;
  char c;
  in >> c; 
  if (c != ',') {
    in.setstate(std::ios_base::failbit);
    return in;
  }
  in >> y;
  point= Cartesian_moving_point_2<Coordinate>(x,y);
  return in;
}

CGAL_KDS_END_INTERNAL_NAMESPACE

CGAL_KDS_BEGIN_NAMESPACE

/*template <>
template <class Coord, class SK>
class To_static< internal::Cartesian_moving_point_2<Coord>, SK>:
  public To_static_base<typename Coord::NT, 
			typename internal::Cartesian_moving_point_2<Coord>,
			typename SK::Point_2> {
  typedef To_static_base<typename Coord::NT, 
			 typename internal::Cartesian_moving_point_2<Coord>,
			 typename SK::Point_2>  P;
public:
  To_static(){}
  typename P::result_type operator()(const typename P::argument_type &arg) const {
    return typename P::result_type(arg.x()(P::time()),
				   arg.y()(P::time()));
  }
  };*/
CGAL_KDS_END_NAMESPACE
#endif
