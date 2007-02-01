/* Copyright 2004
Stanford University

This file is part of the DSR PDB Library.

The DSR PDB Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The DSR PDB Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#ifndef DSR_PDB_POINT_H
#define DSR_PDB_POINT_H
#include <iostream>
#include <dsrpdb/config.h>

#ifdef PDB_USE_CGAL 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_3.h>
namespace dsrpdb {
  // An alternative point class.
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 Point;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Vector_3 Vector;

  struct Squared_distance{
    double operator()(const Point &a, const Point &b) const {
      return CGAL::squared_distance(a,b);
    }
  };
}

#else 
namespace dsrpdb {
  //! A class representing a geometric point.
  /*!  Feel free to replace this with anything you want to use your
    own point. I just construct and call x(), y(), z().
  */
  class Point {
  public:
    //! Default.
    Point(): x_(0), y_(0), z_(0){}
    //! Construct a point.
    Point(double x, double y, double z): x_(x), y_(y), z_(z){}

    //! x
    double x() const {return x_;}
    //! y
    double y() const {return y_;}

    //! y
    double z() const {return z_;}

    //! Get one of the coordinates by index
    double operator[](unsigned int i) const {
      switch(i) {
      case 0: return x_;
      case 1: return y_;
      case 2: return z_;
      default:
	assert(i < 3);
	return -1;
      }
    }

    //! difference between two points
    Point operator-(const Point &o) const {
      return Point(x_-o.x_, y_-o.y_, z_-o.z_);
    }

    //! Sum of two points
    Point operator+(const Point &o) const {
      return Point(x_+o.x_, y_+o.y_, z_+o.z_);
    }

    //! Divide by a scalar
    Point operator/(double d) const {
      return Point(x_/d, y_/d, z_/d);
    }

    //! Dot product
    double operator*(const Point &o) const {
      return x_*o.x_ + y_*o.y_ + z_*o.z_;
    }
  protected:
    double x_, y_, z_;
  };
  
  typedef Point Vector;

  inline std::ostream &operator<<(std::ostream &o, const Point &p) {
    o << p.x() << " " << p.y() << " " << p.z();
    return o;
  }

  inline std::istream &operator>>(std::istream &i, Point &p) {
    double x,y,z;
    i >> x >> y >> z;
    if (!i) return i;
    p= Point(x,y,z);
    return i;
  }

  //! A functor that computes the squared distance between two points.
  struct Squared_distance{
    double operator()(const Point &a, const Point &b) const {
      return (a.x()-b.x())*(a.x()-b.x()) 
	+ (a.y()-b.y())*(a.y()-b.y())
	+ (a.z()-b.z())*(a.z()-b.z());
    }
  };
  
};
#endif

#endif
