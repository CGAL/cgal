// Copyright (c) 1999,2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_BBOX_3_H
#define CGAL_BBOX_3_H

#include <CGAL/basic.h>
#include <CGAL/IO/io.h>
#include <CGAL/Dimension.h>
#include <CGAL/array.h>

namespace CGAL {

template < typename T >
struct Simple_cartesian;

class Bbox_3
{
  cpp11::array<double, 6>   rep;

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<3>  Feature_dimension;

  typedef Simple_cartesian<double>  R;

        Bbox_3() {}

        Bbox_3(double x_min, double y_min, double z_min,
               double x_max, double y_max, double z_max)
	  : rep(CGAL::make_array(x_min, y_min, z_min, x_max, y_max, z_max)) {}

  inline bool operator==(const Bbox_3 &b) const;
  inline bool operator!=(const Bbox_3 &b) const;

  inline int dimension() const;
  double  xmin() const;
  double  ymin() const;
  double  zmin() const;
  double  xmax() const;
  double  ymax() const;
  double  zmax() const;

  inline double min BOOST_PREVENT_MACRO_SUBSTITUTION (int i) const;
  inline double max BOOST_PREVENT_MACRO_SUBSTITUTION (int i) const;

  Bbox_3  operator+(const Bbox_3& b) const;
  Bbox_3& operator+=(const Bbox_3& b);
};

inline
double
Bbox_3::xmin() const
{ return rep[0]; }

inline
double
Bbox_3::ymin() const
{ return rep[1]; }

inline
double
Bbox_3::zmin() const
{ return rep[2]; }

inline
double
Bbox_3::xmax() const
{ return rep[3]; }

inline
double
Bbox_3::ymax() const
{ return rep[4]; }

inline
double
Bbox_3::zmax() const
{ return rep[5]; }

inline
bool
Bbox_3::operator==(const Bbox_3 &b) const
{
  return xmin() == b.xmin() && xmax() == b.xmax()
      && ymin() == b.ymin() && ymax() == b.ymax()
      && zmin() == b.zmin() && zmax() == b.zmax();
}

inline
bool
Bbox_3::operator!=(const Bbox_3 &b) const
{
  return ! (b == *this);
}

inline
int
Bbox_3::dimension() const
{ return 3; }

inline
double
Bbox_3::min BOOST_PREVENT_MACRO_SUBSTITUTION (int i) const
{
  CGAL_kernel_precondition( (i == 0 ) || ( i == 1 ) || ( i == 2) );
  if (i == 0) { return xmin(); }
  if (i == 1) { return ymin(); }
  return zmin();
}

inline
double
Bbox_3::max BOOST_PREVENT_MACRO_SUBSTITUTION (int i) const
{
  CGAL_kernel_precondition( (i == 0 ) || ( i == 1 ) || ( i == 2 ) );
  if (i == 0) { return xmax(); }
  if (i == 1) { return ymax(); }
  return zmax();
}

inline
Bbox_3
Bbox_3::operator+(const Bbox_3& b) const
{
  return Bbox_3((std::min)(xmin(), b.xmin()),
                (std::min)(ymin(), b.ymin()),
                (std::min)(zmin(), b.zmin()),
                (std::max)(xmax(), b.xmax()),
                (std::max)(ymax(), b.ymax()),
                (std::max)(zmax(), b.zmax()));
}

inline
Bbox_3&
Bbox_3::operator+=(const Bbox_3& b)
{
  rep[0] = (std::min)(xmin(), b.xmin());
  rep[1] = (std::min)(ymin(), b.ymin());
  rep[2] = (std::min)(zmin(), b.zmin());
  rep[3] = (std::max)(xmax(), b.xmax());
  rep[4] = (std::max)(ymax(), b.ymax());
  rep[5] = (std::max)(zmax(), b.zmax());
  return *this;
}

inline
bool
do_overlap(const Bbox_3& bb1, const Bbox_3& bb2)
{
    // check for emptiness ??
    if (bb1.xmax() < bb2.xmin() || bb2.xmax() < bb1.xmin())
        return false;
    if (bb1.ymax() < bb2.ymin() || bb2.ymax() < bb1.ymin())
        return false;
    if (bb1.zmax() < bb2.zmin() || bb2.zmax() < bb1.zmin())
        return false;
    return true;
}


inline
std::ostream&
operator<<(std::ostream &os, const Bbox_3& b)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << b.xmin() << ' ' << b.ymin() << ' ' << b.zmin()
		  << ' ' << b.xmax() << ' ' << b.ymax() << ' ' << b.zmax();
    case IO::BINARY :
        write(os, b.xmin());
        write(os, b.ymin());
        write(os, b.zmin());
        write(os, b.xmax());
        write(os, b.ymax());
        write(os, b.zmax());
        return os;
    default:
        os << "Bbox_3((" << b.xmin()
           << ", "       << b.ymin()
           << ", "       << b.zmin() << "), (";
        os <<               b.xmax()
           << ", "       << b.ymax()
           << ", "       << b.zmax() << "))";
        return os;
  }
}

inline
std::istream&
operator>>(std::istream &is, Bbox_3& b)
{
  double xmin, ymin, zmin, xmax, ymax, zmax;

  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
      is >> iformat(xmin) >> iformat(ymin) >> iformat(zmin)
         >> iformat(xmax) >> iformat(ymax) >> iformat(zmax);
        break;
    case IO::BINARY :
        read(is, xmin);
        read(is, ymin);
        read(is, zmin);
        read(is, xmax);
        read(is, ymax);
        read(is, zmax);
        break;
  }
  if (is)
    b = Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
  return is;
}

template <class Input_iterator, class Traits>
Bbox_3 bbox_3(Input_iterator begin, Input_iterator end, const Traits& traits)
{
  if (begin==end) return Bbox_3();
  typename Traits::Construct_bbox_3 get_bbox = traits.construct_bbox_3_object();
  Bbox_3 res = get_bbox( *begin );
  for (++begin; begin!=end; ++begin)
    res += get_bbox( *begin );
  return res;
}

template <class Input_iterator>
Bbox_3 bbox_3(Input_iterator begin, Input_iterator end)
{
  if (begin==end) return Bbox_3();
  Bbox_3 res = begin->bbox();
  for (++begin; begin!=end; ++begin)
    res += begin->bbox();
  return res;
}


} //namespace CGAL

#endif // CGAL_BBOX_3_H
