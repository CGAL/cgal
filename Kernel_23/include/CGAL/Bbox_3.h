// Copyright (c) 1999,2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BBOX_3_H
#define CGAL_BBOX_3_H

#include <CGAL/config.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/result_of.h>
#include <CGAL/IO/io.h>
#include <CGAL/Dimension.h>
#include <CGAL/array.h>
#include <boost/math/special_functions/next.hpp>

namespace CGAL {

template < typename T >
struct Simple_cartesian;

class Bbox_3
{
  std::array<double, 6>   rep;

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<3>  Feature_dimension;

  typedef Simple_cartesian<double>  R;

  Bbox_3()
    : rep(CGAL::make_array( std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity(),
                            - std::numeric_limits<double>::infinity(),
                            - std::numeric_limits<double>::infinity(),
                            - std::numeric_limits<double>::infinity() ))
  {}

  Bbox_3(double x_min, double y_min, double z_min,
         double x_max, double y_max, double z_max)
    : rep(CGAL::make_array(x_min, y_min, z_min, x_max, y_max, z_max))
  {}

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

  inline double min_coord(int i) const { return (min)(i); }
  inline double max_coord(int i) const { return (max)(i); }

  Bbox_3  operator+(const Bbox_3& b) const;
  Bbox_3& operator+=(const Bbox_3& b);

  void dilate(int dist);
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
void
Bbox_3::dilate(int dist)
{
  using boost::math::float_advance;
  rep[0] = float_advance(rep[0],-dist);
  rep[1] = float_advance(rep[1],-dist);
  rep[2] = float_advance(rep[2],-dist);
  rep[3] = float_advance(rep[3],dist);
  rep[4] = float_advance(rep[4],dist);
  rep[5] = float_advance(rep[5],dist);
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
  switch(get_mode(os))
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
    case IO::PRETTY :
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
    double xmin = 0;
    double ymin = 0;
    double zmin = 0;
    double xmax = 0;
    double ymax = 0;
    double zmax = 0;

  switch(get_mode(is))
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
    case IO::PRETTY :
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
