// Copyright (c) 2026 GeometryFactory
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_COLOR_RAMP_H
#define CGAL_IO_COLOR_RAMP_H

#include <list>
#include <CGAL/IO/Color.h>

namespace CGAL {
namespace IO {



class Color_component
{
  typedef std::list<std::pair<double,double> > Values;

public:
  Color_component()
  {
    add(0,0);
    add(1,1);
  }

  Color_component(const double c0, const double c1)
  {
    add(0,c0);
    add(1,c1);
  }


  ~Color_component() {}

  double interpolate(const double v) const
  {
    Values::const_iterator next = next_it(v);

    // next is just after v
    Values::const_iterator prev = --next;
    ++next;

    if ( v>=1 && next != values_.end())
      std::cerr << ".";

    if ( next == values_.end() )
    {
      return prev->second;
    }

    const double& a = prev->first;
    const double& b = next->first;
    return (b-v)/(b-a) * prev->second + (v-a)/(b-a) * next->second;
}

  void add(const double v, double color)
  {
    if ( color > 1 ) { color = 1; }
    if ( color < 0 ) { color = 0; }

    Values::iterator next = next_it(v);
    values_.insert(next, std::make_pair(v,color));
  }


  void rebuild(const double c0, const double c1)
  {
    values_.clear();
    add(1,c1);
    add(0,c0);
  }

  void print() const
  {
    for (Values::const_iterator it = values_.begin(),
         end = values_.end() ; it != end ; ++it )
    {
      std::cout << "<" << it->first << "," << it->second << "> ";
    }

    std::cout << std::endl;
  }

private:
  inline Values::const_iterator next_it(const double v) const;
  inline Values::iterator next_it(const double v);

private:
  Values values_;
};

/*!
  \ingroup PkgStreamSupportRef

  An object of the class `Color_ramp` can map a scalar between 0 and 1 to a `Color`.
 */
class Color_ramp
{
public :
  Color_ramp()
  : r_(), g_(), b_()
  {}

  Color_ramp(const double r0, const double r1,
             const double g0, const double g1,
             const double b0, const double b1)
  {
    r_.rebuild(r0,r1);
    g_.rebuild(g0,g1);
    b_.rebuild(b0,b1);

    r_.add(0.5,(std::max)(r0,r1));
    g_.add(0.5,(std::max)(g0,g1));
    b_.add(0.5,(std::max)(b0,b1));
  }


public :
  inline double r(double v) const;
  inline double g(double v) const;
  inline double b(double v) const;

  /*!
  the interpolation function.*/
  Color interpolate(double v) const
  {
    return Color(255*r(v), 255*g(v), 255*b(v));
  }

  /*!
  A red color ramp
  */
  void build_red()
  {
    r_.rebuild(1,0.5);
    g_.rebuild(1,0);
    b_.rebuild(1,0);

    r_.add(0.3,1);
    r_.add(0.05,1);
    g_.add(0.05,0.95);
    g_.add(0.3,0.5);
    b_.add(0.05,0);
  }

  /*!
  A blue color ramp
  */
  void build_blue()
  {
    r_.rebuild(1,0);
    g_.rebuild(1,0);
    b_.rebuild(1,0.5);

    b_.add(0.1,0.8);
    g_.add(0.1,0.4);
    r_.add(0.1,0.4);
  }


  /*!
  A thermal color ramp
  */
  void build_thermal()
  {
    r_.rebuild(1,0.5);
    g_.rebuild(1,0);
    b_.rebuild(1,0);

    r_.add(0.0,0);
    r_.add(1,1);
    g_.add(0.05,0.8);
    g_.add(0.3,0.5);
    b_.add(0.0,1.0);
    b_.add(1.0,0.0);
  }


  /*!
  A rainbow color ramp
  */
  void build_rainbow()
  {
  r_.rebuild(1,0.5);
  g_.rebuild(1,0);
  b_.rebuild(1,0);

  r_.add(0,0.75);
  g_.add(0,0.75);
  b_.add(0,1);

  r_.add(0.2,0);
  g_.add(0.2,0);
  b_.add(0.2,1);

  r_.add(0.4,0);
  g_.add(0.4,1);
  b_.add(0.4,0);

  r_.add(0.6,1);
  g_.add(0.6,1);
  b_.add(0.6,0);

  r_.add(0.8,1);
  g_.add(0.8,0);
  b_.add(0.8,0);
  }


  void print() const
  {
  std::cout << "r: ";
  r_.print();
  std::cout << "g: ";
  g_.print();
  std::cout << "b: ";
  b_.print();
}

private :
  Color_component r_;
  Color_component g_;
  Color_component b_;
};


inline
Color_component::Values::const_iterator
Color_component::
next_it(const double v) const
{
  Values::const_iterator next = values_.begin();
  while ( next != values_.end() && v >= next->first ) { ++next; }
  return next;
}

inline
Color_component::Values::iterator
Color_component::
next_it(const double v)
{
  Values::iterator next = values_.begin();
  while ( next != values_.end() && v >= next->first ) { ++next; }
  return next;
}

inline
double
Color_ramp::r(double v) const
{ return r_.interpolate(v); }

inline
double
Color_ramp::g(double v) const
{ return g_.interpolate(v);  }

inline
double
Color_ramp::b(double v) const
{ return b_.interpolate(v);  }

} // namespace IO
} // namespace CGAL

#endif // CGAL_IO_COLOR_RAMP
