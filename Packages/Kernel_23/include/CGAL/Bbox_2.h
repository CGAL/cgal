// Copyright (c) 1999,2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Andreas Fabri
 

#ifndef CGAL_BBOX_2_H
#define CGAL_BBOX_2_H

#include <CGAL/basic.h>
#include <CGAL/IO/io.h>
#include <CGAL/Fourtuple.h>

CGAL_BEGIN_NAMESPACE

class Bbox_2
{
  typedef Fourtuple<double>            BBox_rep_2;

  BBox_rep_2 rep;

public:
             Bbox_2() {}

             Bbox_2(double x_min, double y_min,
                    double x_max, double y_max)
		 : rep(x_min, y_min, x_max, y_max) {}

  inline bool       operator==(const Bbox_2 &b) const;
  inline bool       operator!=(const Bbox_2 &b) const;

  inline int        dimension() const;
  inline double     xmin() const;
  inline double     ymin() const;
  inline double     xmax() const;
  inline double     ymax() const;

  inline double     max(int i) const;
  inline double     min(int i) const;

  inline Bbox_2     operator+(const Bbox_2 &b) const;

};

inline
double
Bbox_2::xmin() const
{ return rep.e0; }

inline
double
Bbox_2::ymin() const
{ return rep.e1; }

inline
double
Bbox_2::xmax() const
{ return rep.e2; }

inline
double
Bbox_2::ymax() const
{ return rep.e3; }

inline
bool
Bbox_2::operator==(const Bbox_2 &b) const
{
  return    xmin() == b.xmin() && xmax() == b.xmax()
         && ymin() == b.ymin() && ymax() == b.ymax();
}

inline
bool
Bbox_2::operator!=(const Bbox_2 &b) const
{
  return ! (b == *this);
}

inline
int
Bbox_2::dimension() const
{ return 2; }

inline
double
Bbox_2::min(int i) const
{
  CGAL_kernel_precondition( (i == 0 ) || ( i == 1 ) );
  if(i == 0) { return xmin(); }
  return ymin();
}

inline
double
Bbox_2::max(int i) const
{
  CGAL_kernel_precondition( (i == 0 ) || ( i == 1 ) );
  if(i == 0) { return xmax(); }
  return ymax();
}

inline
Bbox_2
Bbox_2::operator+(const Bbox_2 &b) const
{
  return Bbox_2(std::min(xmin(), b.xmin()),
                std::min(ymin(), b.ymin()),
                std::max(xmax(), b.xmax()),
                std::max(ymax(), b.ymax()));
}

inline
bool
do_overlap(const Bbox_2 &bb1, const Bbox_2 &bb2)
{
    // check for emptiness ??
    if (bb1.xmax() < bb2.xmin() || bb2.xmax() < bb1.xmin())
        return false;
    if (bb1.ymax() < bb2.ymin() || bb2.ymax() < bb1.ymin())
        return false;
    return true;
}

#ifndef CGAL_NO_OSTREAM_INSERT_BBOX_2
inline
std::ostream&
operator<<(std::ostream &os, const Bbox_2 &b)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        os << b.xmin() << ' ' << b.ymin() << ' '
           << b.xmax() << ' ' << b.ymax();
        break;
    case IO::BINARY :
        write(os, b.xmin());
        write(os, b.ymin());
        write(os, b.xmax());
        write(os, b.ymax());
        break;
    default:
        os << "Bbox_2(" << b.xmin() << ", " << b.ymin() << ", "
                        << b.xmax() << ", " << b.ymax() << ")";
        break;
    }
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_BBOX_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_BBOX_2
inline
std::istream&
operator>>(std::istream &is, Bbox_2 &b)
{
    double xmin, ymin, xmax, ymax;

    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> xmin >> ymin >> xmax >> ymax;
        break;
    case IO::BINARY :
        read(is, xmin);
        read(is, ymin);
        read(is, xmax);
        read(is, ymax);
        break;
    }
    b = Bbox_2(xmin, ymin, xmax, ymax);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_BBOX_2

CGAL_END_NAMESPACE

#endif // CGAL_BBOX_2_H
