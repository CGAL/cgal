// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : Bbox_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_BBOX_3_H
#define CGAL_BBOX_3_H

#include <CGAL/basic.h>
#include <CGAL/cartesian_classes.h>
#include <CGAL/Sixtuple.h>

CGAL_BEGIN_NAMESPACE

class Bbox_3 : public Handle_for< Sixtuple<double> >
{
public:
                         Bbox_3();
                         Bbox_3(double x_min, double y_min, double zmin,
                                double x_max, double y_max, double z_max);

  double                 xmin() const;
  double                 ymin() const;
  double                 zmin() const;
  double                 xmax() const;
  double                 ymax() const;
  double                 zmax() const;

  Bbox_3             operator+(const Bbox_3& b) const;
};

inline
Bbox_3::Bbox_3()
{ new ( static_cast< void*>(ptr)) Sixtuple<double>(); }

inline
Bbox_3::Bbox_3(double x_min, double y_min, double z_min,
               double x_max, double y_max, double z_max)
{
  new ( static_cast< void*>(ptr)) Sixtuple<double>(x_min, y_min, z_min,
                                                   x_max, y_max, z_max);
}

inline
double
Bbox_3::xmin() const
{ return ptr->e0; }

inline
double
Bbox_3::ymin() const
{ return ptr->e1; }

inline
double
Bbox_3::zmin() const
{ return ptr->e2; }

inline
double
Bbox_3::xmax() const
{ return ptr->e3; }

inline
double
Bbox_3::ymax() const
{ return ptr->e4; }

inline
double
Bbox_3::zmax() const
{ return ptr->e5; }

inline Bbox_3 Bbox_3::operator+(const Bbox_3& b) const
{
  return Bbox_3(std::min(xmin(), b.xmin()),
                std::min(ymin(), b.ymin()),
                std::min(zmin(), b.zmin()),
                std::max(xmax(), b.xmax()),
                std::max(ymax(), b.ymax()),
                std::max(zmax(), b.zmax()));
}

inline bool do_overlap(const Bbox_3& bb1, const Bbox_3& bb2)
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

#ifndef CGAL_NO_OSTREAM_INSERT_BBOX_3
inline
std::ostream&
operator<<(std::ostream &os, const Bbox_3& b)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << b.xmin() << ' ' << b.ymin() << ' ' << b.zmin();
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
#endif // CGAL_NO_OSTREAM_INSERT_BBOX_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_BBOX_3
inline
std::istream&
operator>>(std::istream &is, Bbox_3& b)
{
  double xmin, ymin, zmin, xmax, ymax, zmax;

  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> xmin >> ymin >> xmax >> ymax;
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
  b = Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
  return is;
}

#endif // CGAL_NO_ISTREAM_EXTRACT_BBOX_3

CGAL_END_NAMESPACE

#endif // CGAL_BBOX_3_H
