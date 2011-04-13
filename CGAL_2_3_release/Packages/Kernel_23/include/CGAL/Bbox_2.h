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
// file          : Bbox_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_BBOX_2_H
#define CGAL_BBOX_2_H

#include <CGAL/basic.h>
#include <CGAL/cartesian_classes.h>
#include <CGAL/Fourtuple.h>

CGAL_BEGIN_NAMESPACE

class Bbox_2 : public Handle_for< Fourtuple<double> >
{
public:
             Bbox_2();
             Bbox_2(double x_min, double y_min,
                    double x_max, double y_max);

  bool       operator==(const Bbox_2 &b) const;
  bool       operator!=(const Bbox_2 &b) const;

  int        dimension() const;
  double     xmin() const;
  double     ymin() const;
  double     xmax() const;
  double     ymax() const;

  double     max(int i) const;
  double     min(int i) const;

  Bbox_2     operator+(const Bbox_2 &b) const;

};


inline
int
Bbox_2::dimension() const
{ return 2; }

inline
double
Bbox_2::xmin() const
{ return ptr->e0; }

inline
double
Bbox_2::ymin() const
{ return ptr->e1; }

inline
double
Bbox_2::xmax() const
{ return ptr->e2; }

inline
double
Bbox_2::ymax() const
{ return ptr->e3; }

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
inline Bbox_2 Bbox_2::operator+(const Bbox_2 &b) const
{
  return Bbox_2(std::min(xmin(), b.xmin()),
                std::min(ymin(), b.ymin()),
                std::max(xmax(), b.xmax()),
                std::max(ymax(), b.ymax()));
}
inline bool do_overlap(const Bbox_2 &bb1, const Bbox_2 &bb2)
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
