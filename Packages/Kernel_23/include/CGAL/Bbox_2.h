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
#include <CGAL/Fourtuple.h>
#include <CGAL/Simple_Handle_for.h>

CGAL_BEGIN_NAMESPACE

class Bbox_2 : public Simple_Handle_for< Fourtuple<double> >
{
  typedef Simple_Handle_for< Fourtuple<double> > BBox_handle_2;
  typedef BBox_handle_2::element_type                   BBox_ref_2;

public:
             Bbox_2()
		 : BBox_handle_2(BBox_ref_2()) {}

             Bbox_2(double x_min, double y_min,
                    double x_max, double y_max)
		 : BBox_handle_2(BBox_ref_2(x_min, y_min,
                                            x_max, y_max)) {}

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
{ return Ptr()->e0; }

inline
double
Bbox_2::ymin() const
{ return Ptr()->e1; }

inline
double
Bbox_2::xmax() const
{ return Ptr()->e2; }

inline
double
Bbox_2::ymax() const
{ return Ptr()->e3; }

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
