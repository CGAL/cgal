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
// file          : Point_d.h
// package       : _d
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schoenherr
//                 Bernd Gaertner
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINT_D_H
#define CGAL_POINT_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/PointHd.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/PointCd.h>
// #include <CGAL/Cartesian/Point_d.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class R_ >
class Point_d : public R_::Point_d_base
{
public:
    typedef  R_   R;
    typedef typename R::RT              RT;
    typedef typename R::FT              FT;
    typedef typename R::Point_d_base    Point;

    Point_d()
    {}

    Point_d(int dim, const Origin &o)
        : Point(dim, o)
    {}

    Point_d(const Point_d<R> &p)
        : Point((Point&)p)
    {}

    Point_d(const Point &p)
        : Point(p)
    {}

    template <class InputIterator>
    Point_d (int dim, InputIterator first, InputIterator last)
        : Point (dim, first, last)
    {}

    Point_d<R>& operator=(const Point_d<R>& p)
    {
          Point::operator=(p);
          return *this;
    }

    bool operator==(const Point_d<R>& p) const
    {
    return Point::operator==(p);
    }

    bool operator!=(const Point_d<R>& p) const
    {
        return !(*this == p);
    }

    int id() const
    {
        return (int) PTR;
    }

    RT homogeneous(int i) const
    {
        return Point::homogeneous(i);
    }

    FT cartesian(int i) const
    {
        return Point::cartesian(i);
    }

    FT operator[](int i) const
    {
        return Point::operator[](i);
    }

    typedef const RT* const_iterator;

    const_iterator begin() const
    {
        return Point::begin();
    }

    const_iterator end() const
    {
        return Point::end();
    }

    int dimension() const
    {
    return Point::dimension();
    }
};

#ifndef CGAL_NO_OSTREAM_INSERT_POINT_D
template < class R >
std::ostream&
operator<<(std::ostream& os, const Point_d<R>& p)
{
  typedef typename  R::Point_d_base    Point;
  return os << (const Point&)p;
}
#endif // CGAL_NO_OSTREAM_INSERT_POINT_D

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINT_D
template < class R >
std::istream&
operator>>(std::istream& is, Point_d<R> &p)
{
  typedef typename  R::Point_d_base    Point;
  return is >> (Point&)p;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINT_D

CGAL_END_NAMESPACE

#endif // CGAL_POINT_D_H
