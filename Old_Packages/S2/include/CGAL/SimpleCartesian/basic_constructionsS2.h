// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/basic_constructionsS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================


#ifndef CGAL_BASIC_CONSTRUCTIONSS2_H
#define CGAL_BASIC_CONSTRUCTIONSS2_H 1

#include <utility>
#include <CGAL/SimpleCartesian/PointS2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
PointS2<FT>
midpoint( PointS2<FT> const& p,
          PointS2<FT> const& q )
{
  FT x,y;
  midpointC2(p.x(),p.y(),q.x(),q.y(),x,y);
  return PointS2<FT>(x,y);
}

template < class FT >
LineS2<FT>
bisector( PointS2<FT> const& p,
          PointS2<FT> const& q )
{
  FT a,b,c;
  bisector_of_pointsC2(p.x(),p.y(),q.x(),q.y(),a,b,c);
  return LineS2<FT>(a,b,c);
}

template < class FT >
PointS2<FT>
circumcenter( PointS2<FT> const& p,
              PointS2<FT> const& q,
              PointS2<FT> const& r)
{
  FT x,y;
  circumcenterC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),x,y);
  return PointS2<FT>(x,y);
}

template < class FT >
FT
squared_circumradius( PointS2<FT> const& p,
                      PointS2<FT> const& q,
                      PointS2<FT> const& r)
{
  return squared_circumradius(p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}

template < class FT >
PointS2<FT>
squared_circumcircle( PointS2<FT> const& p,
                      PointS2<FT> const& q,
                      PointS2<FT> const& r,
                      FT &radius)
{
  FT x,y;
  radius = squared_circumradius(p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),x,y);
  return PointS2<FT>(x,y);
}


template < class FT >
FT
squared_distance( PointS2<FT> const& p,
                  PointS2<FT> const& q)
{
  return squared_distanceC2(p.x(),p.y(),q.x(),q.y());
}

template < class FT >
FT
scaled_distance_to_line( LineS2<FT> const& l,
                         PointS2<FT> const& p)
{
  return squared_distance_to_lineC2(l.a(),l.b(),l.c(),p.x(),p.y());
}

template < class FT >
FT
scaled_distance_to_line( PointS2<FT> const& p,
                         PointS2<FT> const& q,
                         PointS2<FT> const& r)
{
  return squared_distance_to_lineC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}


CGAL_END_NAMESPACE

#endif // CGAL_BASIC_CONSTRUCTIONS_2_H
