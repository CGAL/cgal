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
// release_date  : 2000, October 15
//
// source        : webS3/S3.lw
// file          : include/CGAL/SimpleCartesian/basic_constructionsS3.h
// package       : S3 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.7
// revision_date : 15 Oct 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_BASIC_CONSTRUCTIONSS3_H
#define CGAL_BASIC_CONSTRUCTIONSS3_H

#include <CGAL/SimpleCartesian/PointS3.h>
#include <CGAL/SimpleCartesian/PlaneS3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
PointS3<FT>
midpoint( PointS3<FT> const& p, PointS3<FT> const& q )
{
  FT x,y,z;
  midpointC3(p.x(),p.y(),p.z(),q.x(),q.y(),q.z(),x,y,z);
  return PointS3<FT>(x,y,z);
}

template < class FT >
PointS3<FT>
circumcenter( PointS3<FT> const& p,
              PointS3<FT> const& q,
              PointS3<FT> const& r,
              PointS3<FT> const& s)
{
  FT x,y,z;
  circumcenterC3(p.x(),p.y(),p.z(),
                 q.x(),q.y(),q.z(),
                 r.x(),r.y(),r.z(),
                 s.x(),s.y(),s.z(),
                 x,y,z);
  return PointS3<FT>(x,y,z);
}

template <class FT>
CGAL_KERNEL_LARGE_INLINE
PointS3<FT>
projection(const PointS3<FT>& p, const PlaneS3<FT>& h)
{
  FT x,y,z;
  projection_planeC3(h.a(),h.b(),h.c(),h.d(),
                     p.x(),p.y(),p.z(),
                     x,y,z);
  return PointS3<FT>(x,y,z);
}


template < class FT >
inline
FT
squared_distance(const PointS3<FT>& p, const PointS3<FT>& q)
{
  return squared_distanceC3(p.x(),p.y(),p.z(),q.x(),q.y(),q.z());
}

template < class FT >
inline
FT
scaled_distance_to_plane(const PlaneS3<FT>& h, const PointS3<FT>& p)
{
  return scaled_distance_to_planeC3(h.a(),h.b(),h.c(),h.d(),
                                    p.x(),p.y(),p.z());
}

template < class FT >
inline
FT
scaled_distance_to_plane(const PointS3<FT>& hp,
                         const PointS3<FT>& hq,
                         const PointS3<FT>& hr,
                         const PointS3<FT>& p)
{
  return scaled_distance_to_planeC3(hp.x(),hp.y(),hp.z(),
                                    hq.x(),hq.y(),hq.z(),
                                    hr.x(),hr.y(),hr.z(),
                                    p.x(),p.y(),p.z());
}

template < class FT >
inline
PlaneS3<FT>
bisector(const PointS3<FT>& p, const PointS3<FT>& q)
{ return PlaneS3<FT>( midpoint(p,q), q-p); }

template < class FT >
inline
PointS3<FT>
gp_linear_intersection(const PlaneS3<FT>& f,
                       const PlaneS3<FT>& g,
                       const PlaneS3<FT>& h)
{
  return PointS3<FT>( det3x3_by_formula(-f.d(), f.b(), f.c(),
                                        -g.d(), g.b(), g.c(),
                                        -h.d(), h.b(), h.c()),
                      det3x3_by_formula( f.a(),-f.d(), f.c(),
                                         g.a(),-g.d(), g.c(),
                                         h.a(),-h.d(), h.c()),
                      det3x3_by_formula( f.a(), f.b(),-f.d(),
                                         g.a(), g.b(),-g.d(),
                                         h.a(), h.b(),-h.d()),
                      det3x3_by_formula( f.a(), f.b(), f.c(),
                                         g.a(), g.b(), g.c(),
                                         h.a(), h.b(), h.c()));
}

CGAL_END_NAMESPACE

#endif // CGAL_BASIC_CONSTRUCTIONSS3_H
