// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//

// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/basic_constructions_3.h
// source        : include/CGAL/Cartesian/basic_constructions_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_BASIC_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_BASIC_CONSTRUCTIONS_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#ifndef CGAL_CARTESIAN_POINT_3_H
#include <CGAL/Cartesian/Point_3.h>
#endif // CGAL_CARTESIAN_POINT_3_H
#ifndef CGAL_CARTESIAN_VECTOR_3_H
#include <CGAL/Cartesian/Vector_3.h>
#endif // CGAL_CARTESIAN_VECTOR_3_H
#ifndef CGAL_CARTESIAN_PLANE_3_H
#include <CGAL/Cartesian/Plane_3.h>
#endif // CGAL_CARTESIAN_PLANE_3_H
#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTC3_H
#include <CGAL/constructions/kernel_ftC3.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
VectorC3<R CGAL_CTAG>
cross_product(const VectorC3<R CGAL_CTAG>& v,
              const VectorC3<R CGAL_CTAG>& w)
{
    return VectorC3<R CGAL_CTAG>( v.y() * w.z() - v.z() * w.y() ,
                         v.z() * w.x() - v.x() * w.z() ,
                         v.x() * w.y() - v.y() * w.x() );
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
PointC3<R CGAL_CTAG>
midpoint(PointC3<R CGAL_CTAG> const& p,
         PointC3<R CGAL_CTAG> const& q )
{
  typename R::FT x,y,z;
  midpointC3(p.x(),p.y(),p.z(),q.x(),q.y(),q.z(),x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

template < class R >
PointC3<R CGAL_CTAG>
circumcenter( PointC3<R CGAL_CTAG> const& p,
              PointC3<R CGAL_CTAG> const& q,
              PointC3<R CGAL_CTAG> const& r,
              PointC3<R CGAL_CTAG> const& s)
{
  typename R::FT x,y,z;
  circumcenterC3(p.x(),p.y(),p.z(),
                 q.x(),q.y(),q.z(),
                 r.x(),r.y(),r.z(),
                 s.x(),s.y(),s.z(),
                 x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
projection(const PointC3<R CGAL_CTAG>& p,
           const PlaneC3<R CGAL_CTAG>& h)
{
  typename R::FT x,y,z;
  projectionC3(h.a(),h.b(),h.c(),h.d(),
               p.x(),p.y(),p.z(),
               x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}


template < class R >
inline
typename R::FT
squared_distance(const PointC3<R CGAL_CTAG> &p,
                 const PointC3<R CGAL_CTAG> &q)
{
  return squared_distanceC3(p.x(),p.y(),p.z(),q.x(),q.y(),q.z());
}

template < class R >
inline
typename R::FT
scaled_distance_to_plane(const PlaneC3<R CGAL_CTAG> &h,
                         const PointC3<R CGAL_CTAG> &p)
{
  return scaled_distance_to_planeC3(h.a(),h.b(),h.c(),h.d(),
                                    p.x(),p.y(),p.z());
}

template < class R >
inline
typename R::FT
scaled_distance_to_plane(const PointC3<R CGAL_CTAG> &hp,
                         const PointC3<R CGAL_CTAG> &hq,
                         const PointC3<R CGAL_CTAG> &hr,
                         const PointC3<R CGAL_CTAG> &p)
{
  return scaled_distance_to_planeC3(hp.x(),hp.y(),hp.z(),
                                    hq.x(),hq.y(),hq.z(),
                                    hr.x(),hr.y(),hr.z(),
                                    p.x(),p.y(),p.z());
}


CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_BASIC_CONSTRUCTIONS_3_H
