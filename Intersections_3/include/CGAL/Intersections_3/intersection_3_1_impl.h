// Copyright (c) 1997-2010
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//                 Sebastien Loriot <Sebastien.Loriot@geometryfactory.com>


#include <CGAL/wmult.h>
#include <boost/next_prior.hpp>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {

  template <class K>
  class Plane_3;

  template <class K>
  class Line_3;

  template <class K>
  class Segment_3;

  template <class K>
  class Ray_3;

  template <class K>
  class Sphere_3;

  template <class K>
  class Triangle_3;

  template <class K>
  class Iso_cuboid_3;

// the special plane_3 function
template <class K>
inline 
#if CGAL_INTERSECTION_VERSION < 2
CGAL::Object
#else
typename cpp11::result_of<typename K::Intersect_3(typename K::Plane_3, typename K::Plane_3, typename K::Plane_3)>::type
#endif
intersection(const Plane_3<K> &plane1, const Plane_3<K> &plane2,
             const Plane_3<K> &plane3)
{
  return K().intersect_3_object()(plane1, plane2, plane3);
}

CGAL_INTERSECTION_FUNCTION(Plane_3, Line_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Plane_3, Line_3, 3)

CGAL_INTERSECTION_FUNCTION_SELF(Plane_3, 3)

CGAL_INTERSECTION_FUNCTION_SELF(Line_3, 3)
CGAL_DO_INTERSECT_FUNCTION_SELF(Line_3, 3)

CGAL_INTERSECTION_FUNCTION_SELF(Segment_3, 3)
CGAL_DO_INTERSECT_FUNCTION_SELF(Segment_3, 3)

CGAL_INTERSECTION_FUNCTION(Line_3, Segment_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Line_3, Segment_3, 3)

CGAL_INTERSECTION_FUNCTION(Line_3, Ray_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Line_3, Ray_3, 3)

CGAL_INTERSECTION_FUNCTION(Segment_3, Ray_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Segment_3, Ray_3, 3)

CGAL_INTERSECTION_FUNCTION_SELF(Ray_3, 3)
CGAL_DO_INTERSECT_FUNCTION_SELF(Ray_3, 3)

CGAL_INTERSECTION_FUNCTION(Plane_3, Sphere_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Plane_3, Sphere_3, 3)

CGAL_INTERSECTION_FUNCTION_SELF(Sphere_3, 3)
CGAL_DO_INTERSECT_FUNCTION_SELF(Sphere_3, 3)

CGAL_INTERSECTION_FUNCTION(Plane_3, Ray_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Plane_3, Ray_3, 3)

CGAL_INTERSECTION_FUNCTION(Plane_3, Segment_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Plane_3, Segment_3, 3)

CGAL_INTERSECTION_FUNCTION(Plane_3, Triangle_3, 3)

template <class K>
inline typename
cpp11::result_of<typename K::Intersect_3(typename K::Line_3, Bbox_3)>::type
intersection(const Line_3<K> &a,
	     const Bbox_3 &b) {
  return K().intersect_3_object()(a, b);
}

template <class K>
inline typename
cpp11::result_of<typename K::Intersect_3(typename K::Line_3, Bbox_3)>::type
intersection(const Bbox_3 &a,
             const Line_3<K> &b) {
  return K().intersect_3_object()(a, b);
}

template <class K>
inline typename
cpp11::result_of<typename K::Intersect_3(typename K::Ray_3, Bbox_3)>::type
intersection(const Ray_3<K> &a,
	     const Bbox_3 &b) {
  return K().intersect_3_object()(a, b);
}

template <class K>
inline typename
cpp11::result_of<typename K::Intersect_3(typename K::Ray_3, Bbox_3)>::type
intersection(const Bbox_3 &a,
             const Ray_3<K> &b) {
  return K().intersect_3_object()(a, b);
}

template <class K>
inline typename
cpp11::result_of<typename K::Intersect_3(typename K::Segment_3, Bbox_3)>::type
intersection(const Segment_3<K> &a,
	     const Bbox_3 &b) {
  return K().intersect_3_object()(a, b);
}

template <class K>
inline typename
cpp11::result_of<typename K::Intersect_3(typename K::Segment_3, Bbox_3)>::type
intersection(const Bbox_3 &a,
             const Segment_3<K> &b) {
  return K().intersect_3_object()(a, b);
}

CGAL_INTERSECTION_FUNCTION(Line_3, Iso_cuboid_3, 3)

CGAL_INTERSECTION_FUNCTION(Ray_3, Iso_cuboid_3, 3)

CGAL_INTERSECTION_FUNCTION(Segment_3, Iso_cuboid_3, 3)

CGAL_INTERSECTION_FUNCTION_SELF(Iso_cuboid_3, 3)

CGAL_DO_INTERSECT_FUNCTION_SELF(Plane_3, 3)

template <class R>
inline bool
do_intersect(const Plane_3<R> &plane1, const Plane_3<R> &plane2,
             const Plane_3<R> &plane3) {
  return R().do_intersect_3_object()(plane1, plane2, plane3);
}

CGAL_DO_INTERSECT_FUNCTION_SELF(Iso_cuboid_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Iso_cuboid_3, Line_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Iso_cuboid_3, Ray_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Iso_cuboid_3, Segment_3, 3)

namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
intersection(const typename K::Plane_3  &plane, 
	     const typename K::Line_3 &line, 
	     const K& /*k*/)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Direction_3 Direction_3;
    typedef typename K::RT RT;

    const Point_3 &line_pt = line.point();
    const Direction_3 &line_dir = line.direction();

    RT num = plane.a()*line_pt.hx() + plane.b()*line_pt.hy()
           + plane.c()*line_pt.hz() + wmult_hw((K*)0, plane.d(), line_pt);
    RT den = plane.a()*line_dir.dx() + plane.b()*line_dir.dy()
           + plane.c()*line_dir.dz();
    if (den == 0) {
        if (num == 0) {
            // all line
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Line_3>(line);
        } else {
            // no intersection
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Line_3>();
        }
    }
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Line_3>(Point_3(
        den*line_pt.hx()-num*line_dir.dx(),
        den*line_pt.hy()-num*line_dir.dy(),
        den*line_pt.hz()-num*line_dir.dz(),
        wmult_hw((K*)0, den, line_pt)));
}

template <class K>
inline
typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
intersection(const typename K::Line_3 &line, 
	     const typename K::Plane_3  &plane, 
	     const K& k)
{
  return intersection(plane, line, k);
}

template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Plane_3>::result_type
intersection(const typename K::Plane_3 &plane1, 
	     const typename K::Plane_3 &plane2, 
	     const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Direction_3 Direction_3;
  typedef typename K::Line_3 Line_3;

    typedef typename K::RT RT;
    const RT &a = plane1.a();
    const RT &b = plane1.b();
    const RT &c = plane1.c();
    const RT &d = plane1.d();
    const RT &p = plane2.a();
    const RT &q = plane2.b();
    const RT &r = plane2.c();
    const RT &s = plane2.d();

    RT det = a*q-p*b;
    if (det != 0) {
        Point_3 is_pt = Point_3(b*s-d*q, p*d-a*s, 0, det);
        Direction_3 is_dir = Direction_3(b*r-c*q, p*c-a*r, det);
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(Line_3(is_pt, is_dir));
    }
    det = a*r-p*c;
    if (det != 0) {
        Point_3 is_pt = Point_3(c*s-d*r, 0, p*d-a*s, det);
        Direction_3 is_dir = Direction_3(c*q-b*r, det, p*b-a*q);
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(Line_3(is_pt, is_dir));
    }
    det = b*r-c*q;
    if (det != 0) {
        Point_3 is_pt = Point_3(0, c*s-d*r, d*q-b*s, det);
        Direction_3 is_dir = Direction_3(det, c*p-a*r, a*q-b*p);
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(Line_3(is_pt, is_dir));
    }
// degenerate case
    if (a!=0 || p!=0) {
        if (a*s == p*d)
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(plane1);
        else
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>();
    }
    if (b!=0 || q!=0) {
        if (b*s == q*d)
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(plane1);
        else
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>();
    }
    if (c!=0 || r!=0) {
        if (c*s == r*d)
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(plane1);
        else
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>();
    }
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Plane_3>(plane1);
}

template <class K>
#if CGAL_INTERSECTION_VERSION < 2
CGAL::Object
#else
boost::optional< boost::variant<typename K::Point_3,
                                typename K::Line_3,
                                typename K::Plane_3> >
#endif
intersection(const typename K::Plane_3 &plane1,
	     const typename K::Plane_3 &plane2,
	     const typename K::Plane_3 &plane3,
	     const K& k)
{
    #if CGAL_INTERSECTION_VERSION > 1
    typedef 
      typename boost::optional< 
      boost::variant<typename K::Point_3,
                     typename K::Line_3,
                     typename K::Plane_3> > 
    result_type;
    #endif


    typedef typename K::Point_3      Point_3;
    typedef typename K::Line_3       Line_3;
    typedef typename K::Plane_3      Plane_3;

    // Intersection between plane1 and plane2 can either be
    // a line, a plane, or empty.
    typename Intersection_traits<K, Plane_3, Plane_3>::result_type
      o12 = internal::intersection(plane1, plane2, k);
    
    if(o12) {
      if(const Line_3* l = intersect_get<Line_3>(o12)) {
        // either point or line
        typename Intersection_traits<K, Plane_3, Line_3>::result_type 
          v = internal::intersection(plane3, *l, k);
        if(v) {
          if(const Point_3* p = intersect_get<Point_3>(v))
            #if CGAL_INTERSECTION_VERSION < 2
            return make_object(*p);
            #else
            return result_type(*p);
            #endif
          else if(const Line_3* l = intersect_get<Line_3>(v))
            #if CGAL_INTERSECTION_VERSION < 2
            return make_object(*l);
            #else
            return result_type(*l);
            #endif
        }
      } else if(const Plane_3 *pl = intersect_get<Plane_3>(o12)) {
        // either line or plane
        typename Intersection_traits<K, Plane_3, Plane_3>::result_type 
          v = internal::intersection(plane3, *pl, k);
        if(v) {
          if(const Plane_3* p = intersect_get<Plane_3>(v))
            #if CGAL_INTERSECTION_VERSION < 2
            return make_object(*p);
            #else
            return result_type(*p);
            #endif
          else if(const Line_3* l = intersect_get<Line_3>(v))
            #if CGAL_INTERSECTION_VERSION < 2
            return make_object(*l);
            #else
            return result_type(*l);
            #endif
        }
      }
    }
    
    #if CGAL_INTERSECTION_VERSION < 2
    return Object();
    #else
    return result_type();
    #endif
}


template <class K>
bool
do_intersect(const typename K::Plane_3 &plane, 
	     const typename K::Line_3 &line,
	     const K&)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Direction_3 Direction_3;
    typedef typename K::RT RT;
    const Point_3 &line_pt = line.point();
    const Direction_3 &line_dir = line.direction();

    RT den = plane.a()*line_dir.dx() + plane.b()*line_dir.dy()
           + plane.c()*line_dir.dz();
    if (den != 0)
        return true;
    RT num = plane.a()*line_pt.hx() + plane.b()*line_pt.hy()
           + plane.c()*line_pt.hz() + wmult_hw((K*)0, plane.d(), line_pt);
    if (num == 0) {
        // all line
        return true;
    } else {
        // no intersection
        return false;
    }
}

template <class K>
inline
bool
do_intersect(const typename K::Line_3 &line, 
	     const typename K::Plane_3 &plane, 
	     const K& k)
{
  return do_intersect(plane, line, k);
}

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Line_3>::result_type
intersection(const typename K::Line_3 &l1,
	     const typename K::Line_3 &l2,
	     const K&)
{
  typedef typename K::FT           FT;
  typedef typename K::Point_3      Point_3;
  typedef typename K::Vector_3     Vector_3;

  if(K().has_on_3_object()(l1, l2.point())) {
    const Vector_3& v1 = l1.to_vector();
    const Vector_3& v2 = l2.to_vector();
    if((v1.x() * v2.y() == v1.y() * v2.x()) &&
       (v1.x() * v2.z() == v1.z() * v2.x()) &&
       (v1.y() * v2.z() == v1.z() * v2.y()))
      return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>(l1);
  }
  
  if(K().are_parallel_3_object()(l1,l2)) return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>();
  const Point_3 &p1 = l1.point();
  const Point_3 &p3 = l2.point();
  const Vector_3 &v1 = l1.to_vector();
  const Vector_3 &v2 = l2.to_vector();
  const Point_3 p2 = p1 + v1;
  const Point_3 p4 = p2 + v2;
  if(!K().coplanar_3_object()(p1,p2,p3,p4)) return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>();
  const Vector_3 v3 = p3 - p1;
 const Vector_3 v3v2 = cross_product(v3,v2);
  const Vector_3 v1v2 = cross_product(v1,v2);
  const FT t = ((v3v2.x()*v1v2.x()) + (v3v2.y()*v1v2.y()) + (v3v2.z()*v1v2.z())) /
               (v1v2.squared_length());
  return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>(p1 + (v1 * t));
}

template <class K>
bool
do_intersect(const typename K::Line_3 &l1,
	     const typename K::Line_3 &l2,
	     const K&)
{
  typedef typename K::Point_3      Point_3;
  typedef typename K::Vector_3     Vector_3;

  if(K().has_on_3_object()(l1, l2.point())) return true;
  if(K().are_parallel_3_object()(l1,l2)) return false;
  const Point_3 &p1 = l1.point();
  const Point_3 &p3 = l2.point();
  const Vector_3 &v1 = l1.to_vector();
  const Vector_3 &v2 = l2.to_vector();
  const Point_3 p2 = p1 + v1;
  const Point_3 p4 = p2 + v2;
  return K().coplanar_3_object()(p1,p2,p3,p4);
}

template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Segment_3>::result_type
intersection_collinear_segments(const typename K::Segment_3 &s1,
                                const typename K::Segment_3 &s2,
                                const K& k)
{
  CGAL_precondition(! s1.is_degenerate () && ! s2.is_degenerate () );

  const typename K::Point_3& p=s1[0],q=s1[1],r=s2[0],s=s2[1];
  typename K::Collinear_are_ordered_along_line_3 cln_order=k.collinear_are_ordered_along_line_3_object();
  
  if ( cln_order(p,r,q) ){
    if ( cln_order(p,s,q) )
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(s2);
    if ( cln_order(r,p,s) ){
      if (r!=p) return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>( typename K::Segment_3(r,p) );
      if ( cln_order(r,q,s) ) return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(s1);
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(p);
    }
    return r!=q ? intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>( typename K::Segment_3(r,q) ) 
      : intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(q);
  }

  if ( cln_order(p,s,q) ){
    if ( cln_order(r,p,s) ){
      if (s!=p) return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>( typename K::Segment_3(s,p) );
      if (cln_order(r,q,s)) return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(s1);  
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(p);
    }
   return s!=q ? intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>( typename K::Segment_3(s,q) ) 
     : intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(q);
  }
  
  if ( cln_order(r,p,s) )
    return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(s1); 
  return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>();
}

template<class K>
struct L_p_visitor : public boost::static_visitor< 
  typename Intersection_traits<K, typename K::Segment_3, 
                               typename K::Segment_3>::result_type
  > 
{

  typedef typename Intersection_traits<K, typename K::Segment_3, 
                                         typename K::Segment_3>::result_type result_type;
  L_p_visitor(const typename K::Segment_3& s1, const typename K::Segment_3& s2) :
    s1(s1), s2(s2) { }
  const typename K::Segment_3& s1;
  const typename K::Segment_3& s2;

  result_type
  operator()(const typename K::Point_3& p) const {
    typename K::Collinear_are_ordered_along_line_3 cln_order=K().collinear_are_ordered_along_line_3_object();
    if ( cln_order(s1[0],p,s1[1]) && cln_order(s2[0],p,s2[1]) )
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>(p);
    else
      return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>();
}

  result_type
  operator()(const typename K::Line_3&) const {
    return intersection_collinear_segments(s1,s2,K());
  }
};

template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Segment_3>::result_type
intersection(const typename K::Segment_3 &s1,
	     const typename K::Segment_3 &s2,
	     const K&)
{
  CGAL_precondition(! s1.is_degenerate () && ! s2.is_degenerate () );

  typename Intersection_traits<K, typename K::Line_3, typename K::Line_3>::result_type 
    v = internal::intersection(s1.supporting_line(),s2.supporting_line(), K());

  if(v) {
    #if CGAL_INTERSECTION_VERSION < 2
    // abuse the visitor to do the visitation manually
    L_p_visitor<K> visitor(s1, s2);
    if(const typename K::Point_3* p = object_cast<typename K::Point_3>(&v))
      return visitor(*p);
    if(const typename K::Line_3* l = object_cast<typename K::Line_3>(&v))
      return visitor(*l);
    #else
    return apply_visitor(L_p_visitor<K>(s1, s2) , *v);
    #endif
  }
  return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Segment_3>();
}

template <class K>
inline
bool
do_intersect(const typename K::Segment_3  &s1,
             const typename K::Segment_3  &s2,
             const K & k)
{
  CGAL_precondition(! s1.is_degenerate () && ! s2.is_degenerate () );
  bool b=do_intersect(s1.supporting_line(),s2.supporting_line(),k);
  if (b)
  {
    //supporting_line intersects: points are coplanar
    typename K::Coplanar_orientation_3 cpl_orient=k.coplanar_orientation_3_object();
    ::CGAL::Orientation or1 =  cpl_orient(s1[0],s1[1],s2[0]);
    ::CGAL::Orientation or2 =  cpl_orient(s1[0],s1[1],s2[1]);
    
    if ( or1 == COLLINEAR && or2 ==COLLINEAR )
    {
      //segments are collinear
      typename K::Collinear_are_ordered_along_line_3 cln_order=k.collinear_are_ordered_along_line_3_object();
      return cln_order(s1[0],s2[0],s1[1]) || 
             cln_order(s1[0],s2[1],s1[1]) ||
             cln_order(s2[0],s1[0],s2[1]) ;
    }
    
    if ( or1 != or2 ){
      or1=cpl_orient(s2[0],s2[1],s1[0]);
      return (or1 == COLLINEAR || or1 != cpl_orient(s2[0],s2[1],s1[1]));
    }
  }
  return false;
}

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Segment_3>::result_type
intersection(const typename K::Line_3 &l,
	     const typename K::Segment_3 &s,
	     const K& k)
{
  CGAL_precondition(! l.is_degenerate () && ! s.is_degenerate () );
  
  typename Intersection_traits<K, typename K::Line_3, typename K::Line_3>::result_type 
    v = internal::intersection(l,s.supporting_line(), K());

  if(v) {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3> (v)) {
    typename K::Collinear_are_ordered_along_line_3 cln_order=k.collinear_are_ordered_along_line_3_object();
      if(cln_order(s[0],*p,s[1])) 
        return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Segment_3>(*p);
    } else {
      return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Segment_3>(s);
  }
  }
  
  return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Segment_3>();
}

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Segment_3>::result_type
intersection(const typename K::Segment_3 &s,
	     const typename K::Line_3 &l,
	     const K& k)
{
  return intersection(l,s,k);
}

template <class K>
inline
bool
do_intersect(const typename K::Line_3  &l,
             const typename K::Segment_3  &s,
             const K & k)
{
  CGAL_precondition(! l.is_degenerate () && ! s.is_degenerate () );
  bool b=do_intersect(l,s.supporting_line(),k);
  if (b)
  {
    //supporting_line intersects: points are coplanar
    typename K::Coplanar_orientation_3 cpl_orient=k.coplanar_orientation_3_object();
    typename K::Point_3 p1=l.point(0);
    typename K::Point_3 p2=l.point(1);
    ::CGAL::Orientation or1 =  cpl_orient(p1,p2,s[0]);
       
    if ( or1 == COLLINEAR ) return true;
    
    ::CGAL::Orientation or2 =  cpl_orient(p1,p2,s[1]);
    return or1!=or2;
  }
  return false;
}

template <class K>
inline
bool
do_intersect(const typename K::Segment_3  &s,
             const typename K::Line_3  &l,
             const K & k)
{
  return do_intersect(l,s,k);
}

template <class K>
bool
Ray_3_has_on_collinear_Point_3(
            const typename K::Ray_3 &r,
	    const typename K::Point_3 &p,
	    const K& k)
{
  return
      k.equal_3_object()(r.source(),p)
    ||
      k.equal_3_object() (
        k.construct_direction_3_object()( k.construct_vector_3_object() (r.source(),p) ),
        r.direction()
      );
}

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Ray_3>::result_type
intersection(const typename K::Line_3 &l,
	     const typename K::Ray_3 &r,
	     const K& k)
{
  CGAL_precondition(! l.is_degenerate () && ! r.is_degenerate () );

  typename Intersection_traits<K, typename K::Line_3, typename K::Line_3>::result_type 
    v = internal::intersection(l,r.supporting_line(), k);
  
  if(v) {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v)) {
      if( Ray_3_has_on_collinear_Point_3(r,*p,k) ) 
        return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Ray_3>(*p);
    } else if(intersect_get<typename K::Line_3>(v)) {
      return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Ray_3>(r);
  }
  }
  return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Ray_3>();
}

template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Line_3>::result_type
intersection(const typename K::Ray_3 &r,
	     const typename K::Line_3 &l,
	     const K& k)
{
  return intersection(l,r,k);
}

template <class K>
inline
bool
do_intersect(const typename K::Line_3  &l,
             const typename K::Ray_3  &r,
             const K & k)
{
  CGAL_precondition(! l.is_degenerate () && ! r.is_degenerate () );
  if ( !do_intersect(l,r.supporting_line()) ) return false;
  typename K::Coplanar_orientation_3 pred=k.coplanar_orientation_3_object();
  Orientation p0p1s=pred(l.point(0),l.point(1),r.source());
  if ( p0p1s == COLLINEAR) return true;
  Orientation stp0 =pred(r.source(),r.second_point(),l.point(0));
  if ( stp0 == COLLINEAR )
    return Ray_3_has_on_collinear_Point_3(r,l.point(0),k);
  return p0p1s!=stp0;
}

template <class K>
inline
bool
do_intersect(const typename K::Ray_3  &r,
             const typename K::Line_3  &l,
             const K & k)
{
  return do_intersect(l,r,k);
}

template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Ray_3>::result_type
intersection(const typename K::Segment_3 &s,
	     const typename K::Ray_3 &r,
	     const K& k)
{
  CGAL_precondition(! s.is_degenerate () && ! r.is_degenerate () );

  typename Intersection_traits<K, typename K::Line_3, typename K::Segment_3>::result_type 
    v = internal::intersection(r.supporting_line(),s, K());

  if(v) {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v)) {
      if( Ray_3_has_on_collinear_Point_3(r,*p,k) ) 
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(*p);
    } else if(const typename K::Segment_3* s2 = intersect_get<typename K::Segment_3>(v)) {
      bool has_source=Ray_3_has_on_collinear_Point_3(r,s.source(),k);
      bool has_target=Ray_3_has_on_collinear_Point_3(r,s.target(),k);
      if (has_source){
        if (has_target)
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(*s2);
        else
        {
          if (k.equal_3_object() (r.source(),s.source()))
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(r.source());
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(
              k.construct_segment_3_object()(r.source(),s.source()));
        }
      }
      else{
        if (has_target){
          if (k.equal_3_object() (r.source(),s.target()))
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(r.source());          
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(
              k.construct_segment_3_object()(r.source(),s.target()));
        }
      }
    }
  }

  return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>();
}

template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Segment_3>::result_type
intersection(const typename K::Ray_3 &r,
	     const typename K::Segment_3 &s,
	     const K& k)
{
  return intersection(s,r,k);
}

template <class K>
inline
bool
do_intersect(const typename K::Segment_3  &s,
             const typename K::Ray_3  &r,
             const K & k)
{
  CGAL_precondition(! s.is_degenerate () && ! r.is_degenerate () );
  if ( !do_intersect(s,r.supporting_line()) ) return false;
  typename K::Coplanar_orientation_3 pred=k.coplanar_orientation_3_object();
  Orientation p0p1s=pred(s.point(0),s.point(1),r.source());
  Orientation stp0 =pred(r.source(),r.second_point(),s.point(0));
  if ( p0p1s == COLLINEAR) //s belongs to the supporting line of p0p1
  {
    if ( stp0 == COLLINEAR )//st and p0p1 have the same supporting line
      return Ray_3_has_on_collinear_Point_3(r,s.point(0),k) || Ray_3_has_on_collinear_Point_3(r,s.point(1),k);
    else
      return true;
  }
  if ( stp0 == COLLINEAR )
    return Ray_3_has_on_collinear_Point_3(r,s.point(0),k);
  return p0p1s!=stp0;
}

template <class K>
inline
bool
do_intersect(const typename K::Ray_3  &r,
             const typename K::Segment_3  &s,
             const K & k)
{
  return do_intersect(s,r,k);
}

template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Ray_3>::result_type
intersection(const typename K::Ray_3 &r1,
	     const typename K::Ray_3 &r2,
	     const K& k)
{
  CGAL_precondition(! r1.is_degenerate () && ! r2.is_degenerate () );

  typename Intersection_traits<K, typename K::Line_3, typename K::Ray_3>::result_type 
    v = internal::intersection(r1.supporting_line(),r2, k);

  if(v) {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v)) {
      if(Ray_3_has_on_collinear_Point_3(r1,*p,k)) 
        return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(*p);
    } else if(const typename K::Ray_3* r = intersect_get<typename K::Ray_3>(v)) {
      bool r1_has_s2=Ray_3_has_on_collinear_Point_3(r1,r2.source(),k);
      bool r2_has_s1=Ray_3_has_on_collinear_Point_3(r2,r1.source(),k);
      if (r1_has_s2){
        if (r2_has_s1)
        {
          if (k.equal_3_object()(r1.source(),r2.source()))
            return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(r1.source());
          else {
            return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(
              k.construct_segment_3_object()(r1.source(),r2.source()));          
          }
        }
        else
          return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(*r);
      }
      else{
        if (r2_has_s1)
          return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>(r1);
      }
    }
  }
  return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Ray_3>();

}

template <class K>
inline
bool
do_intersect(const typename K::Ray_3  &r1,
             const typename K::Ray_3  &r2,
             const K & k)
{
  CGAL_precondition(! r1.is_degenerate () && ! r2.is_degenerate () );
  if ( !do_intersect(r1,r2.supporting_line()) ) return false;
  typename K::Coplanar_orientation_3 pred=k.coplanar_orientation_3_object();
  Orientation p0p1s=pred(r1.point(0),r1.point(1),r2.source());
  Orientation stp0 =pred(r2.source(),r2.second_point(),r1.point(0));
  
  if ( p0p1s == COLLINEAR){
    if(stp0 == COLLINEAR ) 
      return  Ray_3_has_on_collinear_Point_3(r2,r1.source(),k) ||
              Ray_3_has_on_collinear_Point_3(r1,r2.source(),k);
    else
      return true;
  }
  if(stp0 == COLLINEAR )
    return Ray_3_has_on_collinear_Point_3(r2,r1.point(0),k);
  return p0p1s!=stp0;
}

template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Sphere_3>::result_type
intersection(const typename K::Plane_3 &p,
             const typename K::Sphere_3 &s,
             const K&)
{
  typedef typename K::Circle_3 Circle_3;
  typedef typename K::Point_3 Point_3;
  typedef typename K::FT FT;
  const FT d2 = CGAL::square(p.a()*s.center().x() + 
                             p.b()*s.center().y() + 
                             p.c()*s.center().z() + p.d()) /
      (square(p.a()) + square(p.b()) + square(p.c()));
  const FT cmp = d2 - s.squared_radius();
  if(CGAL_NTS is_zero(cmp)) { // tangent
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Sphere_3>(p.projection(s.center()));
  } else if(CGAL_NTS is_negative(cmp)) { // intersect
    Point_3 center = p.projection(s.center());
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Sphere_3>(Circle_3(center,s.squared_radius() - d2,p));
  } // do not intersect
  return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Sphere_3>();
}

template <class K>
inline
bool
do_intersect(const typename K::Plane_3 &p,
             const typename K::Sphere_3 &s,
             const K&)
{
  typedef typename K::FT FT;
  const FT d2 = CGAL::square(p.a()*s.center().x() + 
                             p.b()*s.center().y() + 
                             p.c()*s.center().z() + p.d()) /
      (square(p.a()) + square(p.b()) + square(p.c()));
  return d2 <= s.squared_radius();
}

template <class K>
inline
bool
do_intersect(const typename K::Sphere_3 &s,
             const typename K::Plane_3 &p,
             const K&)
{
  return do_intersect(p,s);
}


template <class K>
inline
typename Intersection_traits<K, typename K::Sphere_3, typename K::Plane_3>::result_type
intersection(const typename K::Sphere_3 &s,
             const typename K::Plane_3 &p,
             const K& k)
{
  return intersection(p, s, k);
}

template <class K>
inline
typename Intersection_traits<K, typename K::Sphere_3, typename K::Sphere_3>::result_type
intersection(const typename K::Sphere_3 &s1,
             const typename K::Sphere_3 &s2,
             const K& k)
{
  typedef typename K::Plane_3 Plane_3;
  if(s1.center() == s2.center()) {
    if(s1.squared_radius() == s2.squared_radius()) {
      if(is_zero(s1.squared_radius())) return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>(s1.center());
      else return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>(s1);
    } else return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>();  // cocentrics
  }
  Plane_3 p = K().construct_radical_plane_3_object()(s1,s2);
  
  
  typename Intersection_traits<K, typename K::Sphere_3, typename K::Plane_3>::result_type 
    v = intersection(p, s1, k);


  if(v) {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
      return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>(*p);
    else if(const typename K::Circle_3* c = intersect_get<typename K::Circle_3>(v))
      return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>(*c);
  }

  return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>();
}

template <class K>
inline
bool
do_intersect(const typename K::Sphere_3 &s1,
             const typename K::Sphere_3 &s2,
             const K& k)
{
  typedef typename K::Plane_3 Plane_3;
  if(s1.center() == s2.center()) {
    return s1.squared_radius() == s2.squared_radius();
  }
  Plane_3 p = K().construct_radical_plane_3_object()(s1,s2);
  return do_intersect(p, s1, k);
}

template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Ray_3>::result_type
intersection(const typename K::Plane_3 &plane, 
	     const typename K::Ray_3 &ray, 
	     const K& k)
{
    typedef typename K::Point_3 Point_3;

    typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type 
      v = internal::intersection(plane, ray.supporting_line(), k);

    if(v) {
      if(const Point_3* p = intersect_get<Point_3>(v)) {
        if (ray.collinear_has_on(*p))
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Ray_3>(*p);
        else
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Ray_3>();
    }
    } else {
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Ray_3>();
    }

    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Ray_3>(ray);
}


template <class K>
inline
typename Intersection_traits<K, typename K::Ray_3, typename K::Plane_3>::result_type
intersection(const typename K::Ray_3 &ray, 
	     const typename K::Plane_3 &plane, 
	     const K& k)
{
  return intersection(plane, ray, k);
}



template <class K>
bool
do_intersect(const typename K::Plane_3 &plane, 
	     const typename K::Ray_3 &ray, 
	     const K& k)
{
    typedef typename K::Point_3 Point_3;
    
    typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>
      ::result_type 
      line_intersection = internal::intersection(plane, ray.supporting_line(), k);

    if(!line_intersection)
        return false;
    if(const Point_3 *isp = intersect_get<Point_3>(line_intersection))
        return ray.collinear_has_on(*isp);
    
    return true;
}


template <class K>
inline
bool
do_intersect(const typename K::Ray_3 &ray, 
	     const typename K::Plane_3 &plane, 
	     const K& k)
{
  return do_intersect(plane, ray, k);
}


template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Segment_3>::result_type
intersection(const typename K::Plane_3 &plane, 
	     const typename K::Segment_3 &seg, 
	     const K& k)
{
    typedef typename K::Point_3 Point_3;
    const Point_3 &source = seg.source();
    const Point_3 &target = seg.target();

    Oriented_side source_side = plane.oriented_side(source);
    Oriented_side target_side = plane.oriented_side(target);

    switch (source_side) {
    case ON_ORIENTED_BOUNDARY:
        if (target_side == ON_ORIENTED_BOUNDARY)
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(seg);
        else
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(source);
    case ON_POSITIVE_SIDE:
        switch (target_side) {
        case ON_ORIENTED_BOUNDARY:
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(target);
        case ON_POSITIVE_SIDE:
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>();
        case ON_NEGATIVE_SIDE:
          { 
            // intersection object should be a point, but rounding errors 
            // could lead to a line. In such case, return seg.
            typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
              v = internal::intersection(plane, seg.supporting_line(), k);
            if(v) {
              if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
                return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(*p);
            else
                return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(seg);
            }
          }
        }
    case ON_NEGATIVE_SIDE:
        switch (target_side) {
        case ON_ORIENTED_BOUNDARY:
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(target);
        case ON_POSITIVE_SIDE:
          { 
            // intersection object should be a point, but rounding errors 
            // could lead to a line. In such case, return seg.
            typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
              v = internal::intersection(plane, seg.supporting_line(), k);
            if(v) {
              if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
                return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(*p);
            else 
                return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(seg);
            }
          }
        case ON_NEGATIVE_SIDE:
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>();
        }
    }
    CGAL_kernel_assertion_msg(false, "Supposedly unreachable code.");
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>();
}


template <class K>
inline
typename Intersection_traits<K, typename K::Segment_3, typename K::Plane_3>::result_type
intersection(const typename K::Segment_3 &seg, 
	     const typename K::Plane_3 &plane, 
	     const K& k)
{
  return intersection(plane, seg, k);
}


template <class K>
bool
do_intersect(const typename K::Plane_3  &plane, 
	     const typename K::Segment_3 &seg, 
	     const K&)
{
    typedef typename K::Point_3 Point_3;
    const Point_3 &source = seg.source();
    const Point_3 &target = seg.target();

    Oriented_side source_side = plane.oriented_side(source);
    Oriented_side target_side = plane.oriented_side(target);

    if ( source_side == target_side
       && target_side != ON_ORIENTED_BOUNDARY) {
        return false;
    }
    return true;
}


template <class K>
inline
bool
do_intersect(const typename K::Segment_3 &seg, 
	     const typename K::Plane_3  &plane, 
	     const K& k)
{
  return do_intersect(plane, seg, k);
}

template <class K>
inline
typename Intersection_traits<K, typename K::Plane_3, typename K::Triangle_3>::result_type
intersection(const typename K::Plane_3 &plane, 
	     const typename K::Triangle_3 &tri, 
	     const K& k)
{
  typedef 
  typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type 
  pl_res;

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();
  
  Oriented_side or0=plane.oriented_side(vertex_on(tri,0));
  Oriented_side or1=plane.oriented_side(vertex_on(tri,1));
  Oriented_side or2=plane.oriented_side(vertex_on(tri,2));
  
  if (or0==ON_ORIENTED_BOUNDARY){
    if (or1==ON_ORIENTED_BOUNDARY){
      if (or2==ON_ORIENTED_BOUNDARY) 
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(tri);
      else 
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(k.construct_segment_3_object()
                                                                                   (tri.vertex(0),tri.vertex(1)));
    }
    else{
      if (or2==ON_ORIENTED_BOUNDARY)
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(k.construct_segment_3_object()
                                                                                   (tri.vertex(0),tri.vertex(2)));
      else{
        if (or1==or2)
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(tri.vertex(0));
        else{
          pl_res v = internal::intersection(plane, k.construct_line_3_object()(tri.vertex(1),tri.vertex(2)), k);
          const typename K::Point_3* p = intersect_get<typename K::Point_3>(v);
          CGAL_kernel_assertion(p!=NULL);
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(k.construct_segment_3_object()
                                                                                     (*p,tri.vertex(0)));
        }
      }
    }
  }

  if (or1==ON_ORIENTED_BOUNDARY){
    if (or2==ON_ORIENTED_BOUNDARY)
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(k.construct_segment_3_object()
                                                                                 (tri.vertex(1),tri.vertex(2)));
    if (or2==or0)
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(tri.vertex(1));
    else{
      pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(0),tri.vertex(2)), k);
      const typename K::Point_3* p = intersect_get<typename K::Point_3>(v);
      CGAL_kernel_assertion(p!=NULL);
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(k.construct_segment_3_object()
                                                                                 (*p,tri.vertex(1)));      
    }
  }
  
  if (or2==ON_ORIENTED_BOUNDARY){
    if (or1==or0)
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(tri.vertex(2));
    else{
      pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(0),tri.vertex(1)), k);
      const typename K::Point_3* p = intersect_get<typename K::Point_3>(v);
      CGAL_kernel_assertion(p!=NULL);
      return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>(k.construct_segment_3_object()
                                                                                 (*p,tri.vertex(2)));      
    }
  }
  
  //triangle vertices are not in the plane
  std::vector<typename K::Point_3> pts;
  pts.reserve(2);
  if (or0!=or1){
    pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(0),tri.vertex(1)), k);
    const typename K::Point_3* pt_ptr = intersect_get<typename K::Point_3>(v);
    CGAL_kernel_assertion( pt_ptr!=NULL );    
    pts.push_back( *pt_ptr );
  }
  if (or0!=or2){
    pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(0),tri.vertex(2)), k);
    const typename K::Point_3* pt_ptr = intersect_get<typename K::Point_3>(v);
    CGAL_kernel_assertion( pt_ptr!=NULL );    
    pts.push_back( *pt_ptr );    
  }
  if (or1!=or2){
    pl_res v = intersection(plane, k.construct_line_3_object()(tri.vertex(1),tri.vertex(2)), k);
    const typename K::Point_3* pt_ptr = intersect_get<typename K::Point_3>(v);
    CGAL_kernel_assertion( pt_ptr!=NULL );    
    pts.push_back( *pt_ptr );
  }
  
  if (pts.empty()) return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>();
  
  CGAL_kernel_assertion(pts.size()==2);
  
  return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Triangle_3>( k.construct_segment_3_object()
                                                                              (*pts.begin(),*boost::prior(pts.end())) );
}

template <class K>
inline
typename Intersection_traits<K, typename K::Triangle_3, typename K::Plane_3>::result_type
intersection(const typename K::Triangle_3 &triangle,
	     const typename K::Plane_3  &plane,
	     const K& k)
{
  return intersection(plane, triangle, k);
}

template <class K>
typename Intersection_traits<K, typename K::Line_3, Bbox_3>::result_type
intersection(const typename K::Line_3 &line,
	     const Bbox_3 &box, 
	     const K&)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Direction_3 Direction_3;
    const Point_3 &linepoint = line.point();
    const Direction_3 &linedir = line.direction();
    return intersection_bl<K>(box,
        CGAL::to_double(linepoint.x()),
        CGAL::to_double(linepoint.y()),
        CGAL::to_double(linepoint.z()),
        CGAL::to_double(linedir.dx()),
        CGAL::to_double(linedir.dy()),
        CGAL::to_double(linedir.dz()),
        true, true
        );
}


template <class K>
inline
typename Intersection_traits<K, Bbox_3, typename K::Line_3>::result_type
intersection(const Bbox_3 &box, 
	     const typename K::Line_3 &line, 
	     const K& k)
{
  return intersection(line, box, k);
}


template <class K>
typename Intersection_traits<K, typename K::Ray_3, Bbox_3>::result_type
intersection(const typename K::Ray_3 &ray,
	     const Bbox_3 &box, 
	     const K&)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Direction_3 Direction_3;
    const Point_3 &linepoint = ray.source();
    const Direction_3 &linedir = ray.direction();
    return intersection_bl<K>(box,
        CGAL::to_double(linepoint.x()),
        CGAL::to_double(linepoint.y()),
        CGAL::to_double(linepoint.z()),
        CGAL::to_double(linedir.dx()),
        CGAL::to_double(linedir.dy()),
        CGAL::to_double(linedir.dz()),
        false, true
        );
}


template <class K>
inline
typename Intersection_traits<K, Bbox_3, typename K::Ray_3>::result_type
intersection(const Bbox_3 &box, 
	     const typename K::Ray_3 &ray, 
	     const K& k)
{
  return intersection(ray, box, k);
}



template <class K>
typename Intersection_traits<K, typename K::Segment_3, Bbox_3>::result_type
intersection(const typename K::Segment_3 &seg, 
	     const Bbox_3 &box, 
	     const K&)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Vector_3 Vector_3;
    const Point_3 &linepoint = seg.source();
    const Vector_3 &diffvec = seg.target()-linepoint;
    return intersection_bl<K>(box,
        CGAL::to_double(linepoint.x()),
        CGAL::to_double(linepoint.y()),
        CGAL::to_double(linepoint.z()),
        CGAL::to_double(diffvec.x()),
        CGAL::to_double(diffvec.y()),
        CGAL::to_double(diffvec.z()),
        false, false
        );
}


template <class K>
inline
typename Intersection_traits<K, Bbox_3, typename K::Segment_3>::result_type
intersection(const Bbox_3 &box, 
	     const typename K::Segment_3 &seg, 
	     const K& k)
{
  return intersection(seg, box, k);
}


template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Line_3 &line,
	     const typename K::Iso_cuboid_3 &box, 
	     const K&)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Segment_3 Segment_3;
    typedef typename K::FT FT;
    bool all_values = true;
    FT _min = 0, _max = 0; // initialization to stop compiler warning
    Point_3 const & _ref_point=line.point();
    Vector_3 const & _dir=line.direction().vector();
    Point_3 const & _iso_min=(box.min)();
    Point_3 const & _iso_max=(box.max)();
    for (int i=0; i< _ref_point.dimension(); i++) {
        if (_dir.homogeneous(i) == 0) {
            if (_ref_point.cartesian(i) < _iso_min.cartesian(i)) {
                return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>();
            }
            if (_ref_point.cartesian(i) > _iso_max.cartesian(i)) {
                return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>();
            }
        } else {
            FT newmin, newmax;
            if (_dir.homogeneous(i) > 0) {
                newmin = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            } else {
                newmin = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            }
            if (all_values) {
                _min = newmin;
                _max = newmax;
            } else {
                if (newmin > _min)
                    _min = newmin;
                if (newmax < _max)
                    _max = newmax;
                if (_max < _min) {
                    return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>();
                }
            }
            all_values = false;
        }
    }
    CGAL_kernel_assertion(!all_values);
    if (_max == _min) {
        return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>(Point_3(_ref_point + _dir * _min ));
    }
    return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>(
        Segment_3(_ref_point + _dir*_min, _ref_point + _dir*_max));
}


template <class K>
inline
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Line_3>::result_type
intersection(const typename K::Iso_cuboid_3 &box, 
	     const typename K::Line_3 &line, 
	     const K& k)
{
  return intersection(line, box, k);
}



template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Ray_3 &ray,
	     const typename K::Iso_cuboid_3 &box, 
	     const K&)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Segment_3 Segment_3;
    typedef typename K::FT FT;
    bool all_values = true;
    FT _min = 0, _max = 0; // initialization to prevent compiler warning
    Point_3 const & _ref_point=ray.source();
    Vector_3 const & _dir=ray.direction().vector();
    Point_3 const & _iso_min=(box.min)();
    Point_3 const & _iso_max=(box.max)();

    for (int i=0; i< _ref_point.dimension(); i++) {
        if (_dir.homogeneous(i) == 0) {
            if (_ref_point.cartesian(i) < _iso_min.cartesian(i)) {
                return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>();
            }
            if (_ref_point.cartesian(i) > _iso_max.cartesian(i)) {
                return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>();
            }
        } else {
            FT newmin, newmax;
            if (_dir.homogeneous(i) > 0) {
                newmin = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            } else {
                newmin = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            }
            if (all_values) {
                _max = newmax;
            } else {
                if (newmax < _max)
                    _max = newmax;
            }
            if (newmin > _min)
                 _min = newmin;
            if (_max < _min)
                return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>();
            all_values = false;
        }
    }
    CGAL_kernel_assertion(!all_values);
    if (_max == _min) {
        return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>(Point_3(_ref_point + _dir * _min ));
    }
    return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>(
        Segment_3(_ref_point + _dir*_min, _ref_point + _dir*_max));
}


template <class K>
inline
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Ray_3>::result_type
intersection(const typename K::Iso_cuboid_3 &box, 
	     const typename K::Ray_3 &ray,
	     const K& k)
{
  return intersection(ray, box, k);
}


template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Segment_3 &seg,
	     const typename K::Iso_cuboid_3 &box, 
	     const K&)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Segment_3 Segment_3;
    typedef typename K::FT FT;
    FT _min = 0, _max;

    Point_3 const & _ref_point=seg.source();
    Vector_3 const & _dir=seg.direction().vector();
    Point_3 const & _iso_min=(box.min)();
    Point_3 const & _iso_max=(box.max)();
    int main_dir =
        (CGAL_NTS abs(_dir.x()) > CGAL_NTS abs(_dir.y()) )
            ? (CGAL_NTS abs(_dir.x()) > CGAL_NTS abs(_dir.z()) ? 0 : 2)
            : (CGAL_NTS abs(_dir.y()) > CGAL_NTS abs(_dir.z()) ? 1 : 2);
    _max = (seg.target().cartesian(main_dir)-_ref_point.cartesian(main_dir)) /
            _dir.cartesian(main_dir);

    for (int i=0; i< _ref_point.dimension(); i++) {
        if (_dir.homogeneous(i) == 0) {
            if (_ref_point.cartesian(i) < _iso_min.cartesian(i)) {
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>();
            }
            if (_ref_point.cartesian(i) > _iso_max.cartesian(i)) {
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>();
            }
        } else {
            FT newmin, newmax;
            if (_dir.homogeneous(i) > 0) {
                newmin = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            } else {
                newmin = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            }
            if (newmax < _max)
                _max = newmax;
            if (newmin > _min)
                 _min = newmin;
            if (_max < _min)
                return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>();
        }
    }
    if (_max == _min) {
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>(Point_3(_ref_point + _dir * _min ));
    }
    return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>(
        Segment_3(_ref_point + _dir*_min, _ref_point + _dir*_max));
}


template <class K>
inline
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Segment_3>::result_type
intersection(const typename K::Iso_cuboid_3 &box, 
	     const typename K::Segment_3 &seg,
	     const K& k)
{
  return intersection(seg, box, k);
}


template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Iso_cuboid_3>::result_type
intersection(
    const typename K::Iso_cuboid_3 &icub1,
    const typename K::Iso_cuboid_3 &icub2, 
    const K&)
{
    typedef typename K::Point_3 Point_3;
    typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

    Point_3 min_points[2];
    Point_3 max_points[2];
    min_points[0] = (icub1.min)();
    min_points[1] = (icub2.min)();
    max_points[0] = (icub1.max)();
    max_points[1] = (icub2.max)();
    const int DIM = 3;
    int min_idx[DIM];
    int max_idx[DIM];
    Point_3 newmin;
    Point_3 newmax;
    for (int dim = 0; dim < DIM; ++dim) {
        min_idx[dim] =
          min_points[0].cartesian(dim) >= min_points[1].cartesian(dim) ? 0 : 1;
        max_idx[dim] =
          max_points[0].cartesian(dim) <= max_points[1].cartesian(dim) ? 0 : 1;
        if (min_idx[dim] != max_idx[dim]
                && max_points[max_idx[dim]].cartesian(dim)
                   < min_points[min_idx[dim]].cartesian(dim))
            return intersection_return<typename K::Intersect_3, typename K::Iso_cuboid_3, typename K::Iso_cuboid_3>();
    }
    if (min_idx[0] == min_idx[1] && min_idx[0] == min_idx[2]) {
        newmin = min_points[min_idx[0]];
    } else {
        newmin = Point_3(
            min_idx[0] == 0
                ? wmult_hw((K*)0, min_points[0].hx(), min_points[1])
                : wmult_hw((K*)0, min_points[1].hx(), min_points[0])
            ,
            min_idx[1] == 0
                ? wmult_hw((K*)0, min_points[0].hy(), min_points[1])
                : wmult_hw((K*)0, min_points[1].hy(), min_points[0])
            ,
            min_idx[2] == 0
                ? wmult_hw((K*)0, min_points[0].hz(), min_points[1])
                : wmult_hw((K*)0, min_points[1].hz(), min_points[0])
            ,
            wmult_hw((K*)0, min_points[0].hw(), min_points[1]) );
    }
    if (max_idx[0] == max_idx[1] && max_idx[0] == max_idx[2]) {
        newmax = max_points[max_idx[0]];
    } else {
        newmax = Point_3(
            max_idx[0] == 0
                ? wmult_hw((K*)0, max_points[0].hx(), max_points[1])
                : wmult_hw((K*)0, max_points[1].hx(), max_points[0])
            ,
            max_idx[1] == 0
                ? wmult_hw((K*)0, max_points[0].hy(), max_points[1])
                : wmult_hw((K*)0, max_points[1].hy(), max_points[0])
            ,
            max_idx[2] == 0
                ? wmult_hw((K*)0, max_points[0].hz(), max_points[1])
                : wmult_hw((K*)0, max_points[1].hz(), max_points[0])
            ,
            wmult_hw((K*)0, max_points[0].hw(), max_points[1]) );
    }
    return intersection_return<typename K::Intersect_3, typename K::Iso_cuboid_3, typename K::Iso_cuboid_3>(Iso_cuboid_3(newmin, newmax));
}

template <class R>
inline bool
do_intersect(const Plane_3<R>& plane1, const Plane_3<R>& plane2, const R&)
{
  return bool(intersection(plane1, plane2));
}


template <class R>
inline bool
do_intersect(const Plane_3<R> &plane1, const Plane_3<R> &plane2,
             const Plane_3<R> &plane3, const R&)
{
  return bool(intersection(plane1, plane2, plane3));
}


template <class R>
inline bool
do_intersect(const Iso_cuboid_3<R> &i, const Iso_cuboid_3<R> &j, const R&)
{
  return bool(CGAL::intersection(i, j));
}

template <class R>
inline bool
do_intersect(const Line_3<R> &l, const Iso_cuboid_3<R> &j, const R&)
{
  return bool(CGAL::intersection(l, j));
}

template <class R>
inline bool
do_intersect(const Iso_cuboid_3<R> &j, const Line_3<R> &l, const R&)
{
  return bool(CGAL::intersection(l, j));
}

template <class R>
inline bool
do_intersect(const Ray_3<R> &r, const Iso_cuboid_3<R> &j, const R&)
{
  return bool(CGAL::intersection(r, j));
}

template <class R>
inline bool
do_intersect(const Iso_cuboid_3<R> &j, const Ray_3<R> &r, const R&)
{
  return bool(CGAL::intersection(r, j));
}

template <class R>
inline bool
do_intersect(const Segment_3<R> &s, const Iso_cuboid_3<R> &j, const R&)
{
  return bool(CGAL::intersection(s, j));
}

template <class R>
inline bool
do_intersect(const Iso_cuboid_3<R> &j, const Segment_3<R> &s, const R&)
{
  return bool(CGAL::intersection(s, j));
}

} // namespace internal

} //namespace CGAL
