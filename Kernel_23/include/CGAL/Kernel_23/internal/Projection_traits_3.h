// Copyright (c) 1997-2010, 2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mariette Yvinec, Sebastien Loriot, Mael Rouxel-Labbé

#ifndef CGAL_INTERNAL_PROJECTION_TRAITS_3_H
#define CGAL_INTERNAL_PROJECTION_TRAITS_3_H

#include <CGAL/assertions.h>
#include <CGAL/tags.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>

#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Kernel_23/internal/Has_boolean_tags.h>

namespace CGAL {

namespace internal {

//project Point_3 along coordinate dim
template <class R,int dim>
struct Projector;

//project onto yz
template <class R>
struct Projector<R,0>
{
  typedef typename R::Less_y_3                Less_x_2;
  typedef typename R::Less_z_3                Less_y_2;
  typedef typename R::Compare_y_3             Compare_x_2;
  typedef typename R::Compare_z_3             Compare_y_2;
  typedef typename R::Equal_y_3               Equal_x_2;
  typedef typename R::Equal_z_3               Equal_y_2;

  static typename R::FT x(const typename R::Point_3& p) {return p.y();}
  static typename R::FT y(const typename R::Point_3& p) {return p.z();}
  static typename R::FT x(const typename R::Vector_3& p) {return p.y();}
  static typename R::FT y(const typename R::Vector_3& p) {return p.z();}
  static Bbox_2 bbox(const Bbox_3& bb) { return Bbox_2(bb.ymin(),bb.zmin(),bb.ymax(),bb.zmax()); }
  static const int x_index=1;
  static const int y_index=2;
};
//project onto xz
template <class R>
struct Projector<R,1>
{
  typedef typename R::Less_x_3                Less_x_2;
  typedef typename R::Less_z_3                Less_y_2;
  typedef typename R::Compare_x_3             Compare_x_2;
  typedef typename R::Compare_z_3             Compare_y_2;
  typedef typename R::Equal_x_3               Equal_x_2;
  typedef typename R::Equal_z_3               Equal_y_2;
  static typename R::FT x(const typename R::Point_3& p) {return p.x();}
  static typename R::FT y(const typename R::Point_3& p) {return p.z();}
  static typename R::FT x(const typename R::Vector_3& p) {return p.x();}
  static typename R::FT y(const typename R::Vector_3& p) {return p.z();}
  static Bbox_2 bbox(const Bbox_3& bb) { return Bbox_2(bb.xmin(),bb.zmin(),bb.xmax(),bb.zmax()); }
  static const int x_index=0;
  static const int y_index=2;
};

//project onto xy
template <class R>
struct Projector<R,2>
{
  typedef typename R::Less_x_3                Less_x_2;
  typedef typename R::Less_y_3                Less_y_2;
  typedef typename R::Compare_x_3             Compare_x_2;
  typedef typename R::Compare_y_3             Compare_y_2;
  typedef typename R::Equal_x_3               Equal_x_2;
  typedef typename R::Equal_y_3               Equal_y_2;
  static typename R::FT x(const typename R::Point_3& p) {return p.x();}
  static typename R::FT y(const typename R::Point_3& p) {return p.y();}
  static typename R::FT x(const typename R::Vector_3& p) {return p.x();}
  static typename R::FT y(const typename R::Vector_3& p) {return p.y();}
  static Bbox_2 bbox(const Bbox_3& bb) { return Bbox_2(bb.xmin(),bb.ymin(),bb.xmax(),bb.ymax()); }
  static const int x_index=0;
  static const int y_index=1;
};

template <class R,int dim>
class Construct_bbox_projected_2 {
public:
  typedef typename R::Point_3     Point;
  typedef Bbox_2 result_type;

  Bbox_2 operator()(const Point& p) const { typename R::Construct_bbox_3 bb;  return Projector<R, dim>::bbox(bb(p)); }
};

template <class R,int dim>
class Orientation_projected_3
{
public:
  typedef typename R::Point_3     Point;
  typename R::FT x(const Point &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point &p) const { return Projector<R,dim>::y(p); }

  typename R::Point_2 project(const Point& p) const
  {
    return typename R::Point_2(x(p),y(p));
  }

  CGAL::Orientation operator()(const Point& p,
                               const Point& q,
                               const Point& r) const
    {
      return CGAL::orientation(project(p), project(q), project(r));
    }
};

template <class R, int dim>
class Oriented_side_projected_3
{
public:
  typedef typename R::Segment_2   Segment_2;
  typedef typename R::Triangle_2  Triangle_2;

  typedef typename R::Point_3     Point;
  typedef typename R::Segment_3   Segment_3;
  typedef typename R::Triangle_3  Triangle_3;

  typename R::FT x(const Point& p) const { return Projector<R, dim>::x(p); }
  typename R::FT y(const Point& p) const { return Projector<R, dim>::y(p); }

  typename R::Point_2 project(const Point& p) const
  {
    return typename R::Point_2(x(p), y(p));
  }
  Triangle_2 project(const Triangle_3& t) const
  {
    typename R::Construct_vertex_3 v;
    return Triangle_2(project(v(t, 0)), project(v(t, 1)), project(v(t, 2)));
  }
  Segment_2 project(const Segment_3& s) const
  {
    typename R::Construct_source_3 source;
    typename R::Construct_target_3 target;
    return Segment_2(project(source(s)), project(target(s)));
  }
  CGAL::Oriented_side operator()(const Segment_3& s, const Triangle_3& t) const
  {
    return typename R::Oriented_side_2()(project(s), project(t));
  }
};

template <class R,int dim>
class Side_of_oriented_circle_projected_3
{
public:
  typedef typename R::Point_3     Point;
  typename R::FT x(const Point &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point &p) const { return Projector<R,dim>::y(p); }


  typename R::Point_2 project(const Point& p) const
  {
    return typename R::Point_2(x(p),y(p));
  }
  CGAL::Oriented_side operator() (const Point &p,
                                  const Point &q,
                                  const Point &r,
                                  const Point &s) const
    {
      return CGAL::side_of_oriented_circle(project(p),project(q),project(r),project(s) );
    }
};

template <class R,int dim>
class Side_of_bounded_circle_projected_3
{
public:
  typedef typename R::Point_3     Point;
  typename R::FT x(const Point &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point &p) const { return Projector<R,dim>::y(p); }


  typename R::Point_2 project(const Point& p) const
  {
    return typename R::Point_2(x(p),y(p));
  }
  CGAL::Bounded_side operator() (const Point &p,
                                  const Point &q,
                                  const Point &r,
                                  const Point &s) const
    {
      return CGAL::side_of_bounded_circle(project(p),project(q),project(r),project(s) );
    }

    CGAL::Bounded_side operator() (const Point &p,
                                  const Point &q,
                                  const Point &r) const
    {
      return CGAL::side_of_bounded_circle(project(p),project(q),project(r));
    }
};

template <class R,int dim>
class Compare_distance_projected_3
{
public:
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  Comparison_result operator()(const Point_3& p,const Point_3& q,const Point_3& r) const
  {
    Point_2 p2 = project(p);
    Point_2 q2 = project(q);
    Point_2 r2 = project(r);
    return compare_distance_to_point(p2,q2,r2);
  }
};

template <class R,int dim>
class Collinear_are_ordered_along_line_projected_3
{
public:
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;
  typedef typename R::FT        FT;
  FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  bool operator()(const Point_3& p,const Point_3& q,const Point_3& r) const
  {
    Point_2 p2 = project(p);
    Point_2 q2 = project(q);
    Point_2 r2 = project(r);
    return collinear_are_ordered_along_line(p2,q2,r2);
  }
};

template <class R, int dim>
class Compare_signed_distance_to_line_projected_3
{
public:
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }
  typedef typename R::Comparison_result result_type;

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  result_type operator()(const Point_3& p,
                         const Point_3& q,
                         const Point_3& r,
                         const Point_3& s) const
  {
    return typename R::Compare_signed_distance_to_line_2()
      (  project(p), project(q), project(r), project(s) );
  }
};

template <class R, int dim>
class Less_signed_distance_to_line_projected_3
{
public:
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }
  typedef typename R::Boolean result_type;

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  result_type operator()(const Point_3& p,
                         const Point_3& q,
                         const Point_3& r,
                         const Point_3& s) const
  {
    return typename R::Less_signed_distance_to_line_2()
      (  project(p), project(q), project(r), project(s) );
  }
};

template <class R,int dim>
class Squared_distance_projected_3
{
public:
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;
  typedef typename R::Line_3    Line_3;
  typedef typename R::Line_2    Line_2;
  typedef typename R::Segment_3 Segment_3;
  typedef typename R::Segment_2 Segment_2;
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  RT operator()(const Point_3& p, const Point_3& q) const
  {
          Point_2 p2(project(p));
          Point_2 q2(project(q));
          return squared_distance(p2, q2);
  }

  RT operator()(const Line_3& l, const Point_3& p) const
  {
    Point_2 p2(project(p));
    Line_2 l2(project(l.point(0)), project(l.point(1)));
    return squared_distance(p2, l2);
  }

  RT operator()(const Segment_3& s, const Point_3& p) const
  {
    Point_2 p2(project(p));
    Segment_2 s2(project(s.source()), project(s.target()));
    return squared_distance(p2, s2);
  }
};

template <class R,int dim>
class Construct_centroid_projected_3
{
public:
  typedef typename R::Point_3 Point_3;
  typedef typename R::Point_2 Point_2;
  typedef typename R::FT      RT;
  RT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  RT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p), y(p));
  }

  Point_3 embed(const Point_2& p) const
  {
    RT coords[3];
    coords[Projector<R,dim>::x_index] = p.x();
    coords[Projector<R,dim>::y_index] = p.y();
    coords[dim] = RT(0);

    return Point_3(coords[0], coords[1], coords[2]);
  }

  Point_3 operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
  {
    const Point_2 p2(project(p));
    const Point_2 q2(project(q));
    const Point_2 r2(project(r));
    return embed(CGAL::centroid(p2, q2, r2));
  }
};

template <class R,int dim>
class Compute_determinant_projected_3
{
public:
  typedef typename R::Vector_3 Vector_3;
  typedef typename R::Vector_2 Vector_2;
  typedef typename R::FT       RT;
  RT x(const Vector_3 &v) const { return Projector<R,dim>::x(v); }
  RT y(const Vector_3 &v) const { return Projector<R,dim>::y(v); }

  Vector_2 project(const Vector_3& v) const
  {
    return Vector_2(x(v), y(v));
  }

  RT operator()(const Vector_3& v, const Vector_3& w) const
  {
    const Vector_2 v2(project(v));
    const Vector_2 w2(project(w));
    return CGAL::determinant(v2, w2);
  }
};

template <class R,int dim>
class  Intersect_projected_3
{
public:
  typedef typename R::Point_3   Point_3;
  typedef typename R::Segment_3 Segment_3;
  typedef typename R::Point_2   Point_2;
  typedef typename R::Vector_2  Vector_2;
  typedef typename R::Segment_2 Segment_2;
  typedef typename R::FT        FT;

  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  FT alpha(const Point_2& p, const Point_2& source, const Point_2& target) const
  {
    FT dx = target.x() - source.x();
    FT dy = target.y() - source.y();
    return (CGAL::abs(dx)>CGAL::abs(dy)) ? ( p.x()-source.x() ) / dx : (p.y()-source.y() ) / dy;
  }



  boost::optional< boost::variant<Point_3,Segment_3> >
  operator()(const Segment_3& s1, const Segment_3& s2) const
  {
    typedef  boost::variant<Point_3, Segment_3> variant_type;

    Point_2 s1_source = project(s1.source());
    Point_2 s1_target = project(s1.target());
    Point_2 s2_source = project(s2.source());
    Point_2 s2_target = project(s2.target());
    Segment_2 s1_2(s1_source, s1_target);
    Segment_2 s2_2(s2_source, s2_target);
    CGAL_precondition(!s1_2.is_degenerate());
    CGAL_precondition(!s2_2.is_degenerate());

    //compute intersection points in projected plane
    //We know that none of the segment is degenerate

    auto o = intersection(s1_2,s2_2);
    if(! o){
      return boost::none;
    }

    if(const Segment_2* si = boost::get<Segment_2>(&*o)){
      FT src[3],tgt[3];
      //the third coordinate is the midpoint between the points on s1 and s2
      FT z1 = s1.source()[dim] + ( alpha(si->source(), s1_source, s1_target) * ( s1.target()[dim] - s1.source()[dim] ));
      FT z2 = s2.source()[dim] + ( alpha(si->source(), s2_source, s2_target) * ( s2.target()[dim] - s2.source()[dim] ));
      src[dim] = (z1+z2) / FT(2);


      z1 = s1.source()[dim] + ( alpha(si->target(), s1_source, s1_target) * ( s1.target()[dim] - s1.source()[dim] ));
      z2 = s2.source()[dim] + ( alpha(si->target(), s2_source, s2_target) * ( s2.target()[dim] - s2.source()[dim] ));

      tgt[dim] = (z1+z2) / FT(2);


      src[Projector<R,dim>::x_index] = si->source().x();
      src[Projector<R,dim>::y_index] = si->source().y();
      tgt[Projector<R,dim>::x_index] = si->target().x();
      tgt[Projector<R,dim>::y_index] = si->target().y();
      return boost::make_optional(variant_type(Segment_3( Point_3(src[0],src[1],src[2]),Point_3(tgt[0],tgt[1],tgt[2]) ) ) );
    }


    const Point_2* pi = boost::get<Point_2>(&*o);
    FT coords[3];
    //compute the third coordinate of the projected intersection point onto 3D segments
    FT z1 = s1.source()[dim] + ( alpha(*pi, s1_source, s1_target) * ( s1.target()[dim] - s1.source()[dim] ));
    FT z2 = s2.source()[dim] + ( alpha(*pi, s2_source, s2_target) * ( s2.target()[dim] - s2.source()[dim] ));

    coords[dim] = (z1+z2) / FT(2);
    coords[Projector<R,dim>::x_index] = pi->x();
    coords[Projector<R,dim>::y_index] = pi->y();

    Point_3 res(coords[0],coords[1],coords[2]);
    CGAL_assertion(x(res)==pi->x() && y(res)==pi->y());
    return boost::make_optional(variant_type(res));
  }
};

template <class R, int dim>
class Circumcenter_center_projected
{
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;

  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  Point_3 embed (const Point_2& p) const
  {
    typename R::FT coords[3];
    coords[Projector<R,dim>::x_index]=p.x();
    coords[Projector<R,dim>::y_index]=p.y();
    coords[dim]=typename R::FT(0);
    return Point_3(coords[0],coords[1],coords[2]);
  }

public:
  Point_3 operator() (const Point_3& p1,const Point_3& p2) const
  {
    return embed( CGAL::circumcenter(project(p1),project(p2)) );
  }

  Point_3 operator() (const Point_3& p1,const Point_3& p2,const Point_3& p3) const
  {
    return embed( CGAL::circumcenter(project(p1),project(p2),project(p3)) );
  }
};

template <class R, int dim>
class Compute_area_projected
{
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;

  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }


public:
  typename R::FT operator() (const Point_3& p1,const Point_3& p2,const Point_3& p3) const
  {
    return R().compute_area_2_object() ( project(p1),project(p2),project(p3) );
  }
};

template <class R, int dim>
class Compute_squared_radius_projected
{
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;

  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }


public:
  typename R::FT operator() (const Point_3& p1,const Point_3& p2,const Point_3& p3) const
  {
    return R().compute_squared_radius_2_object() ( project(p1),project(p2),project(p3) );
  }
  typename R::FT operator() (const Point_3& p1,const Point_3& p2) const
  {
    return R().compute_squared_radius_2_object() ( project(p1),project(p2) );
  }

  typename R::FT operator() (const Point_3& p1) const
  {
    return R().compute_squared_radius_2_object() ( project(p1) );
  }
};

template <class R,int dim>
class Compute_scalar_product_projected_3
{
public:
  typedef typename R::Vector_3    Vector_3;
  typedef typename R::FT          FT;
  FT x(const Vector_3 &v) const { return Projector<R,dim>::x(v); }
  FT y(const Vector_3 &v) const { return Projector<R,dim>::y(v); }

  FT operator()(const Vector_3& v1, const Vector_3& v2) const
  {
    return x(v1)*x(v2) + y(v1)*y(v2);
  }
};

template <class R,int dim>
class Compute_squared_length_projected_3
{
  typedef typename R::Vector_3    Vector_3;
  typedef typename R::FT          FT;

  typedef FT result_type;

  FT x(const Vector_3 &v) const { return Projector<R,dim>::x(v); }
  FT y(const Vector_3 &v) const { return Projector<R,dim>::y(v); }

public:
  FT operator()(const Vector_3& v) const
  {
    return CGAL::square(x(v)) + CGAL::square(y(v));
  }
};

template <class R, int dim>
struct Angle_projected_3{
  typedef typename R::Point_3   Point_3;
  typedef typename R::Point_2   Point_2;

  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  CGAL::Angle operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
  {
    return CGAL::angle(project(p), project(q), project(r));
  }
};

// weighted functors
template <class R, int dim>
class Compare_power_distance_projected_3
{
public:
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Weighted_point_2          Weighted_point_2;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Weighted_point_3          Weighted_point_3;
  typedef typename R::FT                        FT;

  FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p), y(p));
  }

  Weighted_point_2 project(const Weighted_point_3& wp) const
  {
    const Point_3& p = R().construct_point_3_object()(wp);
    return Weighted_point_2(Point_2(x(p), y(p)), wp.weight());
  }

  Comparison_result operator()(const Point_3& p,
                               const Weighted_point_3& wq, const Weighted_point_3& wr) const
  {
    Point_2 p2 = project(p);
    Weighted_point_2 wq2 = project(wq);
    Weighted_point_2 wr2 = project(wr);
    return CGAL::compare_power_distance(p2, wq2, wr2);
  }
};

template <class R, int dim>
class Compute_power_product_projected_3
{
public:
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Weighted_point_2          Weighted_point_2;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Weighted_point_3          Weighted_point_3;
  typedef typename R::FT                        FT;

  FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Weighted_point_2 project(const Weighted_point_3& wp) const
  {
    const Point_3& p = R().construct_point_3_object()(wp);
    return Weighted_point_2(Point_2(x(p), y(p)), wp.weight());
  }

  FT operator()(const Weighted_point_3& wp, const Weighted_point_3& wq) const
  {
    Weighted_point_2 wp2(project(wp));
    Weighted_point_2 wq2(project(wq));
    return CGAL::power_product(wp2, wq2);
  }
};

template <class R, int dim>
class Compute_squared_radius_smallest_orthogonal_circle_projected_3
{
public:
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Weighted_point_2          Weighted_point_2;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Weighted_point_3          Weighted_point_3;
  typedef typename R::FT                        FT;

  FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Weighted_point_2 project(const Weighted_point_3& wp) const
  {
    const Point_3& p = R().construct_point_3_object()(wp);
    return Weighted_point_2(Point_2(x(p), y(p)), wp.weight());
  }

  FT operator()(const Weighted_point_3& wp1, const Weighted_point_3& wp2, const Weighted_point_3& wp3) const
  {
    return CGAL::squared_radius_smallest_orthogonal_circle( project(wp1), project(wp2), project(wp3) );
  }

  FT operator()(const Weighted_point_3& wp1, const Weighted_point_3& wp2) const
  {
    return CGAL::squared_radius_smallest_orthogonal_circle( project(wp1), project(wp2) );
  }

  FT operator()(const Weighted_point_3& wp1) const
  {
    return CGAL::squared_radius_smallest_orthogonal_circle( project(wp1) );
  }
};

template <class R, int dim>
class Construct_radical_axis_projected_3
{
public:
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Line_3                    Line_3;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Weighted_point_2          Weighted_point_2;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Weighted_point_3          Weighted_point_3;
  typedef typename R::FT                        FT;

  FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Weighted_point_2 project(const Weighted_point_3& wp) const
  {
    const Point_3& p = R().construct_point_3_object()(wp);
    return Weighted_point_2(Point_2(x(p), y(p)), wp.weight());
  }

  Line_3 embed (const Line_2& l) const
  {
    Point_2 p0 = l.point(0);
    Point_2 p1 = l.point(1);

    FT coords_p0[3];
    coords_p0[Projector<R,dim>::x_index] = p0.x();
    coords_p0[Projector<R,dim>::y_index] = p0.y();
    coords_p0[dim] = FT(0);

    FT coords_p1[3];
    coords_p1[Projector<R,dim>::x_index] = p1.x();
    coords_p1[Projector<R,dim>::y_index] = p1.y();
    coords_p1[dim] = FT(0);

    return Line_3(Point_3(coords_p0[0], coords_p0[1], coords_p0[2]),
                  Point_3(coords_p1[0], coords_p1[1], coords_p1[2]));
  }

public:
  Line_3 operator() (const Weighted_point_3& wp1, const Weighted_point_3& wp2) const
  {
    return embed( CGAL::radical_axis(project(wp1), project(wp2)) );
  }
};

template <class R, int dim>
class Construct_weighted_circumcenter_projected_3
{
public:
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Weighted_point_2          Weighted_point_2;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Weighted_point_3          Weighted_point_3;
  typedef typename R::FT                        FT;

  FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Weighted_point_2 project(const Weighted_point_3& wp) const
  {
    const Point_3& p = R().construct_point_3_object()(wp);
    return Weighted_point_2(Point_2(x(p), y(p)), wp.weight());
  }

  Point_3 embed (const Point_2& p) const
  {
    FT coords[3];
    coords[Projector<R,dim>::x_index] = p.x();
    coords[Projector<R,dim>::y_index] = p.y();
    coords[dim] = FT(0);

    return Point_3(coords[0], coords[1], coords[2]);
  }

public:
  Point_3 operator()(const Weighted_point_3& wp1,
                     const Weighted_point_3& wp2,
                     const Weighted_point_3& wp3) const
  {
    return embed( CGAL::weighted_circumcenter(project(wp1), project(wp2), project(wp3)) );
  }
};

template <class R, int dim>
class Power_side_of_bounded_power_circle_projected_3
{
public:
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Weighted_point_2          Weighted_point_2;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Weighted_point_3          Weighted_point_3;
  typedef typename R::FT                        FT;

  FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Weighted_point_2 project(const Weighted_point_3& wp) const
  {
    const Point_3& p = R().construct_point_3_object()(wp);
    return Weighted_point_2(Point_2(x(p), y(p)), wp.weight());
  }

  CGAL::Bounded_side operator()(const Weighted_point_3 &wp,
                                const Weighted_point_3 &wq,
                                const Weighted_point_3 &wr,
                                const Weighted_point_3 &ws) const
  {
    return CGAL::power_side_of_bounded_power_circle(project(wp), project(wq),
                                                    project(wr), project(ws));
  }

  CGAL::Bounded_side operator()(const Weighted_point_3 &wp,
                                const Weighted_point_3 &wq,
                                const Weighted_point_3 &wr) const
  {
    return CGAL::power_side_of_bounded_power_circle(project(wp), project(wq), project(wr));
  }

  CGAL::Bounded_side operator()(const Weighted_point_3 &wp,
                                const Weighted_point_3 &wq) const
  {
    return CGAL::power_side_of_bounded_power_circle(project(wp), project(wq));
  }
};

template <class R, int dim>
class Power_side_of_oriented_power_circle_projected_3
{
public:
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Weighted_point_2          Weighted_point_2;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Weighted_point_3          Weighted_point_3;
  typedef typename R::FT                        FT;

  FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Weighted_point_2 project(const Weighted_point_3& wp) const
  {
    const Point_3& p = R().construct_point_3_object()(wp);
    return Weighted_point_2(Point_2(x(p), y(p)), wp.weight());
  }

  CGAL::Oriented_side operator()(const Weighted_point_3 &wp,
                                 const Weighted_point_3 &wq,
                                 const Weighted_point_3 &wr,
                                 const Weighted_point_3 &ws) const
  {
    return CGAL::power_side_of_oriented_power_circle(project(wp), project(wq),
                                                     project(wr), project(ws));
  }

  CGAL::Oriented_side operator()(const Weighted_point_3 &wp,
                                 const Weighted_point_3 &wq,
                                 const Weighted_point_3 &wr) const
  {
    return CGAL::power_side_of_oriented_power_circle(project(wp), project(wq), project(wr));
  }

  CGAL::Oriented_side operator()(const Weighted_point_3 &wp,
                                 const Weighted_point_3 &wq) const
  {
    return CGAL::power_side_of_oriented_power_circle(project(wp), project(wq));
  }
};

template <class R, int dim>
class Projection_traits_3;

template <class R, int dim, bool has_filtered_predicates>
struct Projection_traits_base_3 {};

template <class R, int dim>
struct Projection_traits_base_3< R, dim, true> {
  typedef Projection_traits_3<typename R::Exact_kernel, dim>  Exact_kernel;
  Exact_kernel exact_kernel() const { return {}; }
};

// This is for projection traits along a specific canonical plane (xy, yz, xz)
// The generic class for an arbitrary normal is CGAL::Projection_traits_3<K> (not in `internal`)
template <class R, int dim>
class Projection_traits_3
    : public Projection_traits_base_3<
          R, dim, internal::Has_filtered_predicates<R>::value> {
public:
  enum { Has_filtered_predicates = internal::Has_filtered_predicates<R>::value };
  typedef Boolean_tag<Has_filtered_predicates> Has_filtered_predicates_tag;

  typedef Projection_traits_3<R,dim>   Traits;
  typedef R                                                   Rp;
  typedef typename R::FT                                      FT;
  typedef typename Rp::Point_3                                Point_2;
  typedef typename Rp::Weighted_point_3                       Weighted_point_2;
  typedef typename Rp::Segment_3                              Segment_2;
  typedef typename Rp::Vector_3                               Vector_2;
  typedef typename Rp::Triangle_3                             Triangle_2;
  typedef typename Rp::Line_3                                 Line_2;
  typedef typename Rp::Ray_3                                  Ray_2;

  typedef typename Projector<R,dim>::Less_x_2                 Less_x_2;
  typedef typename Projector<R,dim>::Less_y_2                 Less_y_2;
  typedef typename Projector<R,dim>::Compare_x_2              Compare_x_2;
  typedef typename Projector<R,dim>::Compare_y_2              Compare_y_2;
  typedef Orientation_projected_3<Rp,dim>                     Orientation_2;
  typedef Oriented_side_projected_3<Rp,dim>                   Oriented_side_2;
  typedef Angle_projected_3<Rp,dim>                           Angle_2;
  typedef Side_of_oriented_circle_projected_3<Rp,dim>         Side_of_oriented_circle_2;
  typedef Compare_signed_distance_to_line_projected_3<Rp,dim> Compare_signed_distance_to_line_2;
  typedef Less_signed_distance_to_line_projected_3<Rp,dim>    Less_signed_distance_to_line_2;
  typedef Side_of_bounded_circle_projected_3<Rp,dim>          Side_of_bounded_circle_2;
  typedef Compare_distance_projected_3<Rp,dim>                Compare_distance_2;
  typedef Collinear_are_ordered_along_line_projected_3<Rp,dim> Collinear_are_ordered_along_line_2;
  typedef Squared_distance_projected_3<Rp,dim>                Compute_squared_distance_2;
  typedef Intersect_projected_3<Rp,dim>                       Intersect_2;
  typedef Compute_squared_radius_projected<Rp,dim>            Compute_squared_radius_2;
  typedef Compute_scalar_product_projected_3<Rp,dim>          Compute_scalar_product_2;
  typedef Compute_squared_length_projected_3<Rp,dim>          Compute_squared_length_2;

  typedef Compare_power_distance_projected_3<Rp,dim>          Compare_power_distance_2;
  typedef Compute_power_product_projected_3<Rp,dim>           Compute_power_product_2;
  typedef Compute_squared_radius_smallest_orthogonal_circle_projected_3<Rp,dim>
                                                              Compute_squared_radius_smallest_orthogonal_circle_2;
  typedef Construct_radical_axis_projected_3<Rp,dim>          Construct_radical_axis_2;
  typedef Construct_weighted_circumcenter_projected_3<Rp,dim> Construct_weighted_circumcenter_2;
  typedef Power_side_of_bounded_power_circle_projected_3<Rp,dim> Power_side_of_bounded_power_circle_2;
  typedef Power_side_of_oriented_power_circle_projected_3<Rp, dim> Power_side_of_oriented_power_circle_2;
  typedef Construct_bbox_projected_2<Rp,dim>                  Construct_bbox_2;

  typedef Construct_centroid_projected_3<Rp,dim>              Construct_centroid_2;
  typedef Compute_determinant_projected_3<Rp,dim>             Compute_determinant_2;

  typedef typename Rp::Construct_point_3                      Construct_point_2;
  typedef typename Rp::Construct_weighted_point_3             Construct_weighted_point_2;
  typedef typename Rp::Construct_segment_3                    Construct_segment_2;
  typedef typename Rp::Construct_translated_point_3           Construct_translated_point_2;
  typedef typename Rp::Construct_midpoint_3                   Construct_midpoint_2;
  typedef typename Rp::Construct_barycenter_3                 Construct_barycenter_2;
  typedef typename Rp::Construct_vector_3                     Construct_vector_2;
  typedef typename Rp::Construct_scaled_vector_3              Construct_scaled_vector_2;
  typedef typename Rp::Construct_triangle_3                   Construct_triangle_2;
  typedef typename Rp::Construct_line_3                       Construct_line_2;


  struct Less_xy_2 {
    typedef typename R::Boolean result_type;
    bool operator()(const Point_2& p, const Point_2& q) const
    {
      Compare_x_2 cx;
      Comparison_result crx = cx(p,q);
      if(crx == SMALLER){ return true;}
      if(crx == LARGER){return false;}
      Less_y_2 ly;
      return ly(p,q);
    }
  };


  struct Less_yx_2 {
    typedef typename R::Boolean result_type;
    bool operator()(const Point_2& p, const Point_2& q) const
    {
      Compare_y_2 cy;
      Comparison_result cry = cy(p,q);
      if(cry == SMALLER){ return true;}
      if(cry == LARGER){return false;}
      Less_x_2 lx;
      return lx(p,q);
    }
  };

  struct Equal_2 {
    typedef typename R::Boolean result_type;
    bool operator()(const Point_2& p, const Point_2& q) const
    {

      Equal_x_2 eqx;
      Equal_y_2 eqy;
      return eqx(p,q) && eqy(p,q);
    }
  };

  struct Left_turn_2 {
    typedef typename R::Boolean result_type;
    bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {

      Orientation_2 ori;
      return ori(p,q,r) == LEFT_TURN;
    }
  };

  struct Collinear_2 {
    typedef typename R::Boolean result_type;
    bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      Orientation_2 ori;
      return ori(p,q,r) == COLLINEAR;
    }
  };

  //for natural_neighbor_coordinates_2
  typedef typename Projector<R,dim>::Equal_x_2                Equal_x_2;
  typedef typename Projector<R,dim>::Equal_y_2                Equal_y_2;
  typedef Circumcenter_center_projected<Rp,dim>               Construct_circumcenter_2;
  typedef Compute_area_projected<Rp,dim>                      Compute_area_2;
  Construct_circumcenter_2 construct_circumcenter_2_object () const {return Construct_circumcenter_2();}
  Compute_area_2 compute_area_2_object () const {return Compute_area_2();}


  // for compatibility with previous versions
  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;

  Projection_traits_3(){}
  Projection_traits_3(
                   const Projection_traits_3&){}
  Projection_traits_3 &operator=(
            const Projection_traits_3&){return *this;}

  typename Rp::FT x(const Point_2 &p) const { return Projector<R,dim>::x(p); }
  typename Rp::FT y(const Point_2 &p) const { return Projector<R,dim>::y(p); }


 Equal_2
  equal_2_object() const
    { return Equal_2();}

  Left_turn_2
  left_turn_2_object() const
    { return Left_turn_2();}

  Less_x_2
  less_x_2_object() const
    { return Less_x_2();}

  Less_xy_2
  less_xy_2_object() const
    { return Less_xy_2();}

  Less_yx_2
  less_yx_2_object() const
    { return Less_yx_2();}

  Compare_signed_distance_to_line_2
    compare_signed_distance_to_line_2_object() const
    {return Compare_signed_distance_to_line_2();}

  Less_signed_distance_to_line_2
    less_signed_distance_to_line_2_object() const
    {return Less_signed_distance_to_line_2();}

  Less_y_2
  less_y_2_object() const
    { return Less_y_2();}
  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2();}
  Angle_2
  angle_2_object() const {
          return Angle_2();
  }

  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2();}

  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

  Oriented_side_2
  oriented_side_2_object() const
    { return Oriented_side_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
    {return Side_of_oriented_circle_2();}

  Side_of_bounded_circle_2
  side_of_bounded_circle_2_object() const
    {return Side_of_bounded_circle_2();}

  Compare_distance_2
  compare_distance_2_object() const
  {
    return Compare_distance_2();
  }

  Compute_squared_distance_2
  compute_squared_distance_2_object () const
  {
    return Compute_squared_distance_2();
  }

  Compute_squared_radius_2
  compute_squared_radius_2_object () const
  {
    return Compute_squared_radius_2();
  }

  Intersect_2
  intersect_2_object () const
  {
    return Intersect_2();
  }

  Construct_point_2 construct_point_2_object() const
  { return Construct_point_2();}

  Construct_weighted_point_2 construct_weighted_point_2_object() const
  { return Construct_weighted_point_2();}

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_translated_point_2  construct_translated_point_2_object() const
    {return Construct_translated_point_2();}

  Construct_midpoint_2  construct_midpoint_2_object() const
    {return Construct_midpoint_2();}

  Construct_barycenter_2  construct_barycenter_2_object() const
    {return Construct_barycenter_2();}

  Construct_vector_2  construct_vector_2_object() const
    {return Construct_vector_2();}

  Construct_scaled_vector_2  construct_scaled_vector_2_object() const
    {return Construct_scaled_vector_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

  Construct_line_2  construct_line_2_object() const
    {return Construct_line_2();}

  Construct_bbox_2  construct_bbox_2_object() const
    {return Construct_bbox_2();}

  Construct_centroid_2  construct_centroid_2_object() const
    {return Construct_centroid_2();}

  Compute_determinant_2  compute_determinant_2_object() const
    {return Compute_determinant_2();}

  Compute_scalar_product_2 compute_scalar_product_2_object() const
    {return Compute_scalar_product_2();}

  Collinear_2 collinear_2_object() const
    {return Collinear_2();}

  Collinear_are_ordered_along_line_2 collinear_are_ordered_along_line_2_object() const
    {return Collinear_are_ordered_along_line_2();}

  Compute_squared_length_2 compute_squared_length_2_object() const
    {return Compute_squared_length_2();}

  Compare_power_distance_2 compare_power_distance_2_object() const
    {return Compare_power_distance_2();}

  Compute_power_product_2 compute_power_product_2_object() const
    {return Compute_power_product_2();}

  Compute_squared_radius_smallest_orthogonal_circle_2
  compute_squared_radius_smallest_orthogonal_circle_2_object() const
    {return Compute_squared_radius_smallest_orthogonal_circle_2();}

  Construct_radical_axis_2 construct_radical_axis_2_object() const
    {return Construct_radical_axis_2();}

  Construct_weighted_circumcenter_2 construct_weighted_circumcenter_2_object() const
    {return Construct_weighted_circumcenter_2();}

  Power_side_of_bounded_power_circle_2 power_side_of_bounded_power_circle_2_object() const
    {return Power_side_of_bounded_power_circle_2();}

  Power_side_of_oriented_power_circle_2 power_side_of_oriented_power_circle_2_object() const
    {return Power_side_of_oriented_power_circle_2();}
};


} } //namespace CGAL::internal

#endif // CGAL_INTERNAL_PROJECTION_TRAITS_3_H
