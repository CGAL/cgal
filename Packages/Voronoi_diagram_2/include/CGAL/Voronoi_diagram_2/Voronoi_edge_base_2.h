// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>


#ifndef CGAL_VORONOI_DIAGRAM_2_VORONOI_EDGE_BASE_2_H
#define CGAL_VORONOI_DIAGRAM_2_VORONOI_EDGE_BASE_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

template<class DG, class P, class S, class Voronoi_edge,
	 class Voronoi_vertex, class Use_equal_2_obj_tag = Tag_false>
class Voronoi_edge_base_2
{
 private:
  enum Curve_type { SEGMENT = 0, RAY = 1, BISECTOR = 2 };
  
  typedef Voronoi_edge_base_2<DG,P,S,Voronoi_vertex,Voronoi_edge>  Self;

#if 0
  bool equal_2(const Voronoi_edge& o, const Tag_false&) const {
    if ( curve_t_ != o.curve_t_ ) { return false; }

    if ( is_segment() ) {
      return
	north_ == o.north_ && south_ == o.south_ &&
	west_ == o.west_ && east_ == o.east_;
    } else if ( is_ray() ) {
      if ( is_src_ != o.is_src_ ) { return false; }
      return
	north_ == o.north_ && south_ == o.south_ && west_ == o.west_;
    } else {
      CGAL_assertion( is_bisector() );
      return
	north_ == o.north_ && south_ == o.south_;
    }
  }

  bool equal_2(const Voronoi_edge& o, const Tag_true&) const {
    if ( curve_t_ != o.curve_t_ ) { return false; }

    Geom_traits gt;
    if ( is_segment() ) {
      return
	gt.equal_2_object()(north_, o.north_) &&
	gt.equal_2_object()(south_, o.south_) &&
	gt.equal_2_object()(west_, o.west_) &&
	gt.equal_2_object()(east_, o.east_);
    } else if ( is_ray() ) {
      if ( is_src_ != o.is_src_ ) { return false; }
      return
	gt.equal_2_object()(north_, o.north_) &&
	gt.equal_2_object()(south_, o.south_) &&
	gt.equal_2_object()(west_, o.west_);
    } else {
      CGAL_assertion( is_bisector() );
      return
	gt.equal_2_object()(north_, o.north_) &&
	gt.equal_2_object()(south_, o.south_);
    }
  }
#endif

 public:
  typedef DG                                    Delaunay_graph;
  typedef typename Delaunay_graph::Geom_traits  Geom_traits;
  typedef P                                     Point_2;
  typedef S                                     Site_2;
  typedef Voronoi_vertex                        Voronoi_vertex_2;
  typedef Voronoi_edge                          Voronoi_edge_2;

  bool is_bisector() const { return curve_t_ == BISECTOR; }
  bool is_segment() const { return curve_t_ == SEGMENT; }
  bool is_ray() const { return curve_t_ == RAY; }

  bool has_source() const {
    return is_segment() || (is_ray() && is_src_);
  }

  bool has_target() const {
    return is_segment() || (is_ray() && !is_src_);
  }

  Voronoi_vertex_2 source() const {
    CGAL_precondition( has_source() );
    Voronoi_vertex_2 vv;
    vv.set_sites(south_, north_, west_);
    return vv;
  }

  Voronoi_vertex_2 target() const {
    CGAL_precondition( has_target() );
    Voronoi_vertex_2 vv;
    vv.set_sites(north_, south_, east_);
    return vv;
  }

  const Site_2& up() const { return north_; }
  const Site_2& down() const { return south_; }
  const Site_2& left() const {
    CGAL_precondition( has_source() );
    return west_;
  }

  const Site_2& right() const {
    CGAL_precondition( has_target() );
    return east_;
  }

  Voronoi_edge_2 opposite() const {
    Voronoi_edge_2 ve;
    if ( is_segment() ) {
      ve.set_sites(south_, north_, east_, west_);
    } else if ( is_ray() ) {
      ve.set_sites(south_, north_, east_, !is_src_);
    } else {
      CGAL_assertion( is_bisector() );
      ve.set_sites(south_, north_);
    }
    return ve;
  }

#if 0
  bool operator==(const Voronoi_edge_2& o) const {
    return equal_2(o, Use_equal_2_obj_tag());
  }

  bool operator!=(const Voronoi_edge_2& o) const {
    return !(*this == o);
  }
#endif

 protected:
  void set_sites(const Site_2& s0, const Site_2& s1) {
    north_ = s0;
    south_ = s1;
    west_ = north_;
    east_ = south_;
    curve_t_ = BISECTOR;
  }

  void set_sites(const Site_2& s0, const Site_2& s1,
		 const Site_2& s3, bool is_src) {
    north_ = s0;
    south_ = s1;
    west_ = s3;
    east_ = s3;
    curve_t_ = RAY;
    is_src_ = is_src;
  }

  void set_sites(const Site_2& s0, const Site_2& s1,
		 const Site_2& s3, const Site_2& s4) {
    north_ = s0;
    south_ = s1;
    west_ = s3;
    east_ = s4;
    curve_t_ = SEGMENT;
  }

 protected:
  Site_2 north_, south_, west_, east_;
  Curve_type curve_t_;
  bool is_src_;
};


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_VORONOI_EDGE_BASE_2_H
