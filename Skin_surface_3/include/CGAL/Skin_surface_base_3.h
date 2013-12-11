// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SKIN_SURFACE_BASE_3_H
#define CGAL_SKIN_SURFACE_BASE_3_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_3.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/shared_ptr.hpp>

// For the Weighted_converter
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

// Used for the triangulated mixed complex / Voronoi diagram
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <CGAL/Skin_surface_quadratic_surface_3.h>

// Skin surface mesher
#include <CGAL/Marching_tetrahedra_traits_skin_surface_3.h>
#include <CGAL/Skin_surface_marching_tetrahedra_observer_3.h>

// Skin surface subdivider
#include <CGAL/Skin_surface_refinement_policy_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>

namespace CGAL { 

template <class MixedComplexTraits_3> 
class Skin_surface_base_3 {
  typedef MixedComplexTraits_3            Gt;
  typedef Skin_surface_base_3<Gt>         Self;
  
public:
  typedef MixedComplexTraits_3            Geometric_traits;
  typedef typename Gt::Weighted_point     Weighted_point;
  typedef typename Gt::Bare_point         Bare_point;
  typedef typename Gt::FT                 FT;
  // For normal
  typedef typename Gt::Vector_3           Vector;
  
  typedef Regular_triangulation_3<Gt>     Regular;

private:
  typedef Exact_predicates_inexact_constructions_kernel     Filtered_kernel;
public:
  typedef Skin_surface_quadratic_surface_3<Filtered_kernel> Quadratic_surface;

  typedef typename Regular::Vertex_handle                Vertex_handle;
  typedef typename Regular::Edge                         Edge;
  typedef typename Regular::Facet                        Facet;
  typedef typename Regular::Facet_circulator             Facet_circulator;
  typedef typename Regular::Cell_handle                  Cell_handle;
  typedef Triangulation_simplex_3<Regular>               Simplex;

  // pair of a del- and vor-simplex
  typedef std::pair<Simplex,Simplex>                     Anchor_point;

private:
  typedef typename Regular::Finite_vertices_iterator     Finite_vertices_iterator;
  typedef typename Regular::Finite_edges_iterator        Finite_edges_iterator;
  typedef typename Regular::Finite_facets_iterator       Finite_facets_iterator;
  typedef typename Regular::Finite_cells_iterator        Finite_cells_iterator;

public:
  typedef Anchor_point                                  Vertex_info;
  typedef std::pair<Simplex, boost::shared_ptr<Quadratic_surface> >       Cell_info;

private:
  // Triangulated_mixed_complex:
  typedef CGAL::Exact_predicates_exact_constructions_kernel            FK;
  typedef Triangulation_vertex_base_with_info_3<Vertex_info, FK>       Vb;
  typedef Triangulation_cell_base_with_info_3<Cell_info, FK>           Cb;
  typedef Triangulation_data_structure_3<Vb,Cb>                        Tds;
public:
  typedef Triangulation_3<FK, Tds>               TMC;
  typedef typename TMC::Geom_traits              TMC_Geom_traits;
  typedef typename TMC::Finite_vertices_iterator TMC_Vertex_iterator;
  typedef typename TMC::Finite_cells_iterator    TMC_Cell_iterator;
  typedef typename TMC::Vertex_handle            TMC_Vertex_handle;
  typedef typename TMC::Cell_handle              TMC_Cell_handle;
  typedef typename TMC::Point                    TMC_Point;


  // Constructor
  template < class WP_iterator >
  Skin_surface_base_3(WP_iterator begin, WP_iterator end, FT shrink,
                      bool grow_balls = true,
                      Gt gt_ = Gt(), bool _verbose = false);
  
  template <class Polyhedron_3>
  void mesh_surface_3(Polyhedron_3 &p) const;

  template <class Polyhedron_3>
  void subdivide_mesh_3(Polyhedron_3 &p) const;
  


  // Access functions:    
  TMC &triangulated_mixed_complex();
  const FT shrink_factor() const { return shrink; }
  Geometric_traits &geometric_traits() const { return gt; }
  //  TMC &triangulated_mixed_complex() { return _tmc; }
  Regular &regular() { return _regular; }

  // Predicates and constructions
  Sign sign(TMC_Vertex_handle vit) const;
  
  Sign sign(const Bare_point &p, 
            const TMC_Cell_handle start = TMC_Cell_handle()) const;
  
  Sign sign(const Bare_point &p, 
            const Cell_info &info) const;

  // Uses inexact computations to compute the sign
  Sign sign_inexact(const Bare_point &p, 
                    const Cell_info &info) const;
  
  void intersect(TMC_Cell_handle ch, int i, int j, Bare_point &p) const;
  void intersect(Bare_point &p1, Bare_point &p2, 
                 TMC_Cell_handle &s1, TMC_Cell_handle &s2,
                 Bare_point &p) const;
  
  void intersect_with_transversal_segment
    (Bare_point &p,
     const TMC_Cell_handle &start = TMC_Cell_handle()) const;

  
  template <class Gt2>
  static typename Gt2::Bare_point
  get_weighted_circumcenter(const Simplex &s, Gt2 &traits);

  Vector
  normal(const Bare_point &p, 
         TMC_Cell_handle start = TMC_Cell_handle()) const;

  template <class Gt2>
  static typename Gt2::Bare_point
  get_anchor_point(const Anchor_point &anchor, Gt2 &traits);

private:
  void construct_bounding_box(); 

  template< class Traits >
  Skin_surface_quadratic_surface_3<Traits> 
  construct_surface(
    const Simplex &sim, 
    const Traits &traits = typename Geometric_traits::Kernel()) const;

  Sign compare(Cell_info &info, const Bare_point &p1, const Bare_point &p2) const;
  Sign compare(Cell_info &info1, const Bare_point &p1, 
               Cell_info &info2, const Bare_point &p2) const;
  
  TMC_Cell_handle
  locate_in_tmc(const Bare_point &p, 
                TMC_Cell_handle start = TMC_Cell_handle()) const;
private:
  FT shrink;
  Geometric_traits gt;
  Regular _regular;
  // Triangulated mixed complex or Voronoi diagram:
  TMC _tmc;
  
  bool verbose;
};

template <class MixedComplexTraits_3> 
template < class WP_iterator >
Skin_surface_base_3<MixedComplexTraits_3>::
Skin_surface_base_3(WP_iterator begin, WP_iterator end, FT shrink_,
                   bool grow_balls,
                   Gt gt_, bool verbose_) 
: shrink(shrink_), gt(gt_), verbose(verbose_)
{
  gt.set_shrink(shrink);
  CGAL_assertion(begin != end);

  if (grow_balls) {
    for (; begin != end; begin++) {
      regular().insert(Weighted_point(*begin, begin->weight()/shrink_factor()));
    }
  } else {
    regular().insert(begin, end);
  }

  construct_bounding_box();

  if (verbose) {
    std::cerr << "Triangulation ready" << std::endl;
    std::cerr << "Vertices: " << regular().number_of_vertices() << std::endl;
    std::cerr << "Cells:    " << regular().number_of_cells() << std::endl;
  }
}

template <class MixedComplexTraits_3> 
template <class Polyhedron_3>
void
Skin_surface_base_3<MixedComplexTraits_3>::
mesh_surface_3(Polyhedron_3 &p) const {
  typedef Polyhedron_3 Polyhedron;

  typedef Marching_tetrahedra_traits_skin_surface_3<
    Self,
    TMC_Vertex_iterator,
    TMC_Cell_iterator,
    typename Polyhedron::HalfedgeDS>               Traits;
  typedef Skin_surface_marching_tetrahedra_observer_3<
    TMC_Vertex_iterator,
    TMC_Cell_iterator,
    Polyhedron>                                    Observer;

  // Extract the coarse mesh using marching_tetrahedra
  Traits   marching_traits(*this);
  Observer marching_observer;
  marching_tetrahedra_3(_tmc.finite_vertices_begin(), 
                        _tmc.finite_vertices_end(), 
                        _tmc.finite_cells_begin(), 
                        _tmc.finite_cells_end(), 
                        p, 
                        marching_traits,
                        marching_observer);
}

template <class MixedComplexTraits_3> 
template <class Polyhedron_3>
void
Skin_surface_base_3<MixedComplexTraits_3>::
subdivide_mesh_3(Polyhedron_3 &p) const {
  typedef Skin_surface_refinement_policy_3<Self, Polyhedron_3> Policy;
  typedef Skin_surface_sqrt3<Self, Polyhedron_3, Policy>       Subdivider;

  Policy policy(*this);
  Subdivider subdivider(*this, p, policy);
  subdivider.subdivide();
}


template <class MixedComplexTraits_3>
typename Skin_surface_base_3<MixedComplexTraits_3>::TMC &
Skin_surface_base_3<MixedComplexTraits_3>::
triangulated_mixed_complex() {
  return _tmc;
}

template <class MixedComplexTraits_3> 
typename Skin_surface_base_3<MixedComplexTraits_3>::Vector
Skin_surface_base_3<MixedComplexTraits_3>::
normal(const Bare_point &p, 
       TMC_Cell_handle start) const {
  if (start == TMC_Cell_handle()) {
    start = locate_in_tmc(p,start);
  }

  return start->info().second->gradient(p);
}

template <class MixedComplexTraits_3> 
Sign 
Skin_surface_base_3<MixedComplexTraits_3>::
sign(TMC_Vertex_handle vit) const {
  CGAL_assertion(!_tmc.is_infinite(vit));
  TMC_Cell_handle ch = vit->cell();
  if (_tmc.is_infinite(ch)) {
    std::vector<TMC_Cell_handle> nbs;
    _tmc.incident_cells(vit, std::back_inserter(nbs));
    for (typename std::vector<TMC_Cell_handle>::iterator it = nbs.begin();
         _tmc.is_infinite(ch) && (it != nbs.end());
         it++) {
      ch = *it;
    }
  }
  CGAL_assertion(!_tmc.is_infinite(ch));

  // don't use sign, since the point is constructed:
  CGAL_BRANCH_PROFILER(std::string(" NGHK: failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  {
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    Protect_FPU_rounding<true> P;
    try
      {
	Sign result = vit->cell()->info().second->sign(vit->point());
	if (is_certain(result))
	  return result;
      }
    catch (Uncertain_conversion_exception) {}
  }
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
  typedef Exact_predicates_exact_constructions_kernel EK;
  Skin_surface_traits_3<EK> exact_traits(shrink_factor());

  typename Skin_surface_traits_3<EK>::Bare_point p_exact =
    get_anchor_point(vit->info(), exact_traits);
  return construct_surface(vit->cell()->info().first, 
                           EK() ).sign(p_exact);
}

template <class MixedComplexTraits_3> 
Sign 
Skin_surface_base_3<MixedComplexTraits_3>::
sign(const Bare_point &p, 
     const TMC_Cell_handle start) const {
  if (start == TMC_Cell_handle()) {
    return sign(p, locate_in_tmc(p,start)->info());
  } else {
    return sign(p, start->info());
  }
}

template <class MixedComplexTraits_3> 
Sign 
Skin_surface_base_3<MixedComplexTraits_3>::
sign(const Bare_point &p, const Cell_info &info) const {
  CGAL_BRANCH_PROFILER(std::string(" NGHK: failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  {
    Protect_FPU_rounding<true> P;
    try
      {
	Sign result = sign_inexact(p,info);
	if (is_certain(result))
	  return result;
      }
  catch (Uncertain_conversion_exception) {}
  }
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
  return construct_surface
    (info.first, 
     Exact_predicates_exact_constructions_kernel()).sign(p);
}

template <class MixedComplexTraits_3> 
Sign 
Skin_surface_base_3<MixedComplexTraits_3>::
sign_inexact(const Bare_point &p, const Cell_info &info) const {
  return info.second->sign(p);
}

template <class MixedComplexTraits_3> 
void
Skin_surface_base_3<MixedComplexTraits_3>::
intersect(TMC_Cell_handle ch, int i, int j,
          Bare_point &p) const {
  Cartesian_converter<FK, Gt> converter;

  Bare_point p1 = converter(ch->vertex(i)->point());
  Bare_point p2 = converter(ch->vertex(j)->point());

  return intersect(p1, p2, ch, ch, p); 
}

template <class MixedComplexTraits_3> 
void 
Skin_surface_base_3<MixedComplexTraits_3>::
intersect(Bare_point &p1, Bare_point &p2, 
          TMC_Cell_handle &s1, TMC_Cell_handle &s2,
          Bare_point &p) const
{
  typedef typename Bare_point::R  Traits;
  typedef typename Traits::FT FT;
  Cartesian_converter<Traits, 
    typename Geometric_traits::Bare_point::R> converter;

  FT sq_dist = squared_distance(p1,p2);
  // Use value to make the computation robust (endpoints near the surface)
  if (compare(s1->info(), p1, s2->info(), p2) == LARGER) {
    std::swap(p1, p2);
  }
  TMC_Cell_handle sp = s1;

  while ((s1 != s2) && (sq_dist > 1e-8)) {
    p = midpoint(p1, p2);
    sp = locate_in_tmc(converter(p), sp);

    if (sign_inexact(p, sp->info()) == NEGATIVE) { p1 = p; s1 = sp; }
    else { p2 = p; s2 = sp; }

    sq_dist *= .25;
  }
  while (sq_dist > 1e-8) {
    p = midpoint(p1, p2);
    if (sign_inexact(p, s1->info()) == NEGATIVE) { p1 = p; }
    else { p2 = p; }
    sq_dist *= .25;
  }

  p = midpoint(p1, p2);
}

template <class MixedComplexTraits_3> 
void
Skin_surface_base_3<MixedComplexTraits_3>::
intersect_with_transversal_segment(
  Bare_point &p,
  const TMC_Cell_handle &start) const 
{

  typedef typename Geometric_traits::Kernel::Plane_3 Plane;
  typedef typename Geometric_traits::Kernel::Line_3  Line;

  TMC_Cell_handle tet = start;
  if (tet == TMC_Cell_handle()) {
    tet = locate_in_tmc(p, tet);
  }
  
  // get transversal segment:
  Bare_point p1, p2;

  // Compute signs on vertices and sort them:
  int nIn = 0;
  int sortedV[4];
  for (int i=0; i<4; i++) {
    if (sign(tet->vertex(i))==POSITIVE) {
      sortedV[nIn] = i; nIn++;
    } else {
      sortedV[3-i+nIn] = i;
    }
  }

  Cartesian_converter<FK, typename Geometric_traits::Bare_point::R> converter;
  Object obj;
  typename FK::Point_3 tmc_point;
  Bare_point tet_pts[4];
  for (int i=0; i<4; i++) {
    tet_pts[i] = converter(tet->vertex(i)->point());
  }
  if (nIn==1) {
    p1 = tet_pts[sortedV[0]];
    obj = CGAL::intersection(Plane(tet_pts[sortedV[1]],
                                   tet_pts[sortedV[2]],
                                   tet_pts[sortedV[3]]),
                             Line(p1, p));
    if ( !assign(p2, obj) ) {
      CGAL_error_msg("intersection: no intersection.");
    }
  } else if (nIn==2) {
    obj = CGAL::intersection(Plane(tet_pts[sortedV[2]],
                                   tet_pts[sortedV[3]],
                                   p),
                             Line(tet_pts[sortedV[0]],
                                  tet_pts[sortedV[1]]));
    if ( !assign(p1, obj) ) {
      CGAL_error_msg("intersection: no intersection.");
    }
    obj = CGAL::intersection(Plane(tet_pts[sortedV[0]],
                                   tet_pts[sortedV[1]],
                                   p),
                             Line(tet_pts[sortedV[2]],
                                  tet_pts[sortedV[3]]));
    if ( !assign(p2, obj) ) {
      CGAL_error_msg("intersection: no intersection.");
    }
  } else if (nIn==3) {
    p2 = tet_pts[sortedV[3]];
    obj = CGAL::intersection(Plane(tet_pts[sortedV[0]],
                                   tet_pts[sortedV[1]],
                                   tet_pts[sortedV[2]]),
                             Line(p2, p));
    if ( !assign(p1, obj) ) {
      CGAL_error_msg("intersection: no intersection.");
    }
  } else {
    CGAL_error();
  }

  // Find the intersection:
  intersect(p1, p2, tet, tet, p);
}

template <class MixedComplexTraits_3> 
void 
Skin_surface_base_3<MixedComplexTraits_3>::
construct_bounding_box() 
{
  typedef typename Regular::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Regular::Geom_traits     GT;
  typedef typename GT::Point                Weighted_point;
  typedef typename GT::FT                     FT;
  
  Finite_vertices_iterator vit = regular().finite_vertices_begin();
  if (vit != regular().finite_vertices_end()) {
    Bbox_3 bbox = vit->point().bbox();
    FT max_weight=0;
    for (;vit != regular().finite_vertices_end(); ++vit) {
      bbox = bbox + vit->point().bbox();
      if (max_weight < vit->point().weight()) {
        max_weight = vit->point().weight();
      }
    }

    // add a bounding octahedron:
    FT dx = bbox.xmax() - bbox.xmin();
    FT dy = bbox.ymax() - bbox.ymin();
    FT dz = bbox.zmax() - bbox.zmin();
  
    Bare_point mid(bbox.xmin() + dx/2, 
		   bbox.ymin() + dy/2, 
		   bbox.zmin() + dz/2);
    FT dr = 
      (dx+dy+dz+sqrt(CGAL::to_double(max_weight))+.001) / gt.get_shrink();

    Weighted_point wp;
    wp = Weighted_point(Bare_point(mid.x()+dr,
				   mid.y(),
				   mid.z()),-1);
    regular().insert(wp);
    wp = Weighted_point(Bare_point(mid.x()-dr,
				   mid.y(),
				   mid.z()),-1);
    regular().insert(wp);
    wp = Weighted_point(Bare_point(mid.x(),
				   mid.y()+dr,
				   mid.z()),-1);
    regular().insert(wp);
    wp = Weighted_point(Bare_point(mid.x(),
				   mid.y()-dr,
				   mid.z()),-1);
    regular().insert(wp);
    wp = Weighted_point(Bare_point(mid.x(),
				   mid.y(),
				   mid.z()+dr),-1);
    regular().insert(wp);
    wp = Weighted_point(Bare_point(mid.x(),
				   mid.y(),
				   mid.z()-dr),-1);
    regular().insert(wp);
  }
}

template <class MixedComplexTraits_3> 
template <class Gt2>
typename Gt2::Bare_point
Skin_surface_base_3<MixedComplexTraits_3>::
get_anchor_point(const Anchor_point &anchor, Gt2 &traits) {
  typename Gt2::Bare_point p_del, p_vor;
  p_del = get_weighted_circumcenter(anchor.first, traits);
  p_vor = get_weighted_circumcenter(anchor.second, traits);
  return traits.construct_anchor_point_3_object()(p_del,p_vor);
}

template <class MixedComplexTraits_3> 
template< class Traits >
Skin_surface_quadratic_surface_3<Traits> 
Skin_surface_base_3<MixedComplexTraits_3>::
construct_surface(const Simplex &sim, const Traits &) const {
  typedef Skin_surface_quadratic_surface_3<Traits>      Quadratic_surface;
  typedef Weighted_converter_3<Cartesian_converter<
    typename Geometric_traits::Bare_point::R, Traits> > Converter;
  typedef typename Traits::Point_3                      Point;
  typedef typename Traits::FT                           FT;
  typedef CGAL::Weighted_point<Point,FT>                Weighted_point;

  Converter conv;

  switch (sim.dimension()) {
    case 0: {
      Vertex_handle vh = sim;
      return Quadratic_surface(conv(vh->point()), shrink_factor());
    }
    case 1: {
      Edge e = sim;
      Weighted_point p0 = conv(e.first->vertex(e.second)->point());
      Weighted_point p1 = conv(e.first->vertex(e.third)->point());
      return Quadratic_surface(p0, p1, shrink_factor());
    }
    case 2: {
      Facet f = sim;
      Weighted_point p0 = conv(f.first->vertex((f.second+1)&3)->point());
      Weighted_point p1 = conv(f.first->vertex((f.second+2)&3)->point());
      Weighted_point p2 = conv(f.first->vertex((f.second+3)&3)->point());
      return Quadratic_surface(p0,p1,p2, shrink_factor());
    }
    case 3: {
      Cell_handle ch = sim;
      Weighted_point p0 = conv(ch->vertex(0)->point());
      Weighted_point p1 = conv(ch->vertex(1)->point());
      Weighted_point p2 = conv(ch->vertex(2)->point());
      Weighted_point p3 = conv(ch->vertex(3)->point());
      return Quadratic_surface(p0,p1,p2,p3, shrink_factor());
    }
  }
  CGAL_error();
  return Quadratic_surface();
}

template <class MixedComplexTraits_3> 
Sign
Skin_surface_base_3<MixedComplexTraits_3>::
compare(Cell_info &info,
     const Bare_point &p1,
     const Bare_point &p2) const {
  return compare(info, p1, info, p2);
}

template <class MixedComplexTraits_3> 
Sign
Skin_surface_base_3<MixedComplexTraits_3>::
compare(Cell_info &info1,
     const Bare_point &p1,
     Cell_info &info2,
     const Bare_point &p2) const {
  CGAL_BRANCH_PROFILER(std::string(" NGHK: failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  {
    Protect_FPU_rounding<true> P;
    try
      {
	Sign result = CGAL_NTS sign(info1.second->value(p1) -
				    info2.second->value(p2));
	if (is_certain(result))
	  return result;
      }
    catch (Uncertain_conversion_exception) {}
  }
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
    
  return CGAL_NTS sign(
    construct_surface(info1.first,
                      Exact_predicates_exact_constructions_kernel()).value(p1) -
    construct_surface(info2.first,
                      Exact_predicates_exact_constructions_kernel()).value(p2));
}

template <class MixedComplexTraits_3> 
typename Skin_surface_base_3<MixedComplexTraits_3>::TMC_Cell_handle
Skin_surface_base_3<MixedComplexTraits_3>::
locate_in_tmc(const Bare_point &p0, 
              TMC_Cell_handle start) const {
  Cartesian_converter<typename Bare_point::R, TMC_Geom_traits> converter_fk;
  TMC_Point p_inexact = converter_fk(p0);

  // Make sure we continue from here with a finite cell.
  if ( start == TMC_Cell_handle() )
    start = _tmc.infinite_cell();

  int ind_inf;
  if (start->has_vertex(_tmc.infinite_vertex(), ind_inf) )
    start = start->neighbor(ind_inf);

  CGAL_triangulation_precondition(start != TMC_Cell_handle());
  CGAL_triangulation_precondition(!start->has_vertex(_tmc.infinite_vertex()));

  // We implement the remembering visibility/stochastic walk.

  // Remembers the previous cell to avoid useless orientation tests.
  TMC_Cell_handle previous = TMC_Cell_handle();
  TMC_Cell_handle c = start;

  // Now treat the cell c.
  try_next_cell:

  const TMC_Point* pts[4] = { &(c->vertex(0)->point()),
                              &(c->vertex(1)->point()),
                              &(c->vertex(2)->point()),
                              &(c->vertex(3)->point()) };

  // For the remembering stochastic walk,
  // we need to start trying with a random index :
  boost::rand48 rng;  
  boost::uniform_smallint<> four(0, 3);
  boost::variate_generator<boost::rand48&, boost::uniform_smallint<> > die4(rng, four);
  int i = die4();
  // For the remembering visibility walk (Delaunay only), we don't :
  // int i = 0;

  Orientation o;
  for (int j=0; j != 4; ++j, i = (i+1)&3) {
    TMC_Cell_handle next = c->neighbor(i);
    if (previous == next) {
      continue;
    }
    // We temporarily put p at i's place in pts.
    const TMC_Point* backup = pts[i];
    pts[i] = &p_inexact;
    {
      Protect_FPU_rounding<true> P;
      try {
	o = TMC_Geom_traits().orientation_3_object()(*pts[0], *pts[1], *pts[2], *pts[3]);
      } catch (Uncertain_conversion_exception) {
	Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
	typedef Exact_predicates_exact_constructions_kernel EK;
	Cartesian_converter<typename Geometric_traits::Bare_point::R, EK> converter_ek;
	
	Skin_surface_traits_3<EK> exact_traits(shrink_factor());
	
	typename EK::Point_3 e_pts[4];
	
	// We know that the 4 vertices of c are positively oriented.
	// So, in order to test if p is seen outside from one of c's facets,
	// we just replace the corresponding point by p in the orientation
	// test.  We do this using the array below.
	for (int k=0; k<4; k++) {
	  if (k != i) {
	    e_pts[k] = get_anchor_point(c->vertex(k)->info(), exact_traits);
	  } else {
	    e_pts[k] = converter_ek(p0);
	  }
	}
	o = orientation(e_pts[0], e_pts[1], e_pts[2], e_pts[3]);
      }
    }

    if ( o != NEGATIVE ) {
      pts[i] = backup;
      continue;
    }
    if ( next->has_vertex(_tmc.infinite_vertex()) ) {
      std::cout << "We are outside the convex hull." << std::endl;
      return next;
    }
    previous = c;
    c = next;
    goto try_next_cell;
  }
  
  CGAL_assertion(c->vertex(0) != _tmc.infinite_vertex());
  CGAL_assertion(c->vertex(1) != _tmc.infinite_vertex());
  CGAL_assertion(c->vertex(2) != _tmc.infinite_vertex());
  CGAL_assertion(c->vertex(3) != _tmc.infinite_vertex());
  
  return c;
}
  
template <class MixedComplexTraits_3> 
template <class Gt2>
typename Gt2::Bare_point
Skin_surface_base_3<MixedComplexTraits_3>::
get_weighted_circumcenter(const Simplex &s, Gt2 &traits) {
  Vertex_handle vh;
  Edge           e;
  Facet          f;
  Cell_handle   ch;

  Weighted_converter_3<
    Cartesian_converter<typename Gt::Bare_point::R, 
                        typename Gt2::Bare_point::R> > converter;

  typename Gt2::Bare_point result;
  switch(s.dimension()) {
  case 0: 
    {
      vh = s;
      result = converter(vh->point());
      break;
    }
  case 1:
    {
      e = s;
      result = traits.construct_weighted_circumcenter_3_object()
        (converter(e.first->vertex(e.second)->point()),
         converter(e.first->vertex(e.third)->point()));
      break;
    }
  case 2: 
    {
      f = s;
      result = traits.construct_weighted_circumcenter_3_object()
        (converter(f.first->vertex((f.second+1)&3)->point()),
         converter(f.first->vertex((f.second+2)&3)->point()),
         converter(f.first->vertex((f.second+3)&3)->point()));
      break;
    }
  case 3: 
    {
      ch = s;
      result = traits.construct_weighted_circumcenter_3_object()
        (converter(ch->vertex(0)->point()),
         converter(ch->vertex(1)->point()),
         converter(ch->vertex(2)->point()),
         converter(ch->vertex(3)->point()));
      break;
    }
  default:
    {
      CGAL_error();
    }
  }
  return result;
}

} //namespace CGAL 

#endif // CGAL_SKIN_SURFACE_BASE_3_H
