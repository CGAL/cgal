// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// $URL$
// $Id$
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_UNION_OF_BALLS_3_H
#define CGAL_UNION_OF_BALLS_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/triangulate_power_diagram_3.h>

// Contains the cell base of the triangulated mixed complex (and VD)
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Triangulation_simplex_3.h>

CGAL_BEGIN_NAMESPACE 

template <class UnionOfBallsTraits_3> 
class Union_of_balls_3 {
  typedef UnionOfBallsTraits_3            Gt;
  typedef Union_of_balls_3<Gt>            Self;
public:
  typedef UnionOfBallsTraits_3            Geometric_traits;
  typedef typename Gt::Weighted_point     Weighted_point;
  typedef typename Gt::Bare_point         Bare_point;
  typedef typename Gt::FT                 FT;
  // NGHK:: added for the Delaunay mesher
  typedef typename Gt::Sphere_3           Sphere;
  typedef typename Gt::Vector_3           Vector;
  
  typedef Regular_triangulation_3<Gt>     Regular;

private:
  typedef Exact_predicates_inexact_constructions_kernel     Filtered_kernel;
public:
  typedef Skin_surface_quadratic_surface_3<Filtered_kernel> 
                                                         Quadratic_surface;
public:
  typedef typename Regular::Vertex_handle                Vertex_handle;
  typedef typename Regular::Edge                         Edge;
  typedef typename Regular::Facet                        Facet;
  typedef typename Regular::Facet_circulator             Facet_circulator;
  typedef typename Regular::Cell_handle                  Cell_handle;
  typedef Triangulation_simplex_3<Regular>               Simplex;

  // pair of a del- and vor-simplex
  typedef std::pair<Simplex,Simplex>                     Anchor_point;

  //private:
  typedef typename Regular::Finite_vertices_iterator     Finite_vertices_iterator;
  typedef typename Regular::Finite_edges_iterator        Finite_edges_iterator;
  typedef typename Regular::Finite_facets_iterator       Finite_facets_iterator;
  typedef typename Regular::Finite_cells_iterator        Finite_cells_iterator;

public:
  typedef Anchor_point                                  Vertex_info;
  typedef std::pair<Simplex, Quadratic_surface *>       Cell_info;
private:
  // Triangulated_mixed_complex:
  typedef Simple_cartesian<Interval_nt_advanced>                       FK;
  typedef Triangulation_vertex_base_with_info_3<Vertex_info, FK>       Vb;
  typedef Triangulation_cell_base_with_info_3<Cell_info, FK>           Cb;
  typedef Triangulation_data_structure_3<Vb,Cb>                        Tds;
public:
  typedef Triangulation_3<FK, Tds>                                     TMC;
private:
  typedef typename TMC::Finite_vertices_iterator TMC_Vertex_iterator;
  typedef typename TMC::Finite_cells_iterator    TMC_Cell_iterator;
  typedef typename TMC::Vertex_handle            TMC_Vertex_handle;
  typedef typename TMC::Cell_handle              TMC_Cell_handle;
  typedef typename TMC::Point                    TMC_Point;
////public:
////  typedef UnionOfBallsTraits_3           Gt;
////  typedef typename Gt::Weighted_point   Weighted_point;
////  typedef typename Gt::Bare_point       Bare_point;
////  typedef typename Gt::RT               RT;
////  
////  typedef Regular_triangulation_3<Gt>     Regular;
//////   typedef Triangulation_data_structure_3 <
//////     Triangulation_vertex_base_3<GT>,
//////     Triangulated_mixed_complex_cell_3<GT, PolyhedronKernel_3> >
//////                                           Triangulated_mixed_complex_tds;
////
////  // defining the triangulated mixed complex:
////  typedef Exact_predicates_inexact_constructions_kernel    TMC_Traits;
////public:
////
////#ifdef CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE
////  typedef Skin_surface_quadratic_surface_3<TMC_Traits>   Quadratic_surface;
////#else
////  typedef Skin_surface_quadratic_surface_3<Simple_cartesian<double> > 
////                                                         Quadratic_surface;
////#endif // CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE
////
////  typedef Triangulation_3<
////    TMC_Traits,
////    Triangulation_data_structure_3
////    < Triangulated_mixed_complex_vertex_3<TMC_Traits>,
////      Triangulated_mixed_complex_cell_3<TMC_Traits,Quadratic_surface> > 
////  >                                                   Triangulated_mixed_complex;
////  typedef typename Triangulated_mixed_complex::Vertex_handle TMC_Vertex_handle;
////  typedef typename Triangulated_mixed_complex::Cell_handle   TMC_Cell_handle;
////  typedef typename TMC_Traits::Point_3                       TMC_Point;

public:  
  template < class WP_iterator >
  Union_of_balls_3(WP_iterator begin, WP_iterator end, 
		 Gt gt = Gt(),
		 bool verbose = false
		 ) 
    : gt(gt), verbose(verbose) {

    CGAL_assertion(begin != end);

    Regular regular;
    regular.insert(begin, end);
    construct_bounding_box(regular);
  
    if (verbose) {
      std::cerr << "Triangulation ready" << std::endl;
      std::cerr << "Vertices: " << regular.number_of_vertices() << std::endl;
      std::cerr << "Cells:    " << regular.number_of_cells() << std::endl;
    }
    
    // Construct the triangulated mixed complex:
    Triangulated_mixed_complex_observer_3<TMC, Self> observer(1);
    triangulate_power_diagram_3(regular, tmc, observer, verbose);
    
    CGAL_assertion(tmc.is_valid());
    if (verbose) {
      std::cerr << "Triangulated mixed complex ready" << std::endl;
      std::cerr << "Vertices: " << tmc.number_of_vertices() << std::endl;
      std::cerr << "Cells:    " << tmc.number_of_cells() << std::endl;
    }
  }

  FT shrink_factor() const {
    return 1;
  }

  template <class Polyhedron_3>
  void mesh_union_of_balls_3(Polyhedron_3 &p) const;

  template <class Polyhedron_3>
  void subdivide_union_of_balls_mesh_3(Polyhedron_3 &p) const;
  
  const TMC &triangulated_mixed_complex() const {
    return tmc;
  }
  
  TMC_Cell_handle locate(const TMC_Point &p) const{
    return tmc.locate(p);
  }

  Sign sign(const TMC_Cell_handle &ch, const Bare_point &p) const {
    return sign(ch->info(), p);
  }
  Sign sign(Cell_info &info, const Bare_point &p) const {
    try
      {
        CGAL_PROFILER(std::string("NGHK: calls to    : ") + 
                      std::string(CGAL_PRETTY_FUNCTION));
        Protect_FPU_rounding<true> P;
        Sign result = info.second->sign(p);
        if (! is_indeterminate(result))
          return result;
      }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("NGHK: failures of : ") + 
                  std::string(CGAL_PRETTY_FUNCTION));
    Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
    return construct_surface
      (info.first, 
       Exact_predicates_exact_constructions_kernel()).sign(p);
  }
  
  Sign sign(TMC_Vertex_handle vit) const {
    CGAL_assertion(!tmc.is_infinite(vit));
    TMC_Cell_handle ch = vit->cell();
    if (tmc.is_infinite(ch)) {
      std::vector<TMC_Cell_handle> nbs;
      tmc.incident_cells(vit, std::back_inserter(nbs));
      for (typename std::vector<TMC_Cell_handle>::iterator it = nbs.begin();
           tmc.is_infinite(ch) && (it != nbs.end());
           it++) {
        ch = *it;
      }
    }
    CGAL_assertion(!tmc.is_infinite(ch));
  
    // don't use sign, since the point is constructed:
    try
      {
        CGAL_PROFILER(std::string("NGHK: calls to    : ") + 
                      std::string(CGAL_PRETTY_FUNCTION));
        Protect_FPU_rounding<true> P;
        Sign result = vit->cell()->info().second->sign(vit->point());
        if (! is_indeterminate(result))
          return result;
      }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("NGHK: failures of : ") + 
                  std::string(CGAL_PRETTY_FUNCTION));
    Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
    typedef Exact_predicates_exact_constructions_kernel EK;
    Skin_surface_traits_3<EK> exact_traits(1);

    typename Skin_surface_traits_3<EK>::Bare_point p_exact =
      get_anchor_point(vit->info(), exact_traits);
    return construct_surface(vit->cell()->info().first, 
                             EK() ).sign(p_exact);
  }
  
  void intersect_with_transversal_segment
    (Bare_point &p,
     const TMC_Cell_handle &start = TMC_Cell_handle()) const;

  TMC_Cell_handle
  locate_mixed(const Bare_point &p, 
               TMC_Cell_handle start = TMC_Cell_handle()) const;
  

  Quadratic_surface
  construct_surface(const Simplex &sim) const {
    return construct_surface(sim, typename Geometric_traits::Kernel());
  }
  template< class Traits >
  Skin_surface_quadratic_surface_3<Traits> 
  construct_surface(const Simplex &sim, const Traits &traits) const;

  template <class Gt2>
  typename Gt2::Bare_point
  get_anchor_point(const Anchor_point &anchor, Gt2 &traits) const {
    return get_weighted_circumcenter(anchor.first, traits);
  }
  
  template <class Gt2>
  typename Gt2::Bare_point
  get_weighted_circumcenter(const Simplex &s, Gt2 &traits) const {
    Cartesian_converter<typename Gt::Bare_point::R, 
      typename Gt2::Bare_point::R> converter;
    switch(s.dimension()) {
    case 0: 
      {
        Vertex_handle vh = s;
        typename Gt::Weighted_point wp = vh->point();
        return converter(wp.point());
      }
    case 1:
      {
        Edge e = s;
        return traits.construct_weighted_circumcenter_3_object()
          (converter(e.first->vertex(e.second)->point()),
           converter(e.first->vertex(e.third)->point()));
      }
    case 2: 
      {
        Facet f = s;
        return traits.construct_weighted_circumcenter_3_object()
          (converter(f.first->vertex((f.second+1)&3)->point()),
           converter(f.first->vertex((f.second+2)&3)->point()),
           converter(f.first->vertex((f.second+3)&3)->point()));
      }
    case 3: 
      {
        Cell_handle ch = s;
        return traits.construct_weighted_circumcenter_3_object()
          (converter(ch->vertex(0)->point()),
           converter(ch->vertex(1)->point()),
           converter(ch->vertex(2)->point()),
           converter(ch->vertex(3)->point()));
      }
    }
    CGAL_assertion(false);
    return typename Gt2::Point_3();
  }

  void intersect(TMC_Cell_handle ch, int i, int j,
                 //TMC_Vertex_handle p2, 
                 Bare_point &p) const {
//    typedef typename Bare_point::R  Traits;
//    typedef typename Traits::FT FT;
//    Cartesian_converter<FK, 
//      typename Geometric_traits::Bare_point::R> converter;
//
//    Cell_info &cell_info = ch->info();
//    Bare_point p1 = converter(ch->vertex(i)->point());
//    Bare_point p2 = converter(ch->vertex(j)->point());
//
//    FT sq_dist = squared_distance(p1,p2);
//    // Use value to make the computation robust (endpoints near the surface)
//    if (less(cell_info, p2, p1)) std::swap(p1, p2);
//    while (sq_dist > 1e-8) {
//      p = midpoint(p1, p2);
//      if (sign(cell_info, p) == NEGATIVE) { p1 = p; }
//      else { p2 = p; }
//      sq_dist *= .25;
//    }
//
//    p = midpoint(p1, p2);
  }

  FT
  less(Cell_info &info1,
       const Bare_point &p1,
       Cell_info &info2,
       const Bare_point &p2) const {
    try
      {
        CGAL_PROFILER(std::string("NGHK: calls to    : ") + 
                      std::string(CGAL_PRETTY_FUNCTION));
        Protect_FPU_rounding<true> P;
        Sign result = CGAL_NTS sign(info2.second->value(p2) -
                                    info1.second->value(p1));
        if (! is_indeterminate(result))
          return result==POSITIVE;
      }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("NGHK: failures of : ") + 
                  std::string(CGAL_PRETTY_FUNCTION));
    Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
      
    return 
      CGAL_NTS sign(construct_surface(info2.first,
                                      Exact_predicates_exact_constructions_kernel()
                                      ).value(p2) -
                    construct_surface(info1.first,
                                      Exact_predicates_exact_constructions_kernel()
                                      ).value(p1));
  }
  void intersect(Bare_point &p1, Bare_point &p2, 
                 TMC_Cell_handle &s1, TMC_Cell_handle &s2,
                 Bare_point &p) const {
    typedef typename Bare_point::R  Traits;
    typedef typename Traits::FT FT;
    Cartesian_converter<Traits, 
      typename Geometric_traits::Bare_point::R> converter;

    FT sq_dist = squared_distance(p1,p2);
    // Use value to make the computation robust (endpoints near the surface)
    if (less(s2->info(), p2, s1->info(), p1)) std::swap(p1, p2);
    TMC_Cell_handle sp = s1;

    while ((s1 != s2) && (sq_dist > 1e-8)) {
      p = midpoint(p1, p2);
      sp = locate_mixed(converter(p), sp);

      if (sign(sp, p) == NEGATIVE) { p1 = p; s1 = sp; }
      else { p2 = p; s2 = sp; }

      sq_dist *= .25;
    }
    while (sq_dist > 1e-8) {
      p = midpoint(p1, p2);
      if (sign(s1, p) == NEGATIVE) { p1 = p; }
      else { p2 = p; }
      sq_dist *= .25;
    }

    p = midpoint(p1, p2);
  }

private:
  void construct_bounding_box(Regular &regular);

  Gt &gt;
  TMC tmc;
  bool verbose;
  mutable Random rng;
};

template <class UnionOfBallsTraits_3> 
void 
Union_of_balls_3<UnionOfBallsTraits_3>::
construct_bounding_box(Regular &regular) 
{
  typedef typename Regular::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Regular::Geom_traits     GT;
  typedef typename GT::Bare_point             Point;
  typedef typename GT::Point                Weighted_point;
  typedef typename GT::RT                     RT;
  
  Finite_vertices_iterator vit = regular.finite_vertices_begin();
  if (vit != regular.finite_vertices_end()) {
    Bbox_3 bbox = vit->point().bbox();
    RT max_weight=vit->point().weight();
    while (++vit != regular.finite_vertices_end()) {
      bbox = bbox + vit->point().bbox();
      if (max_weight < vit->point().weight())
	max_weight = vit->point().weight();
    }

    // add a bounding octahedron:
    RT dx = bbox.xmax() - bbox.xmin();
    RT dy = bbox.ymax() - bbox.ymin();
    RT dz = bbox.zmax() - bbox.zmin();
  
    Bare_point mid(bbox.xmin() + dx/2, bbox.ymin() + dy/2, bbox.zmin() + dz/2);
    RT dr = sqrt(CGAL::to_double(max_weight)) + 1;
  
    regular.insert(Weighted_point(
      Bare_point(bbox.xmax()+(dy+dz+dr),mid.y(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(bbox.xmin()-(dy+dz+dr),mid.y(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),bbox.ymax()+(dx+dz+dr),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),bbox.ymin()-(dx+dz+dr),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),mid.y(),bbox.zmax()+(dx+dy+dr)),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),mid.y(),bbox.zmin()-(dx+dy+dr)),-1));
  }
}

template <class MixedComplexTraits_3> 
template <class Polyhedron_3>
void
Union_of_balls_3<MixedComplexTraits_3>::mesh_union_of_balls_3(Polyhedron_3 &p) const {
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
  marching_tetrahedra_3(triangulated_mixed_complex().finite_vertices_begin(), 
			triangulated_mixed_complex().finite_vertices_end(), 
			triangulated_mixed_complex().finite_cells_begin(), 
			triangulated_mixed_complex().finite_cells_end(), 
			p, 
			marching_traits,
			marching_observer);
}

template <class MixedComplexTraits_3> 
template< class Traits >
Skin_surface_quadratic_surface_3<Traits> 
Union_of_balls_3<MixedComplexTraits_3>::
construct_surface(const Simplex &sim, const Traits &traits) const {
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
      return Quadratic_surface(conv(vh->point()), 1);
      break;
    }
    case 1: {
      Edge e = sim;
      Weighted_point p0 = conv(e.first->vertex(e.second)->point());
      Weighted_point p1 = conv(e.first->vertex(e.third)->point());
      return Quadratic_surface(p0, p1, 1);
      break;
    }
    case 2: {
      Facet f = sim;
      Weighted_point p0 = conv(f.first->vertex((f.second+1)&3)->point());
      Weighted_point p1 = conv(f.first->vertex((f.second+2)&3)->point());
      Weighted_point p2 = conv(f.first->vertex((f.second+3)&3)->point());
      return Quadratic_surface(p0,p1,p2, 1);
      break;
    }
    case 3: {
      Cell_handle ch = sim;
      Weighted_point p0 = conv(ch->vertex(0)->point());
      Weighted_point p1 = conv(ch->vertex(1)->point());
      Weighted_point p2 = conv(ch->vertex(2)->point());
      Weighted_point p3 = conv(ch->vertex(3)->point());
      return Quadratic_surface(p0,p1,p2,p3, 1);
      break;
    }
  }
  CGAL_assertion(false);
  return Quadratic_surface();
}

template <class MixedComplexTraits_3> 
template <class Polyhedron_3>
void
Union_of_balls_3<MixedComplexTraits_3>::subdivide_union_of_balls_mesh_3(Polyhedron_3 &p) const {
  std::cout << "Skin_Surface_3.subdivide_skin_surface_mesh_3(p)" << std::endl;

  typedef Skin_surface_refinement_policy_3<Self, Polyhedron_3> Policy;
  typedef Skin_surface_sqrt3<Self, Polyhedron_3, Policy>       Subdivider;

  Policy policy(*this);
  Subdivider subdivider(*this, p, policy);
  subdivider.subdivide();
}

template <class MixedComplexTraits_3> 
void
Union_of_balls_3<MixedComplexTraits_3>::
intersect_with_transversal_segment(
  Bare_point &p,
  const TMC_Cell_handle &start) const 
{

  typedef typename Geometric_traits::Kernel::Plane_3 Plane;
  typedef typename Geometric_traits::Kernel::Line_3  Line;

  TMC_Cell_handle tet = locate_mixed(p, start);
  
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
      CGAL_assertion_msg(false,"intersection: no intersection.");
    }
  } else if (nIn==2) {
    obj = CGAL::intersection(Plane(tet_pts[sortedV[2]],
                                   tet_pts[sortedV[3]],
                                   p),
                             Line(tet_pts[sortedV[0]],
                                  tet_pts[sortedV[1]]));
    if ( !assign(p1, obj) ) {
      CGAL_assertion_msg(false,"intersection: no intersection.");
    }
    obj = CGAL::intersection(Plane(tet_pts[sortedV[0]],
                                   tet_pts[sortedV[1]],
                                   p),
                             Line(tet_pts[sortedV[2]],
                                  tet_pts[sortedV[3]]));
    if ( !assign(p2, obj) ) {
      CGAL_assertion_msg(false,"intersection: no intersection.");
    }
  } else if (nIn==3) {
    p2 = tet_pts[sortedV[3]];
    obj = CGAL::intersection(Plane(tet_pts[sortedV[0]],
                                   tet_pts[sortedV[1]],
                                   tet_pts[sortedV[2]]),
                             Line(p2, p));
    if ( !assign(p1, obj) ) {
      CGAL_assertion_msg(false,"intersection: no intersection.");
    }
  } else {
    CGAL_assertion(false);
  }

  // Find the intersection:
  intersect(p1, p2, tet, tet, p);
}

template <class MixedComplexTraits_3> 
typename Union_of_balls_3<MixedComplexTraits_3>::TMC_Cell_handle
Union_of_balls_3<MixedComplexTraits_3>::
locate_mixed(const Bare_point &p0, 
             TMC_Cell_handle start) const {
  Cartesian_converter<typename Geometric_traits::Bare_point::R, FK> converter_fk;
  typename FK::Point_3 p_inexact = converter_fk(p0);

  Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);

  // Make sure we continue from here with a finite cell.
  if ( start == TMC_Cell_handle() )
    start = tmc.infinite_cell();

  int ind_inf;
  if (start->has_vertex(tmc.infinite_vertex(), ind_inf) )
    start = start->neighbor(ind_inf);

  CGAL_triangulation_precondition(start != TMC_Cell_handle());
  CGAL_triangulation_precondition(!start->has_vertex(tmc.infinite_vertex()));

  // We implement the remembering visibility/stochastic walk.

  // Remembers the previous cell to avoid useless orientation tests.
  TMC_Cell_handle previous = TMC_Cell_handle();
  TMC_Cell_handle c = start;

  // Now treat the cell c.
             try_next_cell:

  const typename FK::Point_3* pts[4] = { &(c->vertex(0)->point()),
                                         &(c->vertex(1)->point()),
                                         &(c->vertex(2)->point()),
                                         &(c->vertex(3)->point()) };

  // For the remembering stochastic walk,
  // we need to start trying with a random index :
  int i = rng.template get_bits<2>();
  // For the remembering visibility walk (Delaunay only), we don't :
  // int i = 0;

  Orientation o;
  for (int j=0; j != 4; ++j, i = (i+1)&3) {
    TMC_Cell_handle next = c->neighbor(i);
    if (previous == next) {
      continue;
    }
    // We temporarily put p at i's place in pts.
    const typename FK::Point_3* backup = pts[i];
    pts[i] = &p_inexact;
    try {
      o = orientation(*pts[0], *pts[1], *pts[2], *pts[3]);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      typedef Exact_predicates_exact_constructions_kernel EK;
      Cartesian_converter<typename Geometric_traits::Bare_point::R, EK> converter_ek;

      Skin_surface_traits_3<EK> exact_traits(shrink_factor());
   
      typename EK::Point_3 pts[4];

      // We know that the 4 vertices of c are positively oriented.
      // So, in order to test if p is seen outside from one of c's facets,
      // we just replace the corresponding point by p in the orientation
      // test.  We do this using the array below.
      for (int j=0; j<4; j++) {
        if (j != i) {
          pts[j] = get_anchor_point(c->vertex(j)->info(), exact_traits);
        } else {
          pts[j] = converter_ek(p0);
        }
      }
 

    }
    if ( o != NEGATIVE ) {
      pts[i] = backup;
      continue;
    }
    if ( next->has_vertex(tmc.infinite_vertex()) ) {
      std::cout << "We are outside the convex hull." << std::endl;
      return next;
    }
    previous = c;
    c = next;
    goto try_next_cell;
  }
  
  CGAL_assertion(c->vertex(0) != tmc.infinite_vertex());
  CGAL_assertion(c->vertex(1) != tmc.infinite_vertex());
  CGAL_assertion(c->vertex(2) != tmc.infinite_vertex());
  CGAL_assertion(c->vertex(3) != tmc.infinite_vertex());

  return c;
}

CGAL_END_NAMESPACE

#endif // CGAL_UNION_OF_BALLS_3_H
