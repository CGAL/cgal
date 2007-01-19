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

#ifndef CGAL_SKIN_SURFACE_3_H
#define CGAL_SKIN_SURFACE_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/triangulate_mixed_complex_3.h>

// Needed for the (Delaunay) surface mesher
#include <CGAL/Skin_surface_mesher_oracle_3.h>
#include <CGAL/Triangulation_simplex_3.h>

// For point location in the mixed complex
#include <CGAL/Random.h>
#include <CGAL/Skin_surface_traits_3.h>

#include <CGAL/Skin_surface_marching_tetrahedra_observer_3.h>
#include <CGAL/Skin_surface_refinement_policy_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>

CGAL_BEGIN_NAMESPACE 

template <class MixedComplexTraits_3> 
class Skin_surface_3 {
  typedef MixedComplexTraits_3            Gt;
  typedef Skin_surface_3<Gt>              Self;
public:
  typedef MixedComplexTraits_3            Geometric_traits;
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

private:
  // NGHK: added for the (Delaunay) surface mesher, document
  typedef Exact_predicates_inexact_constructions_kernel  Mesher_Gt;
  typedef Skin_surface_mesher_oracle_3<Mesher_Gt,Self> Surface_mesher_traits_3;

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
public:
  template < class WP_iterator >
  Skin_surface_3(WP_iterator begin, WP_iterator end, 
		 FT shrink,
		 bool grow_balls = true,
		 Gt gt_ = Gt(),
		 bool _verbose = false
		 ) 
    : gt(gt_), verbose(_verbose) {

    gt.set_shrink(shrink);
    CGAL_assertion(begin != end);

    if (grow_balls) {
      for (; begin != end; begin++) {
	regular.insert(Weighted_point(*begin, begin->weight()/shrink_factor()));
      }
    } else {
      regular.insert(begin, end);
    }
    construct_bounding_box(regular);
  
    if (verbose) {
      std::cerr << "Triangulation ready" << std::endl;
      std::cerr << "Vertices: " << regular.number_of_vertices() << std::endl;
      std::cerr << "Cells:    " << regular.number_of_cells() << std::endl;
    }
    

    // Construct the Triangulated_mixed_complex:
    Triangulated_mixed_complex_observer_3<TMC, Self> observer(shrink_factor());
    triangulate_mixed_complex_3(regular, shrink_factor(), tmc, observer, true);

    CGAL_assertion(tmc.dimension() == 3);
    { // NGHK: debug code:
      CGAL_assertion(tmc.is_valid());
      std::vector<TMC_Vertex_handle> ch_vertices;
      tmc.incident_vertices(tmc.infinite_vertex(), 
			    std::back_inserter(ch_vertices));
      for (typename std::vector<TMC_Vertex_handle>::iterator
	     vit = ch_vertices.begin(); vit != ch_vertices.end(); vit++) {
	CGAL_assertion(sign(*vit) == POSITIVE);
      }
    }
				
//     mc_triangulator = new CMCT(regular, verbose);

  }

  // This class has to be a friend:
  template <class SkinSurface_3, 
	    class Vertex_iterator, class Cell_iterator, 
	    class HalfedgeDS>
  friend class Marching_tetrahedra_traits_skin_surface_3;

  template <class Polyhedron_3>
  void mesh_skin_surface_3(Polyhedron_3 &p) const;
  
  // This class has to be a friend:
  template <class SkinSurface_3, class Polyhedron_3>
  friend class Skin_surface_refinement_policy_3;

  template <class Polyhedron_3>
  void subdivide_skin_surface_mesh_3(Polyhedron_3 &p) const;
  
  Sign sign(const Bare_point &p, 
	    const TMC_Cell_handle start = TMC_Cell_handle()) const {
    if (start == TMC_Cell_handle()) {
      return sign(locate_mixed(p,start), p);
    } else {
      return sign(start, p);
    }
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
    Skin_surface_traits_3<EK> exact_traits(shrink_factor());

    typename Skin_surface_traits_3<EK>::Bare_point p_exact =
      get_anchor_point(vit->info(), exact_traits);
    return construct_surface(vit->cell()->info().first, 
			     EK() ).sign(p_exact);
  }
  Vector
  normal(const Bare_point &p, 
	 const TMC_Cell_handle start = TMC_Cell_handle()) const {
    if (start == TMC_Cell_handle()) {
      return get_normal(locate_mixed(p,start), p);
    } else {
      return get_normal(start, p);
    }
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

  template <class Gt2>
  typename Gt2::Bare_point
  get_anchor_point(const Anchor_point &anchor, Gt2 &traits) const {
    typename Gt2::Bare_point p_del, p_vor;
    p_del = get_weighted_circumcenter(anchor.first, traits);
    p_vor = get_weighted_circumcenter(anchor.second, traits);
    return traits.construct_anchor_point_3_object()(p_del,p_vor);
  }
private:
  bool is_infinite_mixed_cell(const Simplex &s) const {
    switch (s.dimension()) {
    case 0:
      {
	Vertex_handle vh = s;
	std::list<Vertex_handle> nbs;
	regular.incident_vertices(vh, std::back_inserter(nbs));
	for (typename std::list<Vertex_handle>::iterator it = nbs.begin();
	     it != nbs.end(); it ++) {
	  if (regular.is_infinite(*it)) return true;
	}
	return false;
      }
    case 1:
      {
	Edge e = s;
	Facet_circulator fcir, fstart;
	fcir = fstart = regular.incident_facets(e);
	do {
	  if (regular.is_infinite(*fcir)) return true;
	} while (++fcir != fstart);
	return false;
      }
    case 2:
      {
	Facet f = s;
	return (regular.is_infinite(f.first) ||
		regular.is_infinite(f.first->neighbor(f.second)));
      }
    case 3:
      {
	return false;
      }
    default:
      {
	CGAL_assertion(false);
      }
    }
    CGAL_assertion(false);
    return false;
  }
//   CMCT_Cell locate_tet(const Bare_point &p, const Simplex &s) const {
//     Skin_surface_traits_3<Exact_predicates_exact_constructions_kernel> 
//       traits(shrink_factor());
    
//     std::vector<CMCT_Cell> cells;

//     switch (s.dimension()) {
//     case 0:
//       {
// 	Vertex_handle vh = s;
// 	CGAL_assertion(!regular.is_infinite(vh));
// 	mc_triangulator->construct_0_cell(vh, std::back_inserter(cells));
// 	break;
//       }
//     case 1:
//       {
// 	Edge e = s;
// 	CGAL_assertion(!regular.is_infinite(e));
// 	mc_triangulator->construct_1_cell(e, std::back_inserter(cells));
// 	break;
//       }
//     case 2:
//       {
// 	Facet f = s;
// 	CGAL_assertion(!regular.is_infinite(f));
// 	mc_triangulator->construct_2_cell(f, std::back_inserter(cells));
// 	break;
//       }
//     case 3:
//       {
// 	Cell_handle ch = s;
// 	CGAL_assertion(!regular.is_infinite(ch));
// 	mc_triangulator->construct_3_cell(ch, std::back_inserter(cells));
// 	break;
//       }
//     default:
//       {
// 	CGAL_assertion(false);
//       }
//     }


//     bool found = false;

//     for (typename std::vector<CMCT_Cell>::iterator it = cells.begin();
// 	 ((!found)&&(it != cells.end())); it++) {
//       if (mc_triangulator->bounded_side(p, *it, traits) != ON_UNBOUNDED_SIDE) {
// 	return *it;
//       }
//     }

//     return CMCT_Cell();
//   }
public:
  TMC_Cell_handle
  locate_mixed(const Bare_point &p, 
	       TMC_Cell_handle start = TMC_Cell_handle()) const;
  
  // exact computation of the sign on a vertex of the TMC
//   Sign sign(const CMCT_Vertex_handle vh) const {
//     typedef Exact_predicates_exact_constructions_kernel K;
//     Skin_surface_traits_3<K> traits(shrink_factor());
    
//     typename K::Point_3 p = mc_triangulator->location(vh, traits);

//     return construct_surface(vh->first, K()).sign(p);
//   }

  // Trivial caching: check wether the surface is the same as the previous:
//   mutable Skin_surface_quadratic_surface_3<
//     Simple_cartesian<Interval_nt_advanced> > previous_sign_surface;
//   mutable Simplex                            previous_sign_simplex;

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
  FT
  less(Cell_info &info,
       const Bare_point &p1,
       const Bare_point &p2) const {
    try
    {
      CGAL_PROFILER(std::string("NGHK: calls to    : ") + 
		    std::string(CGAL_PRETTY_FUNCTION));
      Protect_FPU_rounding<true> P;
      Sign result = CGAL_NTS sign(info.second->value(p2) -
				  info.second->value(p1));
      if (! is_indeterminate(result))
        return result==POSITIVE;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("NGHK: failures of : ") + 
		  std::string(CGAL_PRETTY_FUNCTION));
    Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
      
    return 
      CGAL_NTS sign(construct_surface(info.first,
				      Exact_predicates_exact_constructions_kernel()
				      ).value(p2) -
		    construct_surface(info.first,
				      Exact_predicates_exact_constructions_kernel()
				      ).value(p1));
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
  FT
  value(const Bare_point &p) const {
    TMC_Cell_handle ch = locate_mixed(p);
    return value(Simplex(ch),p);
  }

  FT
  value(const Simplex &sim, const Bare_point &p) const {
    return 
      construct_surface(sim, typename Geometric_traits::Kernel()).value(p);
  }
  FT
  value(TMC_Cell_handle ch, const Bare_point &p) const {
    CGAL_assertion(!tmc.is_infinite(ch));
    
    return value(ch->info(), p);
  }
  Vector
  get_normal(TMC_Cell_handle ch, const Bare_point &p) const {
    CGAL_assertion(!tmc.is_infinite(ch));
    
    return ch->info().second->gradient(p);
  }

  // Move the point in the direction of the gradient
  void to_surface(Bare_point &p,
		  const TMC_Cell_handle &start = TMC_Cell_handle()) const {
    Bare_point p1 = p;
    TMC_Cell_handle ch1 = start;
    if (start != TMC_Cell_handle()) {
      ch1 = locate_mixed(p,ch1);
    }
    Sign sign1 = sign(ch1, p1);

    Vector n = get_normal(ch1,p);
    if (sign1 == POSITIVE) n = -value(ch1,p)*n;

    int k=2;
    Bare_point p2 = p+k*n;
    TMC_Cell_handle ch2 = locate_mixed(p2, ch1);
    while (sign(ch2,p2) == sign1) {
      k++;
      p1 = p2;
      ch1 = ch2;
      p2 = p+k*n;
      ch2 = locate_mixed(p2, ch2);
    }
    intersect(p1,p2, ch1,ch2, p);
  }
//   void intersect(const CMCT_Vertex_handle vh1,
// 		 const CMCT_Vertex_handle vh2,
// 		 Bare_point &p) const {
//     Bare_point p1 = mc_triangulator->location(vh1, gt);
//     Bare_point p2 = mc_triangulator->location(vh2, gt);
//     Simplex s1 = vh1->first;
//     Simplex s2 = vh2->first;
//     intersect(p1,p2, s1,s2, p);
//   }

//   void intersect(const CMCT_Vertex_handle vh1,
// 		 const CMCT_Vertex_handle vh2,
// 		 const Simplex &s,
// 		 Bare_point &p) const {
//     Bare_point p1 = mc_triangulator->location(vh1, gt);
//     Bare_point p2 = mc_triangulator->location(vh2, gt);
//     Simplex sp = s;
//     intersect(p1,p2, sp,sp, p);
//   }

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

  void intersect(TMC_Cell_handle ch, int i, int j,
		 //TMC_Vertex_handle p2, 
		 Bare_point &p) const {
    typedef typename Bare_point::R  Traits;
    typedef typename Traits::FT FT;
    Cartesian_converter<FK, 
                        typename Geometric_traits::Bare_point::R> converter;

    Cell_info &cell_info = ch->info();
    Bare_point p1 = converter(ch->vertex(i)->point());
    Bare_point p2 = converter(ch->vertex(j)->point());

    FT sq_dist = squared_distance(p1,p2);
    // Use value to make the computation robust (endpoints near the surface)
    if (less(cell_info, p2, p1)) std::swap(p1, p2);
    while (sq_dist > 1e-8) {
      p = midpoint(p1, p2);
      if (sign(cell_info, p) == NEGATIVE) { p1 = p; }
      else { p2 = p; }
      sq_dist *= .25;
    }

    p = midpoint(p1, p2);
  }

  void intersect_with_transversal_segment
  (Bare_point &p,
   const TMC_Cell_handle &start = TMC_Cell_handle()) const 
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

  Quadratic_surface
  construct_surface(const Simplex &sim) const {
    return construct_surface(sim, typename Geometric_traits::Kernel());
  }
  template< class Traits >
  Skin_surface_quadratic_surface_3<Traits> 
  construct_surface(const Simplex &sim, const Traits &traits) const {
    typedef Skin_surface_quadratic_surface_3<Traits>      Quadratic_surface;
    typedef Weighted_converter_3<Cartesian_converter<
      typename Geometric_traits::Bare_point::R, Traits> > Converter;
    typedef typename Traits::Point_3                      Point;
    typedef typename Traits::FT                           FT;
    typedef CGAL::Weighted_point<Point,FT>                Weighted_point;

    Converter conv;

    switch (sim.dimension()) {
    case 0:
      {
	Vertex_handle vh = sim;
	return Quadratic_surface(conv(vh->point()), shrink_factor());
	break;
      }
    case 1:
      {
	Edge e = sim;
	Weighted_point p0 = conv(e.first->vertex(e.second)->point());
	Weighted_point p1 = conv(e.first->vertex(e.third)->point());
	return Quadratic_surface(p0, p1, shrink_factor());
	break;
      }
    case 2:
      {
	Facet f = sim;
	Weighted_point p0 = conv(f.first->vertex((f.second+1)&3)->point());
	Weighted_point p1 = conv(f.first->vertex((f.second+2)&3)->point());
	Weighted_point p2 = conv(f.first->vertex((f.second+3)&3)->point());
	return Quadratic_surface(p0,p1,p2, shrink_factor());
	break;
      }
    case 3:
      {
	Cell_handle ch = sim;
	Weighted_point p0 = conv(ch->vertex(0)->point());
	Weighted_point p1 = conv(ch->vertex(1)->point());
	Weighted_point p2 = conv(ch->vertex(2)->point());
	Weighted_point p3 = conv(ch->vertex(3)->point());
	return Quadratic_surface(p0,p1,p2,p3, shrink_factor());
	break;
      }
    }
    CGAL_assertion(false);
    return Quadratic_surface();
  }

  // Access to the implicit triangulated mixed complex:
  TMC_Vertex_iterator tmc_vertices_begin() const 
  { return tmc.finite_vertices_begin(); }
  TMC_Vertex_iterator tmc_vertices_end() const 
  { return tmc.finite_vertices_end(); }
  TMC_Cell_iterator tmc_cells_begin() const 
  { return tmc.finite_cells_begin(); }
  TMC_Cell_iterator tmc_cells_end() const
  { return tmc.finite_cells_end(); }

  // NGHK: added for the (Delaunay) surface mesher, document
  Sphere bounding_sphere() const {
    return _bounding_sphere;
  }
  FT squared_error_bound() const {
    return .01;
  }

  typename Mesher_Gt::FT 
  get_density(const typename Mesher_Gt::Point_3 &p) const {
    // NGHK: Make adaptive
    return 1;
  }
public:
  const Regular &get_regular_triangulation() const {
    return regular;
  }
  FT shrink_factor() const {
    return gt.get_shrink();
  }
  const TMC &triangulated_mixed_complex() const {
    return tmc;
  }

private:
  void construct_bounding_box(Regular &regular);

  Regular regular;
  Gt gt;
  bool verbose;
  Sphere _bounding_sphere;
  mutable Random rng;

  // Triangulated mixed complex:
  TMC tmc;
  // We want to construct this object later (the pointer):
//   CMCT *mc_triangulator;
};

template <class MixedComplexTraits_3> 
void 
Skin_surface_3<MixedComplexTraits_3>::
construct_bounding_box(Regular &regular) 
{
  typedef typename Regular::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Regular::Geom_traits     GT;
  typedef typename GT::Bare_point             Point;
  typedef typename GT::Point                Weighted_point;
  typedef typename GT::FT                     FT;
  
  Finite_vertices_iterator vit = regular.finite_vertices_begin();
  if (vit != regular.finite_vertices_end()) {
    Bbox_3 bbox = vit->point().bbox();
    FT max_weight=vit->point().weight();
    while (++vit != regular.finite_vertices_end()) {
      bbox = bbox + vit->point().bbox();
      if (max_weight < vit->point().weight())
	max_weight = vit->point().weight();
    }

    // add a bounding octahedron:
    FT dx = bbox.xmax() - bbox.xmin();
    FT dy = bbox.ymax() - bbox.ymin();
    FT dz = bbox.zmax() - bbox.zmin();
  
    Bare_point mid(bbox.xmin() + dx/2, bbox.ymin() + dy/2, bbox.zmin() + dz/2);
    FT dr = sqrt(CGAL::to_double(max_weight)) + .001;
  
    regular.insert(Weighted_point(
      Bare_point(bbox.xmax()+(dy+dz+dr)/gt.get_shrink(),mid.y(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(bbox.xmin()-(dy+dz+dr)/gt.get_shrink(),mid.y(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),bbox.ymax()+(dx+dz+dr)/gt.get_shrink(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),bbox.ymin()-(dx+dz+dr)/gt.get_shrink(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),mid.y(),bbox.zmax()+(dx+dy+dr)/gt.get_shrink()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),mid.y(),bbox.zmin()-(dx+dy+dr)/gt.get_shrink()),-1));

    // Set the bounding sphere for the Delaunay mesher
    _bounding_sphere = Sphere(mid, dr*dr+1);
  }
}

template <class MixedComplexTraits_3> 
typename Skin_surface_3<MixedComplexTraits_3>::TMC_Cell_handle
Skin_surface_3<MixedComplexTraits_3>::
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
  
//   Simplex prev, s;
//   if (start == Simplex()) {
//     Cell_handle ch = regular.locate(p);
//     if (regular.is_infinite(ch->vertex(0))) { s = ch->vertex(1); }
//     else { s = ch->vertex(0); }
//   } else {
//     s = start;
//   }
//   CGAL_assertion(s != Simplex());
//   // random walk, start with vh:
//   CGAL_assertion(regular.dimension() == 3);

//   // For storing a simplex
//   Cell_handle ch; int i1,i2;

//   // Traits class object:
//   typename Gt::Side_of_mixed_cell_3 
//     side_tester = gt.side_of_mixed_cell_3_object();
// //   std::cout << "[";
//  try_next_cell:

// //   std::cout << s.dimension();
//   switch (s.dimension()) {
//   case 0:
//     {
//       Vertex_handle vh = s;
//       std::vector<Vertex_handle> nbs;
//       nbs.reserve(64);
//       regular.incident_vertices(vh, std::back_inserter(nbs));
//       int nrNbs = nbs.size();
      
//       int index = rng.get_int(0,nrNbs);

//       for (int i=0; i<nrNbs; i++, index = (index+1)%nrNbs) {
// 	if (!regular.is_infinite(nbs[index])) {
// 	  if (prev != Simplex(nbs[index])) {
// 	    if (side_tester(vh->point(), nbs[index]->point(), p) == POSITIVE) {
// 	      prev = s;
// 	      regular.is_edge(vh, nbs[index], ch, i1, i2);
// 	      s = Edge(ch,i1,i2);
// 	      goto try_next_cell;
// 	    }
// 	    }
// 	}
//       }
//       break;
//     }
//   case 1:
//     {
//       Edge e = s;
//       Vertex_handle vh1 = e.first->vertex(e.second);
//       Vertex_handle vh2 = e.first->vertex(e.third);

//       Facet_circulator fcir;
//       fcir = regular.incident_facets(e);
//       int nrFacets = circulator_size(fcir);
      
//       // 2 additional neighbors for vertices
//       int index = rng.get_int(0,nrFacets+2);
//       if (index < nrFacets-1) for (int i=0; i<index; i++) fcir++;

//       for (int i=0; i<nrFacets+2; i++, index = (index+1)%(nrFacets+2)) {
// 	if (index < nrFacets) {
// 	  // Check incident facets:
// 	  if (!regular.is_infinite(*fcir)) {
// 	    if (prev != Simplex(fcir)) {
// 	      i1 = (*fcir).first->index(vh1);
// 	      i2 = (*fcir).first->index(vh2);
// 	      Vertex_handle vh3 = (*fcir).first->vertex(6-(*fcir).second-i1-i2);

// 	      if (side_tester(vh1->point(), vh2->point(), vh3->point(),
// 			      p) == POSITIVE) {
// 		prev = s;
// 		s = fcir;
// 		goto try_next_cell;
// 	      }
// 	      }
// 	  }
// 	  fcir++;
// 	} else {
// 	  // Check incident vertices:
// 	  if (index==nrFacets) {
// 	    if (prev != Simplex(vh1)) {
// 	      if (side_tester(vh1->point(), vh2->point(), p) == NEGATIVE) {
// 		prev = s;
// 		s = vh1;
// 		goto try_next_cell;
// 	      }
// 	      }
// 	  } else {
// 	    if (prev != Simplex(vh2)) {
// 	      if (side_tester(vh2->point(), vh1->point(), p) == NEGATIVE) {
// 		prev = s;
// 		s = vh2;
// 		goto try_next_cell;
// 	      }
// 	      }
// 	  }
// 	}
//       }
//       break;
//     }
//   case 2:
//     {
//       Facet f = s;
//       // 3x towards edge, 2x towards cell
//       int index = rng.get_int(0,5); 
//       for (int i=0; i<5; i++, index = (index+1)%5) {
// 	if (index > 2) {
// 	  // Check incident cells
// 	  ch = f.first;
// 	  i1 = f.second;
// 	  if (index == 3) {
// 	    ch = ch->neighbor(i1);
// 	    i1 = ch->index(f.first);
// 	  }
// 	  CGAL_assertion(!regular.has_vertex(f, ch->vertex(i1)));
// 	  if (!regular.is_infinite(ch->vertex(i1))) {
// 	    if (prev != Simplex(ch)) {
// 	      if (side_tester(ch->vertex((i1+1)&3)->point(),
// 			      ch->vertex((i1+2)&3)->point(),
// 			      ch->vertex((i1+3)&3)->point(),
// 			      ch->vertex(i1)->point(), p) == POSITIVE) {
// 		prev = s;
// 		s = ch;
// 		goto try_next_cell;
// 	      }
// 	      }
// 	  }
// 	} else {
// 	  // Check incident edges (index = 0,1,2)
// 	  i1 = (f.second+1)&3;
// 	  i2 = (f.second+2)&3;
// 	  int i3 = (f.second+3)&3;
// 	  if (index == 1) std::swap(i1,i3);
// 	  if (index == 2) std::swap(i2,i3);
// 	  Vertex_handle vh1 = f.first->vertex(i1);
// 	  Vertex_handle vh2 = f.first->vertex(i2);
// 	  Vertex_handle vh3 = f.first->vertex(i3);
	  
// 	  if (prev != Simplex(Edge(f.first,i1,i2))) {
// 	    if (side_tester(vh1->point(), vh2->point(), vh3->point(), 
// 			    p) == NEGATIVE) {
// 	      prev = s;
// 	      s = Edge(f.first,i1,i2);
// 	      goto try_next_cell;
// 	    }
// 	    }
// 	}
//       }
//       break;
//     }
//   case 3:
//     {
//       Cell_handle ch = s;
//       int index = rng.get_int(0,4); 
//       for (int i=0; i<4; i++, index = (index+1)&3) {
// 	if (prev != Simplex(Facet(ch, index))) {
// 	  if (side_tester(ch->vertex((index+1)&3)->point(),
// 			  ch->vertex((index+2)&3)->point(),
// 			  ch->vertex((index+3)&3)->point(),
// 			  ch->vertex(index)->point(), p) == NEGATIVE) {
// 	    prev = s;
// 	    s = Facet(ch, index);
// 	    goto try_next_cell;
// 	  }
// 	  }
//       }
//       break;
//     }
//   default:
//     {
//       CGAL_assertion(false);
//     }
//   }
// //   std::cout << "]";

//   return s;
// }

template <class MixedComplexTraits_3> 
template <class Polyhedron_3>
void
Skin_surface_3<MixedComplexTraits_3>::mesh_skin_surface_3(Polyhedron_3 &p) const {
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
  marching_tetrahedra_3(tmc_vertices_begin(), 
			tmc_vertices_end(), 
			tmc_cells_begin(), 
			tmc_cells_end(), 
			p, 
			marching_traits,
			marching_observer);
}

template <class MixedComplexTraits_3> 
template <class Polyhedron_3>
void
Skin_surface_3<MixedComplexTraits_3>::subdivide_skin_surface_mesh_3(Polyhedron_3 &p) const {
  std::cout << "Skin_Surface_3.subdivide_skin_surface_mesh_3(p)" << std::endl;

  typedef Skin_surface_refinement_policy_3<Self, Polyhedron_3> Policy;
  typedef Skin_surface_sqrt3<Self, Polyhedron_3, Policy>       Subdivider;

  Policy policy(*this);
  Subdivider subdivider(*this, p, policy);
  subdivider.subdivide();
}


  

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_3_H
