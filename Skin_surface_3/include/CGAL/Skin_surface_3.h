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

#include <CGAL/triangulate_mixed_complex_3.h>

// Needed for the (Delaunay) surface mesher
#include <CGAL/Skin_surface_mesher_oracle_3.h>
#include <CGAL/Triangulation_simplex_3.h>

// For point location in the mixed complex
#include <CGAL/Random.h>
#include <CGAL/Skin_surface_traits_3.h>

CGAL_BEGIN_NAMESPACE 

template <class MixedComplexTraits_3> 
class Skin_surface_3 {
  typedef MixedComplexTraits_3            Gt;
  typedef Skin_surface_3<Gt>              Self;
public:
  typedef MixedComplexTraits_3            Geometric_traits;
  typedef typename Gt::Weighted_point     Weighted_point;
  typedef typename Weighted_point::Weight RT;
  // NGHK:: added for the Delaunay mesher
  typedef typename Gt::Sphere_3           Sphere;
  typedef typename Weighted_point::Point  Bare_point;
  typedef typename Gt::Vector_3           Vector;
  
  typedef Regular_triangulation_3<Gt>     Regular;

  typedef Exact_predicates_inexact_constructions_kernel     Filtered_kernel;
  typedef Skin_surface_quadratic_surface_3<Filtered_kernel> 
                                                         Quadratic_surface;
public:
  typedef typename Regular::Vertex_handle                Vertex_handle;
  typedef typename Regular::Edge                         Edge;
  typedef typename Regular::Facet                        Facet;
  typedef typename Regular::Facet_circulator             Facet_circulator;
  typedef typename Regular::Cell_handle                  Cell_handle;
  typedef Triangulation_simplex_3<Regular>               Simplex;

  typedef typename Regular::Finite_vertices_iterator     Finite_vertices_iterator;
  typedef typename Regular::Finite_edges_iterator        Finite_edges_iterator;
  typedef typename Regular::Finite_facets_iterator       Finite_facets_iterator;
  typedef typename Regular::Finite_cells_iterator        Finite_cells_iterator;

  typedef Combinatorial_mixed_complex_triangulator_3<Regular>  CMCT;
  typedef typename CMCT::Vertex_handle                    CMCT_Vertex_handle;
  typedef typename CMCT::Vertex_iterator                  CMCT_Vertex_iterator;
  typedef typename CMCT::Cell                             CMCT_Cell;
  typedef typename CMCT::Cell_iterator                    CMCT_Cell_iterator;


  // NGHK: added for the (Delaunay) surface mesher, document
  typedef Exact_predicates_inexact_constructions_kernel  Mesher_Gt;
  typedef Skin_surface_mesher_oracle_3<Mesher_Gt,Self> Surface_mesher_traits_3;

private:
  
public:
  template < class WP_iterator >
  Skin_surface_3(WP_iterator begin, WP_iterator end, 
		 RT shrink_factor,
		 bool grow_balls = true,
		 Gt gt_ = Gt(),
		 bool _verbose = false
		 ) 
    : gt(gt_), verbose(_verbose) {

    gt.set_shrink(shrink_factor);
    CGAL_assertion(begin != end);

    if (grow_balls) {
      for (; begin != end; begin++) {
	regular.insert(Weighted_point(*begin, begin->weight()/gt.get_shrink()));
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
    

    mc_triangulator = new CMCT(regular, verbose);

  }

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
  CMCT_Cell locate_tet(const Bare_point &p, const Simplex &s) const {
    Skin_surface_traits_3<Exact_predicates_exact_constructions_kernel> 
      traits(gt.get_shrink());
    
    std::vector<CMCT_Cell> cells;

    switch (s.dimension()) {
    case 0:
      {
	Vertex_handle vh = s;
	CGAL_assertion(!regular.is_infinite(vh));
	mc_triangulator->construct_0_cell(vh, std::back_inserter(cells));
	break;
      }
    case 1:
      {
	Edge e = s;
	CGAL_assertion(!regular.is_infinite(e));
	mc_triangulator->construct_1_cell(e, std::back_inserter(cells));
	break;
      }
    case 2:
      {
	Facet f = s;
	CGAL_assertion(!regular.is_infinite(f));
	mc_triangulator->construct_2_cell(f, std::back_inserter(cells));
	break;
      }
    case 3:
      {
	Cell_handle ch = s;
	CGAL_assertion(!regular.is_infinite(ch));
	mc_triangulator->construct_3_cell(ch, std::back_inserter(cells));
	break;
      }
    default:
      {
	CGAL_assertion(false);
      }
    }


    bool found = false;

    for (typename std::vector<CMCT_Cell>::iterator it = cells.begin();
	 ((!found)&&(it != cells.end())); it++) {
      if (mc_triangulator->bounded_side(p, *it, traits) != ON_UNBOUNDED_SIDE) {
	return *it;
      }
    }

    return CMCT_Cell();
  }
  Simplex locate_mixed(const Bare_point &p, 
		       const Simplex &start = Simplex()) const;

  // exact computation of the sign on a vertex of the TMC
  Sign sign(const CMCT_Vertex_handle vh) const {
    typedef Exact_predicates_exact_constructions_kernel K;
    Skin_surface_traits_3<K> traits(gt.get_shrink());
    
    typename K::Point_3 p = mc_triangulator->location(vh, traits);

    return construct_surface(vh->first, K()).sign(p);
  }
  Sign sign(const Bare_point &p, const Simplex &start = Simplex()) const {
    return get_sign(locate_mixed(p,start), p);
  }

  // Trivial caching: check wether the surface is the same as the previous:
  mutable Skin_surface_quadratic_surface_3<
    Simple_cartesian<Interval_nt_advanced> > previous_sign_surface;
  mutable Simplex                            previous_sign_simplex;
  Sign get_sign(const Simplex &sim, const Bare_point &p) const {
    if (previous_sign_simplex != sim) {
      previous_sign_simplex = sim;
      previous_sign_surface = 
	construct_surface(sim, 
			  Simple_cartesian<Interval_nt_advanced>());
    }
    try
    {
      CGAL_PROFILER(std::string("NGHK: calls to    : ") + 
		    std::string(CGAL_PRETTY_FUNCTION));
      Protect_FPU_rounding<true> P;
      Sign result = previous_sign_surface.sign(p);
      if (! is_indeterminate(result))
        return result;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("NGHK: failures of : ") + 
		  std::string(CGAL_PRETTY_FUNCTION));
    Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
    return construct_surface
      (sim, 
       Exact_predicates_exact_constructions_kernel()).sign(p);
  }
  RT
  value(const Bare_point &p) const {
    Simplex sim = locate_mixed(p);
    return value(sim,p);
  }

  RT
  value(const Simplex &sim, const Bare_point &p) const {
    return 
      construct_surface(sim, typename Geometric_traits::Kernel()).value(p);
  }
  Vector
  normal(const Bare_point &p, const Simplex &start = Simplex()) const {
    return get_normal(locate_mixed(p,start), p);
  }
  Vector
  get_normal(const Simplex &mc, const Bare_point &p) const {
    return construct_surface(mc).gradient(p);
  }

  // Move the point in the direction of the gradient
  void to_surface(Bare_point &p,
		  const Simplex &start = Simplex()) const {
    Bare_point p1 = p;
    Simplex s1 = locate_mixed(p,start);
    Sign sign1 = get_sign(s1, p1);

    Vector n = get_normal(s1,p);
    if (sign1 == POSITIVE) n = -value(s1,p)*n;

    int k=2;
    Bare_point p2 = p+k*n;
    Simplex s2 = locate_mixed(p2, s1);
    while (get_sign(s2,p2) == sign1) {
      k++;
      p1 = p2;
      s1 = s2;
      p2 = p+k*n;
      s2 = locate_mixed(p2, s2);
    }
    intersect(p1,p2, s1,s2, p);
  }
  void intersect(const CMCT_Vertex_handle vh1,
		 const CMCT_Vertex_handle vh2,
		 Bare_point &p) const {
    Bare_point p1 = mc_triangulator->location(vh1, gt);
    Bare_point p2 = mc_triangulator->location(vh2, gt);
    Simplex s1 = vh1->first;
    Simplex s2 = vh2->first;
    intersect(p1,p2, s1,s2, p);
  }
  void intersect(const CMCT_Vertex_handle vh1,
		 const CMCT_Vertex_handle vh2,
		 const Simplex &s,
		 Bare_point &p) const {
    Bare_point p1 = mc_triangulator->location(vh1, gt);
    Bare_point p2 = mc_triangulator->location(vh2, gt);
    Simplex sp = s;
    intersect(p1,p2, sp,sp, p);
  }

  void intersect(Bare_point &p1, Bare_point &p2, 
		 Simplex &s1, Simplex &s2,
		 Bare_point &p) const {
    typedef typename Bare_point::R  Traits;
    typedef typename Traits::RT RT;
    Cartesian_converter<Traits, 
                        typename Geometric_traits::Bare_point::R> converter;

    RT sq_dist = squared_distance(p1,p2);
    // Use value to make the computation robust (endpoints near the surface)
    if (value(s1, p1) > value(s2, p2)) std::swap(p1, p2);
    Simplex sp = s1;

    while ((s1 != s2) && (sq_dist > 1e-8)) {
      p = midpoint(p1, p2);
      sp = locate_mixed(converter(p), sp);

      if (get_sign(sp, p) == NEGATIVE) { p1 = p; s1 = sp; }
      else { p2 = p; s2 = sp; }

      sq_dist *= .25;
    }
    while (sq_dist > 1e-8) {
      p = midpoint(p1, p2);
      if (get_sign(s1, p) == NEGATIVE) { p1 = p; }
      else { p2 = p; }
      sq_dist *= .25;
    }

    p = midpoint(p1, p2);
  }

  void intersect_with_transversal_segment
  (Bare_point &p,
   const Simplex &start = Simplex()) const 
  {

    typedef typename Geometric_traits::Kernel::Plane_3 Plane;
    typedef typename Geometric_traits::Kernel::Line_3  Line;

    Simplex sim = locate_mixed(p, start);
    CMCT_Cell tet = locate_tet(p, sim);
    
    // get transversal segment:
    Bare_point p1, p2;

    // Compute signs on vertices and sort them:
    int nIn = 0;
    int sortedV[4];
    for (int i=0; i<4; i++) {
      if (sign(tet.vertex(i))==POSITIVE) {
        sortedV[nIn] = i; nIn++;
      } else {
        sortedV[3-i+nIn] = i;
      }
    }

    Object obj;
    if (nIn==1) {
      p1 = mc_triangulator->location(tet.vertex(sortedV[0]), gt);
      obj = CGAL::intersection(
        Plane
	(mc_triangulator->location(tet.vertex(sortedV[1]), gt),
	 mc_triangulator->location(tet.vertex(sortedV[2]), gt),
	 mc_triangulator->location(tet.vertex(sortedV[3]), gt)),
        Line(p1, p));
      if ( !assign(p2, obj) ) {
        CGAL_assertion_msg(false,"intersection: no intersection.");
      }
    } else if (nIn==2) {
      obj = CGAL::intersection(
        Plane
	(mc_triangulator->location(tet.vertex(sortedV[2]), gt),
	 mc_triangulator->location(tet.vertex(sortedV[3]), gt),
          p),
        Line(
          mc_triangulator->location(tet.vertex(sortedV[0]), gt),
	  mc_triangulator->location(tet.vertex(sortedV[1]), gt)));
      if ( !assign(p1, obj) ) {
        CGAL_assertion_msg(false,"intersection: no intersection.");
      }
      obj = CGAL::intersection(
       Plane
	(mc_triangulator->location(tet.vertex(sortedV[0]), gt),
	 mc_triangulator->location(tet.vertex(sortedV[1]), gt),
          p),
        Line(
          mc_triangulator->location(tet.vertex(sortedV[2]), gt),
	  mc_triangulator->location(tet.vertex(sortedV[3]), gt)));
      if ( !assign(p2, obj) ) {
        CGAL_assertion_msg(false,"intersection: no intersection.");
      }
    } else if (nIn==3) {
      p2 = mc_triangulator->location(tet.vertex(sortedV[3]), gt);
      obj = CGAL::intersection(
        Plane(
          mc_triangulator->location(tet.vertex(sortedV[0]), gt),
          mc_triangulator->location(tet.vertex(sortedV[1]), gt),
          mc_triangulator->location(tet.vertex(sortedV[2]), gt)),
        Line(p2, p));
      if ( !assign(p1, obj) ) {
        CGAL_assertion_msg(false,"intersection: no intersection.");
      }
    } else {
      CGAL_assertion(false);
    }

    // Find the intersection:
    intersect(p1, p2, sim, sim, p);
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
    typedef typename Traits::RT                           RT;
    typedef CGAL::Weighted_point<Point,RT>                Weighted_point;

    Converter conv;

    switch (sim.dimension()) {
    case 0:
      {
	Vertex_handle vh = sim;
	return Quadratic_surface(conv(vh->point()), gt.get_shrink());
	break;
      }
    case 1:
      {
	Edge e = sim;
	Weighted_point p0 = conv(e.first->vertex(e.second)->point());
	Weighted_point p1 = conv(e.first->vertex(e.third)->point());
	return Quadratic_surface(p0, p1, gt.get_shrink());
	break;
      }
    case 2:
      {
	Facet f = sim;
	Weighted_point p0 = conv(f.first->vertex((f.second+1)&3)->point());
	Weighted_point p1 = conv(f.first->vertex((f.second+2)&3)->point());
	Weighted_point p2 = conv(f.first->vertex((f.second+3)&3)->point());
	return Quadratic_surface(p0,p1,p2, gt.get_shrink());
	break;
      }
    case 3:
      {
	Cell_handle ch = sim;
	Weighted_point p0 = conv(ch->vertex(0)->point());
	Weighted_point p1 = conv(ch->vertex(1)->point());
	Weighted_point p2 = conv(ch->vertex(2)->point());
	Weighted_point p3 = conv(ch->vertex(3)->point());
	return Quadratic_surface(p0,p1,p2,p3, gt.get_shrink());
	break;
      }
    }
    CGAL_assertion(false);
    return Quadratic_surface();
  }

  // Access to the implicit triangulated mixed complex:
  CMCT_Vertex_iterator cmct_vertices_begin() const
  { return mc_triangulator->vertices_begin(); }
  CMCT_Vertex_iterator cmct_vertices_end() const
  { return mc_triangulator->vertices_end(); }
  CMCT_Cell_iterator cmct_cells_begin() const
  { return mc_triangulator->cells_begin(); }
  CMCT_Cell_iterator cmct_cells_end() const
  { return mc_triangulator->cells_end(); }
  

  // NGHK: added for the (Delaunay) surface mesher, document
  Sphere bounding_sphere() const {
    return _bounding_sphere;
  }
  RT squared_error_bound() const {
    return .01;
  }

  typename Mesher_Gt::RT 
  get_density(const typename Mesher_Gt::Point_3 &p) const {
    // NGHK: Make adaptive
    return 1;
  }
  const Regular &get_regular_triangulation() const {
    return regular;
  }
  RT get_shrink_factor() const {
    return gt.get_shrink();
  }

private:
  void construct_bounding_box(Regular &regular);

  Regular regular;
  Gt gt;
  bool verbose;
  Sphere _bounding_sphere;
  mutable Random rng;

  // We want to construct this object later (the pointer):
  CMCT *mc_triangulator;
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
    RT dr = sqrt(CGAL::to_double(max_weight)) + .001;
  
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
typename Skin_surface_3<MixedComplexTraits_3>::Simplex 
Skin_surface_3<MixedComplexTraits_3>::
locate_mixed(const Bare_point &p, const Simplex &start) const {
  Simplex /*prev,*/ s;
  if (start == Simplex()) {
    Cell_handle ch = regular.locate(p);
    if (regular.is_infinite(ch->vertex(0))) { s = ch->vertex(1); }
    else { s = ch->vertex(0); }
  } else {
    s = start;
  }
  CGAL_assertion(s != Simplex());
  // random walk, start with vh:
  CGAL_assertion(regular.dimension() == 3);

  // For storing a simplex
  Cell_handle ch; int i1,i2;

  // Traits class object:
  typename Gt::Side_of_mixed_cell_3 
    side_tester = gt.side_of_mixed_cell_3_object();
//   std::cout << "[";
 try_next_cell:

//   std::cout << s.dimension();
  switch (s.dimension()) {
  case 0:
    {
      Vertex_handle vh = s;
      std::vector<Vertex_handle> nbs;
      nbs.reserve(64);
      regular.incident_vertices(vh, std::back_inserter(nbs));
      int nrNbs = nbs.size();
      
      int index = rng.get_int(0,nrNbs);

      for (int i=0; i<nrNbs; i++, index = (index+1)%nrNbs) {
	if (!regular.is_infinite(nbs[index])) {
	  //if (prev != Simplex(nbs[index])) {
	    if (side_tester(vh->point(), nbs[index]->point(), p) == POSITIVE) {
	      //prev = s;
	      bool b = regular.is_edge(vh, nbs[index], ch, i1, i2);
	      CGAL_assertion(b);
	      s = Edge(ch,i1,i2);
	      goto try_next_cell;
	    }
	    //}
	}
      }
      break;
    }
  case 1:
    {
      Edge e = s;
      Vertex_handle vh1 = e.first->vertex(e.second);
      Vertex_handle vh2 = e.first->vertex(e.third);

      Facet_circulator fcir;
      fcir = regular.incident_facets(e);
      int nrFacets = circulator_size(fcir);
      
      // 2 additional neighbors for vertices
      int index = rng.get_int(0,nrFacets+2);
      if (index < nrFacets-1) for (int i=0; i<index; i++) fcir++;

      for (int i=0; i<nrFacets+2; i++, index = (index+1)%(nrFacets+2)) {
	if (index < nrFacets) {
	  // Check incident facets:
	  if (!regular.is_infinite(*fcir)) {
	    //if (prev != Simplex(fcir)) {
	      i1 = (*fcir).first->index(vh1);
	      i2 = (*fcir).first->index(vh2);
	      Vertex_handle vh3 = (*fcir).first->vertex(6-(*fcir).second-i1-i2);

	      if (side_tester(vh1->point(), vh2->point(), vh3->point(),
			      p) == POSITIVE) {
		//prev = s;
		s = fcir;
		goto try_next_cell;
	      }
	      //}
	  }
	  fcir++;
	} else {
	  // Check incident vertices:
	  if (index==nrFacets) {
	    //if (prev != Simplex(vh1)) {
	      if (side_tester(vh1->point(), vh2->point(), p) == NEGATIVE) {
		//prev = s;
		s = vh1;
		goto try_next_cell;
	      }
	      //}
	  } else {
	    //if (prev != Simplex(vh2)) {
	      if (side_tester(vh2->point(), vh1->point(), p) == NEGATIVE) {
		//prev = s;
		s = vh2;
		goto try_next_cell;
	      }
	      //}
	  }
	}
      }
      break;
    }
  case 2:
    {
      Facet f = s;
      // 3x towards edge, 2x towards cell
      int index = rng.get_int(0,5); 
      for (int i=0; i<5; i++, index = (index+1)%5) {
	if (index > 2) {
	  // Check incident cells
	  ch = f.first;
	  i1 = f.second;
	  if (index == 3) {
	    ch = ch->neighbor(i1);
	    i1 = ch->index(f.first);
	  }
	  CGAL_assertion(!regular.has_vertex(f, ch->vertex(i1)));
	  if (!regular.is_infinite(ch->vertex(i1))) {
	    //if (prev != Simplex(ch)) {
	      if (side_tester(ch->vertex((i1+1)&3)->point(),
			      ch->vertex((i1+2)&3)->point(),
			      ch->vertex((i1+3)&3)->point(),
			      ch->vertex(i1)->point(), p) == POSITIVE) {
		//prev = s;
		s = ch;
		goto try_next_cell;
	      }
	      //}
	  }
	} else {
	  // Check incident edges (index = 0,1,2)
	  i1 = (f.second+1)&3;
	  i2 = (f.second+2)&3;
	  int i3 = (f.second+3)&3;
	  if (index == 1) std::swap(i1,i3);
	  if (index == 2) std::swap(i2,i3);
	  Vertex_handle vh1 = f.first->vertex(i1);
	  Vertex_handle vh2 = f.first->vertex(i2);
	  Vertex_handle vh3 = f.first->vertex(i3);
	  
	  //if (prev != Simplex(Edge(f.first,i1,i2))) {
	    if (side_tester(vh1->point(), vh2->point(), vh3->point(), 
			    p) == NEGATIVE) {
	      //prev = s;
	      s = Edge(f.first,i1,i2);
	      goto try_next_cell;
	    }
	    //}
	}
      }
      break;
    }
  case 3:
    {
      Cell_handle ch = s;
      int index = rng.get_int(0,4); 
      for (int i=0; i<4; i++, index = (index+1)&3) {
	//if (prev != Simplex(Facet(ch, index))) {
	  if (side_tester(ch->vertex((index+1)&3)->point(),
			  ch->vertex((index+2)&3)->point(),
			  ch->vertex((index+3)&3)->point(),
			  ch->vertex(index)->point(), p) == NEGATIVE) {
	    //prev = s;
	    s = Facet(ch, index);
	    goto try_next_cell;
	  }
	  //}
      }
      break;
    }
  default:
    {
      CGAL_assertion(false);
    }
  }
//   std::cout << "]";

  return s;
}

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_3_H
