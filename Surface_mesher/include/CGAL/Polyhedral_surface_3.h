// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POLYHEDRAL_SURFACE_3_H
#define CGAL_POLYHEDRAL_SURFACE_3_H

#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
#  define CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
#endif

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Data_structure_using_octree_3.h>
#include <CGAL/Surface_mesher/Polyhedral_oracle.h>
#include <CGAL/Surface_mesher/Has_edges.h>
#include <iostream>
#include <vector>
#include <algorithm> // random_shuffle

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
#include <boost/format.hpp>
#include <CGAL/Timer.h>
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION

namespace CGAL {

  namespace Surface_mesher {
    template <
      class Surface,
      class Point_creator,
      class Visitor,
      class Tag
      >
    class Polyhedral_oracle;
  } // end namespace Surface_mesher

template <class GT, class Has_edges_tag_ = Surface_mesher::Has_no_edges >
class Polyhedral_surface_3 : public Has_edges_tag_
{
public:
  typedef GT Geom_traits;
  typedef Has_edges_tag_ Has_edges_tag;

  class Normalized_geom_traits : public Geom_traits 
  {
  public:
    typedef typename 
    Kernel_traits<typename Geom_traits::Point_3>::Kernel::Point_3 Point_3;
  };

  typedef Data_structure_using_octree_3<Normalized_geom_traits> Subfacets_octree;
  typedef Data_structure_using_octree_3<Normalized_geom_traits> Subsegments_octree;
  typedef typename GT::Point_3 Point_3;
  typedef typename GT::Vector_3 Vector_3;
  typedef typename GT::FT FT;

  typedef Polyhedral_surface_3<GT, Has_edges_tag> Self;

  template <
    class Surface,
    class Point_creator,
    class Visitor,
    class Tag
    >
  friend class Surface_mesher::Polyhedral_oracle;

  typedef Surface_mesher::Polyhedral_oracle<Self> Surface_mesher_traits_3;

  typedef typename Subfacets_octree::Bbox Bbox;

  Polyhedral_surface_3(std::istream& input_file,
		       const FT cosine_bound = FT(1)/FT(2))
    : subfacets_octree(), subsegments_octree(false, true, false),
      input_points()
  {
    GT gt = GT();
    typename GT::Construct_orthogonal_vector_3 orthogonal_vector = 
      gt.construct_orthogonal_vector_3_object();
    typename GT::Compute_squared_length_3 squared_length =
      gt.compute_squared_length_3_object();
    typename GT::Compute_scalar_product_3 scalar_product = 
      gt.compute_scalar_product_3_object();

    typedef CGAL::Polyhedron_3<GT> Polyhedron_3;

    class Facet_ortho_vector {
      GT gt;
    public:
      Facet_ortho_vector(GT gt) : gt(gt) {}

      typename GT::Vector_3 
      operator()(const typename Polyhedron_3::Facet& f) {
	typename Polyhedron_3::Halfedge_around_facet_const_circulator 
	  edges_circ = f.facet_begin();
	const Point_3& p1 = edges_circ++->vertex()->point();
	const Point_3& p2 = edges_circ++->vertex()->point();
	const Point_3& p3 = edges_circ++->vertex()->point();
	return gt.construct_orthogonal_vector_3_object()(p1, p2, p3);
      }
    };
    Facet_ortho_vector facet_ortho_vector(gt);

    const FT cosine_squared_bound = cosine_bound * cosine_bound;

    Polyhedron_3 polyhedron;

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    CGAL::Timer timer;
    std::cerr << "Creating polyhedron... ";
    timer.start();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    CGAL::scan_OFF(input_file, polyhedron, true);
    CGAL_assertion(input_file);
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    {
      std::ofstream dump("input_dump.off");
      dump << polyhedron;
    }
    timer.stop();
    std::cerr << 
      ::boost::format("done (%1%s)\n"
		      "  number of vertices: %2%\n"
		      "  number of edges:    %3%\n"
		      "  number of facets:   %4%\n")
      % timer.time()
      % polyhedron.size_of_vertices()
      % ( polyhedron.size_of_halfedges() / 2 )
      % polyhedron.size_of_facets();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION

    for(typename Polyhedron_3::Vertex_const_iterator vit = 
          polyhedron.vertices_begin();
        vit != polyhedron.vertices_end();
        ++vit)
    {
      subfacets_octree.add_constrained_vertex(vit->point());
      input_points.insert(vit->point());
    }

    typename Polyhedron_3::size_type facet_index = 0;
    for(typename Polyhedron_3::Facet_const_iterator fit = 
          polyhedron.facets_begin();
        fit != polyhedron.facets_end();
        ++fit, ++facet_index)
    {
      CGAL_assertion(fit->is_triangle());
      typename Polyhedron_3::Halfedge_around_facet_const_circulator 
        edges_circ = fit->facet_begin();

      const Point_3& p1 = edges_circ++->vertex()->point();
      const Point_3& p2 = edges_circ++->vertex()->point();
      const Point_3& p3 = edges_circ++->vertex()->point();

      subfacets_octree.add_constrained_facet(p1, p2, p3);
#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	std::cerr << ::boost::format("new facet in subfacets_octree: #%4% (%1%, %2%, %3%)\n")
	  % p1 %  p2 % p3 % facet_index;
#endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
    }

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    std::cerr << "Creating subfacets_octree... ";
    timer.reset();
    timer.start();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    subfacets_octree.create_data_structure();
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    timer.stop();
    std::cerr <<
      ::boost::format("done (%1%s)\n"
		      "  number of facets:      %3%\n"
		      "  number of constraints: %4% (in subfacets_octree)\n")
      % timer.time()
      % subfacets_octree.number_of_vertices()
      % subfacets_octree.number_of_facets()
      % subfacets_octree.number_of_constraints();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    if( this->has_edges() ) {

      typedef std::map<Point_3, int> Edges_vertex_counter;
      Edges_vertex_counter edges_vertex_counter;

      for(typename Polyhedron_3::Edge_const_iterator eit = 
	    polyhedron.edges_begin();
	  eit != polyhedron.edges_end();
	  ++eit)
      {
	typename Polyhedron_3::Halfedge_const_handle opposite = eit->opposite();

	bool insert_that_edge = false;
	if( eit->is_border_edge() ) 
	  insert_that_edge = true;
	else	
	{
	  //       CGAL_assertion(eit->is_triangle());
	  //       CGAL_assertion(opposite->is_triangle());
	  if(eit->facet()->facet_degree() != 3)
	    std::cerr << "warning: degree=" << eit->facet_degree() << "\n";
	  if(opposite->facet()->facet_degree() != 3)
	    std::cerr << "warning: degree(opposite)=" << opposite->facet_degree() << "\n";

	  const Vector_3 v1 = facet_ortho_vector(*eit->facet());
	  const Vector_3 v2 = facet_ortho_vector(*eit->opposite()->facet());

	  const FT product = scalar_product(v1, v2);

	  insert_that_edge = product < FT(0) ||
	    product * product < 
	    cosine_squared_bound * squared_length(v1) * squared_length(v2);
	}
	if(insert_that_edge)
	{ 
	  const Point_3 pa = eit->vertex()->point();
	  const Point_3 pb = opposite->vertex()->point();
	  ++edges_vertex_counter[pa];
	  ++edges_vertex_counter[pb];
	  subsegments_octree.add_constrained_edge(pa, pb);
#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  std::cerr << 
	    ::boost::format("new edge in subsegments_octree: (%1%, %2%)")
	    % pa %  pb;
	  if(eit->is_border_edge())
	    std:: cerr << " (on border)\n";
	  else
	  {
	    struct Triangle {
	      typename GT::Triangle_3 operator()(typename Polyhedron_3::Facet facet) {
		CGAL_assertion(facet.is_triangle());
		typename Polyhedron_3::Halfedge_around_facet_const_circulator 
		  edges_circ = facet.facet_begin();
		const Point_3& p1 = edges_circ++->vertex()->point();
		const Point_3& p2 = edges_circ++->vertex()->point();
		const Point_3& p3 = edges_circ++->vertex()->point();
		return typename GT::Triangle_3(p1, p2, p3);
	      }
	    };
	    std:: cerr << " (not on border) ";
	    std::cerr << 
	      ::boost::format("facets=(%3%, %4%) normals=(%1%, %2%)\n")
	      % facet_ortho_vector(*eit->facet())
	      % facet_ortho_vector(*eit->opposite()->facet())
	      % Triangle()(*(eit->facet()))
	      % Triangle()(*(eit->opposite()->facet()));
	  }
#endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	}
      }

      for(typename Edges_vertex_counter::const_iterator it = 
	    edges_vertex_counter.begin();
	  it != edges_vertex_counter.end();
	  ++it)
      {
	input_points.erase(it->first);
	if(it->second != 2)
	{
#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  std::cerr << ::boost::format("corner point: (%1%)\n")
	    % it->first;
#endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  corner_points.push_back(it->first);
	}
	else
	{
#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  std::cerr << ::boost::format("edge point: (%1%)\n")
	    % it->first;
#endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  edges_points.push_back(it->first);
	}
      }
      if(!corner_points.empty() && edges_points.empty())
      {
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
	std::cerr << "Incorrect input data. "
		  << "Swap corner vertices and edges vertices...\n";
#endif
	std::swap(corner_points, edges_points);
      }

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      std::cerr << "Shuffle edges vertices... ";
#endif
      std::random_shuffle(edges_points.begin(), edges_points.end());
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      std::cerr << "done\n";
#endif

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      std::cerr <<
	::boost::format("number of corner vertices: %1%\n"
			"number of edges vertices:  %2%\n")
	% corner_points.size()
	% edges_points.size();

      std::cerr << "Creating subsegments_octree... ";
      timer.reset();
      timer.start();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION

      subsegments_octree.create_data_structure();
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      timer.stop();

    std::cerr <<
      ::boost::format("done (%1%s)\n"
		      "  number of edges:       %3%\n"
		      "  number of constraints: %4% (in subsegments_octree)\n")
      % timer.time()
      % subsegments_octree.number_of_vertices()
      % subsegments_octree.number_of_edges()
      % subsegments_octree.number_of_constraints();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    } // end "if(this->has_edges())"
    
//     subfacets_octree.input(input_file,
//                            std::back_inserter(input_points));
    bounding_box = subfacets_octree.bbox();
    bounding_box = bounding_box + subsegments_octree.bbox();
    bounding_box_sq_radius = bounding_box.xmax()-bounding_box.xmin();
    bounding_box_sq_radius =
      CGAL_NTS max BOOST_PREVENT_MACRO_SUBSTITUTION 
      (bounding_box_sq_radius,
       FT(bounding_box.ymax()-bounding_box.ymin()));
    bounding_box_sq_radius =
      CGAL_NTS max BOOST_PREVENT_MACRO_SUBSTITUTION
      (bounding_box_sq_radius,
       FT(bounding_box.zmax()-bounding_box.zmin()));
    bounding_box_sq_radius /= 2;
    bounding_box_sq_radius *= bounding_box_sq_radius;
  } // end of Polyhedral_surface_3 constructor

  const Bbox& bbox() const
  {
    return bounding_box;
  }

  const FT& bounding_sphere_squared_radius() const
  {
    return bounding_box_sq_radius;
  }

public:
  Subfacets_octree subfacets_octree;
  Subsegments_octree subsegments_octree;
  std::set<Point_3> input_points;
  std::vector<Point_3> corner_points;
  std::vector<Point_3> edges_points;

  Bbox bounding_box;
  FT bounding_box_sq_radius;
};

} // end namespace CGAL

#endif // CGAL_POLYHEDRAL_SURFACE_3_H
