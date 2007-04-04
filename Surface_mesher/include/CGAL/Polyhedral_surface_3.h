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

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Data_structure_using_octree_3.h>
#include <CGAL/Surface_mesher/Polyhedral_oracle.h>
#include <iostream>
#include <vector>

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
#include <boost/format.hpp>
#include <CGAL/Timer.h>
#endif

namespace CGAL {

template <class GT>
class Polyhedral_surface_3
{
public:
  typedef GT Geom_traits;

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

  typedef Polyhedral_surface_3<GT> Self;

  typedef Surface_mesher::Polyhedral_oracle<Self> Surface_mesher_traits_3;

  typedef typename Subfacets_octree::Bbox Bbox;

  Polyhedral_surface_3(std::istream& input_file,
		       const FT cosine_bound = FT(1)/FT(2))
    : subfacets_octree(), subsegments_octree(false, true, false),
      input_points()
  {
    const FT cosine_squared_bound = cosine_bound * cosine_bound;

    typedef CGAL::Polyhedron_3<GT> Polyhedron_3;
    Polyhedron_3 polyhedron;

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    CGAL::Timer timer;
    std::cerr << "Creating polyhedron... ";
    timer.start();
#endif
    CGAL::scan_OFF(input_file, polyhedron, true);
    CGAL_assertion(input_file);
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
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
#endif

    input_points.reserve(polyhedron.size_of_vertices());

    for(typename Polyhedron_3::Vertex_const_iterator vit = 
          polyhedron.vertices_begin();
        vit != polyhedron.vertices_end();
        ++vit)
    {
      subfacets_octree.add_constrained_vertex(vit->point());
      input_points.push_back(vit->point());
    }

    for(typename Polyhedron_3::Facet_const_iterator fit = 
          polyhedron.facets_begin();
        fit != polyhedron.facets_end();
        ++fit)
    {
      CGAL_assertion(fit->is_triangle());
      typename Polyhedron_3::Halfedge_around_facet_const_circulator 
        edges_circ = fit->facet_begin();

      const Point_3& p1 = edges_circ++->vertex()->point();
      const Point_3& p2 = edges_circ++->vertex()->point();
      const Point_3& p3 = edges_circ++->vertex()->point();

      subfacets_octree.add_constrained_facet(p1, p2, p3);
    }

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    std::cerr << "Creating subfacets_octree... ";
    timer.reset();
    timer.start();
#endif
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
#endif

    for(typename Polyhedron_3::Edge_const_iterator eit = 
          polyhedron.edges_begin();
        eit != polyhedron.edges_end();
        ++eit)
    {
      typename Polyhedron_3::Halfedge_const_handle opposite = eit->opposite();

      if( eit->is_border_edge() ) 
      {
	  subsegments_octree.add_constrained_edge(eit->vertex()->point(),
						  opposite->vertex()->point());
      }
      else	
      {
	//       CGAL_assertion(eit->is_triangle());
	//       CGAL_assertion(opposite->is_triangle());
	if(eit->facet()->facet_degree() != 3) std::cerr << "degree=" << eit->facet_degree() << "\n";
	if(opposite->facet()->facet_degree() != 3) std::cerr << "degree=" << opposite->facet_degree() << "\n";

	//       typename Polyhedron_3::Halfedge_around_facet_const_circulator 
	//         edges_circ = eit->facet_begin();

	//       const Point_3& p1 = edges_circ++->vertex()->point();
	//       const Point_3& p2 = edges_circ++->vertex()->point();
	//       const Point_3& p3 = edges_circ++->vertex()->point();

	//       edges_circ = opposite->facet_begin();
	//       const Point_3& p4 = edges_circ++->vertex()->point();
	//       const Point_3& p5 = edges_circ++->vertex()->point();
	//       const Point_3& p6 = edges_circ++->vertex()->point();

	typename GT::Construct_orthogonal_vector_3 orthogonal_vector = 
	  GT().construct_orthogonal_vector_3_object();

	typename GT::Compute_squared_length_3 squared_length =
	  GT().compute_squared_length_3_object();

	typename GT::Compute_scalar_product_3 scalar_product = 
	  GT().compute_scalar_product_3_object();

	const Vector_3 v1 = orthogonal_vector(eit->facet()->plane());
	const Vector_3 v2 = orthogonal_vector(opposite->facet()->plane());
	// (p4, p6, p5) in that order, because 'opposite' is in opposite
	// orientation.

	const FT product = scalar_product(v1, v2);

	if(product < FT(0) ||
	   product * product < 
	   cosine_squared_bound * squared_length(v1) * squared_length(v2))
	{
	  subsegments_octree.add_constrained_edge(eit->vertex()->point(),
						  opposite->vertex()->point());
	}
      }
    }

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    std::cerr << "Creating subsegments_octree... ";
    timer.reset();
    timer.start();
#endif
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
#endif
    
//     subfacets_octree.input(input_file,
//                            std::back_inserter(input_points));
  }

  Bbox bbox() const
  {
    return subfacets_octree.bbox();
  }

public:
  Subfacets_octree subfacets_octree;
  Subsegments_octree subsegments_octree;
  std::vector<Point_3> input_points;
};

} // end namespace CGAL

#endif // CGAL_POLYHEDRAL_SURFACE_3_H
