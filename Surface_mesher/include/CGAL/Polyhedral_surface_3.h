// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POLYHEDRAL_SURFACE_3_H
#define CGAL_POLYHEDRAL_SURFACE_3_H

#include <CGAL/config.h>
#ifdef CGAL_DONT_SUBMIT

#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
#  define CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
#endif

#ifndef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
#  define CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE 1
#endif

//#include <CGAL/Polyhedron_3.h>
#include <CGAL/enriched_polyhedron.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Inverse_index.h>
#include <CGAL/circulator.h>
#include <CGAL/iterator.h>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <CGAL/make_surface_mesh.h>
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
#  include <CGAL/Data_structure_using_octree_3.h>
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
#  include <CGAL/Surface_mesher/Intersection_data_structure_3.h>
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
#  include <pinpolyhedron.h>
#endif

#include <CGAL/Surface_mesher/Polyhedral_oracle.h>
#include <CGAL/Surface_mesher/Has_edges.h>
#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <algorithm> // random_shuffle, distance, copy, set_intersection

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
      class Tag,
      bool
      >
    class Polyhedral_oracle;
    
  } // end namespace Surface_mesher
  
  template <class Polyhedron>
  struct Polyhedral_surface_3_type_helper {
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef Halfedge_handle Edge_handle;

    typedef typename Polyhedron::Traits::Point_3 Point_3;

    struct Vertex_node {
      typedef std::set<int> Incident_edges; // set of edge id
      Incident_edges incident_edges;
    };
    
    struct Edge_node {
      typedef std::vector<Edge_handle> Sub_edges;
      typedef std::set<int> Incident_facets; // set of facet id
      Sub_edges sub_edges;
      Incident_facets incident_facets;
    };
    
    struct Facet_node {
      typedef std::vector<Facet_handle> Sub_facets;
      Sub_facets sub_facets;
    };

    struct Incidence_graph {
      typedef std::vector<Vertex_node> Vertices_nodes;
      typedef std::vector<Edge_node> Edges_nodes;
      typedef std::vector<Facet_node> Facet_nodes;
      typedef typename Vertices_nodes::iterator Vertex_iterator;
      typedef typename Edges_nodes::iterator Edge_iterator;
      typedef typename Facet_nodes::iterator Facet_iterator;
      Vertices_nodes vertices;
      Edges_nodes edges;
      Facet_nodes facets;

      void clear()
      {
        vertices.clear();
        edges.clear();
        facets.clear();
      }

      int vertex_id(Vertex_iterator vit)
      {
        return std::distance(vertices.begin(), vit);
      }

      int edge_id(Edge_iterator eit)
      {
        return std::distance(edges.begin(), eit);
      }

      int facet_id(Facet_iterator fit)
      {
        return std::distance(facets.begin(), fit);
      }

      void gl_draw_edge(const int edge_index)
      {
        ::glBegin(GL_LINES);
        typedef typename Edge_node::Sub_edges Sub_edges;
        for(typename Sub_edges::iterator 
              eit = edges[edge_index].sub_edges.begin(),
              end = edges[edge_index].sub_edges.end();
            eit != end;
            ++eit)
        {
          const Point_3& p1 = (*eit)->opposite()->vertex()->point();
          const Point_3& p2 = (*eit)->vertex()->point();
          ::glVertex3f(p1[0],p1[1],p1[2]);
          ::glVertex3f(p2[0],p2[1],p2[2]);       
        }
        ::glEnd();
      } // end gl_draw_edge(...)

      void gl_draw_facet(const int facet_index,
                         bool smooth_shading,
                         bool use_normals,
                         bool inverse_normals = false)
      {
        // draw triangles
        ::glBegin(GL_TRIANGLES);
        typedef typename Facet_node::Sub_facets Sub_facets;
        for(typename Sub_facets::iterator 
              fit = facets[facet_index].sub_facets.begin(),
              end_fit = facets[facet_index].sub_facets.end();
            fit != end_fit;
            ++fit)
        {
          // one normal per face
          if(use_normals && !smooth_shading)
          {
            typedef typename Polyhedron::Facet Facet;
            const typename Facet::Normal_3& normal = (*fit)->normal();
            if(inverse_normals)
              ::glNormal3f(-normal[0],-normal[1],-normal[2]);
            else
              ::glNormal3f(normal[0],normal[1],normal[2]);
          }
          typename Polyhedron::Halfedge_around_facet_circulator 
            he = (*fit)->facet_begin(),
            end_he = he;
          do
          {
            // one normal per vertex
            if(use_normals && smooth_shading)
            {
              const typename Polyhedron::Facet::Normal_3& normal = he->vertex()->normal();
              if(inverse_normals)
                ::glNormal3f(-normal[0],-normal[1],-normal[2]);
              else      
                ::glNormal3f(normal[0],normal[1],normal[2]);
            }
            // polygon assembly is performed per vertex
            const Point_3& point  = he->vertex()->point();
            ::glVertex3d(point[0],point[1],point[2]);
          }
          while(++he != end_he);
        }
        ::glEnd(); // end polygon assembly
      } // end gl_draw_facet(...)

    }; // end Incidence_graph
  }; // end Polyhedral_surface_3_type_helper

template <class GT,
          class Has_edges_tag_ = Surface_mesher::Has_no_edges,
          class Polyhedron_3 = Enriched_polyhedron<GT,
                                                   Enriched_items> >
class Polyhedral_surface_3 : 
  public Has_edges_tag_,
  public Polyhedron_3
{
public:
  typedef GT Geom_traits;
  typedef Has_edges_tag_ Has_edges_tag;
  typedef Polyhedron_3 Polyhedron;

  typedef Polyhedral_surface_3<Geom_traits,
                               Has_edges_tag,
                               Polyhedron> Self;

  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef Halfedge_handle Edge_handle;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator
                                  Halfedge_around_vertex_circulator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator
                                  Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;

  struct Compare_vertex_iterators {
    bool operator()(const typename Polyhedron_3::Vertex_iterator& va,
                    const typename Polyhedron_3::Vertex_iterator& vb) const
    {
      return (&*va)<(&*vb);
    }
  };

  using Polyhedron::halfedges_begin;
  using Polyhedron::halfedges_end;
  using Polyhedron::edges_begin;
  using Polyhedron::edges_end;
  using Polyhedron::facets_begin;
  using Polyhedron::facets_end;
  

  typedef Polyhedral_surface_3_type_helper<Polyhedron> Type_helper;
  typedef typename Type_helper::Vertex_node Graph_vertex_node;
  typedef typename Type_helper::Edge_node Graph_edge_node;
  typedef typename Type_helper::Facet_node Graph_facet_node;
  typedef typename Type_helper::Incidence_graph Incidence_graph;

  class Normalized_geom_traits : public Geom_traits 
  {
  public:
    typedef typename 
    Kernel_traits<typename Geom_traits::Point_3>::Kernel::Point_3 Point_3;
  };

  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename Geom_traits::Triangle_3 Triangle_3;
  typedef typename Geom_traits::Vector_3 Vector_3;
  typedef Bbox_3 Bbox;


#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
  typedef Data_structure_using_octree_3<Normalized_geom_traits> Subfacets_tree;
  typedef Data_structure_using_octree_3<Normalized_geom_traits> Subsegments_tree;
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
  typedef Intersection_data_structure_3<Normalized_geom_traits, 
					Triangle_3> Subfacets_tree;
  typedef Intersection_data_structure_3<Normalized_geom_traits,
					Segment_3> Subsegments_tree;
#endif
  typedef boost::shared_ptr<Subfacets_tree> Subfacets_tree_ptr;
  typedef boost::shared_ptr<Subsegments_tree> Subsegments_tree_ptr;
//   typedef typename Subsegments_tree::Point_with_index
//   Intersection_point;

  template <
    class Surface,
    class Point_creator,
    class Visitor,
    class Tag
    >
  friend class Surface_mesher::Polyhedral_oracle;

  typedef Surface_mesher::Polyhedral_oracle<Self> Surface_mesher_traits_3;

  Polyhedral_surface_3(const double sharp_edges_angle_lower_bound = 60.,
		       const double sharp_edges_angle_upper_bound = 180.,
                       const double sharp_vertices_angle_lower_bound = 30.,
                       const double sharp_vertices_angle_upper_bound = 180.) :
    subfacets_tree_ptr(new Subfacets_tree()),
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
    subsegments_tree_ptr(new Subsegments_tree(false, true, false)),
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
    subsegments_tree_ptr(new Subsegments_tree()),
#endif
    input_vertices_ptr(new Input_vertices()),
    corner_vertices_ptr(new Corner_vertices()),
    edges_vertices_ptr(new Edges_vertices()),
    sharp_edges_angle_lower_bound(sharp_edges_angle_lower_bound),
    sharp_edges_angle_upper_bound(sharp_edges_angle_upper_bound),
    sharp_vertices_angle_lower_bound(sharp_vertices_angle_lower_bound),
    sharp_vertices_angle_upper_bound(sharp_vertices_angle_upper_bound)
  {
  }

  Polyhedral_surface_3(std::istream& input_file,
		       const double sharp_edges_angle_lower_bound = 60.,
		       const double sharp_edges_angle_upper_bound = 180.,
                       bool auto_construct_octree = true) :
    subfacets_tree_ptr(new Subfacets_tree()),
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
    subsegments_tree_ptr(new Subsegments_tree(false, true, false)),
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
    subsegments_tree_ptr(new Subsegments_tree()),
#endif
    input_vertices_ptr(new Input_vertices()),
    corner_vertices_ptr(new Corner_vertices()),
    edges_vertices_ptr(new Edges_vertices()),
    sharp_edges_angle_lower_bound(sharp_edges_angle_lower_bound),
    sharp_edges_angle_upper_bound(sharp_edges_angle_upper_bound)
  {

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    CGAL::Timer timer;
    std::cerr << "Creating polyhedron... ";
    timer.start();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    CGAL::scan_OFF(input_file, *this, true);
    if(!input_file)
      return;
    this->compute_bounding_box();
    this->compute_normals();

    bounding_box = Bbox(this->xmin(),
                        this->ymin(),
                        this->zmin(),
                        this->xmax(),
                        this->ymax(),
                        this->zmax());

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
    bounding_box_sq_radius *= 3;

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    timer.stop();
    std::cerr << 
      ::boost::format("done (%1%s)\n"
		      "  number of vertices: %2%\n"
		      "  number of edges:    %3%\n"
		      "  number of facets:   %4%\n")
      % timer.time()
      % this->size_of_vertices()
      % ( this->size_of_halfedges() / 2 )
      % this->size_of_facets();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    if(auto_construct_octree)
    {
      compute_sharp_edges_incidence_graph();
      construct_octree();
    }
  } // end of Polyhedral_surface_3 constructor

  void construct_octree()
  {
#if 0
    GT gt = GT();
    typename GT::Construct_orthogonal_vector_3 orthogonal_vector = 
      gt.construct_orthogonal_vector_3_object();
    typename GT::Compute_squared_length_3 squared_length =
      gt.compute_squared_length_3_object();
    typename GT::Compute_scalar_product_3 scalar_product = 
      gt.compute_scalar_product_3_object();

//    typedef CGAL::Polyhedron_3<GT> Polyhedron_3;
//    typedef Enriched_polyhedron<GT,Enriched_items> Polyhedron_3;

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

    const double cosine_bound = 
      std::cos(CGAL_PI * sharp_edges_angle_lower_bound / 360. );
    
    const FT cosine_squared_bound = cosine_bound * cosine_bound;
#endif // if 0
    typedef typename Polyhedron_3::Vertex_const_iterator Vertex_const_it;

#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
    double (*vertices_array)[3] =
      (double (*)[3])new double [3*this->size_of_vertices()];
//       static_cast<double(*)[3]>(new double [3*this->size_of_vertices()]);
    int (*facets_array)[3] = 
      (int (*)[3])new int [3*this->size_of_facets()];
    int i = 0;
    typedef Inverse_index<Vertex_const_it> Index;
    Index index( this->vertices_begin(), this->vertices_end());
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
    for(typename Polyhedron_3::Vertex_const_iterator
          vit = this->vertices_begin(),
          end = this->vertices_end();
        vit != end;
        ++vit)
    {
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
      subfacets_tree_ptr->add_constrained_vertex(vit->point());
#endif
//       input_vertices_ptr->insert(vit);
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
      vertices_array[i][0] = vit->point().x();
      vertices_array[i][1] = vit->point().y();
      vertices_array[i][2] = vit->point().z();
      ++i;
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
    }

    typename Polyhedron_3::size_type facet_index = 0;
    for(typename Polyhedron_3::Facet_const_iterator fit = 
          this->facets_begin();
        fit != this->facets_end();
        ++fit, ++facet_index)
    {
      CGAL_assertion(fit->is_triangle());
      typename Polyhedron_3::Halfedge_around_facet_const_circulator 
        edges_circ = fit->facet_begin();

      const Point_3& p1 = edges_circ++->vertex()->point();
      const Point_3& p2 = edges_circ++->vertex()->point();
      const Point_3& p3 = edges_circ++->vertex()->point();

#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
      edges_circ = fit->facet_begin();
      facets_array[facet_index][0] = index[edges_circ++->vertex()];
      facets_array[facet_index][1] = index[edges_circ++->vertex()];
      facets_array[facet_index][2] = index[edges_circ++->vertex()];
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON


#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
      CGAL_assertion(fit->tag() >= 0);
      subfacets_tree_ptr->add_constrained_facet(p1, p2, p3, fit->tag());
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
      subfacets_tree_ptr->add_element(Triangle_3(p1, p2, p3));
#endif     
#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	std::cerr << ::boost::format("new facet: #%4% (%1%, %2%, %3%)\n")
	  % p1 %  p2 % p3 % facet_index;
#endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
    }

#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
# ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    std::cerr << "Creating the PointInPolyhedron data structure... ";
    timer.reset();
    timer.start();
# endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    pinpolyhedron_ptr = 
      PointInPolyhedron_ptr(new PointInPolyhedron(vertices_array,
						  this->size_of_vertices(),
						  facets_array,
						  this->size_of_facets()));
    delete [] vertices_array;
    delete [] facets_array;
# ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    timer.stop();
    std::cerr << ::boost::format("done (%1%s)\n") % timer.time();
# endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    std::cerr << "Creating subfacets_tree... ";
    CGAL::Timer timer;
    timer.start();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    subfacets_tree_ptr->create_data_structure();
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    timer.stop();
# ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
    std::cerr <<
      ::boost::format("done (%1%s)\n"
		      "  number of facets:      %3%\n"
		      "  number of constraints: %4% (in subfacets_tree)\n")
      % timer.time()
      % subfacets_tree_ptr->number_of_vertices()
      % subfacets_tree_ptr->number_of_facets()
      % subfacets_tree_ptr->number_of_constraints();
# endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
# ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
    std::cerr <<
      ::boost::format("done (%1%s)\n"
		      "  number of facets:      %2%\n")
      % timer.time()
      % subfacets_tree_ptr->number_of_elements();
# endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    if( this->has_edges() ) {

      for(typename Polyhedron_3::Edge_const_iterator eit = 
	    this->edges_begin();
	  eit != this->edges_end();
	  ++eit)
      {
        if(eit->tag() >= 0)
        {
	  const Point_3 pa = eit->vertex()->point();
	  const Point_3 pb = eit->opposite()->vertex()->point();
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
	  subsegments_tree_ptr->add_constrained_edge(pa, pb, eit->tag());
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
	  subsegments_tree_ptr->add_element(Segment_3(pa,  pb));
# ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  std::cerr << 
	    ::boost::format("new edge: (%1%, %2%)")
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
	      ::boost::format("facets=(%1%, %2%)")
	      % Triangle()(*(eit->facet()))
	      % Triangle()(*(eit->opposite()->facet()));
	  }
# endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
          
        }
      }

#if 0
      typedef std::map<Point_3, int> Edges_vertex_counter;
      Edges_vertex_counter edges_vertex_counter;

      for(typename Polyhedron_3::Edge_const_iterator eit = 
	    this->edges_begin();
	  eit != this->edges_end();
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
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
	  subsegments_tree_ptr->add_constrained_edge(pa, pb);
#endif
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
	  subsegments_tree_ptr->add_element(Segment_3(pa,  pb));
# ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  std::cerr << 
	    ::boost::format("new edge: (%1%, %2%)")
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
# endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
#endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
	}
      }

      for(typename Edges_vertex_counter::const_iterator it = 
	    edges_vertex_counter.begin();
	  it != edges_vertex_counter.end();
	  ++it)
      {
	input_vertices_ptr->erase(it->first);
	if(it->second != 2)
	{
#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  std::cerr << ::boost::format("corner point: (%1%)\n")
	    % it->first;
#endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  corner_vertices_ptr->push_back(it->first);
	}
	else
	{
#ifdef CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  std::cerr << ::boost::format("edge point: (%1%)\n")
	    % it->first;
#endif // CGAL_POLYHEDRAL_SURFACE_VERBOSE_CONSTRUCTION
	  edges_vertices_ptr->push_back(it->first);
	}
      }
//       if(!corner_vertices_ptr->empty() && edges_vertices_ptr->empty())
//       {
// #ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
// 	std::cerr << "Incorrect input data. "
// 		  << "Swap corner vertices and edges vertices...\n";
// #endif
// 	std::swap(corner_vertices, edges_vertices);
//       }

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      std::cerr << "Shuffle edges vertices... ";
#endif
      std::random_shuffle(edges_vertices_ptr->begin(), edges_vertices_ptr->end());
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      std::cerr << "done\n";
#endif

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      std::cerr <<
	::boost::format("number of corner vertices: %1%\n"
			"number of edges vertices:  %2%\n")
	% corner_vertices_ptr->size()
	% edges_vertices_ptr->size();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION


#endif // if 0

#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      std::cerr << "Creating subsegments_tree... ";
      timer.reset();
      timer.start();
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION

      subsegments_tree_ptr->create_data_structure();
#ifdef CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
      timer.stop();

# ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
      std::cerr <<
	::boost::format("done (%1%s)\n"
			"  number of edges:       %3%\n"
			"  number of constraints: %4% (in subsegments_tree)\n")
	% timer.time()
	% subsegments_tree_ptr->number_of_vertices()
	% subsegments_tree_ptr->number_of_edges()
	% subsegments_tree_ptr->number_of_constraints();
# endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
# ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
      std::cerr <<
	::boost::format("done (%1%s)\n"
			"  number of edges:       %2%\n")
	% timer.time()
	% subsegments_tree_ptr->number_of_elements();
# endif // CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_INTERSECTION_DATA_STRUCTURE
#endif // CGAL_SURFACE_MESHER_DEBUG_POLYHEDRAL_SURFACE_CONSTRUCTION
    } // end "if(this->has_edges())"
//     subfacets_tree_ptr->input(input_file,
//                            std::back_inserter(input_vertices));
  } // end construct_octree()

  void set_sharp_edges_angle_bounds(double lower_bound,
                                    double upper_bound) {
    sharp_edges_angle_lower_bound = lower_bound;
    sharp_edges_angle_upper_bound = upper_bound;
  }

  void set_sharp_vertices_angle_bounds(double lower_bound,
                                       double upper_bound) {
    sharp_vertices_angle_lower_bound = lower_bound;
    sharp_vertices_angle_upper_bound = upper_bound;
  }

  unsigned int tag_border_edges()
  {
    unsigned int nb = 0;
    for(Halfedge_iterator he = edges_begin();
        he != edges_end();
        he++)
    {
      const bool tag = ( he->sharp() ||
                         he->is_border() ||
                         he->opposite()->is_border() );
      he->sharp() = tag;
      he->opposite()->sharp() = tag;
      nb += tag ? 1 : 0;
    }
    return nb;
  }

  void compute_sharp_edges_incidence_graph()
  {
    this->tag_sharp_edges(sharp_edges_angle_lower_bound);
    this->tag_border_edges();
    construct_incidence_graph();
  }

  const Bbox& bbox() const
  {
    return bounding_box;
  }

  const FT& bounding_sphere_squared_radius() const
  {
    return bounding_box_sq_radius;
  }

  static void new_sub_edge(Graph_edge_node& edge_node,
                           const int edge_index,
                           Halfedge_handle& he)
  {
    he->tag(edge_index);
    he->opposite()->tag(edge_index);
    edge_node.sub_edges.push_back(he);
  }

  static void new_sub_facet(Graph_facet_node& facet_node,
                            const int facet_index,
                            Facet_handle& fh)
  {
    fh->tag(facet_index);
    facet_node.sub_facets.push_back(fh);
  }

  Vertex_handle follow_the_edge(Graph_edge_node& edge_node,
                                const int edge_index,
                                Halfedge_handle he) const
  {
    using boost::tie;

    bool can_continue;
    Halfedge_handle next_he;

    tie(can_continue, next_he) = can_follow_sharp_edges(he);
    while(can_continue && next_he->tag() < 0)
    {
      CGAL_assertion(he->vertex()->tag()<0);
      he->vertex()->tag(edge_index);
      new_sub_edge(edge_node, edge_index, next_he);
      he = next_he;
      tie(can_continue, next_he) = Self::can_follow_sharp_edges(he);
    }
    if(can_continue)
    {
      he->vertex()->tag(edge_index);
      return Vertex_handle();
    }
    else
      return he->vertex();
  }

  void setup_incident_vertex_to_an_edge(const int edge_index,
                                        const Vertex_handle& v)
  {
    if(v->tag() < 0)
    {
      Graph_vertex_node vertex_node;
      vertex_node.incident_edges.insert(edge_index);
      v->tag(incidence_graph.vertices.size());
      incidence_graph.vertices.push_back(vertex_node);
    }
    else
      incidence_graph.vertices[v->tag()].incident_edges.insert(edge_index);
  }

  void construct_incidence_graph()
  {
    incidence_graph.clear();
    input_vertices_ptr->clear();
    corner_vertices_ptr->clear();
    edges_vertices_ptr->clear();
    this->tag_vertices(-1);
    this->tag_halfedges(-1);
    this->tag_facets(-1);
    for(Halfedge_iterator he = halfedges_begin();
        he != halfedges_end();
        he++)
    {
      if(he->sharp() && he->tag() < 0)
      {
        const int edge_index = incidence_graph.edges.size();
        Graph_edge_node edge_node;

        new_sub_edge(edge_node, edge_index, he);
        const Vertex_handle& v1 = follow_the_edge(edge_node,
                                                  edge_index,
                                                  he);
        if(v1 != Vertex_handle())
          setup_incident_vertex_to_an_edge(edge_index, v1);
        const Vertex_handle& v2 = follow_the_edge(edge_node,
                                                  edge_index,
                                                  he->opposite());
        if(v2 != Vertex_handle())
          setup_incident_vertex_to_an_edge(edge_index, v2);
        incidence_graph.edges.push_back(edge_node);
      } // end "if halfedge he is sharp and not tagged"
    } // end "for all halfedges"

    for(Facet_iterator fit = facets_begin(), end = facets_end();
        fit != end;
        ++fit)
    {
      if(fit->tag() < 0)
      {
        const int facet_index = incidence_graph.facets.size();
        Graph_facet_node facet_node;
        construct_graph_facet(fit,
                              facet_index,
                              facet_node,
                              incidence_graph);
        incidence_graph.facets.push_back(facet_node);
      }
    }

    // display to std::cerr, for debugging
    std::cerr << "INCIDENCE GRAPH SUMMARY\n";
    std::cerr << "Number of vertices: " << incidence_graph.vertices.size()
              << "\nNumber of edges: " << incidence_graph.edges.size()
              << "\nNumber of facets: " << incidence_graph.facets.size()
              << std::endl;
    std::cerr << "VERTICES\n";
    for(int i = 0, end_i = incidence_graph.vertices.size(); i < end_i; ++i)
    {
      std::cerr << boost::format("Vertex #%1% incident edges: ") % i;
      for(typename Graph_vertex_node::Incident_edges::iterator
            sub_edge_id_it = incidence_graph.vertices[i].incident_edges.begin(),
            end = incidence_graph.vertices[i].incident_edges.end();
          sub_edge_id_it != end; ++sub_edge_id_it)
      {
        std::cerr << 
          boost::format("%1% ") % *sub_edge_id_it;
      }
      std::cerr << std::endl;      
    }
    std::cerr << "EDGES\n";
    for(int i = 0, end_i = incidence_graph.edges.size(); i < end_i; ++i)
    {
      std::cerr << 
        boost::format("Edge #%1% number of sub-edges: %2%\n")
        % i % incidence_graph.edges[i].sub_edges.size();
      std::cerr << 
        boost::format("Edge #%1% incident_facets: ") % i;
      for(typename Graph_edge_node::Incident_facets::iterator
            facet_id_it = incidence_graph.edges[i].incident_facets.begin(),
            end = incidence_graph.edges[i].incident_facets.end();
          facet_id_it != end; ++facet_id_it)
      {
        std::cerr << 
          boost::format("%1% ") % *facet_id_it;
      }
      std::cerr << std::endl;      
    }
    std::cerr << "FACETS\n";
    for(int i = 0, end_i = incidence_graph.facets.size(); i < end_i; ++i)
    {
      std::cerr << 
        boost::format("Facet #%1% number of sub-facets: %2%\n")
        % i % incidence_graph.facets[i].sub_facets.size();
    }


    for(typename Polyhedron_3::Vertex_iterator 
          vit = this->vertices_begin(),
          end = this->vertices_end();
        vit != end;
        ++vit)
    {
      switch(this->type(vit)) {
      case Polyhedron::SMOOTH:
        vit->tag(vit->halfedge()->facet()->tag());
        CGAL_assertion(vit->tag()>=0);
        input_vertices_ptr->insert(vit);
        break;
      case Polyhedron::DART:
      case Polyhedron::CORNER:
        CGAL_assertion(vit->tag()>=0);
        corner_vertices_ptr->push_back(vit);
      case Polyhedron::CREASE_REGULAR:
      case Polyhedron::CREASE_IRREGULAR:
        edges_vertices_ptr->push_back(vit);
      }
    }
  }

  void construct_graph_facet(Facet_handle seed_facet,
                             const int facet_index,
                             Graph_facet_node& facet_node,
                             Incidence_graph& graph)
  {
    typedef std::queue<Facet_handle> Facets_queue;
    Facets_queue facets_queue;
    
    facets_queue.push(seed_facet);

    while(!facets_queue.empty())
    {
      Facet_handle fh = facets_queue.front();
      facets_queue.pop();
      if(fh->tag() < 0)
      {
        fh->tag(facet_index);
        facet_node.sub_facets.push_back(fh);
          
        Halfedge_around_facet_circulator fit = fh->facet_begin();
        Halfedge_around_facet_circulator end = fit;
        CGAL_For_all(fit,end)
        {
          const int edge_index = fit->tag();
          const Halfedge_handle& neighbor_facet_he = fit->opposite();
          if(!neighbor_facet_he->is_border())
          {
            if(edge_index < 0) // if not yet handled...
            {
              facets_queue.push(neighbor_facet_he->facet());
            }
            else
            {
              graph.edges[edge_index].incident_facets.insert(facet_index);
            }
          }
        } // end for all halfedge around the facet fh
      } // end if fh not yet handled
    }
  } // end of construct_graph_facet(...)

  std::pair<bool, Halfedge_handle>
  can_follow_sharp_edges(const Halfedge_handle& he) const
  {
    typename Geom_traits::Construct_vector_3 vector
      = Geom_traits().construct_vector_3_object();

    const Vertex_handle& v = he->vertex();

    Halfedge_handle result = Halfedge_handle();
    unsigned int nb_sharp_edges = 0;

    Halfedge_around_vertex_circulator other_he = v->vertex_begin();
    Halfedge_around_vertex_circulator end = other_he;
    CGAL_For_all(other_he,end)
    {
      if(other_he->sharp())
      {
        ++nb_sharp_edges;
        if(he != other_he) 
          result = other_he;
      }
    }
    CGAL_assertion(this->nb_sharp_edges(v)==nb_sharp_edges);

    if(nb_sharp_edges == 2)
    {
      CGAL_assertion(result->vertex() == he ->vertex());
      const Vector_3 vector_he = vector(he->opposite()->vertex()->point(),
                                        he->vertex()->point());
      const Vector_3 vector_result_opposite =
        vector(result->vertex()->point(),
               result->opposite()->vertex()->point());
      return std::make_pair(true, result->opposite());
    }
    else 
      return std::make_pair(false, Halfedge_handle());
  }

#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_OCTREE
  void gl_draw_facet_octree()
  {
    subfacets_tree_ptr->gl_draw_nodes();
  }

  void gl_draw_edges_octree()
  {
    subsegments_tree_ptr->gl_draw_nodes();
  }
#endif

  // draw edges
  void gl_draw_sharp_edges_with_names(const float line_width,
                                      unsigned char r,
                                      unsigned char g,
                                      unsigned char b)
  {
    ::glLineWidth(line_width);
    ::glColor3ub(r,g,b);

    for(Halfedge_iterator he = edges_begin();
        he != edges_end();
        he++)
    {
      if(he->sharp())
      {
        const Point_3& a =  he->opposite()->vertex()->point();
        const Point_3& b =  he->vertex()->point();
        ::glPushName(he->tag()+incidence_graph.vertices.size());
        ::glBegin(GL_LINES);
        ::glVertex3d(a[0],a[1],a[2]);
        ::glVertex3d(b[0],b[1],b[2]);
        ::glEnd();
        ::glPopName();
      }
    }
  }

  void gl_draw_direct_triangles_with_name(bool smooth_shading,
                                          bool use_normals,
                                          bool inverse_normals = false)
  {
    // draw triangles
    Facet_iterator pFacet = facets_begin();
    for(;pFacet != facets_end();pFacet++)
    {
      ::glPushName(pFacet->tag()
                   + incidence_graph.vertices.size()
                   + incidence_graph.edges.size());
      ::glBegin(GL_TRIANGLES);
      gl_draw_facet(pFacet,smooth_shading,use_normals,inverse_normals);
      ::glEnd(); // end polygon assembly
      ::glPopName();
    }
  }

  void gl_draw_almost_all_triangles(int facet_index_not_drawed,
                                    bool smooth_shading,
                                    bool use_normals,
                                    bool inverse_normals = false)
  {
    // draw triangles
    ::glBegin(GL_TRIANGLES);
    Facet_iterator pFacet = facets_begin();
    for(;pFacet != facets_end();pFacet++)
      if(pFacet->tag() != facet_index_not_drawed)
        gl_draw_facet(pFacet,smooth_shading,use_normals,inverse_normals);
    ::glEnd(); // end polygon assembly
  }

  template <typename Vertex_handle>
  bool vertices_not_on_same_curve(const Vertex_handle& v1,
                                  const Vertex_handle& v2) const
  {
    struct Display_curves_indices {
      void operator()(const std::set<int>& container) const {
        std::cerr << "(";
        for(std::set<int>::const_iterator 
              it = container.begin(),
              end = container.end();
            it != end;)
        {
          std::cerr << *it;
          if(++it!=end)
            std::cerr << " ";
        }
        std::cerr << ")";
      }
    };
//     if(v1->point().dimension() < 0) {
//       if(v2->point().dimension() < 0)
//         return false; // both in volume
//       else
//         return true; // v1 in volume and v2 on a surface
//     }
//     else
//       if(v2->point().dimension() < 0)
//         return true; // v2 in volume, and v1 on a surface

    if(v1->point().dimension() < 0 || v2->point().dimension() < 0)
      return true;
    std::set<int> incident_edges_v1, incident_edges_v2, intersection;
    incident_edges(v1, CGAL::inserter(incident_edges_v1));
    incident_edges(v2, CGAL::inserter(incident_edges_v2));
#ifdef CGAL_SURFACE_MESHER_DEBUG_INCIDES
    Display_curves_indices()(incident_edges_v1);
    Display_curves_indices()(incident_edges_v2);
#endif // CGAL_SURFACE_MESHER_DEBUG_INCIDES
    std::set_intersection(incident_edges_v1.begin(), incident_edges_v1.end(),
                          incident_edges_v2.begin(), incident_edges_v2.end(),
                          CGAL::inserter(intersection));
    return intersection.empty();
  }

  template <typename Vertex_handle>
  bool vertices_not_on_same_surface_patch(const Vertex_handle& v1,
                                          const Vertex_handle& v2,
                                          const Vertex_handle& v3) const
  {
//     if(v1->point().dimension() < 0) {
//       if(v2->point().dimension() < 0 && v3->point().dimension() < 0)
//         return false; // all three vertices in volume
//       else
//         return true; // v1 in volume, and (v2 or v3) on a surface
//     }
//     else
//       if(v2->point().dimension() < 0 || v3->point().dimension() < 0
    if(v1->point().dimension() < 0 ||
       v2->point().dimension() < 0 ||
       v3->point().dimension() < 0)
      return true;
    std::set<int> incident_facets_v1, incident_facets_v2, incident_facets_v3, intersection;
    incident_facets(v1, CGAL::inserter(incident_facets_v1));
    incident_facets(v2, CGAL::inserter(incident_facets_v2));
    incident_facets(v3, CGAL::inserter(incident_facets_v3));

    std::set_intersection(incident_facets_v1.begin(), incident_facets_v1.end(),
                          incident_facets_v2.begin(), incident_facets_v2.end(),
                          CGAL::inserter(intersection));
    if(intersection.empty())
      return true;
    std::set<int> intersection2;
    std::set_intersection(intersection.begin(), intersection.end(),
                          incident_facets_v3.begin(), incident_facets_v3.end(),
                          CGAL::inserter(intersection2));
    return intersection2.empty();
  }

  template <typename Vertex_handle, typename OutputIterator>
  void incident_edges(const Vertex_handle& v,
                      OutputIterator out_it) const
  {
    CGAL_assertion(v->point().element_index() >= 0);
    switch(v->point().dimension())
    {
    case 0: {
      const int index = v->point().element_index();
      CGAL_assertion(index <  incidence_graph.vertices.size() );
      CGAL_assertion(incidence_graph.vertices[index].incident_edges.size() >= 1 );
      std::copy(incidence_graph.vertices[index].incident_edges.begin(),
                incidence_graph.vertices[index].incident_edges.end(),
                out_it);
      break;
    }
    case 1:
      *out_it++ = v->point().element_index();
      break;
    default: 
      // if dimension=2, nothing to output
      // if dimension=-1 (point in volume, nothing to output)
      break;
    }
  }

  template <typename Vertex_handle, typename OutputIterator>
  void incident_facets(const Vertex_handle& v,
                       OutputIterator out_it) const
  {
    CGAL_assertion(v->point().element_index() >= 0);
    switch(v->point().dimension())
    {
    case 0: {
      std::vector<int> incident_edges;
      const int index = v->point().element_index();
      CGAL_assertion(index <  incidence_graph.vertices.size() );
      std::copy(incidence_graph.vertices[index].incident_edges.begin(),
                incidence_graph.vertices[index].incident_edges.end(),
                std::back_inserter(incident_edges));
      for(std::vector<int>::const_iterator
            edge_index_it = incident_edges.begin(),
            end = incident_edges.end();
          edge_index_it != end;
          ++edge_index_it)
      {
        std::copy(incidence_graph.edges[*edge_index_it].incident_facets.begin(),
                  incidence_graph.edges[*edge_index_it].incident_facets.end(),
                  out_it);
      }
      break;
    }
    case 1: {
      const int index = v->point().element_index();
      CGAL_assertion(index <  incidence_graph.edges.size() );
      std::copy(incidence_graph.edges[index].incident_facets.begin(),
                incidence_graph.edges[index].incident_facets.end(),
                out_it);
      break;
    }
    default: // dimension=2
      *out_it++ = v->point().element_index();
      break;
    }
  }

public:
  Subfacets_tree_ptr subfacets_tree_ptr;
  Subsegments_tree_ptr subsegments_tree_ptr;
  typedef std::set<Vertex_handle, Compare_vertex_iterators> Input_vertices;
  boost::shared_ptr<Input_vertices> input_vertices_ptr;
  typedef std::vector<Vertex_handle> Corner_vertices;
  boost::shared_ptr<Corner_vertices> corner_vertices_ptr;
  typedef std::vector<Vertex_handle> Edges_vertices;
  boost::shared_ptr<Edges_vertices> edges_vertices_ptr;
#ifdef CGAL_SURFACE_MESHER_POLYHEDRAL_SURFACE_USE_PINPOLYHEDRON
  typedef boost::shared_ptr<PointInPolyhedron> PointInPolyhedron_ptr;
  PointInPolyhedron_ptr pinpolyhedron_ptr;
#endif

  Bbox bounding_box;
  FT bounding_box_sq_radius;
  double sharp_edges_angle_lower_bound;
  double sharp_edges_angle_upper_bound;
  double sharp_vertices_angle_lower_bound;
  double sharp_vertices_angle_upper_bound;
  Incidence_graph incidence_graph;
};

} // end namespace CGAL

#endif

#endif // CGAL_POLYHEDRAL_SURFACE_3_H
