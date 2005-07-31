// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_CONSTRAINED_REGULAR_TRIANGULATION_3_H
#define CGAL_CONSTRAINED_REGULAR_TRIANGULATION_3_H

#include <iostream>
#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/utility.h>

#include <list>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <CGAL/Constraint_hierarchy_2.h>
#include <CGAL/Delaunay_mesher_2.h>

CGAL_BEGIN_NAMESPACE

template <class CT_3, class Traits_2_3>
struct Special_mesh_traits_2 : public Traits_2_3
{
  typedef typename CT_3::Vertex_handle Vertex_handle_3;
  typedef typename CT_3::Cell_handle Cell_handle_3;
  typedef typename CT_3::Weighted_point Weighted_point;
  typedef typename Traits_2_3::Point_2 Point_2;

  CT_3& t;
  
  Special_mesh_traits_2<CT_3, Traits_2_3>(CT_3& tr) : t(tr) {};
  
  typedef int Quality;
  
  struct Is_bad
  {
    CT_3& t;
    
    Is_bad(CT_3& tr) : t(tr) {};
    
    bool operator()(const Point_2& pa,
		    const Point_2& pb,
		    const Point_2& pc,
		    Quality& q) const 
    {
      q = 0; // constant, do not sort faces
      
      const Vertex_handle_3& va = t.insert(Weighted_point(pa));
      const Vertex_handle_3& vb = t.insert(Weighted_point(pb));
      const Vertex_handle_3& vc = t.insert(Weighted_point(pc));
      
      Cell_handle_3 c;
      int i, j, k;
      return ! t.is_facet(va, vb, vc, c, i, j, k);
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(t); }
};

/** PLC_loader<CT_3. CT_2>.
    CT_3: Contrained_triangulation_3,
    CT_2: Contrained_Delaunay_triangulation_2,
 
    \todo{TODO:
    Il faudrait propifier ca.
    Avec un namespace PLC, peut-etre, contenant les definitions de types.}
*/
template <typename CT_3, typename CT_2>
struct PLC_loader
{
  typedef PLC_loader<CT_3, CT_2 > Self;

  typedef typename CT_2::Point Point_3;
  typedef typename CT_2::Vertex_handle Vertex_handle_2;
  typedef typename CT_2::Face_handle Face_handle_2;
  typedef typename CT_2::Geom_traits Geom_traits_2;
  typedef typename CT_3::Weighted_point Weighted_point;
  typedef typename CT_3::Vertex_handle Vertex_handle_3;
  typedef typename CT_3::Cell_handle Cell_handle_3;
  typedef typename CT_3::Locate_type Locate_type_3;
  typedef typename CT_3::Geom_traits Traits_3;
  typedef typename Traits_3::Vector_3 Vector_3;
  typedef typename Traits_3::FT FT;

  typedef typename CT_3::Sharp_vertices Sharp_vertices;
  typedef typename CT_3::Sharp_vertices_iterator Sharp_vertices_iterator;
  typedef typename CT_3::Finite_facets_iterator Finite_facets_iterator;
  
  struct Conform_extras {
    CT_3& t;

    explicit Conform_extras(CT_3& tr) : t(tr) {};
    
    bool is_bad(const CT_2& tr2,
		const Face_handle_2& fh, const int index) const {
      const Point_3& a = fh->vertex(tr2.cw(index) )->point();
      const Point_3& b = fh->vertex(tr2.ccw(index))->point();

      const Vertex_handle_3 va = t.insert(Weighted_point(a));
      const Vertex_handle_3 vb = t.insert(Weighted_point(b));

      Cell_handle_3 c;
      int i, j;
      return not t.is_edge(va, vb, c, i, j);
    }

    // for signal_inserted_before_vertex_in_edge
    Vertex_handle_2 va, vb;

    void signal_before_inserted_vertex_in_edge(const CT_2 &,
					       const Face_handle_2& fh, 
					       const int edge_index,
					       const Point_3&)
    {
      va = fh->vertex(CT_2::cw(edge_index) );
      vb = fh->vertex(CT_2::ccw(edge_index));
    }

    void signal_after_inserted_vertex_in_edge(const CT_2&,
					      const Vertex_handle_2& v) const
    {
      const Point_3& p = v->point();

      const Point_3& a = va->point();
      const Point_3& b = vb->point();

      const Vertex_handle_3 v3a = t.insert(Weighted_point(a));
      const Vertex_handle_3 v3b = t.insert(Weighted_point(b));

      Cell_handle_3 c;
      int i, j;
      if( t.is_edge(v3a, v3b, c, i, j) )
	t.insert(Weighted_point(p), CT_3::EDGE, c, i, j);
      else
	t.insert(Weighted_point(p), v3a->cell());
      // MY_DEBUG
      std::cerr << "inserted in edge: " << p << std::endl
		<< "(edge " << a << ", " << b << ")" << std::endl
		<< "t.number_of_vertices()=" << t.number_of_vertices()
		<< std::endl;
    }

    template<class EdgeIt, class FaceIt>
    void signal_before_inserted_vertex_in_face(const CT_2&,
					       const Face_handle_2& fh,
					       EdgeIt,
					       EdgeIt,
					       FaceIt,
					       FaceIt,
					       const Point_3& p)
    {
      // MY_DEBUG
      std::cerr << "inserted in faces: " << p << std::endl
		<< "(face "
		<< fh->vertex(0)->point() << ", "
		<< fh->vertex(1)->point() << ", "
		<< fh->vertex(2)->point() << ")" << std::endl;
      
    }

    void signal_after_inserted_vertex_in_face(const CT_2&,
					      const Vertex_handle_2& v) const
    {
      t.insert(Weighted_point(v->point()));
      // MY_DEBUG
      std::cerr << "t.number_of_vertices()=" << t.number_of_vertices()
		<< std::endl;
    }

  }; //end of class Conform_extras
  
  typedef std::pair<Point_3, Point_3> Constraint_2;

  /** \todo{ TODO: faire un peu mieux qu'hard-coder la std::list...} */
  typedef std::list<Constraint_2> Constraints_2;
  typedef typename Constraints_2::const_iterator Constraints_2_iterator;

  typedef std::list<Point_3> Seeds;
  typedef typename Seeds::const_iterator Seeds_iterator;

  struct Constrained_face {
    Constraints_2 edges;
    Vector_3 normal;
    Seeds seeds;
    bool seeds_in_face;
  };

  typedef typename std::list<Constrained_face>::const_iterator CF_it;

  typedef std::list<const Constrained_face*> Vertex_3d_context;

  struct Vertex_2d_context_in_edge {
    const Constrained_face* face;
    const Constraint_2* edge;
    int position; // 0 or 1;
//     bool operator<(const Vertex_2d_context_in_edge& other) const
//     {
//       return this < &other;
//     }
  };

  typedef std::list<Vertex_2d_context_in_edge> Vertex_2d_context;


  typedef typename Vertex_3d_context::const_iterator
                                              Vertex_3d_context_iterator;
  typedef typename Vertex_2d_context::const_iterator
                                              Vertex_2d_context_iterator;
  
  struct Vertex_context
  {
    Vertex_2d_context context_2;
    Vertex_3d_context context_3;
  };

  typedef std::map<Vertex_handle_3, Vertex_context>
                                                Vertices_context;
  typedef typename Vertices_context::const_iterator
                                                Vertices_context_iterator;

  // datas
  Vertices_context vertices_context;
  CT_3& t;
  Traits_3 traits_3;
  
  Cell_handle_3 cell_hint;

  PLC_loader(CT_3& tr, Traits_3 traits = Traits_3())
    : t(tr),
      traits_3(traits), 
      cell_hint(Cell_handle_3())
  {};

  void add_vertex_context(const Point_3& p,
			  const CF_it& cf_it,
			  const Constraints_2_iterator& c_it,
			  const int pos)
  {
    // MY_DEBUG
    std::cerr << "add_context(" << p << ")\n";
    std::cerr << "  " << c_it->first << ", "
	      << c_it->second
	      << std::endl
	      << "t.size()=" << t.number_of_vertices() << std::endl;
    for(typename CT_3::Finite_vertices_iterator
	  it = t.finite_vertices_begin();
	it != t.finite_vertices_end();
	++it)
      std::cerr << it->point() << std::endl;

    Locate_type_3 lt;
    int i, j;
    Cell_handle_3 c = t.locate(Weighted_point(p), lt, i, j);
    Vertex_handle_3 v;
    if( lt != CT_3::VERTEX)
      v = t.insert(p);
    else
      v = c->vertex(i);

    Vertex_context& context = vertices_context[v];

    // MY_DEBUG
    std::cerr << &*v << " " << &context << " " << lt << std::endl;
    
    Vertex_2d_context_in_edge c_in_edge;

    c_in_edge.face = &*cf_it;
    c_in_edge.edge = &*c_it;
    c_in_edge.position = pos;
    
    context.context_2.push_back(c_in_edge);
    context.context_3.push_back(&*cf_it);
  }

  void load_triangulation(std::istream& is)
  {
    int number_of_faces;
    is >> number_of_faces;

    std::list<Constrained_face> list_of_faces;

    for(int i = 0; i < number_of_faces; ++i)
      {
	int number_of_segments;
	int number_of_seeds;

	Constrained_face face;

	is >> number_of_segments >> number_of_seeds;

	for(int j = 0; j < number_of_segments; ++j)
	  {
	    Point_3 a;
	    Point_3 b;

	    is >> a >> b;
	    face.edges.push_back(std::make_pair(a, b));
	  }

	Vector_3 normal;

	is >> normal;

	/** \todo{TODO:
	    Comment calculer la normal à une face ?
	    } */

	face.normal = normal;

	for(int k = 0; k < number_of_seeds; ++k)
	  {
	    Point_3 seed;

	    is >> seed;
	    face.seeds.push_back(seed);
	  }
	
	face.seeds_in_face = false;

	list_of_faces.push_back(face);
      }

    create_triangulation(list_of_faces.begin(),
			 list_of_faces.end());
  }

  /** CD_it is an forward iterator of value_type Constrained_face. */
  void create_triangulation( const CF_it& constrained_faces_begin,
			     const CF_it& constrained_faces_end )
  {
    for(CF_it cf_it = constrained_faces_begin;
	cf_it != constrained_faces_end; ++cf_it)
      {
// 	typedef std::set<Point_3> Points_of_face;
// 	typedef typename Points_of_face::const_iterator
// 	                                        Points_of_face_iterator;

// 	Points_of_face points_of_face; // set of points of the face *cf_it

	for(Constraints_2_iterator c_it = cf_it->edges.begin();
	    c_it != cf_it->edges.end();
	    ++c_it) // *c_it is a Contraint_2, that is, a pair of Point_3
	  {
	    std::cerr << "edge(" << c_it->first // MY_DEBUG
		      << ", " << c_it->second << ")" << std::endl;
	    add_vertex_context(c_it->first, cf_it, c_it, 0);
	    add_vertex_context(c_it->second, cf_it, c_it, 1);
	  }
      }

    for(Vertices_context_iterator v_context_it = vertices_context.begin();
	v_context_it != vertices_context.end();
	++v_context_it)
      {
	const Vertex_handle_3& vh = v_context_it->first;

	// 3d context
	if(v_context_it->second.context_3.size()>1)
	  {
	    // double for:
	    //   for i=0 to end
	    //     for j=0 to i-1
	    for(Vertex_3d_context_iterator v_3d_context_it1 = 
		  v_context_it->second.context_3.begin();
		v_3d_context_it1 != v_context_it->second.context_3.end();
		++v_3d_context_it1)
	      for(Vertex_3d_context_iterator v_3d_context_it2 =
		    v_context_it->second.context_3.begin();
		  v_3d_context_it2 != v_3d_context_it1;
		  ++ v_3d_context_it2)
		{
		  typename Traits_3::Compute_scalar_product_3 scalar = 
		    traits_3.compute_scalar_product_3_object();
		  /** \todo{TODO, WARNING, ERROR
		      Pourquoi angle_3 ne prend-il pas des
		      Vector_3 ??}
		  */
		  if( scalar((*v_3d_context_it1)->normal,
			     (*v_3d_context_it2)->normal)
		      > typename Traits_3::RT(0) )
		    t.sharp_vertices[vh]=FT(0);
		}
	  }

	// 2d context
	if(v_context_it->second.context_2.size()>1)
	  {
	    // double for:
	    //   for i=0 to end
	    //     for j=0 to i-1
	    for(Vertex_2d_context_iterator v_2d_context_it1 = 
		  v_context_it->second.context_2.begin();
		v_2d_context_it1 != v_context_it->second.context_2.end();
		++v_2d_context_it1)
	      for(Vertex_2d_context_iterator v_2d_context_it2 =
		    v_context_it->second.context_2.begin();
		  v_2d_context_it2 != v_2d_context_it1;
		  ++ v_2d_context_it2)
		{
		  const Point_3& a = vh->point().point();

		  Point_3 b1;
		  Point_3 b2;

		  if( v_2d_context_it1->position == 0 )
		    b1 = v_2d_context_it1->edge->first;
		  else
		    b1 = v_2d_context_it1->edge->second;

		  if( v_2d_context_it2->position == 0 )
		    b2 = v_2d_context_it2->edge->first;
		  else
		    b2 = v_2d_context_it2->edge->second;

		  if ( b1 != b2 )
		    {
		      typename Traits_3::Angle_3 angle_3 = 
			traits_3.angle_3_object();
		      
		      if( angle_3(b1, a, b2) == CGAL::ACUTE )
			t.sharp_vertices[vh]=FT(0);
		    }
		}
	  }
      } // here, t.sharp_vertices has been fully filled.

    // sharp vertices protection!
    for(Sharp_vertices_iterator sh_v_it = t.sharp_vertices.begin();
	sh_v_it != t.sharp_vertices.end();
	++sh_v_it)
      {
	const Vertex_handle_3& vh = sh_v_it->first;
	Vertex_context& vertex_context = vertices_context[vh];
	Vertex_2d_context vertex_2d_context = vertex_context.context_2;
	Vertex_3d_context vertex_3d_context = vertex_context.context_3;

	typename Traits_3::Compute_squared_distance_3 distance_3 = 
	  traits_3.compute_squared_distance_3_object();

	std::set<Constraint_2> simplified_2d_context;
	// pointers to iterators, because iterators of std::list do not have
	// operator < 
	// TODO, WARNING: je ne suis pas sur que ca soit legal, de se fier
	// aux adresses des iterateurs.
	// ERROR: D'ailleurs, seuls les random access iterators sont censes
	// avoir  un operator <. :-(

	// MY_DEBUG
	std::cerr << "point " << vh->point().point() << std::endl;
	std::cerr << "size=" << vertex_2d_context.size() << std::endl;

	for(Vertex_2d_context_iterator v_2d_context_it = 
	      vertex_2d_context.begin();
	    v_2d_context_it != vertex_2d_context.end();
	    ++v_2d_context_it)
	  {	// MY_DEBUG
	    std::cerr << v_2d_context_it->edge << ":"
		      << v_2d_context_it->edge->first << " "
		      << v_2d_context_it->edge->second << std::endl;

	    const Point_3& a = v_2d_context_it->edge->first;
	    const Point_3& b = v_2d_context_it->edge->second;
	    simplified_2d_context.insert(std::make_pair(a,b));
	  }
	std::cerr << "simplified.size()=" << simplified_2d_context.size()
		  << std::endl;
	// -- brute force LFS for sharp vertices --

	// initialization
	FT squared_lfs =
	  distance_3(vertex_2d_context.begin()->edge->first,
		     vertex_2d_context.begin()->edge->second);
	// MY_DEBUG
	std::cerr << "i. lfs " << squared_lfs << std::endl;
	
	for(CF_it cf_it = constrained_faces_begin;
	    cf_it != constrained_faces_end; ++cf_it)
	  {
	    // vertex/face distance
	    if( std::find(vertex_3d_context.begin(), vertex_3d_context.end(),
			  &*cf_it) == vertex_3d_context.end() )
	      {
		typename Traits_3::Construct_plane_3 construct_plane_3 =
		  traits_3.construct_plane_3_object();
		
		typename Traits_3::Plane_3 plane = 
		  construct_plane_3(cf_it->edges.begin()->first,
				    cf_it->normal);
		squared_lfs = CGAL::min(distance_3(vh->point().point(),
						   plane),
					squared_lfs);
		/** \todo{
		    TODO, ERROR:
		    C'est faux. Il faudrait vraiment calculer la distance
		    entre la face et le sharp_vertex. La c'est la distance
		    entre le plan porteur de la face et le sharp vertex.
		    } */


		// MY_DEBUG
		std::cerr << "v/f lfs " << squared_lfs << std::endl;
	      }

	    // vertex/edge distance
	    for(Constraints_2_iterator c_it = cf_it->edges.begin();
		c_it != cf_it->edges.end();
		++c_it)
	      {
		const Point_3& a = c_it->first;
		const Point_3& b = c_it->second;
		if( simplified_2d_context.find(std::make_pair(a, b)) == 
		    simplified_2d_context.end() && 
		    simplified_2d_context.find(std::make_pair(b, a)) == 
		    simplified_2d_context.end() )
		  {
		    CGAL_assertion(vh->point().point() != c_it->first);
		    CGAL_assertion(vh->point().point() != c_it->second);
		    typename Traits_3::Construct_segment_3
		      construct_segment_3 =
		        traits_3.construct_segment_3_object();
		    typename Traits_3::Segment_3 segment_3 =
		      construct_segment_3(c_it->first, c_it->second);
		    
		    squared_lfs = 
		      CGAL::min(distance_3(vh->point().point(), segment_3),
				squared_lfs);

		    // MY_DEBUG
		    std::cerr << "v/e lfs " << squared_lfs << std::endl;
		  }
	      }
	  }
	
	std::cerr << "lfs(" << vh->point().point() << ") = "
		  << squared_lfs << std::endl;

	t.sharp_vertices[vh] = squared_lfs;

      }

    for(CF_it cf_it = constrained_faces_begin;
	cf_it != constrained_faces_end; ++cf_it)
      {
	typedef Delaunay_mesh_2<CT_2, Conform_extras> Mesh_2;

	// MY_DEBUG
	std::cerr << "face: ";

	Conform_extras extra = Conform_extras(t);

	Geom_traits_2 gt = Geom_traits_2(t);

	Mesh_2* tr2 = new Mesh_2(gt, extra);

	typedef std::queue<CGAL::Triple<Vertex_handle_3,
                                        Vertex_handle_3,
                                        Vertex_handle_3> > Facet_queue;

	Facet_queue queue;

	// MY_DEBUG
	std::cerr << cf_it->edges.size() << " edges." << std::endl;

	for(Constraints_2_iterator c_it = cf_it->edges.begin();
	    c_it != cf_it->edges.end();
	    ++c_it) // *c_it is a Contraint_2, that is, a pair of Point_3
	  tr2->insert(c_it->first, c_it->second); // insert the
						 // constrained edge in tr2

	tr2->set_seeds(cf_it->seeds.begin(),
		       cf_it->seeds.end(),
		       cf_it->seeds_in_face);

	tr2->refine_mesh();

	for(typename CT_2::Finite_edges_iterator e_it_2 = 
	      tr2->finite_edges_begin();
	    e_it_2 != tr2->finite_edges_end();
	    ++e_it_2)
	  {
	    const Point_3& a = e_it_2->first->
	      vertex(tr2->cw(e_it_2->second) )->point();
	    const Point_3& b = e_it_2->first->
	      vertex(tr2->ccw(e_it_2->second))->point();

	    const Vertex_handle_3& va = t.insert(Weighted_point(a));
	    const Vertex_handle_3& vb = t.insert(Weighted_point(b));

	    std::cerr << "e: " << t.insert_constrained_edge(va, vb)
		      << std::endl;
	  }

	// MY_DEBUG
	std::cerr << "t.number_of_vertices()=" << t.number_of_vertices()
		  << std::endl;

	for(typename CT_2::Finite_faces_iterator f_it_2 =
	      tr2->finite_faces_begin();
	    f_it_2 != tr2->finite_faces_end();
	    ++f_it_2)
	  {
	    if( ! f_it_2->is_marked() ) continue;

	    Vertex_handle_3 v0 = t.insert(f_it_2->vertex(0)->point());
	    Vertex_handle_3 v1 = t.insert(f_it_2->vertex(1)->point());
	    Vertex_handle_3 v2 = t.insert(f_it_2->vertex(2)->point());

	    std::cerr << "f: " << t.insert_constrained_facet(v0, v1, v2)
		      << std::endl
		      << "(" << v0->point()
		      << ", " << v1->point() 
		      << ", " << v2->point()
		      << ")" << std::endl;
	    int n = 0;
	    for(Finite_facets_iterator it = t.finite_facets_begin();
		it != t.finite_facets_end();
		++it)
	      if(t.is_constrained(*it))
		++n;
	    std::cerr << n << "facets" << std::endl;
	  }
	delete tr2;
      }
  } // end of PLC_loader::create_triangulation(begin, end)

	
}; // end of class PLC_loader

template <class CT_3>
struct Conformer_2D
{
  typedef typename CT_3::Bare_point Point_3;
  typedef typename CT_3::Vertex_handle Vertex_handle_3;
  typedef typename CT_3::Cell_handle Cell_handle_3;
  typedef typename CT_3::Locate_type Locate_type_3;
  typedef typename CT_3::Geom_traits Traits_3;
  typedef typename Traits_3::Vector_3 Vector_3;
  typedef typename Traits_3::FT FT;
  
  typedef typename CT_3::Sharp_vertices Sharp_vertices;
  typedef typename CT_3::Sharp_vertices_iterator Sharp_vertices_iterator;

  typedef std::pair<Vertex_handle_3, Vertex_handle_3> Contraint_2;
  typedef std::queue<Contraint_2> Queue_2D;

  // data members
  CT_3& t;
  Traits_3 traits;
  Queue_2D queue_2d;

  Conformer_2D(CT_3& tr) : t(tr), traits() {};

  void push_2d_constraint(const Contraint_2& contraint)
    {
      queue_2d.push(contraint);
    }

  bool is_edge(const Contraint_2& contraint)
    {
      Cell_handle_3 c;
      int i,j;
      return t.is_edge(contraint->first, contraint->second,
		       c, i, j);
    }

  Point_3 middle(const Contraint_2& contraint) const 
    {
      typename Traits_3::Construct_midpoint_3 midpoint = 
	traits.construct_midpoint_3_object();

      return midpoint(contraint->first->point().point(),
		      contraint->second->point().point());
    }

  void conform()
    {
      while( ! queue_2d.empty() )
	{
	  Contraint_2 contraint = queue_2d.front();
	  queue_2d.pop();

	  if( is_edge(contraint) ) continue; // next

	  const Point_3 p = middle(contraint);

	}
    }
};	  
	    
/** Create a CT from an OFF file, in an std::istream */
template <typename CT_3, typename CT_2>
struct Off_loader
{
  typedef typename CT_2::Point Point_3;
  typedef typename CT_2::Vertex_handle Vertex_handle_2;
  typedef typename CT_2::Face_handle Face_handle_2;
  typedef typename CT_2::Geom_traits Geom_traits_2;
  typedef typename CT_3::Weighted_point Weighted_point;
  typedef typename CT_3::Vertex_handle Vertex_handle_3;
  typedef typename CT_3::Cell_handle Cell_handle_3;
  typedef typename CT_3::Locate_type Locate_type_3;
  typedef typename CT_3::Geom_traits Traits_3;
  typedef typename Traits_3::Vector_3 Vector_3;
  typedef typename Traits_3::FT FT;

  typedef typename CT_3::Sharp_vertices Sharp_vertices;
  typedef typename CT_3::Sharp_vertices_iterator Sharp_vertices_iterator;

  typedef PLC_loader<CT_3, CT_2> PLC_load;
  typedef typename PLC_load::Constrained_face Constrained_face;

  // datas
  CT_3& t;
  std::list<Constrained_face> faces;

  Off_loader(CT_3& tr)
    : t(tr), faces()
  {};

  void load_triangulation(std::istream& is, bool verbose = false )
  {
    t.clear();
    off_file_input(is, verbose);
    PLC_load(t).create_triangulation(faces.begin(), faces.end());
  }

  bool
  off_file_input(std::istream& is, bool verbose)
  {
    faces.clear();

    File_scanner_OFF scanner(is, verbose);
    if (! is) {
      if (scanner.verbose()) {
	std::cerr 
	  << " " << std::endl
	  << "Constrained_regular_triangulation_3::off_file_input"
	  << std::endl
	  << " input error: file format is not OFF." << std::endl;
      }
      return false;
    }
    
    std::vector<Point_3> vp(scanner.size_of_vertices());
    
    // insert points
    int i;
    for ( i = 0; i < scanner.size_of_vertices(); i++) {
      Point_3 p;
      file_scan_vertex( scanner, p);
      vp[i] = p;
      scanner.skip_to_next_vertex( i);
    }
    
    if ( ! is ) {
      is.clear( std::ios::badbit);
      return false;
    }
    
    // inserts constrained edges and facets
    for ( i = 0; i < scanner.size_of_facets(); i++) {
      Constrained_face face;

      Integer32 no;
      scanner.scan_facet( no, i);
      if( ! is ) {
	is.clear( std::ios::badbit);
	return false;
      }

      std::vector<Point_3> points(no);
      Integer32 index0;
      Integer32 before;
      for(Integer32 k = 0; k < no; ++k)
	{
	  Integer32 index;
	  scanner.scan_facet_vertex_index( index, i);

	  if( k == 0 )
	    index0 = index;
	  else
	    {
	      face.edges.push_back(std::make_pair(vp[before], vp[index]));
	      if( k == (no - 1) )
		face.edges.push_back(std::make_pair(vp[index], vp[index0]));
	    }
	  before = index;
	}
      /** \todo{ Use kernel } */
      face.normal = CGAL::cross_product( points[1] - points[0],
					 points[no - 1] - points[0] );
      faces.push_back(face);
    }
    return true;
  }

};
  

/** TODO
    \todo{Ne pas deriver d'une Regular_triangulation...} */
template <class Tr, class CT_2>
class Constrained_regular_triangulation_3 : public Tr
{
public:
  typedef Constrained_regular_triangulation_3 Self;
  typedef Tr Triangulation;
  typedef typename Triangulation::Geom_traits Geom_traits;

  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Cell_handle Cell_handle;
  typedef typename Geom_traits::Bare_point Bare_point;
  typedef typename Triangulation::Weighted_point Weighted_point;
  typedef Weighted_point Point_3;
  typedef typename Triangulation::Locate_type Locate_type;

  typedef typename Triangulation::Facet Facet;
  typedef typename Triangulation::Edge Edge;

  typedef typename Triangulation::Finite_vertices_iterator
                                                 Finite_vertices_iterator;
  typedef typename Triangulation::Finite_edges_iterator
                                                 Finite_edges_iterator;
  typedef typename Triangulation::Finite_facets_iterator
                                                 Finite_facets_iterator;

  typedef typename Triangulation::Facet_circulator Facet_circulator;

  typedef typename Geom_traits::FT FT;

  typedef Constraint_hierarchy_2<Vertex_handle, bool>
                                                 Constraints_2_hierarchy;

  // TODO
  /** \todo{clear(), copy(), swap(), and in Mesh_2 too. */

  void clear()
  {
    Triangulation::clear();
    sharp_vertices.clear();
    hierarchy_2.clear();
  }

  // --- PUBLIC DATA MEMBERS ---

  typedef std::map<Vertex_handle, FT>            Sharp_vertices;
  typedef typename Sharp_vertices::const_iterator 
                                                 Sharp_vertices_iterator;

  Sharp_vertices sharp_vertices;
  Constraints_2_hierarchy hierarchy_2;

  // --- PRIVATE DATA MEMBERS

private:
  typedef std::queue<std::pair<Vertex_handle,
			       Vertex_handle> > Edges_to_be_conformed;
  typedef std::queue<Triple<Vertex_handle,
			    Vertex_handle,
			    Vertex_handle> > Facets_to_be_conformed;
  
  Facets_to_be_conformed facets_to_be_conformed;
  Edges_to_be_conformed edges_to_be_conformed;
public:

  // --- IO FUNCTIONS ---

  Vertex_handle off_file_input( std::istream& is, bool verbose = false);
  void off_file_output( std::ostream& os);

  // --- CONSTUCTORS ---
  Constrained_regular_triangulation_3(const Geom_traits& gt = Geom_traits())
    : Triangulation(gt)
  {}

  // --- CONSTRAINTS HANDLING ---

  bool is_constrained(const Facet& f) const
  {
    return f.first->is_constrained(f.second);
  }

  bool is_constrained(Edge e) const 
  {
    const Cell_handle& c = e.first;
    return is_constrained(c->vertex(e.second), c->vertex(e.third));
  }

  bool is_constrained(const Vertex_handle& va,
                      const Vertex_handle& vb) const
  {
    return va->is_adjacent_by_constraint(vb);
  }

  // --- INSERTIONS ---

  void insert_bounding_box()
  {
    FT xmin, xmax;
    FT ymin, ymax;
    FT zmin, zmax;

    Finite_vertices_iterator vi = this->finite_vertices_begin();

    xmin=xmax=vi->point().x();
    ymin=ymax=vi->point().y();
    zmin=zmax=vi->point().z();

    while(vi != this->finite_vertices_end())
      {
	if(vi->point().x() < xmin) xmin=vi->point().x();
	if(vi->point().x() > xmax) xmax=vi->point().x();
	if(vi->point().y() < ymin) ymin=vi->point().y();
	if(vi->point().y() > ymax) ymax=vi->point().y();
	if(vi->point().z() < ymin) ymin=vi->point().z();
	if(vi->point().z() > ymax) ymax=vi->point().z();
	vi++;
      }

    FT xcenter=(xmin+xmax)/2;
    FT ycenter=(ymin+ymax)/2;
    FT zcenter=(zmin+zmax)/2;
    FT xspan = (xmax-xmin)/2;
    FT yspan = (ymax-ymin)/2;
    FT zspan = (zmax-zmin)/2;

    xmax+=xspan/2;
    xmin-=xspan/2;
    ymax+=yspan/2;
    ymin-=yspan/2;
    zmax+=zspan/2;
    zmin-=zspan/2;

    Vertex_handle va = insert(Bare_point(-xmin, -ymin, -zmin));
    Vertex_handle vb = insert(Bare_point( xmin, -ymin, -zmin));
    Vertex_handle vc = insert(Bare_point( xmin,  ymin, -zmin));
    Vertex_handle vd = insert(Bare_point(-xmin,  ymin, -zmin));
    Vertex_handle ve = insert(Bare_point(-xmin, -ymin,  zmin));
    Vertex_handle vf = insert(Bare_point( xmin, -ymin,  zmin));
    Vertex_handle vg = insert(Bare_point( xmin,  ymin,  zmin));
    Vertex_handle vh = insert(Bare_point(-xmin,  ymin,  zmin));
  }

  // -- points insertions --
  Vertex_handle insert(const Bare_point & p, Cell_handle start = NULL);
  Vertex_handle insert(const Weighted_point & p, Cell_handle start = NULL);

  Vertex_handle insert(const Weighted_point & p, Locate_type lt,
	               Cell_handle c, int li, int);

  bool insert_constrained_edge(const Vertex_handle& va,
			       const Vertex_handle& vb,
			       bool update_hierachy = true)
  {
    Cell_handle ch;
    int i, j;

    if( !is_edge(va, vb, ch, i, j) )
      return false; /** \todo Conform the edge. */
    else
      {
	va->set_is_adjacent_by_constraint(vb, true);
	vb->set_is_adjacent_by_constraint(va, true);
	if( update_hierachy ) hierarchy_2.insert_constraint(va, vb);
      }
    return true;
  }

  void remove_constrained_edge(Edge e)
  {
    const Cell_handle& c = e.first;
    c->vertex(e.second)->set_is_adjacent_by_constraint(c->vertex(e.third),
						       false);
    c->vertex(e.third)->set_is_adjacent_by_constraint(c->vertex(e.second),
						      false);
  }

  bool insert_constrained_facet(const Vertex_handle& va,
				const Vertex_handle& vb,
				const Vertex_handle& vc)
  {
    Cell_handle c;
    int i, j, k;
    if( !is_facet(va, vb, vc,
		  c, i, j, k) )
      return false; /** \todo Force the facet into the triangulation. */
    else
      {
	const int l = 6-i-j-k;
	const Cell_handle& n = c->neighbor(l);
	c->set_constrained(l, true);
	n->set_constrained(n->index(c), true);
	return true;
      }
  }
private:

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q) const
  {
      CGAL_precondition(equal(p, q));
      return geom_traits().power_test_3_object()(p, q);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r) const
  {
      CGAL_precondition(collinear(p, q, r));
      return geom_traits().power_test_3_object()(p, q, r);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r, const Weighted_point &s) const
  {
      CGAL_precondition(coplanar(p, q, r, s));
      return geom_traits().power_test_3_object()(p, q, r, s);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r, const Weighted_point &s,
	     const Weighted_point &t) const
  {
      return geom_traits().power_test_3_object()(p, q, r, s, t);
  }

  bool in_conflict_3(const Weighted_point &p, const Cell_handle c) const
  {
      return side_of_power_sphere(c, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_2(const Weighted_point &p, const Cell_handle c, int i) const
  {
      return side_of_power_circle(c, i, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_1(const Weighted_point &p, const Cell_handle c) const
  {
      return side_of_power_segment(c, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_0(const Weighted_point &p, const Cell_handle c) const
  {
      return power_test(c->vertex(0)->point(), p) == ON_POSITIVE_SIDE;
  }

  class Conflict_tester_3
  {
      const Weighted_point &p;
      const Self *t;
      mutable std::vector<Vertex_handle> cv;

  public:

      Conflict_tester_3(const Weighted_point &pt, const Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const Cell_handle c) const
      {
	  // We mark the vertices so that we can find the deleted ones easily.
	  if (t->in_conflict_3(p, c))
	  {
	      for (int i=0; i<4; i++)
	      {
		  Vertex_handle v = c->vertex(i);
		  if (v->cell() != NULL)
		  {
		      cv.push_back(v);
		      v->set_cell(NULL);
		  }
	      }
	      return true;
	  }
	  return false;
      }

      std::vector<Vertex_handle> & conflict_vector()
      {
	  return cv;
      }
  };

  class Conflict_tester_2
  {
      const Weighted_point &p;
      const Self *t;
      mutable std::vector<Vertex_handle> cv;

  public:

      Conflict_tester_2(const Weighted_point &pt, const Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const Cell_handle c) const
      {
	  if (t->in_conflict_2(p, c, 3))
	  {
	      for (int i=0; i<3; i++)
	      {
		  Vertex_handle v = c->vertex(i);
		  if (v->cell() != NULL)
		  {
		      cv.push_back(v);
		      v->set_cell(NULL);
		  }
	      }
	      return true;
	  }
	  return false;
      }

      std::vector<Vertex_handle> & conflict_vector()
      {
	  return cv;
      }
  };

  friend class Conflict_tester_3;
  friend class Conflict_tester_2;


  template <class Conflict_test,
            class OutputIteratorBoundaryFacets,
            class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets,
         OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts_2(Cell_handle c, const Conflict_test &tester,
	           Triple<OutputIteratorBoundaryFacets,
                          OutputIteratorCells,
		          OutputIteratorInternalFacets> it) const
  {
    CGAL_triangulation_precondition( dimension()==2 );
    CGAL_triangulation_precondition( tester(c) );

    c->set_in_conflict_flag(1);
    *it.second++ = c;

    for (int i=0; i<3; ++i) {
      Cell_handle test = c->neighbor(i);
      if (test->get_in_conflict_flag() == 1) {
	  if (c < test)
	      *it.third++ = Facet(c, i); // Internal facet.
          continue; // test was already in conflict.
      }
      if (test->get_in_conflict_flag() == 0) {
	  if (tester(test)) {
	      if (c < test)
		  *it.third++ = Facet(c, i); // Internal facet.
              it = find_conflicts_2(test, tester, it);
	      continue;
	  }
	  test->set_in_conflict_flag(2); // test is on the boundary.
      }
      *it.first++ = Facet(c, i);
    }
    return it;
  }

  // Note: the code duplication between _2 and _3 should be avoided one day.
  template <class Conflict_test,
            class OutputIteratorBoundaryFacets,
            class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets,
         OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts_3(Cell_handle c, const Conflict_test &tester,
	           Triple<OutputIteratorBoundaryFacets,
                          OutputIteratorCells,
		          OutputIteratorInternalFacets> it) const
  {
    CGAL_triangulation_precondition( dimension()==3 );
    CGAL_triangulation_precondition( tester(c) );

    c->set_in_conflict_flag(1);
    *it.second++ = c;

    for (int i=0; i<4; ++i) {
      Cell_handle test = c->neighbor(i);
      if (test->get_in_conflict_flag() == 1) { // test was already in conflict.
	  if (c < test)
	      *it.third++ = Facet(c, i); // Internal facet.
          continue;
      }
      if ( !c->is_constrained(i) && test->get_in_conflict_flag() == 0) {
	  if (tester(test)) {
	      if (c < test)
		  *it.third++ = Facet(c, i); // Internal facet.
              it = find_conflicts_3(test, tester, it);
	      continue;
	  }
	  test->set_in_conflict_flag(2); // test is on the boundary.
      }
      *it.first++ = Facet(c, i);
    }
    return it;
  }

  // This one takes a function object to recursively determine the cells in
  // conflict, then calls _tds._insert_in_hole().
  template < class Conflict_test >
  Vertex_handle
  insert_conflict_2(Cell_handle c, const Conflict_test &tester)
  {
    CGAL_triangulation_precondition( dimension() == 2 );
    CGAL_triangulation_precondition( c != NULL );
    CGAL_triangulation_precondition( tester(c) );

    std::vector<Cell_handle> cells;
    cells.reserve(32);

    Facet facet;

    // Find the cells in conflict
    find_conflicts_2(c, tester, make_triple(Oneset_iterator<Facet>(facet),
		                            std::back_inserter(cells),
				            Emptyset_iterator()));

    // Create the new cells and delete the old.
    Vertex_handle vh = this->_tds._insert_in_hole(cells.begin(), cells.end(),
						  facet.first, facet.second);
    // Restore contraint status.
    /** \todo{WARNING, TODO Est-ce nécessaire en dimension 2 ?} */
    cells.clear();
    incident_cells(vh, std::back_inserter(cells));
    for(typename std::vector<Cell_handle>::iterator cit = cells.begin();
	cit != cells.end();
	++cit)
      {
	const int index = (*cit)->index(vh);
	for(int i = 0; i<4; ++i)
	  if( i != index)
	    (*cit)->set_constrained(i, false);
	Cell_handle ch = (*cit)->neighbor(index);
	const int mirror = (*cit)->mirror_index(index);
	(*cit)->set_constrained(index, ch->is_constrained(mirror));
      }
    return vh;
  }

  // This one takes a function object to recursively determine the cells in
  // conflict, then calls _tds._insert_in_hole().
  template < class Conflict_test >
  Vertex_handle
  insert_conflict_3(Cell_handle c, const Conflict_test &tester)
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    CGAL_triangulation_precondition( c != NULL );
    CGAL_triangulation_precondition( tester(c) );

    std::set<Cell_handle> cells;
    typedef typename std::set<Cell_handle>::iterator Cells_iterator;

    std::set<Edge> edges;
    typedef typename std::set<Edge>::iterator Edges_iterator;

    Facet facet;

    // Find the cells in conflict
    find_conflicts_3(c, tester, make_triple(Oneset_iterator<Facet>(facet),
		                            std::inserter(cells, 
							  cells.begin()),
				            Emptyset_iterator()));

    //    bool go_on = true; // If go_on is false, the point is not inserted.

    /** \todo{ TODO: optimiser ca, si possible} */
    for(Cells_iterator c_it = cells.begin();
	c_it != cells.end();
	++c_it)
      for( int k = 0; k < 4; k++)
	if( cells.find((*c_it)->neighbor(k)) == cells.end() )
	  // (c_it, k) is a sub-facet of the zone's boundary 
	  for( int i = 0; i < 4; ++i)
	    for( int j = 0; j < i ; ++j)
	      if( k != i && k != j )
		{
		  const int facet = 6-i-j-k; //sub-facet (i, j, k)
		  // if 'facet' is in the interior of the zone and
		  // contrained, the edge (i, j) is on the boundary of the
		  // zone (because k is on the boundary) and incident to
		  // 'facet' that is _in_ the zone and constrained.
		  if( (*c_it)->is_constrained(facet) &&
		      cells.find((*c_it)->neighbor(facet)) != cells.end() )
		    {
// 		      go_on = go_on && ! is_encroached(Edge(*c_it, i, j), 
// 						       tester.p.point());
		      edges.insert(Edge(*c_it, i, j));
		    }
		}

    // Create the new cells and delete the old.
    Vertex_handle vh = this->_tds._insert_in_hole(cells.begin(), cells.end(),
						  facet.first, facet.second);
    // Restore contraint status.
    cells.clear();
    incident_cells(vh, std::inserter(cells, cells.begin()));
    for(Cells_iterator cit = cells.begin();
	cit != cells.end();
	++cit)
      {
	const int index = (*cit)->index(vh);
	for(int i = 0; i<4; ++i)
	  if( i != index)
	    (*cit)->set_constrained(i, false);
	Cell_handle ch = (*cit)->neighbor(index);
	const int mirror = (*cit)->mirror_index(index);
	(*cit)->set_constrained(index, ch->is_constrained(mirror));
      }
    for(Edges_iterator e_it = edges.begin();
	e_it != edges.end();
	++e_it)
      {
	const Vertex_handle& va = e_it->first->vertex(e_it->second);
	const Vertex_handle& vb = e_it->first->vertex(e_it->third);

	/** \todo{ TODO: Optimiser, car insert_constrained_facet fait un
	    is_facet} */
	CGAL_assertion_code(bool b = ) insert_constrained_facet(va, vb, vh);
	CGAL_assertion(b);
      }

    return vh;
  }

public:
  template <class OutputIteratorBoundaryFacets,
            class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets,
         OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts(const Point_3 &p, Cell_handle c,
	         OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit,
		 OutputIteratorInternalFacets ifit) const
  {
      CGAL_triangulation_precondition(dimension() >= 2);

      std::vector<Cell_handle> cells;
      cells.reserve(32);
      std::vector<Facet> facets;
      facets.reserve(64);

      if (dimension() == 2) {
          Conflict_tester_2 tester(p, this);
	  ifit = find_conflicts_2(c, tester,
                                  make_triple(std::back_inserter(facets),
		                              std::back_inserter(cells),
                                              ifit)).third;
      }
      else {
          Conflict_tester_3 tester(p, this);
	  ifit = find_conflicts_3(c, tester,
                                  make_triple(std::back_inserter(facets),
		                              std::back_inserter(cells),
                                              ifit)).third;
      }

      // Reset the conflict flag on the boundary.
      for(typename std::vector<Facet>::iterator fit=facets.begin();
          fit != facets.end(); ++fit) {
        fit->first->neighbor(fit->second)->set_in_conflict_flag(0);
	*bfit++ = *fit;
      }

      // Reset the conflict flag in the conflict cells.
      for(typename std::vector<Cell_handle>::iterator ccit=cells.begin();
        ccit != cells.end(); ++ccit) {
        (*ccit)->set_in_conflict_flag(0);
	*cit++ = *ccit;
      }
      return make_triple(bfit, cit, ifit);
  }
};



template <class Tr, class CT_2>
inline 
typename Constrained_regular_triangulation_3<Tr, CT_2>::Vertex_handle
Constrained_regular_triangulation_3<Tr, CT_2>::
insert(const Bare_point & p, Cell_handle start) 
{
  return this->insert(Weighted_point(p), start);
}

template <class Tr, class CT_2>
typename Constrained_regular_triangulation_3<Tr, CT_2>::Vertex_handle
Constrained_regular_triangulation_3<Tr, CT_2>::
insert(const Weighted_point & p, Cell_handle start) 
{
    Locate_type lt;
    int li, lj;
    Cell_handle c = locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
}

template <class Tr, class CT_2>
typename Constrained_regular_triangulation_3<Tr, CT_2>::Vertex_handle
Constrained_regular_triangulation_3<Tr, CT_2>::
insert(const Weighted_point & p, Locate_type lt,
       Cell_handle c, int li, int lj)
{
  Vertex_handle v1, v2;
  bool insert_in_constrained_edge = false;

  if ( lt == Triangulation::EDGE &&
       c->vertex(li)->is_adjacent_by_constraint(c->vertex(lj)) )
    {
      remove_constrained_edge(Edge(c, li, lj));
      insert_in_constrained_edge = true;
      v1=c->vertex(li); //endpoint of the constraint
      v2=c->vertex(lj); // endpoint of the constraint
    }

  switch (dimension()) {
  case 3:
    {
      // TODO :
      // In case the point is completely equal (including weight), then we need
      // to discard it (don't update the triangulation, nor hide it), right ?
      if (! in_conflict_3(p, c)) {  // new point is hidden
          if (lt == Tr::VERTEX)
              return c->vertex(li); // by coinciding point
          else
              return NULL;          // by cell
      }

      // Should I mark c's vertices too ?
      Conflict_tester_3 tester(p, this);
      Vertex_handle v = insert_conflict_3(c, tester);
      v->set_point(p);
      for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
      {
        if ((*it)->cell() == NULL)
	{
          // vertex has to be deleted
          tds().delete_vertex(*it);
	}
      }

      if (insert_in_constrained_edge)
	{
	  // hierarchy_2.split_constraint(v1,v2,v);

	  CGAL_assertion_code(bool b1 =) insert_constrained_edge(v1, v, false);
	  CGAL_assertion_code(bool b2 =) insert_constrained_edge(v2, v, false);
	  CGAL_assertion( b1 && b2);

          //	  test_if_encroached(v1, v);
          //	  test_if_encroached(v2, v);
	}

      // TODO : manage the hidden points.
      return v;
    }
  case 2:
    {
      switch (lt) {
      case Tr::OUTSIDE_CONVEX_HULL:
      case Tr::CELL:
      case Tr::FACET:
      case Tr::EDGE:
      case Tr::VERTEX:
	{
          if (! in_conflict_2(p, c, 3)) {  // new point is hidden
              if (lt == Tr::VERTEX)
                  return c->vertex(li); // by coinciding point
              else
                  return NULL;          // by face
          }

	  Conflict_tester_2 tester(p, this);
	  Vertex_handle v = insert_conflict_2(c, tester);
	  v->set_point(p);
          for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
	  {
            if ((*it)->cell() == NULL)
	    {
              // vertex has to be deleted
              tds().delete_vertex(*it);
	    }
	  }
	  if (insert_in_constrained_edge)
	    hierarchy_2.split_constraint(v1,v2,v);

	  return v;
	}
      case Tr::OUTSIDE_AFFINE_HULL:
	{
	  // if the 2d triangulation is Regular, the 3d
	  // triangulation will be Regular
	  return Tr::insert_outside_affine_hull(p);
	}
      }
    }//dim 2
  case 1:
    {
      switch (lt) {
      case Tr::OUTSIDE_CONVEX_HULL:
      case Tr::EDGE:
      case Tr::VERTEX:
	{
          if (! in_conflict_1(p, c)) {  // new point is hidden
              if (lt == Tr::VERTEX)
                  return c->vertex(li); // by coinciding point
              else
                  return NULL;          // by edge
          }

	  Cell_handle bound[2];
          // corresponding index: bound[j]->neighbor(1-j) is in conflict.
	  std::vector<Vertex_handle>  hidden_vertices;
	  std::vector<Cell_handle>    conflicts;
          conflicts.push_back(c);

          // We get all cells in conflict,
          // and remember the 2 external boundaries.

	  for (int j = 0; j<2; ++j) {
	    Cell_handle n = c->neighbor(j);
	    while ( in_conflict_1( p, n) ) {
	      conflicts.push_back(n);
              hidden_vertices.push_back(n->vertex(j));
	      n = n->neighbor(j);
	    }
	    bound[j] = n;
	  }

          // We preserve the order (like the orientation in 2D-3D).

	  Vertex_handle v = tds().create_vertex();
	  v->set_point(p);
          Cell_handle c0 = tds().create_face(v, bound[0]->vertex(0), NULL);
          Cell_handle c1 = tds().create_face(bound[1]->vertex(1), v, NULL);
          tds().set_adjacency(c0, 1, c1, 0);
          tds().set_adjacency(bound[0], 1, c0, 0);
          tds().set_adjacency(c1, 1, bound[1], 0);
          bound[0]->vertex(0)->set_cell(bound[0]);
          bound[1]->vertex(1)->set_cell(bound[1]);
          v->set_cell(c0);

	  tds().delete_cells(conflicts.begin(), conflicts.end());
	  tds().delete_vertices(hidden_vertices.begin(), hidden_vertices.end());
	  if (insert_in_constrained_edge)
	    hierarchy_2.split_constraint(v1,v2,v);

	  return v;
	}
      case Tr::OUTSIDE_AFFINE_HULL:
	return Tr::insert_outside_affine_hull(p);
      case Tr::FACET:
      case Tr::CELL:
	// impossible in dimension 1
        CGAL_assertion(false);
	return NULL;
      }
    }
  case 0:
    {
        // We need to compare the weights when the points are equal.
        if (lt == Tr::VERTEX && in_conflict_0(p, c)) {
            CGAL_assertion(li == 0);
            c->vertex(li)->set_point(p); // replace by heavier point
        }
        else
            return Tr::insert(p, c);
    }
  default :
    {
      return Tr::insert(p, c);
    }
  }
}

template <class Tr, class CT_2>
typename Tr::Vertex_handle
Constrained_regular_triangulation_3<Tr, CT_2>::
off_file_input(std::istream& is, bool verbose)
{
  Vertex_handle vinf(0);
  File_scanner_OFF scanner(is, verbose);
  if (! is) {
    if (scanner.verbose()) {
         std::cerr 
	   << " " << std::endl
	   << "Constrained_regular_triangulation_3::off_file_input"
	   << std::endl
	   << " input error: file format is not OFF." << std::endl;
    }
    return vinf;
  }
  
  clear();
  
  std::vector<Vertex_handle> vvh(scanner.size_of_vertices());
  //  std::map<Vh_pair, Edge> edge_map;

  // insert points
  int i;
  for ( i = 0; i < scanner.size_of_vertices(); i++) {
    Bare_point p;
    file_scan_vertex( scanner, p);
    vvh[i] = insert(p); // insert the point in the triangulation
    scanner.skip_to_next_vertex( i);
  }

  if ( ! is ) {
    is.clear( std::ios::badbit);
    return vinf;
  }

  // inserts constrained edges and facets
  for ( i = 0; i < scanner.size_of_facets(); i++) {
    Integer32 no;
    scanner.scan_facet( no, i);
    if( ! is || (no != 3 && no!= 2) ) {
      if ( scanner.verbose()) {
	std::cerr 
	  << " " << std::endl
	  << "Constrained_regular_triangulation_3::off_file_input"
	  << "edge or facet " << i << "does not have 2 or 3 vertices." 
	  << std::endl;
      }
      is.clear( std::ios::badbit);
      return vinf;
    }

    Integer32 index0;
    Integer32 index1;
    scanner.scan_facet_vertex_index( index0, i);
    scanner.scan_facet_vertex_index( index1, i);
    if( no == 3 ) // facet
      {
	Integer32 index2;
	scanner.scan_facet_vertex_index( index2, i);
	/*bool r = */insert_constrained_facet(vvh[index0], vvh[index1],
					  vvh[index2]);
	//	CGAL_assertion(r);
      }
    else // edge
      {
	/*bool r = */insert_constrained_edge(vvh[index0], vvh[index1]);
	//	CGAL_assertion(r);
      }
  }
  return vinf;
}

template <class Tr, class CT_2>
void
Constrained_regular_triangulation_3<Tr, CT_2>::
off_file_output(std::ostream& os)
{
  std::set<Facet> constrained_facets;

  for(Finite_facets_iterator it = this->finite_facets_begin();
      it != this->finite_facets_end();
      ++it)
    if(is_constrained(*it))
      constrained_facets.insert(*it);

  os << "OFF" << std::endl
     << this->number_of_vertices() << " "
     << constrained_facets.size() << " " 
     << "0" << std::endl;

  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int counter = 0;
  for(Finite_vertices_iterator vit = this->finite_vertices_begin();
      vit != this->finite_vertices_end();
      ++vit)
    {
      V[vit] = counter++;
      os << vit->point().point() << std::endl;
    }

  for(typename std::set<Facet>::const_iterator fit =
	constrained_facets.begin();
      fit != constrained_facets.end();
      ++fit)
    {
      os << "3" << std::endl;
      for(int i = 0; i < 4; ++i)
	if(i != fit->second)
	  os << V[fit->first->vertex(i)] << " ";
      os << std::endl;
    }
}

CGAL_END_NAMESPACE

#endif // CGAL_CONSTRAINED_REGULAR_TRIANGULATION_3_H
