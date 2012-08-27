// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_H
#define CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_H

#include <map>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Apollonius_graph_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Apollonius_graph_hierarchy_vertex_base_2.h>

#include <CGAL/Apollonius_graph_2/Traits_wrapper_2.h>

namespace CGAL {


// parameterization of the  hierarchy
const unsigned int ag_hierarchy_2__ratio    = 30;
const unsigned int ag_hierarchy_2__minsize  = 20;
const unsigned int ag_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

template < class Gt,
	   class Agds = Triangulation_data_structure_2<
             Apollonius_graph_hierarchy_vertex_base_2<
               Apollonius_graph_vertex_base_2<Gt,true> >,
               Triangulation_face_base_2<Gt> >,
	   class LTag = Tag_false >
class Apollonius_graph_hierarchy_2
  : public Apollonius_graph_2< Gt, Agds, LTag >
{
private:
  typedef Apollonius_graph_2<Gt,Agds,LTag>  Apollonius_graph_base;
  typedef Apollonius_graph_base             Ag_base;

public:
  typedef Gt                                Geom_traits;
  typedef typename Gt::Site_2               Site_2;
  typedef typename Gt::Point_2              Point_2;

  typedef typename Ag_base::Vertex_handle    Vertex_handle;
  typedef typename Ag_base::Face_handle      Face_handle;
  typedef typename Ag_base::Edge             Edge;

  typedef typename Ag_base::Face_circulator       Face_circulator;
  typedef typename Ag_base::Edge_circulator       Edge_circulator;
  typedef typename Ag_base::Vertex_circulator     Vertex_circulator;

  typedef typename Ag_base::All_faces_iterator    All_faces_iterator;
  typedef typename Ag_base::Finite_faces_iterator Finite_faces_iterator;

  typedef typename Ag_base::All_vertices_iterator All_vertices_iterator;
  typedef typename Ag_base::Finite_vertices_iterator 
                                                     Finite_vertices_iterator;

  typedef typename Ag_base::All_edges_iterator    All_edges_iterator;
  typedef typename Ag_base::Finite_edges_iterator Finite_edges_iterator;

  typedef typename Ag_base::Sites_iterator         Sites_iterator;
  typedef typename Ag_base::Visible_sites_iterator Visible_sites_iterator;
  typedef typename Ag_base::Hidden_sites_iterator  Hidden_sites_iterator;

  typedef typename Ag_base::size_type              size_type;

  using Ag_base::insert_first;
  using Ag_base::insert_second;
  using Ag_base::insert_third;
  using Ag_base::is_hidden;
  using Ag_base::incircle;
  using Ag_base::edge_interior;
  using Ag_base::insert_degree_2;
  using Ag_base::initialize_conflict_region;
  using Ag_base::expand_conflict_region;
  using Ag_base::retriangulate_conflict_region;
  using Ag_base::is_infinite;

public:
  // CREATION
  //---------
  Apollonius_graph_hierarchy_2
  (const Geom_traits& gt = Geom_traits());

  template<class Input_iterator>
  Apollonius_graph_hierarchy_2(Input_iterator first,
			       Input_iterator beyond,
			       const Geom_traits& gt = Geom_traits())
    : Apollonius_graph_base(gt)
  {
    init_hierarchy(gt);
    insert(first, beyond);
  }

  Apollonius_graph_hierarchy_2
  (const Apollonius_graph_hierarchy_2& agh);

  Apollonius_graph_hierarchy_2&
  operator=(const Apollonius_graph_hierarchy_2& agh);

  ~Apollonius_graph_hierarchy_2();

protected:
  // used to initialize the hierarchy at construction time
  void init_hierarchy(const Geom_traits& gt);

public:
  // INSERTION
  //----------
  template < class Input_iterator >
  size_type insert(Input_iterator first, Input_iterator beyond)
  {
    // copy the sites to a local container
    typename Apollonius_graph_base::Site_list wp_list;
    //    Site_list wp_list;
    for (Input_iterator it = first; it != beyond; ++it) {
      wp_list.push_back(*it);
    }

    // sort by decreasing weight
    typename Apollonius_graph_base::Site_less_than_comparator
      less_than(this->geom_traits());
    std::sort(wp_list.begin(), wp_list.end(), less_than);

    // now insert
    typename Apollonius_graph_base::Site_list_iterator lit;
    for (lit = wp_list.begin(); lit != wp_list.end(); ++lit) {
      insert(*lit);
    }


    // store how many sites where in the range
    std::size_t num = wp_list.size();

    // clear the local container
    wp_list.clear();

    // return the number of sites in range
    return num;
  }

  Vertex_handle insert(const Site_2& p);
  Vertex_handle insert(const Site_2& p,
		       Vertex_handle vnear) {
    // the following statement has been added in order to avoid
    // a g++3.2.1_FreeBSD-RELEASE warning
    vnear = Vertex_handle();
    return insert(p);
  }

public:
  // REMOVAL
  //--------
  void remove(Vertex_handle v);

public:
  // NEAREST NEIGHBOR LOCATION
  //--------------------------
public:
  Vertex_handle nearest_neighbor(const Point_2& p) const;
  inline Vertex_handle nearest_neighbor(const Point_2& p,
					Vertex_handle /* vnear */) const {
    return nearest_neighbor(p);
  }

public:
  // VALIDITY CHECK
  //---------------
  bool is_valid(bool verbose = false, int level = 1) const;

public:
  // MISCELLANEOUS
  //--------------
  void clear();
  void swap(Apollonius_graph_hierarchy_2& agh);

  // I/O
  //----
  void file_input(std::istream& is);
  void file_output(std::ostream& os) const;

private:
  // private methods
  void
  nearest_neighbor(const Point_2& p,
		   Vertex_handle vnear[ag_hierarchy_2__maxlevel]) const; 

  int random_level();

  void copy(const Apollonius_graph_hierarchy_2 &agh);

private:
  // class variables
  // here is the stack of graphs which form the hierarchy
  Apollonius_graph_base*   hierarchy[ag_hierarchy_2__maxlevel];
  boost::rand48  random; // random generator

public:
  template<class OutputItFaces, class OutputItBoundaryEdges,
	   class OutputItHiddenVertices>
  boost::tuples::tuple<OutputItFaces, OutputItBoundaryEdges,
		       OutputItHiddenVertices>
  get_conflicts_and_boundary_and_hidden_vertices(const Site_2& p,
						 OutputItFaces fit,
						 OutputItBoundaryEdges eit,
						 OutputItHiddenVertices vit,
						 Vertex_handle start =
						 Vertex_handle()) const
  {
    Vertex_handle vnearest = nearest_neighbor(p.point(), start);
    return this->get_all(p, fit, eit, vit, vnearest, false);
  }

  template<class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces, OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Site_2& p,
			     OutputItFaces fit,
			     OutputItBoundaryEdges eit,
			     Vertex_handle start =
			     Vertex_handle()) const {
    boost::tuples::tuple<OutputItFaces,OutputItBoundaryEdges,Emptyset_iterator>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     fit,
						     eit,
						     Emptyset_iterator(),
						     start);
    return std::make_pair( boost::tuples::get<0>(tup),
			   boost::tuples::get<1>(tup) );
  }


  template<class OutputItBoundaryEdges, class OutputItHiddenVertices>
  std::pair<OutputItBoundaryEdges, OutputItHiddenVertices>
  get_boundary_of_conflicts_and_hidden_vertices(const Site_2& p,
						OutputItBoundaryEdges eit,
						OutputItHiddenVertices vit,
						Vertex_handle start =
						Vertex_handle()) const {
    boost::tuples::tuple<Emptyset_iterator,OutputItBoundaryEdges,
      OutputItHiddenVertices>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     Emptyset_iterator(),
						     eit,
						     vit,
						     start);
    return std::make_pair( boost::tuples::get<1>(tup),
			   boost::tuples::get<2>(tup) );
  }


  template<class OutputItFaces, class OutputItHiddenVertices>
  std::pair<OutputItFaces, OutputItHiddenVertices>
  get_conflicts_and_hidden_vertices(const Site_2& p,
				    OutputItFaces fit,
				    OutputItHiddenVertices vit,
				    Vertex_handle start =
				    Vertex_handle()) const {
    boost::tuples::tuple<OutputItFaces,Emptyset_iterator,
      OutputItHiddenVertices>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     fit,
						     Emptyset_iterator(),
						     vit,
						     start);
    return std::make_pair( boost::tuples::get<0>(tup),
			   boost::tuples::get<2>(tup) );
  }

  template<class OutputItFaces>
  OutputItFaces get_conflicts(const Site_2& p,
			      OutputItFaces fit,
			      Vertex_handle start = Vertex_handle()) const {
    boost::tuples::tuple<OutputItFaces,Emptyset_iterator,Emptyset_iterator>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     fit,
						     Emptyset_iterator(),
						     Emptyset_iterator(),
						     start);
    return boost::tuples::get<0>(tup);
  }

  template<class OutputItBoundaryEdges>
  OutputItBoundaryEdges
  get_boundary_of_conflicts(const Site_2& p,
			    OutputItBoundaryEdges eit,
			    Vertex_handle start = Vertex_handle()) const {
    boost::tuples::tuple<Emptyset_iterator,OutputItBoundaryEdges,
      Emptyset_iterator>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     Emptyset_iterator(),
						     eit,
						     Emptyset_iterator(),
						     start);
    return boost::tuples::get<1>(tup);
  }

  template<class OutputItHiddenVertices>
  OutputItHiddenVertices
  get_hidden_vertices(const Site_2& p,
		      OutputItHiddenVertices vit,
		      Vertex_handle start = Vertex_handle()) const {
    boost::tuples::tuple<Emptyset_iterator,Emptyset_iterator,
      OutputItHiddenVertices>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     Emptyset_iterator(),
						     Emptyset_iterator(),
						     vit,
						     start);
    return boost::tuples::get<2>(tup);
  }
};


template<class Gt, class Agds, class LTag>
std::ostream& operator<<(std::ostream& os,
			 const Apollonius_graph_hierarchy_2<Gt,Agds,LTag>& agh)
{
  agh.file_output(os);
  return os;
}

template<class Gt, class Agds, class LTag>
std::istream& operator>>(std::istream& is,
			 Apollonius_graph_hierarchy_2<Gt,Agds,LTag>& agh)
{
  agh.file_input(is);
  return is;
}

} //namespace CGAL


#include <CGAL/Apollonius_graph_2/Apollonius_graph_hierarchy_2_impl.h>


#endif // CGAL_APOLLONIUS_GRAPH_HIERARCHY_2_H
