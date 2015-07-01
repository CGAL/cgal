// Copyright (c) 2013 INRIA Sophia-Antipolis (France),
//               2014-2015 GeometryFactory (France).
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
// Author(s) : Lakulish Antani, Christophe Delage, Jane Tournois, Pierre Alliez
//

#ifndef CGAL_LIPSCHITZ_SIZING_FIELD_2_H
#define CGAL_LIPSCHITZ_SIZING_FIELD_2_H

#include <list>

#include <CGAL/basic.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

#include <CGAL/Mesh_2/Sizing_field_2.h>

namespace CGAL
{

template <typename Tr>
class Lipschitz_sizing_field_2
  : public virtual Sizing_field_2<Tr>
{
public:
  typedef typename Tr::Geom_traits      Geom_traits;
  typedef typename Geom_traits::Point_2 Point;
    
  typedef Delaunay_triangulation_2<Geom_traits> Delaunay_triangulation;
  typedef typename Delaunay_triangulation::All_faces_iterator Face_iterator;
  typedef typename Delaunay_triangulation::Finite_vertices_iterator Vertex_iterator;
  typedef typename Delaunay_triangulation::Face_circulator Face_circulator;
    
  typedef Search_traits_2<Geom_traits> Tree_traits;
  typedef Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Search_tree;
    
  typedef Apollonius_graph_traits_2<Geom_traits> Apollonius_traits;
  typedef Apollonius_graph_2<Apollonius_traits> Apollonius_graph;
  typedef typename Apollonius_traits::Site_2 Site;

public:    
  typedef std::list<Site> Site_set_2;

public:

  // comparison functor for sorting sites
  class Compare_site_2
  {
  public:
    bool operator ()(Site s1, Site s2)
    {
      return (s1.weight() > s2.weight());
    }
  };

public:
  // default constructor: empty point set and K = 1.0
  Lipschitz_sizing_field_2()
    : K(1.0)
  {
    sites = Site_set_2();

    generate_delaunay();
    extract_poles();
    generate_sites();
    generate_apollonius();
  }
    
  // constructor with point set and K as parameters
  template <typename InputIterator>
  Lipschitz_sizing_field_2(InputIterator first,
			   InputIterator beyond,
			   const double k = 1.0)
    : K(k)
  {
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  std::cout << "Building sizing field..." << std::flush;
#endif
    points = std::list<Point>(first, beyond);
    generate_delaunay();
    extract_poles();
    generate_sites();
    generate_apollonius();
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  std::cout << "done." << std::endl;
#endif
  }

  // constructor from a triangulation
  Lipschitz_sizing_field_2(Tr& tr, const double k = 1.0)
    : K(k)
  {
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
    std::cout << "Building sizing field..." << std::flush;
#endif
    points = std::list<Point>(tr.points_begin(), tr.points_end());
    generate_delaunay();
    extract_poles();
    generate_sites();
    generate_apollonius();
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
    std::cout << "done." << std::endl;
#endif
  }

  // assignment operator, copies point set and K
  Lipschitz_sizing_field_2& operator =(const Lipschitz_sizing_field_2& f)
  {
    points.clear();
    if(!poles.empty())
      poles.clear();

    if(!sites.empty())
      sites.clear();

    dt.clear();
    ag.clear();

    points = f.points;
    K = f.K;

    generate_delaunay();
    extract_poles();
    generate_sites();
    generate_apollonius();

    return *this;
  }

  double operator()(const Point& p) const
  {
    if(points.empty() || points.size() == 1)
      return K;
    Site ns = (*ag.nearest_neighbor(p)).site();
    return K * weighted_distance(p, ns);
  }

  void set_K(double k)
  {
    K = k;
  }

  double get_K()
  {
    return K;
  }

  std::list<Point>& get_poles()
  {
    return poles;
  }

protected:

  double K;
  Site_set_2 sites;
  std::list<Point> points;
  std::list<Point> poles;
  Apollonius_graph ag;
  Delaunay_triangulation dt;

  // generate the delaunay triangulation (for voronoi diagram)
  void generate_delaunay()
  {
    if(points.empty())
      return;

    for(typename std::list<Point>::iterator ppi = points.begin();
	ppi != points.end();ppi++)
      dt.insert(*ppi);
  }
    

  // pole extraction
  void extract_poles()
  {
    if (points.empty())
      return;
    Vertex_iterator vi;
    for (vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
      {
	Face_circulator fc = dt.incident_faces(vi);
	Face_circulator c = fc;
	std::list<Point> vv;
			    
	if (fc != NULL)
	  {
	    do
	      {
		if (!dt.is_infinite(c) && dt.triangle(c).area() != 0)
		  vv.push_back(dt.dual(c));
	      } while (++c != fc);
	  }
			    
	// find the farthest voronoi vertex from this point
	typename std::list<Point>::iterator maxp = vv.begin();
	for (typename std::list<Point>::iterator pi = vv.begin(); pi != vv.end(); pi++)
	  {
	    if (squared_distance(*pi, (*vi).point()) > squared_distance(*maxp, (*vi).point()))
	      {
		maxp = pi;
	      }
	  }
	poles.push_back(*maxp);
			    
	// find the farthest voronoi vertex from this point in the other half-plane
	typename std::list<Point>::iterator maxp2 = vv.begin();
	for (typename std::list<Point>::iterator pi = vv.begin(); pi != vv.end(); pi++)
	  {
	    if (angle((Point) *pi, (Point) (*vi).point(), (Point) *maxp) == OBTUSE &&
		squared_distance(*pi, (*vi).point()) > squared_distance(*maxp2, (*vi).point()))
	      {
		maxp2 = pi;
	      }
	  }
	poles.push_back(*maxp2);
      }
  }

  void generate_sites()
  {
    if(poles.empty())
      return;

    Search_tree tree_poles(poles.begin(), poles.end());
    Search_tree tree_points(points.begin(), points.end());

    for(typename std::list<Point>::iterator pi = points.begin();
	pi != points.end();
	pi++)
      {

	// estimate lfs (distance to medial axis estimate)
	Neighbor_search search_poles(tree_poles, *pi, 1);
	double lfs = std::sqrt((search_poles.begin())->second);

	// compute distance to nearest input point
	Neighbor_search search_points(tree_points, *pi, 2); // notice 2
	typename Neighbor_search::iterator it = search_points.begin();
	it++; // the first one is...itself
	double d = std::sqrt(it->second);
					
	// take min(lfs,distance to nearest point) as sizing
	double sizing = (std::min)(lfs,d);

	sites.push_back(Site(*pi, -sizing / K));
      }
  }

  // generate apollonius graph
  // first sorting the sites in decreasing order of weights
  void generate_apollonius()
  {
    if(sites.empty())
      return;
    sites.sort(Compare_site_2());
    ag.insert(sites.begin(),sites.end());
  }
			    
  double weighted_distance(const Point p,
			   const Site s) const
  {
    return std::sqrt(squared_distance(p, s.point())) - s.weight();
  }

};

}//namespace CGAL

#endif //CGAL_LIPSCHITZ_SIZING_FIELD_2_H
