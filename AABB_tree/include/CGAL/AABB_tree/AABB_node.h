// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
// Author(s)     : 

#ifndef CGAL_AABB_NODE_H
#define CGAL_AABB_NODE_H

#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <vector>

#include "Ray_3_Bbox_3_do_intersect.h"
#include "Bbox_3_Bbox_3_do_intersect.h"
#include "Segment_3_Bbox_3_do_intersect.h"
#include "Plane_3_Bbox_3_do_intersect.h"
#include "Triangle_3_Bbox_3_do_intersect.h"
#include "Line_3_Bbox_3_do_intersect.h"
#include "Plucker_ray_3_Bbox_3_do_intersect.h"
#include "Sphere_3_Bbox_do_intersect.h"

CGAL_BEGIN_NAMESPACE

template <class Kernel, class Input, class PSC>
class AABB_node
{
public:

  // type with fixed (double) floating point arithmetic
  typedef CGAL::Bbox_3 Bbox;

  // basic kernel object types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Ray_3 Ray;
  typedef typename Kernel::Line_3 Line;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Plane_3 Plane;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Sphere_3 Sphere;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::Triangle_3 Triangle;

  typedef AABB_node<Kernel,Input,PSC> Node;
  typedef typename std::vector<Input>::iterator Iterator;
  typedef std::pair<Point, Input> Point_with_input;

  typedef typename PSC::Traits PSC_kernel;
  typedef typename PSC_kernel::Point_3 PSC_point;
  typedef typename PSC_kernel::Vector_3 PSC_vector;
  typedef typename PSC_kernel::Triangle_3 PSC_triangle;
  typedef CGAL::Cartesian_converter<PSC_kernel, Kernel > Converter;

private:

  // bounding box
  Bbox m_bbox;

  // children nodes
  // either pointing towards children (if the children are not leaves)
  // or pointing toward Input primitives (if the children are leaves)
  void *m_left_child;
  void *m_right_child;

public:

  // life cycle
  AABB_node()
  {
    m_left_child = m_right_child = NULL;
  }

  ~AABB_node() {}

public:

  // deletes the whole subtree rooted at this node (except this node, of course).
  // node = number of primitives contained in this node.
  void cleanup_recurse(int nb_primitives)
  {
    switch(nb_primitives)
    {
    case 2:
      break;
    case 3:
      delete static_cast<Node*>(m_right_child);
      break;
    default:
      static_cast<Node*>(m_left_child)->cleanup_recurse(nb_primitives/2);
      static_cast<Node*>(m_right_child)->cleanup_recurse(nb_primitives - nb_primitives/2);
      delete static_cast<Node*>(m_left_child);
      delete static_cast<Node*>(m_right_child);
    }
  }

private:

  Bbox compute_bbox(Input f)
  {
    const PSC_point a = f->halfedge()->vertex()->point();
    const PSC_point b = f->halfedge()->next()->vertex()->point();
    const PSC_point c = f->halfedge()->next()->next()->vertex()->point();
    return a.bbox() + b.bbox() + c.bbox();
  }

  // compute bbox for iterator range of input primitives
  Bbox bbox(const PSC& psc, Iterator a, Iterator b)
  {
    Bbox bbox = compute_bbox(*a);
    for(++a; a != b; ++a)
      bbox = bbox + compute_bbox(*a);
    return bbox;
  }

public:
  // builds the tree by recursive expansion.
  // [a,b[ is the range of primitives to be added to the tree.
  // node is the length of this range.
  void expand(const PSC& psc, Iterator a, Iterator b, int range)
  {
    m_bbox = bbox(psc, a, b);

    // sort primitives along longest axis aabb
    sort_primitives(a, b);

    switch(range)
    {
    case 2:
      m_left_child = &(*a);
      m_right_child = &(*(++a));
      break;
    case 3:
      m_left_child = &(*a);
      m_right_child = static_cast<Node*>(this)+1;
      static_cast<Node*>(m_right_child)->expand(psc, a+1, b, 2);	
      break;
    default:
      m_left_child = static_cast<Node*>(this)+1;
      m_right_child = (static_cast<Node*>(this))+(range/2);
      static_cast<Node*>(m_left_child)->expand(psc, a, a+range/2, range/2);
      static_cast<Node*>(m_right_child)->expand(psc, a+range/2, b, range - range/2);				
    }
  }

private:

  static Point centroid(Input f)
  {
    const PSC_point a = f->halfedge()->vertex()->point();
    const PSC_point b = f->halfedge()->next()->vertex()->point();
    const PSC_point c = f->halfedge()->next()->next()->vertex()->point();
    // somehow CGAL::centroid does not compile
    PSC_vector u = a - CGAL::ORIGIN;
    PSC_vector v = b - CGAL::ORIGIN;
    PSC_vector w = c - CGAL::ORIGIN;
    Converter convert;
    return convert(CGAL::ORIGIN + (u + v + w) / 3.0);
  }

  template<typename Input_>
  static bool lower_x(const Input_& i1, const Input_& i2)
  {
    return lower_x(i1,i2);
  }

  template<typename Input_>
  static bool lower_y(const Input_& i1, const Input_& i2)
  {
    return lower_y(i1,i2);
  }

  template<typename Input_>
  static bool lower_z(const Input_& i1, const Input_& i2)
  {
    return lower_z(i1,i2);
  }

  static bool lower_x(const Input& f1, const Input& f2)
  {
    return centroid(f1).x() < centroid(f2).x();
  }
  static bool lower_y(const Input& f1, const Input& f2)
  {
    return centroid(f1).y() < centroid(f2).y();
  }
  static bool lower_z(const Input& f1, const Input& f2)
  {
    return centroid(f1).z() < centroid(f2).z();
  }

  void sort_primitives(Iterator a, Iterator b)
  {
    switch(longest_axis())
    {
    case 0: // sort along x
      std::sort(a,b,lower_x<Input>);
      break;
    case 1: // sort along y
      std::sort(a,b,lower_y<Input>);
      break;
    default: // sort along z
      std::sort(a,b,lower_z<Input>);
    }
  }

  int longest_axis()
  {
    FT max_size = std::max(xsize(),std::max(ysize(),zsize()));
    if(max_size == xsize())
      return 0; // axis along x
    if(max_size == ysize())
      return 1; // axis along y
    return 2; // axis along z
  }

  // size of bounding box along each axis
  FT xsize() { return m_bbox.xmax() - m_bbox.xmin(); }
  FT ysize() { return m_bbox.ymax() - m_bbox.ymin(); }
  FT zsize() { return m_bbox.zmax() - m_bbox.zmin(); }

public:


  // compute the first intersection encountered.
  // nb_primitives = number of primitives contained in this node.
  template<class QueryType, class ResultType>
  void first_intersection(const QueryType& q,
    ResultType& result, 
    int nb_primitives)
  {
    traversal<First_intersection_traits<QueryType, ResultType> >(q, result, nb_primitives);
  }

  template<class QueryType, class ResultType>
  class First_intersection_traits
  {
  private:
    ResultType& r;
  public:
    bool go_further()
    {
      return !r.first;
    }
    First_intersection_traits(ResultType& result) : r(result) {}

    bool intersection(const QueryType& q, const Input& i)
    {
      ResultType result;
      if(Node::intersection(q, i, result.second))
      {
	r.first = true;
	r.second = result.second;
	return true;
      }
      return false;
    }
    bool do_intersect(const QueryType& q, const Node& node)
    {
      return Node::do_intersect(q, node);
    }
  };

  // -----------------------------------------------------------//
  // -----------------------LINE ORACLES-------------------------//
  // -----------------------------------------------------------//

  static bool intersection(const Line& line,
    Input f,
    Point& result)
  {
    Triangle t = triangle(f);
    if(CGAL::do_intersect(t,line)) 
    {
      CGAL::Object inter = CGAL::intersection(t.supporting_plane(),line);
      if(CGAL::assign(result, inter))
	return true;
    }
    return false;
  }

  static bool intersection(const Line& line, 
    const typename PSC::Facet_handle& f, 
    Point_with_input& p)
  {
    Point p_alone;
    if(Node::intersection(line, f, p_alone))
    {
      p = Point_with_input(p_alone, f);
      return true;
    }
    return false;
  }

  static bool do_intersect(const Line& line, 
    const typename PSC::Facet_handle& f)
  {
    return PSC::do_intersect(line, f);
  }

  static bool do_intersect(const Line& line,
    const Node& node)
  {
    return CGAL::do_intersect(line, node.m_bbox);
  }



  // -----------------------------------------------------------//
  // -----------------------SEGMENT ORACLES-------------------------//
  // -----------------------------------------------------------//

  static bool intersection(const Segment& segment,
    Input f,
    Point& result)
  {
    Triangle t = triangle(f);
    if(CGAL::do_intersect(t,segment)) 
    {
      CGAL::Object inter = CGAL::intersection(t.supporting_plane(),segment);
      if(CGAL::assign(result, inter))
	return true;
    }
    return false;
  }

  static bool intersection(const Segment& segment, 
    const typename PSC::Facet_handle& f, 
    Point_with_input& p)
  {
    Point p_alone;
    if(Node::intersection(segment, f, p_alone))
    {
      p = Point_with_input(p_alone, f);
      return true;
    }
    return false;
  }

  static Triangle triangle(Input f)
  {
    Converter convert;
    const PSC_point a = f->halfedge()->vertex()->point();
    const PSC_point b = f->halfedge()->next()->vertex()->point();
    const PSC_point c = f->halfedge()->next()->next()->vertex()->point();
    return convert(PSC_triangle(a,b,c));
  }

  static bool do_intersect(const Segment& segment, Input f)
  {
    return CGAL::do_intersect(triangle(f),segment);
  }

  static bool do_intersect(const Segment& segment,
    const Node& node)
  {
    return CGAL::do_intersect(segment, node.m_bbox);
  }

  // -----------------------------------------------------------//
  // -----------------------RAY ORACLES-------------------------//
  // -----------------------------------------------------------//

  static bool intersection(const Ray& ray,
    Input f,
    Point& result)
  {
    Triangle t = triangle(f);
    if(CGAL::do_intersect(t,ray)) 
    {
      CGAL::Object inter = CGAL::intersection(t.supporting_plane(),ray);
      if(CGAL::assign(result, inter))
	return true;
    }
    return false;
  }

  static bool intersection(const Ray& ray, 
    const typename PSC::Facet_handle& f, 
    Point_with_input& p)
  {
    Point p_alone;
    if(Node::intersection(ray, f, p_alone))
    {
      p = Point_with_input(p_alone, f);
      return true;
    }
    return false;
  }

  static bool do_intersect(const Ray& ray, 
    const typename PSC::Facet_handle& f)
  {
    return PSC::do_intersect(ray, f);
  }

  static bool do_intersect(const Ray& ray,
    const Node& node)
  {
    return CGAL::do_intersect(ray, node.m_bbox);
  }



  // -----------------------------------------------------------//
  // ----------------------GENERAL QUERY------------------------//
  // -----------------------------------------------------------//

  // general traversal query, the traits class allows to use it for the various 
  // traversal methods we need: listing, counting, detecting intersections, drawing the boxes.

  template<class Traits, class QueryType, class ResultType>
  void traversal(const QueryType& query,
    ResultType& result,
    int nb_primitives)
  {
    Traits traits(result);
    bool left_test;
    switch(nb_primitives)
    {
    case 2:
      left_test = traits.intersection(query, *static_cast<Input*>(m_left_child));
      if((! left_test) || (traits.go_further()))
	traits.intersection(query, *static_cast<Input*>(m_right_child));
      break;
    case 3:
      left_test = traits.intersection(query, *static_cast<Input*>(m_left_child));
      if((! left_test) || (traits.go_further()))
	if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
	  static_cast<Node*>(m_right_child)->traversal<Traits>(query, result, 2);
      break;
    default:
      if(traits.do_intersect(query, *static_cast<Node*>(m_left_child)))
      {
	static_cast<Node*>(m_left_child)->traversal<Traits>(query, result, nb_primitives/2);
	if(traits.go_further())
	  if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
	    static_cast<Node*>(m_right_child)->traversal<Traits>(query, result, nb_primitives - nb_primitives/2);
      }
      else
	if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
	  static_cast<Node*>(m_right_child)->traversal<Traits>(query, result, nb_primitives - nb_primitives/2);
    }
  }


};

CGAL_END_NAMESPACE

#endif // CGAL_AABB_NODE_

