// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETHZ (Suisse).
// Copyrigth (c) 2009  GeometryFactory (France)
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
// Author(s)     :  Camille Wormser, Jane Tournois, Pierre Alliez, Laurent Rineau

#ifndef CGAL_AABB_NODE_H
#define CGAL_AABB_NODE_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <vector>

#include <boost/mpl/vector.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/contains.hpp>

#include <CGAL/AABB_tree/Ray_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Bbox_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Segment_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Plane_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Triangle_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Line_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Sphere_3_Bbox_do_intersect.h>
#include <CGAL/AABB_tree/Triangle_3_segment_3_intersection.h>
#include <CGAL/AABB_tree/Triangle_3_ray_3_intersection.h>

namespace CGAL {

template <class Kernel, class Input, class PSC>
class AABB_node
{
public:

  // Bounding box with fixed (double) floating point arithmetic
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
  typedef typename PSC_kernel::FT PSC_FT;
  typedef typename PSC_kernel::Point_3 PSC_point;
  typedef typename PSC_kernel::Vector_3 PSC_vector;
  typedef typename PSC_kernel::Triangle_3 PSC_triangle;

  // Converter of number types, points and vectors
  typedef CGAL::Cartesian_converter<PSC_kernel, Kernel > Converter;

private:

  // bounding box
  Bbox m_bbox;

  // children nodes:
  // either pointing towards children (if the children are not leaves)
  // or pointing toward Input primitives (if the children are leaves).
  void *m_left_child;
  void *m_right_child;

public:

  // life cycle
  AABB_node()
  {
    m_left_child = m_right_child = NULL;
  }

  ~AABB_node() {}

  Bbox bbox() const { return m_bbox; }
  // -----------------------------------------------------------//
  // -----------------RECURSIVE MEMBER FUNCTIONS----------------//
  // -----------------------------------------------------------//

  // builds the tree by recursive expansion.
  // [a,b[ is the range of primitives to be added to the tree.
  // 'range' is the length of this range.
  void expand(Iterator a, Iterator b, int range)
  {
    m_bbox = bbox(a, b);

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
      static_cast<Node*>(m_right_child)->expand(a+1, b, 2);	
      break;
    default:
      m_left_child = static_cast<Node*>(this)+1;
      m_right_child = (static_cast<Node*>(this))+(range/2);
      static_cast<Node*>(m_left_child)->expand(a, a+range/2, range/2);
      static_cast<Node*>(m_right_child)->expand(a+range/2, b, range - range/2);				
    }
  }

  // deletes the whole subtree rooted at this node (except this node, of course).
  // nb_primitives = number of primitives contained in this node.
  void cleanup_recurse(int nb_primitives)
  {
    switch(nb_primitives)
    {
    case 2:
      break;
    case 3:
      delete static_cast<Node*>(m_right_child); m_right_child = NULL;
      break;
    default:
      static_cast<Node*>(m_left_child)->cleanup_recurse(nb_primitives/2);
      static_cast<Node*>(m_right_child)->cleanup_recurse(nb_primitives - nb_primitives/2);
      delete static_cast<Node*>(m_left_child); m_left_child = NULL;
      delete static_cast<Node*>(m_right_child); m_right_child = NULL;
    }
  }

  // general traversal query, the traits class allows to use it for the various 
  // traversal methods we need: listing, counting, detecting intersections, drawing the boxes.

  template<class Traits, class QueryType>
  void traversal(const QueryType& query,
    Traits& traits,
    const int nb_primitives) const
  {
    switch(nb_primitives)
    {
    case 2:
      {
        const bool left_test =
          traits.intersection(query, *static_cast<Input*>(m_left_child));
        if((! left_test) || (traits.go_further()))
          traits.intersection(query, *static_cast<Input*>(m_right_child));
      }
      break;
    case 3:
      {
        const bool left_test =
          traits.intersection(query, *static_cast<Input*>(m_left_child));
        if((! left_test) || (traits.go_further()))
          if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
            static_cast<Node*>(m_right_child)->traversal<Traits>(query, traits, 2);
        break;
      }
    default:
      if(traits.do_intersect(query, *static_cast<Node*>(m_left_child)))
      {
	static_cast<Node*>(m_left_child)->traversal<Traits>(query, traits, nb_primitives/2);
	if(traits.go_further())
	  if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
	    static_cast<Node*>(m_right_child)->traversal<Traits>(query, traits, nb_primitives - nb_primitives/2);
      }
      else
	if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
	  static_cast<Node*>(m_right_child)->traversal<Traits>(query, traits, nb_primitives - nb_primitives/2);
    }
  }


private:

  // DEPENDENCIES OF 'expand'

  // compute bbox for one input primitive
  Bbox compute_bbox(Input f)
  {
    CGAL_PROFILER("Call to AABB_node::compute_bbox(Input)");
    const PSC_point a = f->halfedge()->vertex()->point();
    const PSC_point b = f->halfedge()->next()->vertex()->point();
    const PSC_point c = f->halfedge()->next()->next()->vertex()->point();
    return a.bbox() + b.bbox() + c.bbox();
  }

  // compute bbox for iterator range of input primitives
  Bbox bbox(Iterator a, Iterator b)
  {
    Bbox bbox = compute_bbox(*a);
    for(++a; a != b; ++a)
      bbox = bbox + compute_bbox(*a);
    return bbox;
  }

  static Point centroid(Input f)
  {
    CGAL_PROFILER("Call to AABB_node::centroid(Input)");
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
    CGAL_PROFILER("Call to AABB_node::less_x(Input, Input)");
    const PSC_FT& ax1 = f1->halfedge()->vertex()->point().x();
    const PSC_FT& bx1 = f1->halfedge()->next()->vertex()->point().x();
    const PSC_FT& cx1 = f1->halfedge()->next()->next()->vertex()->point().x();

    const PSC_FT& ax2 = f2->halfedge()->vertex()->point().x();
    const PSC_FT& bx2 = f2->halfedge()->next()->vertex()->point().x();
    const PSC_FT& cx2 = f2->halfedge()->next()->next()->vertex()->point().x();

    return (ax1+bx1+cx1) < (ax2+bx2+cx2);
  }
  static bool lower_y(const Input& f1, const Input& f2)
  {
    CGAL_PROFILER("Call to AABB_node::less_y(Input, Input)");
    const PSC_FT& ay1 = f1->halfedge()->vertex()->point().y();
    const PSC_FT& by1 = f1->halfedge()->next()->vertex()->point().y();
    const PSC_FT& cy1 = f1->halfedge()->next()->next()->vertex()->point().y();

    const PSC_FT& ay2 = f2->halfedge()->vertex()->point().y();
    const PSC_FT& by2 = f2->halfedge()->next()->vertex()->point().y();
    const PSC_FT& cy2 = f2->halfedge()->next()->next()->vertex()->point().y();

    return (ay1+by1+cy1) < (ay2+by2+cy2);
  }
  static bool lower_z(const Input& f1, const Input& f2)
  {
    CGAL_PROFILER("Call to AABB_node::less_z(Input, Input)");
    const PSC_FT& az1 = f1->halfedge()->vertex()->point().z();
    const PSC_FT& bz1 = f1->halfedge()->next()->vertex()->point().z();
    const PSC_FT& cz1 = f1->halfedge()->next()->next()->vertex()->point().z();

    const PSC_FT& az2 = f2->halfedge()->vertex()->point().z();
    const PSC_FT& bz2 = f2->halfedge()->next()->vertex()->point().z();
    const PSC_FT& cz2 = f2->halfedge()->next()->next()->vertex()->point().z();

    return (az1+bz1+cz1) < (az2+bz2+cz2);
  }

  enum Axis { CGAL_AXIS_X = 0, CGAL_AXIS_Y, CGAL_AXIS_Z};

  void sort_primitives(Iterator a, Iterator b)
  {
    Iterator middle = a + (b - a)/2;
    switch(longest_axis())
    {
    case CGAL_AXIS_X: // sort along x
      std::nth_element(a,middle,b,lower_x<Input>);
      break;
    case CGAL_AXIS_Y: // sort along y
      std::nth_element(a,middle,b,lower_y<Input>);
      break;
    default: // sort along z
      std::nth_element(a,middle,b,lower_z<Input>);
    }
  }

  // Return the node's longest axis as 0 (X), 1 (Y) or 2 (Z)
  Axis longest_axis()
  {
    const double dx = xsize();
    const double dy = xsize();
    const double dz = xsize();
    if(dx>=dy) {
      if(dx>=dz) {
        return CGAL_AXIS_X; // axis along x
      }
      else {//  dz>dx and dx>=dy
        return CGAL_AXIS_Z;
      }
    }
    else { // dy>dx
      if(dy>=dz) {
        return CGAL_AXIS_Y;
      }
      else { // dz>dy and dy>dx 
        return CGAL_AXIS_Z;
      }
    }
  }

  // size of bounding box along each axis
  double xsize() const { return m_bbox.xmax() - m_bbox.xmin(); }
  double ysize() const { return m_bbox.ymax() - m_bbox.ymin(); }
  double zsize() const { return m_bbox.zmax() - m_bbox.zmin(); }

public:
  double max_length() const { 
    return (std::max)(xsize(),(std::max)(ysize(),zsize()));
  }

  // HELPER FUNCTION of all predicates below
  static Triangle triangle(Input f)
  {
    Converter convert;
    const PSC_point a = f->halfedge()->vertex()->point();
    const PSC_point b = f->halfedge()->next()->vertex()->point();
    const PSC_point c = f->halfedge()->next()->next()->vertex()->point();
    return convert(PSC_triangle(a,b,c));
  }

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
    const Input& f, 
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
    const Node& node)
  {
    return CGAL::do_intersect(line, node.m_bbox);
  }

  // -----------------------------------------------------------//
  // -----------------------GENERIC ORACLES---------------------//
  // -----------------------------------------------------------//
  // In the following, Query can be any of Plane, Ray, Line, or Segment.
  
  typedef boost::mpl::vector<Plane, Ray, Line, Segment> Allowed_query_types;

  // The following function template is restricted to that Query can only be in
  // {Ray, Line, Segment}. It return type is bool.
  // The trick uses enable_if and the Boost MPL.
  template <class Query, class Result>
  static
  typename boost::enable_if<
    typename boost::mpl::contains<Allowed_query_types,
                                  Query>::type,
    bool>::type
  intersection(const Query& query,
               Input f,
               Result& result)
  {
    Triangle t = triangle(f);
    if(CGAL::do_intersect(t, query)) 
    {
      CGAL::Object inter = CGAL::intersection(t, query);
      if(CGAL::assign(result, inter))
	return true;
    }
    return false;
  }

  // The following function template is restricted to that Query can only be in
  // {Ray, Line, Segment}. It return type is bool.
  // The trick uses enable_if and the Boost MPL.
  template <class Query>
  static
  typename boost::enable_if<
    typename boost::mpl::contains<Allowed_query_types,
                                  Query>::type,
    bool>::type
  intersection(const Query& query, 
    const Input& f, 
    Point_with_input& p)
  {
    Point p_alone;
    if(Node::intersection(query, f, p_alone))
    {
      p = Point_with_input(p_alone, f);
      return true;
    }
    return false;
  }

  // The following function template is restricted to that Query can only be in
  // {Ray, Line, Segment}. It return type is bool.
  // The trick uses enable_if and the Boost MPL.
  template <class Query>
  static
  typename boost::enable_if<
    typename boost::mpl::contains<Allowed_query_types,
                                  Query>::type,
    bool>::type
  do_intersect(const Query& query, Input f)
  {
    return CGAL::do_intersect(triangle(f), query);
  }

  // The following function template is restricted to that Query can only be in
  // {Ray, Line, Segment}. It return type is bool.
  // The trick uses enable_if and the Boost MPL.
  template <class Query>
  static
  typename boost::enable_if<
    typename boost::mpl::contains<Allowed_query_types,
                                  Query>::type,
    bool>::type
  do_intersect(const Query& query,
    const Node& node)
  {
    return CGAL::do_intersect(query, node.m_bbox);
  }

}; // end class template AABB_node

} // end namespace CGAL

#endif // CGAL_AABB_NODE_H

