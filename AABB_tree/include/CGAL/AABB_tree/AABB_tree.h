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
// Author(s)     :  Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_AABB_TREE_H
#define CGAL_AABB_TREE_H

#include <vector>
#include <list>
#include <stack>
#include <CGAL/AABB_tree/AABB_node.h>
#include <boost/mpl/vector.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/contains.hpp>

namespace CGAL {

template <class Kernel, class Input, class PSC>
class AABB_tree
{
public:

  // basic kernel object types
  typedef typename Kernel::FT FT;
  typedef typename CGAL::Bbox_3 Bbox;
  typedef typename Kernel::Ray_3 Ray;
  typedef typename Kernel::Line_3 Line;
  typedef typename Kernel::Plane_3 Plane;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid;

  // Nodes of the tree
  typedef AABB_node<Kernel,Input,PSC> Node;
  typedef typename Node::Point_with_input Point_with_input;

private:

  // set of input primitives (halfedge or face handles)
  std::vector<Input> m_data;

  // single root node
  Node *m_root;

public:
  // life cycle
  AABB_tree() : m_root(NULL) {}
  ~AABB_tree()
  {
    cleanup();
  }

  void cleanup()
  {
    m_data.clear();
    delete [] m_root; m_root = NULL;
  }

  // build tree when input = face_handle
  bool build_faces(PSC& psc)
  {
    cleanup();
    set_face_data(psc);
    if(!empty())
    {
      m_root = new Node[m_data.size()-1]();
      m_root->expand(m_data.begin(),m_data.end(), m_data.size());
      return true;
    }
    return false;
  }

private:

  void set_face_data(PSC& psc)
  {
    unsigned int nbf = psc.size_of_facets();
    m_data.reserve(nbf);
    typename PSC::Facet_iterator f;
    for(f = psc.facets_begin(); f != psc.facets_end(); f++)
      m_data.push_back(f);
  }

public:

  bool empty()
  {
    return m_data.size() < 2; // TODO: change this requirement to < 1
  }

  // --------------------RAY/SEGMENT ORACLES----------------------//

  typedef boost::mpl::vector<Ray, Line, Segment> Allowed_query_types;

  // The following function template is restricted to that T can only be in
  // {Ray, Line, Segment}. It return type is bool.
  // The trick uses enable_if and the Boost MPL.
  template <class T>
  typename boost::enable_if<
    typename boost::mpl::contains<Allowed_query_types,
                                  T>::type,
    bool>::type
  first_intersection(const T& x,
                     Point_with_input& pwh)
  {
    std::pair<bool,Point_with_input> result;
    m_root->first_intersection(x, result, m_data.size());
    if(result.first)
    {
      pwh = result.second;
      return true;
    }
    return false;
  }

}; // end class AABB_tree

} // end namespace CGAL

#endif // CGAL_AABB_TREE_H
