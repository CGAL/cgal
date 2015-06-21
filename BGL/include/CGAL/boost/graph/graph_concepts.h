// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Philipp Moeller

#ifndef CGAL_GRAPH_CONCEPTS_H
#define CGAL_GRAPH_CONCEPTS_H

#include <boost/graph/graph_concepts.hpp>
#include <boost/concept/detail/concept_def.hpp>

namespace CGAL {
namespace concepts {

BOOST_concept(HalfedgeGraph,(G))
  : boost::concepts::Graph<G>
{
  typedef typename boost::graph_traits<G>::halfedge_descriptor             halfedge_descriptor;

  BOOST_CONCEPT_USAGE(HalfedgeGraph)
  {
    BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<halfedge_descriptor>));
    BOOST_CONCEPT_ASSERT((boost::EqualityComparable<halfedge_descriptor>));
    BOOST_CONCEPT_ASSERT((boost::Assignable<halfedge_descriptor>));


    e = edge(h, g);
    h = halfedge(e, g);
    h = halfedge(v, g);
    hp = halfedge(u, v, g);
    h = opposite(h, g);
    v = source(h, g);
    v = target(h, g);
    h = next(h, g);
    h = prev(h, g);
    const_constraints(g);
  }

  void const_constraints(const G& cg)
  {
    e = edge(h, cg);
    h = halfedge(e, cg);
    h = halfedge(v, cg);
    hp = halfedge(u, v, cg);
    h = opposite(h, cg);
    v = source(h, cg);
    v = target(h, cg);
    h = next(h, cg);
    h = prev(h, cg);
  }
  
  G g;
  
  typename boost::graph_traits<G>::vertex_descriptor v, u;
  typename boost::graph_traits<G>::edge_descriptor e;
  typename boost::graph_traits<G>::halfedge_descriptor h;
  std::pair<halfedge_descriptor, bool> hp;
};

BOOST_concept(HalfedgeListGraph,(G))
  : HalfedgeGraph<G>
{
  typedef typename boost::graph_traits<G>::halfedge_iterator   halfedge_iterator;
  typedef typename boost::graph_traits<G>::halfedges_size_type halfedges_size_type;
  
  BOOST_CONCEPT_USAGE(HalfedgeListGraph)
  {
    // BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<halfedge_iterator>));
    h_num = num_halfedges(g);
    p = halfedges(g);
    this->h = *p.first;
    const_constraints(g);
  }

  void const_constraints(const G& cg)
  {
    h_num = num_halfedges(cg);
    p = halfedges(cg);
    this->h = *p.first;
  }
  
  G g;
  halfedges_size_type h_num;
  std::pair<halfedge_iterator, halfedge_iterator> p;
};

BOOST_concept(FaceGraph,(G))
  : HalfedgeGraph<G>
{
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;

  BOOST_CONCEPT_USAGE(FaceGraph)
  {
    BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<face_descriptor>));
    BOOST_CONCEPT_ASSERT((boost::EqualityComparable<face_descriptor>));
    BOOST_CONCEPT_ASSERT((boost::Assignable<face_descriptor>));

    f = face(h, g);
    h = halfedge(f, g);
    boost::graph_traits<G>::null_face();
    const_constraints(g);
  }

  void const_constraints(const G& cg)
  {
    f = face(h, cg);
    h = halfedge(f, cg);
  }

  G g;
  face_descriptor f;
  typename boost::graph_traits<G>::halfedge_descriptor h;
};

BOOST_concept(FaceListGraph,(G))
  : FaceGraph<G>
{
  typedef typename boost::graph_traits<G>::face_iterator face_iterator;
  typedef typename boost::graph_traits<G>::faces_size_type faces_size_type;

  BOOST_CONCEPT_USAGE(FaceListGraph)
  {
    // BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<face_iterator>));
    p = faces(g);
    nf = num_faces(g);
    this->f = *p.first;
    const_constraints(g);
  }

  void const_constraints(const G& cg)
  {
    p = faces(cg);
    nf = num_faces(cg);
    this->f = *p.first;
  }

  G g;
  std::pair<face_iterator, face_iterator> p;
  typename boost::graph_traits<G>::faces_size_type nf;
};

BOOST_concept(MutableFaceGraph,(G))
  : FaceGraph<G>
{
  BOOST_CONCEPT_USAGE(MutableFaceGraph)
  {
    f = add_face(g);
    remove_face(f, g);
    set_face(h, f, g);
    set_halfedge(f, h, g);
  }
  G g;
  typename boost::graph_traits<G>::face_descriptor f;
  typename boost::graph_traits<G>::halfedge_descriptor h;
};

BOOST_concept(MutableHalfedgeGraph,(G))
  : HalfedgeGraph<G>
{
  BOOST_CONCEPT_USAGE(MutableHalfedgeGraph)
  {
    v = add_vertex(g);
    remove_vertex(v, g);
    e = add_edge(g);
    remove_edge(e, g);
    set_target(h1, v, g);
    set_next(h1, h2, g);
    set_halfedge(v, h1, g);
  }

  G g;
  typename boost::graph_traits<G>::edge_descriptor e;
  typename boost::graph_traits<G>::vertex_descriptor v;
  typename boost::graph_traits<G>::halfedge_descriptor h1, h2;
};

} // concepts

using CGAL::concepts::HalfedgeGraphConcept;
using CGAL::concepts::HalfedgeListGraphConcept;
using CGAL::concepts::FaceGraphConcept;
using CGAL::concepts::FaceListGraphConcept;
using CGAL::concepts::MutableFaceGraphConcept;
using CGAL::concepts::MutableHalfedgeGraphConcept;
} // CGAL

#include <boost/concept/detail/concept_undef.hpp>

#endif /* CGAL_GRAPH_CONCEPTS_H */
