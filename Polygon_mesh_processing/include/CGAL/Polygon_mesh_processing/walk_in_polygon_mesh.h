// Copyright (c) 2014 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_WALK_IN_POLYGON_MESH_H
#define CGAL_WALK_IN_POLYGON_MESH_H

#include <CGAL/assertions.h>
#include <CGAL/Polyhedron_simplex_type.h>

#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <CGAL/boost/graph/helpers.h>

#include <set>
#include <vector>
#include <stdexcept>

#ifdef CGAL_DEBUG_WALK_IN_POLYGON_MESH
#include <iostream>
#define CGAL_WALKER_TRACE(X) std::cout << X;
#else
#define CGAL_WALKER_TRACE(X)
#endif

namespace CGAL{

struct Walk_is_cycling_exception : std::logic_error
{
  Walk_is_cycling_exception()
    : std::logic_error("Loop detected during the walk")
  {}
};

namespace walker_impl{

  /// \todo I think this was used in the clipping-snapping branch
  struct Walk_on_polyhedron_border_crossing_exception : std::logic_error
  {
    Walk_on_polyhedron_border_crossing_exception()
      : std::logic_error("Trying to walking outside of the polyhedron")
    {}
  };

  /// get the set of faces incident to a given simplex
  template <class PM>
  void get_incident_faces(Polyhedron_simplex_type simplex_type,
                          typename boost::graph_traits<PM>::halfedge_descriptor h,
                          std::set<typename boost::graph_traits<PM>::face_descriptor>& incident_faces,
                          const PM& pm)
  {
    switch(simplex_type)
    {
      case POLYHEDRON_FACET:
        incident_faces.insert(face(h, pm));
      break;
      case POLYHEDRON_EDGE:
        if (!is_border(h, pm)) incident_faces.insert(face(h, pm));
        if (!is_border(opposite(h, pm), pm)) incident_faces.insert(face(opposite(h, pm), pm));
      break;
      case POLYHEDRON_VERTEX:
      {
        BOOST_FOREACH(typename boost::graph_traits<PM>::halfedge_descriptor hd,
                      halfedges_around_target(h, pm))
        {
          if (!is_border(hd, pm))
            incident_faces.insert(face(hd, pm));
        }
      }
      break;
      default:
      break;
    }
  }

  template <class PM>
  typename boost::graph_traits<PM>::halfedge_descriptor
  canonical_hedge(typename boost::graph_traits<PM>::halfedge_descriptor h,
                 const PM& pm)
  {
    return h<opposite(h, pm)?h:opposite(h, pm);
  }

  template<class PM>
  typename boost::graph_traits<PM>::halfedge_descriptor
  common_edge( typename boost::graph_traits<PM>::face_descriptor f1,
               typename boost::graph_traits<PM>::face_descriptor f2,
               const PM& pm)
  {
    CGAL_assertion(f1!=boost::graph_traits<PM>::null_face());
    CGAL_assertion(f2!=boost::graph_traits<PM>::null_face());
    BOOST_FOREACH(typename boost::graph_traits<PM>::halfedge_descriptor h,
                  halfedges_around_face(halfedge(f1, pm), pm))
    {
      if ( face(opposite(h, pm), pm)==f2 )
        return h;
    }
    return boost::graph_traits<PM>::null_halfedge();
  }
} // end of internal namespace

/// \todo document me
/// We walk from (src, source_type) to (tgt, target_type) intersecting
/// all edges encountered. The predicate inside_edge_pred is responsible to indicate the correct
/// direction to move.
/// \param source_type is different from POLYHEDRON_NONE
template< class PolygonMesh,
          class EdgeIntersectionPredicate,
          class EdgeIntersectionVisitor >
bool walk_in_polygon_mesh(const PolygonMesh& pm,
                          typename boost::graph_traits<PolygonMesh>::halfedge_descriptor src,
                          Polyhedron_simplex_type source_type,
                          typename boost::graph_traits<PolygonMesh>::halfedge_descriptor tgt,
                          Polyhedron_simplex_type target_type,
                          const EdgeIntersectionPredicate& inside_edge_pred,
                          EdgeIntersectionVisitor& visitor)
{
  typedef boost::graph_traits<PolygonMesh> graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;

  using namespace walker_impl;

  /// info returned by the predicate when an edge of the input triangle mesh
  /// is intersected in its interior
  typedef typename EdgeIntersectionPredicate::Barycentric_NT Barycentric_NT;

/// First check that the projected vertices do not fall in a common simplex
  halfedge_descriptor common_face = graph_traits::null_halfedge();
  bool new_edge_to_add = true;
  switch( source_type )
  {
    case POLYHEDRON_FACET:
      CGAL_WALKER_TRACE("\n1-POLYHEDRON_FACET\n")
      switch(target_type)
      {
        case POLYHEDRON_FACET:
          CGAL_WALKER_TRACE("2-POLYHEDRON_FACET\n")
          if (face(src, pm) == face(tgt, pm))
            common_face=src;
        break;
        case POLYHEDRON_EDGE:
          CGAL_WALKER_TRACE("2-POLYHEDRON_EDGE\n")
          if ( face(src, pm)==face(tgt, pm) ||
               face(src, pm)==face(opposite(tgt, pm), pm) )
            common_face=src;
        break;
        case POLYHEDRON_VERTEX:
        {
          CGAL_WALKER_TRACE("2-POLYHEDRON_VERTEX\n")
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(tgt, pm))
            if (face(hd, pm)==face(src, pm))
            {
              common_face=src;
              break;
            }
        }
        break;
        case POLYHEDRON_NONE:
          CGAL_WALKER_TRACE("2-POLYHEDRON_NONE\n")
        break;
      }
    break;
    case POLYHEDRON_EDGE:
      CGAL_WALKER_TRACE("\n1-POLYHEDRON_EDGE\n")
      switch(target_type)
      {
        case POLYHEDRON_FACET:
          CGAL_WALKER_TRACE("2-POLYHEDRON_FACET\n")
          if ( face(tgt, pm)==face(src, pm) ||
               face(tgt, pm)==face(opposite(src, pm), pm) )
            common_face=tgt;
        break;
        case POLYHEDRON_EDGE:
          CGAL_WALKER_TRACE("2-POLYHEDRON_EDGE\n")
          if (canonical_hedge(src, pm)==canonical_hedge(tgt, pm))
          {
            /// we explicitly do not update common_face
            /// so that the visitor knows src and tgt are completly
            /// inside an edge
            new_edge_to_add=false;
          }
          else
          {
            if ( !is_border(src, pm) &&
                  ( face(src, pm) == face(tgt, pm) ||
                    face(src, pm) == face(opposite(tgt, pm), pm) ) )
              common_face=src;
            else{
              if (!is_border(opposite(src, pm), pm)){
                if ( face(tgt, pm) == face(opposite(src, pm), pm) )
                  common_face=tgt;
                else
                  if ( face(opposite(tgt, pm), pm) == face(opposite(src, pm), pm) )
                    common_face=opposite(tgt, pm);
              }
            }
          }
        break;
        case POLYHEDRON_VERTEX:
          CGAL_WALKER_TRACE("2-POLYHEDRON_VERTEX\n")
          if ( target(src, pm)==target(tgt, pm) ||
               source(src, pm)==target(tgt, pm) )
          {
            common_face = target(src, pm)==target(tgt, pm)?src:opposite(src, pm);
            new_edge_to_add=false;
          }
          else
          {
            if (!is_border(src, pm) &&
                 target(next(src, pm), pm) == target(tgt, pm))
              common_face=src;
            else
              if (!is_border(opposite(src, pm), pm) &&
                   target(next(opposite(src, pm), pm), pm) == target(tgt, pm))
                common_face=opposite(src, pm);
          }
        break;
        case POLYHEDRON_NONE:
          CGAL_WALKER_TRACE("2-POLYHEDRON_NONE\n")
        break;
      }
    break;
    case POLYHEDRON_VERTEX:
      CGAL_WALKER_TRACE("\n1-POLYHEDRON_VERTEX\n")
      switch(target_type)
      {
        case POLYHEDRON_FACET:
        {
          CGAL_WALKER_TRACE("2-POLYHEDRON_FACET\n")
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(src, pm))
            if (face(hd, pm)==face(tgt, pm))
            {
              common_face=tgt;
              break;
            }
        }
        break;
        case POLYHEDRON_EDGE:
          CGAL_WALKER_TRACE("2-POLYHEDRON_EDGE\n")
          if ( target(tgt, pm)==target(src, pm) ||
               source(tgt, pm)==target(src, pm) )
          {
            common_face=target(tgt, pm)==target(src, pm)?
                         opposite(tgt, pm):tgt;
            new_edge_to_add=false;
          }
          else
          {
            if (!is_border(tgt, pm) &&
                 target(next(tgt, pm), pm) == target(src, pm))
              common_face=tgt;
            else
              if (!is_border(opposite(tgt, pm), pm) &&
                   target(next(opposite(tgt, pm), pm), pm) == target(src, pm))
                common_face=opposite(tgt, pm);
          }
        break;
        case POLYHEDRON_VERTEX:
        {
          /// @TODO Handle the case of when src and tgt are the same vertices
          // in case of identical vertices, it might indicate that we want to walk along a cycle.
          // If not, the predicate should return false for all candidates
          if (target(src, pm) == target(tgt, pm) ) break;

          CGAL_WALKER_TRACE("2-POLYHEDRON_VERTEX\n")
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(src, pm))
            if (source(hd, pm)==target(tgt, pm))
            {
              common_face=opposite(hd, pm);
              new_edge_to_add=false;
              break;
            }
        }
        break;
        case POLYHEDRON_NONE:
          CGAL_WALKER_TRACE("2-POLYHEDRON_NONE\n")
        break;
      }
    break;
    case POLYHEDRON_NONE:
    CGAL_assertion("Should never be here\n");
  }

  if (!new_edge_to_add)
  {
    // In this case, the edge to add is already on an edge of the input mesh
    visitor.on_walk_end_on_edge(common_face);
    return true;
  }
  if ( common_face != graph_traits::null_halfedge() )
  {
    CGAL_assertion( !is_border(common_face, pm) );
    visitor.on_walk_end(face(common_face, pm));
    return true;
  }


/// We now do the walk in the triangulation and find edge-edge intersection points
  std::set <face_descriptor> target_faces; /// \todo try using a boost::flat_set instead
  get_incident_faces(target_type, tgt, target_faces, pm);
  std::set <face_descriptor> faces_already_visited;
//  halfedge_descriptor first_halfedge_intersected;

  while( true )
  {
    std::vector<halfedge_descriptor> hedges_to_test;
    switch( source_type )
    {
      case POLYHEDRON_FACET:
        CGAL_WALKER_TRACE("   i-POLYHEDRON_FACET\n")
        hedges_to_test.push_back( src );
        hedges_to_test.push_back( next(src, pm) );
        hedges_to_test.push_back( prev(src, pm) );
      break;
      case POLYHEDRON_EDGE:
        CGAL_WALKER_TRACE("   i-POLYHEDRON_EDGE\n")
        // if the v_src is on an edge, we don't know in which direction to go
        // This is only true the first time of the loop
        if ( faces_already_visited.empty() &&
             !is_border(opposite(src, pm), pm) )
        {
          hedges_to_test.push_back( next(opposite(src, pm), pm) );
          hedges_to_test.push_back( prev(opposite(src, pm), pm) );
        }
        /// @TODO what if src is a border edge?
        hedges_to_test.push_back( next(src, pm) );
        hedges_to_test.push_back( prev(src, pm) );
      break;
      case POLYHEDRON_VERTEX:
      {
        CGAL_WALKER_TRACE("   i-POLYHEDRON_VERTEX\n")
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(src, pm))
        {
          if ( faces_already_visited.count(face(hd, pm))==0 //we don't want to go back, so we filter faces
               && !is_border(hd, pm)){
              hedges_to_test.push_back( prev(hd, pm) );
          }
        }
      }
      break;
      case POLYHEDRON_NONE:
        CGAL_assertion("Should never be here\n");
    }

    halfedge_descriptor intersected_edge = graph_traits::null_halfedge();
    boost::optional< std::pair<Barycentric_NT, boost::variant<bool, Barycentric_NT> > > inter_res;

    CGAL_WALKER_TRACE("   Looking for intersections hedges_to_test.size()="<< hedges_to_test.size() << "\n")
    std::size_t k=0;
    for (;k<hedges_to_test.size();++k)
    {
      CGAL_WALKER_TRACE("     k="<< k << ", testing "
                                << get(boost::vertex_point, pm, source(hedges_to_test[k], pm))
                                << " "
                                << get(boost::vertex_point, pm, target(hedges_to_test[k], pm))
                                << "\n")
      intersected_edge=hedges_to_test[k];
      inter_res =  inside_edge_pred( canonical_hedge(intersected_edge, pm) );
      if ( inter_res != boost::none){
        visitor.set_beta(inter_res->first);
        CGAL_WALKER_TRACE("     intersection found\n")
//        if (first_halfedge_intersected==graph_traits::null_halfedge())
//          first_halfedge_intersected=canonical_hedge(intersected_edge, pm);
//        else
//          if (first_halfedge_intersected==canonical_hedge(intersected_edge, pm))
//          {
//            CGAL_WALKER_TRACE("     ERROR: loop detected. Abording.\n");
//            throw Walk_is_cycling_exception();
//          }
#if 0
        if ( is_border_edge(intersected_edge, pm) )
          throw Walk_on_polyhedron_border_crossing_exception();
          // \todo use the boolean in the API about incomplete overlay
#endif
        break;
      }
    }

    if ( k==hedges_to_test.size() )
    {
       // this handles the case when the src is on the border
       // and the walk directly goes out
      if (target_type==POLYHEDRON_NONE) break;
      // this handle the case when the src is on the border
      // and the tgt is inside or on the boundary of the supporting mesh
      // but the path is going out
      CGAL_WALKER_TRACE("WARNING: Halting but tgt not reached!\n")
      return false;


      CGAL_WALKER_TRACE("ERROR DID NOT FIND THE INTERSECTION!!!!!!!!!!!!!!!\n")
      CGAL_assertion(!"NO INTERSECTION FOUND");
      break;
    }

    std::set <face_descriptor> current_incident_facets;
    get_incident_faces(source_type, src, current_incident_facets, pm);
    faces_already_visited.clear();

    const Barycentric_NT* barycentric_coord =
      boost::get<Barycentric_NT>( &(inter_res->second) );
    if (!barycentric_coord)
    {
      // we reached a vertex
      halfedge_descriptor hedge=canonical_hedge(intersected_edge, pm);
      bool is_target = boost::get<bool>( inter_res->second );
      if (!is_target) hedge=opposite(hedge, pm);
      vertex_descriptor vh=target(hedge, pm);

      Polyhedron_simplex_type prev_type=source_type;
      halfedge_descriptor prev_source=src;
      source_type=POLYHEDRON_VERTEX;
      src=hedge;

      std::vector<halfedge_descriptor> common_faces;

      bool did_break=false;

      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(hedge, pm))
      {
        if ( !is_border(hd, pm) ){ //note that target_faces and current_incident_facets do not contain any NULL face
          if ( !did_break && target_faces.count(face(hd, pm))!=0 )
            did_break=true;
          if ( current_incident_facets.count( face(hd, pm) )!=0 )
          {
            common_faces.push_back(hd);
            faces_already_visited.insert(face(hd, pm));
          }
        }
      }

      //insert current edge
      if (common_faces.size()==1)
      {
        // vertex-vertex case (on the mesh border)
        if (prev_type==POLYHEDRON_VERTEX){
          CGAL_assertion(target(common_faces[0], pm)==vh);
          visitor.on_input_edge(
            source(common_faces[0], pm)==target(prev_source, pm)
            ? common_faces[0]
            : opposite(next(common_faces[0], pm), pm));
        }
        else{
          // edge->vertex case (on the mesh border)
          if (prev_type==POLYHEDRON_EDGE && (target(prev_source, pm)==vh || source(prev_source, pm)==vh))
            visitor.on_input_edge(target(prev_source, pm)==vh
              ? prev_source:opposite(prev_source, pm));
          else
            // we reached a vertex opposite to the edge or in a face
            visitor.found_vertex_intersection(common_faces[0]);
        }
      }
      else
      {
        CGAL_assertion(common_faces.size()==2);
        //this is a vertex-vertex case or a edge-vertex
        halfedge_descriptor edge =
          common_edge(face(common_faces[0], pm), face(common_faces[1], pm), pm);
        visitor.on_input_edge(target(edge, pm)==vh?edge:opposite(edge, pm));
      }

      if ( did_break )
      {
        CGAL_assertion((hedge==canonical_hedge(intersected_edge, pm)) == is_target);
        common_faces.clear();

        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(hedge, pm))
        {
          if ( !is_border(hd, pm) && target_faces.count( face(hd, pm) )!=0 )
            common_faces.push_back(hd);
        }

        //insert the last edge
        if (common_faces.size()==1){
          // the last visited edge is a complete border edge
          if (target_type==POLYHEDRON_VERTEX){
            CGAL_assertion(target(common_faces[0], pm)==vh);
            visitor.on_walk_end_on_edge(
              source(common_faces[0], pm)==target(tgt, pm)
              ? opposite(common_faces[0], pm)
              : next(common_faces[0], pm) );
          }
          else
          {
            // the walk ends on a border edge
            if(target_type==POLYHEDRON_EDGE &&
               (target(tgt, pm)==vh || source(tgt, pm)==vh) )
            {
              visitor.on_walk_end_on_edge(target(tgt, pm)==vh
                ? opposite(tgt, pm):tgt);
            }
            else
              visitor.on_walk_end(face(common_faces[0], pm));
          }
        }
        else
        {
          CGAL_assertion(common_faces.size()==2);
          //this is a vertex-vertex case or a edge-vertex
          halfedge_descriptor edge =
            common_edge(face(common_faces[0], pm), face(common_faces[1], pm), pm);
          visitor.on_walk_end_on_edge(target(edge, pm)==vh?opposite(edge, pm):edge);
        }
        break;
      }
    }
    else{
      // we crossed the interior of an edge
      faces_already_visited.insert(face(src, pm));
      source_type=POLYHEDRON_EDGE;
      src=opposite(intersected_edge, pm);
      // barycentric_coord must have been computed using `canonical_hedge(intersected_edge, pm)`
      visitor.found_edge_intersection( intersected_edge, *barycentric_coord );

      if ( target_faces.count( face(src, pm) )!=0 )
      {
        CGAL_assertion ( target_type != POLYHEDRON_VERTEX ||
                         ( target(src, pm) != target(tgt, pm) &&
                           source(src, pm) != target(tgt, pm)) );
        //insert the last edge
        if ( !is_border(opposite(intersected_edge, pm), pm) )
          visitor.on_walk_end(face(opposite(intersected_edge, pm), pm));
        break;
      }
    }
    if ( is_border_edge(intersected_edge, pm) ){
      if (target_type!=POLYHEDRON_NONE){
        CGAL_WALKER_TRACE("WARNING: Halting but tgt not reached!\n")
        return false;
      }
      break; // the walk reached the boundary of the domain
    }

    if(visitor.do_break())
      break;

  }
  return true;
}

}

#undef CGAL_WALKER_TRACE

#endif // CGAL_WALK_IN_POLYGON_MESH_H
