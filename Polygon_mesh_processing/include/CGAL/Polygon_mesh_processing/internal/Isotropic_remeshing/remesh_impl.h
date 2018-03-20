// Copyright (c) 2015 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jane Tournois


#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/property_map.h>
#include <CGAL/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <CGAL/tags.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/range.hpp>
#include <boost/range/join.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/container/flat_set.hpp>

#include <map>
#include <list>
#include <vector>
#include <iterator>
#include <fstream>

#ifdef CGAL_PMP_REMESHING_DEBUG
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#define CGAL_DUMP_REMESHING_STEPS
#define CGAL_PMP_REMESHING_VERBOSE_PROGRESS
#endif

#ifdef CGAL_PMP_REMESHING_VERY_VERBOSE
#define CGAL_PMP_REMESHING_VERBOSE
#define CGAL_PMP_REMESHING_VERBOSE_PROGRESS
#endif

#ifdef CGAL_PMP_REMESHING_VERBOSE_PROGRESS
#define CGAL_PMP_REMESHING_VERBOSE
#endif


namespace CGAL {

namespace PMP = Polygon_mesh_processing;

namespace Polygon_mesh_processing {
namespace internal {

  enum Halfedge_status {
    PATCH,       //h and hopp belong to the patch to be remeshed
    PATCH_BORDER,//h belongs to the patch, hopp is MESH
    MESH,        //h and hopp belong to the mesh, not the patch
    MESH_BORDER  //h belongs to the mesh, face(hopp, pmesh) == null_face()
  };

  // A property map
  template<typename Descriptor>
  struct No_constraint_pmap
  {
  public:
    typedef Descriptor                          key_type;
    typedef bool                                value_type;
    typedef value_type&                         reference;
    typedef boost::read_write_property_map_tag  category;

    friend bool get(const No_constraint_pmap& , const key_type& ) {
      return false;
    }
    friend void put(No_constraint_pmap& , const key_type& , const bool ) {}
  };

  template <typename PM, typename FaceRange, typename FaceIndexMap>
  struct Border_constraint_pmap
  {
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PM>::edge_descriptor edge_descriptor;
    typedef FaceIndexMap FIMap;

    boost::shared_ptr< std::set<edge_descriptor> > border_edges_ptr;
    const PM* pmesh_ptr_;

  public:
    typedef edge_descriptor                     key_type;
    typedef bool                                value_type;
    typedef value_type&                         reference;
    typedef boost::read_write_property_map_tag  category;

    Border_constraint_pmap()
      : border_edges_ptr(new std::set<edge_descriptor>() )
      , pmesh_ptr_(NULL)
    {}
    Border_constraint_pmap(const PM& pmesh
                         , const FaceRange& faces
                         , const FIMap& fimap)
      : border_edges_ptr(new std::set<edge_descriptor>() )
      , pmesh_ptr_(&pmesh)
    {
      std::vector<halfedge_descriptor> border;
      PMP::border_halfedges(faces, *pmesh_ptr_, std::back_inserter(border)
        , PMP::parameters::face_index_map(fimap));

      BOOST_FOREACH(halfedge_descriptor h, border)
        border_edges_ptr->insert(edge(h, *pmesh_ptr_));
    }

    friend bool get(const Border_constraint_pmap<PM, FaceRange, FIMap>& map,
                    const edge_descriptor& e)
    {
      CGAL_assertion(map.pmesh_ptr_!=NULL);
      return map.border_edges_ptr->count(e)!=0;
    }

    friend void put(Border_constraint_pmap<PM, FaceRange, FIMap>& map,
                    const edge_descriptor& e,
                    const bool is)
    {
      CGAL_assertion(map.pmesh_ptr_ != NULL);
      if (is)
        map.border_edges_ptr->insert(e);
      else
        map.border_edges_ptr->erase(e);
    }
  };


  template <typename PM,
            typename EdgeIsConstrainedMap,
            typename FaceIndexMap>
  struct Connected_components_pmap
  {
    typedef typename boost::graph_traits<PM>::face_descriptor   face_descriptor;
    typedef std::size_t                                         Patch_id;
    typedef FaceIndexMap                                        FIMap;
    typedef EdgeIsConstrainedMap                                ECMap;
    typedef Connected_components_pmap<PM, ECMap, FIMap>         CCMap;

    typedef CGAL::dynamic_face_property_t<Patch_id> Face_property_tag;
    typedef typename boost::property_map<PM, Face_property_tag >::type Patch_ids_map;
    Patch_ids_map patch_ids_map;
    std::size_t nb_cc;

  public:
    typedef face_descriptor                     key_type;
    typedef Patch_id                            value_type;
    typedef Patch_id&                           reference;
    typedef boost::read_write_property_map_tag  category;

    //note pmesh is a non-const ref because properties are added and removed
    //modify the mesh data structure, but not the mesh itself
    Connected_components_pmap(PM& pmesh
                            , EdgeIsConstrainedMap ecmap
                            , FIMap fimap
                            , const bool do_init = true)
    {
      patch_ids_map = get(Face_property_tag(),pmesh);
      if (do_init)
      {
#ifdef CGAL_PMP_REMESHING_VERBOSE
        std::cout << "Compute connected components property map." << std::endl;
#endif
        nb_cc
          = PMP::connected_components(pmesh,
                                      patch_ids_map,
                                      PMP::parameters::edge_is_constrained_map(ecmap)
                                     .face_index_map(fimap));
        if(nb_cc == 1){
          patch_ids_map = Patch_ids_map();
        }
      }
    }


    friend value_type get(const CCMap& m, const key_type& f)
    {
      if (m.nb_cc == 1)
        return 0;
      return get(m.patch_ids_map, f);
    }
    friend void put(CCMap& m, const key_type& f, const value_type i)
    {
      if (m.nb_cc != 1)
        put(m.patch_ids_map, f, i);
    }
  };

  template<typename PM,
           typename EdgeConstraintMap,
           typename VertexPointMap>
  bool constraints_are_short_enough(const PM& pmesh,
                                    EdgeConstraintMap ecmap,
                                    VertexPointMap vpmap,
                                    const double& high)
  {
    double sqh = high*high;
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PM>::edge_descriptor     edge_descriptor;
    BOOST_FOREACH(edge_descriptor e, edges(pmesh))
    {
      if (get(ecmap, e))
      {
        halfedge_descriptor h = halfedge(e, pmesh);
        if (sqh < CGAL::squared_distance(get(vpmap, source(h, pmesh)),
                                         get(vpmap, target(h, pmesh))))
        {
          return false;
        }
      }
    }
    return true;
  }

  template<typename PolygonMesh
         , typename VertexPointMap
         , typename GeomTraits
         , typename EdgeIsConstrainedMap
         , typename VertexIsConstrainedMap
         , typename FacePatchMap
         , typename FaceIndexMap
  >
  class Incremental_remesher
  {
    typedef PolygonMesh PM;
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PM>::edge_descriptor     edge_descriptor;
    typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;
    typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;

    typedef typename GeomTraits::Point_3    Point;
    typedef typename GeomTraits::Vector_3   Vector_3;
    typedef typename GeomTraits::Plane_3    Plane_3;
    typedef typename GeomTraits::Triangle_3 Triangle_3;

    typedef Incremental_remesher<PM, VertexPointMap
                               , GeomTraits
                               , EdgeIsConstrainedMap
                               , VertexIsConstrainedMap
                               , FacePatchMap
                               , FaceIndexMap
                               > Self;

  private:
    typedef typename boost::property_traits<FacePatchMap>::value_type Patch_id;
    typedef std::vector<Triangle_3>                      Triangle_list;
    typedef std::vector<Patch_id>                        Patch_id_list;
    typedef std::map<Patch_id,std::size_t>               Patch_id_to_index_map;

    typedef CGAL::AABB_triangle_primitive<GeomTraits,
                     typename Triangle_list::iterator>    AABB_primitive;
    typedef CGAL::AABB_traits<GeomTraits, AABB_primitive> AABB_traits;
    typedef CGAL::AABB_tree<AABB_traits>                  AABB_tree;

    typedef typename boost::property_map<
      PM, CGAL::dynamic_halfedge_property_t<Halfedge_status> >::type Halfedge_status_pmap;

  public:
    Incremental_remesher(PolygonMesh& pmesh
                       , VertexPointMap& vpmap
                       , const bool protect_constraints
                       , EdgeIsConstrainedMap ecmap
                       , VertexIsConstrainedMap vcmap
                       , FacePatchMap fpmap
                       , FaceIndexMap fimap
                       , const bool build_tree = true)//built by the remesher
      : mesh_(pmesh)
      , vpmap_(vpmap)
      , build_tree_(build_tree)
      , has_border_(false)
      , input_triangles_()
      , input_patch_ids_()
      , protect_constraints_(protect_constraints)
      , patch_ids_map_(fpmap)
      , ecmap_(ecmap)
      , vcmap_(vcmap)
      , fimap_(fimap)
    {
      halfedge_status_pmap_ = get(CGAL::dynamic_halfedge_property_t<Halfedge_status>(),
                                  pmesh);
      BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh_))
        put(halfedge_status_pmap_, h, MESH);

      CGAL_assertion(CGAL::is_triangle_mesh(mesh_));
    }

    ~Incremental_remesher()
    {
      if (build_tree_){
        for(std::size_t i=0; i < trees.size();++i){
          delete trees[i];
        }
      }
    }

    template<typename FaceRange>
    void init_remeshing(const FaceRange& face_range)
    {
      tag_halfedges_status(face_range); //called first

      constrain_patch_corners(face_range);

      BOOST_FOREACH(face_descriptor f, face_range)
      {
        if (is_degenerate_triangle_face(halfedge(f,mesh_),mesh_,vpmap_,GeomTraits())){
          continue;
        }
        Patch_id pid = get_patch_id(f);
        input_triangles_.push_back(triangle(f));
        input_patch_ids_.push_back(pid);
        std::pair<typename Patch_id_to_index_map::iterator, bool> 
          res = patch_id_to_index_map.insert(std::make_pair(pid,0));
        if(res.second){
          res.first->second =  patch_id_to_index_map.size()-1;
        }
      }
      CGAL_assertion(input_triangles_.size() == input_patch_ids_.size());

      if (!build_tree_)
        return;
      trees.resize(patch_id_to_index_map.size());
      for(std::size_t i=0; i < trees.size(); ++i){
        trees[i] = new AABB_tree();
      }
      typename Triangle_list::iterator it;
      typename Patch_id_list::iterator pit;
      for(it = input_triangles_.begin(), pit = input_patch_ids_.begin();
          it != input_triangles_.end();
          ++it, ++pit){
        trees[patch_id_to_index_map[*pit]]->insert(it);
      }
      for(std::size_t i=0; i < trees.size(); ++i){
        trees[i]->build();
        trees[i]->accelerate_distance_queries();
      }
    }


    // split edges of edge_range that have their length > high
    template<typename EdgeRange>
    void split_long_edges(const EdgeRange& edge_range,
                          const double& high)
    {
      typedef boost::bimap<
        boost::bimaps::set_of<halfedge_descriptor>,
        boost::bimaps::multiset_of<double, std::greater<double> > >  Boost_bimap;
      typedef typename Boost_bimap::value_type                       long_edge;

#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "Split long edges (" << high << ")...";
      std::cout.flush();
#endif
      double sq_high = high*high;

      //collect long edges
      Boost_bimap long_edges;
      BOOST_FOREACH(edge_descriptor e, edge_range)
      {
        double sqlen = sqlength(e);
        if (sqlen > sq_high)
        {
          long_edges.insert(long_edge(halfedge(e, mesh_), sqlen));
          put(ecmap_, e, false);
        }
        else
          put(ecmap_, e, true);
      }

      //split long edges
      unsigned int nb_splits = 0;
      while (!long_edges.empty())
      {
        //the edge with longest length
        typename Boost_bimap::right_map::iterator eit = long_edges.right.begin();
        halfedge_descriptor he = eit->second;
        double sqlen = eit->first;
        long_edges.right.erase(eit);

        //split edge
        Point refinement_point = this->midpoint(he);
        halfedge_descriptor hnew = CGAL::Euler::split_edge(he, mesh_);
        CGAL_assertion(he == next(hnew, mesh_));
        ++nb_splits;

        //move refinement point
        vertex_descriptor vnew = target(hnew, mesh_);
        put(vpmap_, vnew, refinement_point);
#ifdef CGAL_PMP_REMESHING_VERY_VERBOSE
        std::cout << "   refinement point : " << refinement_point << std::endl;
#endif

        //check sub-edges
        double sqlen_new = 0.25 * sqlen;
        if (sqlen_new > sq_high)
        {
          //if it was more than twice the "long" threshold, insert them
          long_edges.insert(long_edge(hnew, sqlen_new));
          long_edges.insert(long_edge(next(hnew, mesh_), sqlen_new));
        }
        else
        {
          put(ecmap_, edge(hnew, mesh_), true);
          put(ecmap_, edge(next(hnew, mesh_), mesh_), true);
        }

        //insert new edges to keep triangular faces, and update long_edges
        if (!is_border(hnew, mesh_))
        {
          CGAL::Euler::split_face(hnew, next(next(hnew, mesh_), mesh_), mesh_);
        }

        //do it again on the other side if we're not on boundary
        halfedge_descriptor hnew_opp = opposite(hnew, mesh_);
        if (!is_border(hnew_opp, mesh_))
        {
          CGAL::Euler::split_face(prev(hnew_opp, mesh_), next(hnew_opp, mesh_), mesh_);
        }
      }
#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << " done (" << nb_splits << " splits)." << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("0-border_split.off");
#endif
    }

    // PMP book :
    // "visits all edges of the mesh
    //if an edge is longer than the given threshold `high`, the edge
    //is split at its midpoint and the two adjacent triangles are bisected (2-4 split)"
    void split_long_edges(const double& high)
    {
      typedef boost::bimap<
        boost::bimaps::set_of<halfedge_descriptor>,
        boost::bimaps::multiset_of<double, std::greater<double> > >  Boost_bimap;
      typedef typename Boost_bimap::value_type                       long_edge;

#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "Split long edges (" << high << ")..." << std::endl;
#endif
      double sq_high = high*high;

      //collect long edges
      Boost_bimap long_edges;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        if (!is_split_allowed(e))
          continue;
        double sqlen = sqlength(e);
        if(sqlen > sq_high)
          long_edges.insert(long_edge(halfedge(e, mesh_), sqlen));
      }

      //split long edges
      unsigned int nb_splits = 0;
      while (!long_edges.empty())
      {
        //the edge with longest length
        typename Boost_bimap::right_map::iterator eit = long_edges.right.begin();
        halfedge_descriptor he = eit->second;
        double sqlen = eit->first;
        long_edges.right.erase(eit);

#ifdef CGAL_PMP_REMESHING_VERBOSE_PROGRESS
        std::cout << "\r\t(" << long_edges.left.size() << " long edges, ";
        std::cout << nb_splits << " splits)";
        std::cout.flush();
#endif

        if (protect_constraints_ && !is_longest_on_faces(edge(he, mesh_)))
          continue;

        //collect patch_ids
        Patch_id patch_id = get_patch_id(face(he, mesh_));
        Patch_id patch_id_opp = get_patch_id(face(opposite(he, mesh_), mesh_));

        //split edge
        Point refinement_point = this->midpoint(he);
        halfedge_descriptor hnew = CGAL::Euler::split_edge(he, mesh_);
        CGAL_assertion(he == next(hnew, mesh_));
        ++nb_splits;

        //move refinement point
        vertex_descriptor vnew = target(hnew, mesh_);
        put(vpmap_, vnew, refinement_point);
#ifdef CGAL_PMP_REMESHING_VERY_VERBOSE
        std::cout << "   Refinement point : " << refinement_point << std::endl;
#endif

        //after splitting
        halfedge_descriptor hnew_opp = opposite(hnew, mesh_);
        halfedge_added(hnew, status(he));
        halfedge_added(hnew_opp, status(opposite(he, mesh_)));

        //check sub-edges
        double sqlen_new = 0.25 * sqlen;
        if (sqlen_new > sq_high)
        {
          //if it was more than twice the "long" threshold, insert them
          long_edges.insert(long_edge(hnew,              sqlen_new));
          long_edges.insert(long_edge(next(hnew, mesh_), sqlen_new));
        }

        //insert new edges to keep triangular faces, and update long_edges
        if (!is_on_border(hnew))
        {
          halfedge_descriptor hnew2 = CGAL::Euler::split_face(hnew,
                                                              next(next(hnew, mesh_), mesh_),
                                                              mesh_);
          Halfedge_status snew = (is_on_patch(hnew) || is_on_patch_border(hnew))
            ? PATCH
            : MESH;
          halfedge_added(hnew2,                  snew);
          halfedge_added(opposite(hnew2, mesh_), snew);
          set_patch_id(face(hnew2, mesh_), patch_id);
          set_patch_id(face(opposite(hnew2, mesh_), mesh_), patch_id);

          if (snew == PATCH)
          {
            double sql = sqlength(hnew2);
            if (sql > sq_high)
              long_edges.insert(long_edge(hnew2, sql));
          }
        }

        //do it again on the other side if we're not on boundary
        if (!is_on_border(hnew_opp))
        {
          halfedge_descriptor hnew2 = CGAL::Euler::split_face(prev(hnew_opp, mesh_),
                                                              next(hnew_opp, mesh_),
                                                              mesh_);
          Halfedge_status snew = (is_on_patch(hnew_opp) || is_on_patch_border(hnew_opp))
             ? PATCH
            : MESH;
          halfedge_added(hnew2,                  snew);
          halfedge_added(opposite(hnew2, mesh_), snew);
          set_patch_id(face(hnew2, mesh_), patch_id_opp);
          set_patch_id(face(opposite(hnew2, mesh_), mesh_), patch_id_opp);

          if (snew == PATCH)
          {
            double sql = sqlength(hnew2);
            if (sql > sq_high)
              long_edges.insert(long_edge(hnew2, sql));
          }
        }
      }
#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << " done ("<< nb_splits << " splits)." << std::endl;
#endif

#ifdef CGAL_PMP_REMESHING_DEBUG
      CGAL_expensive_assertion(is_triangle_mesh(mesh_));
      debug_status_map();
      debug_self_intersections();
#endif

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("1-edge_split.off");
#endif
    }

    // PMP book :
    // "collapses and thus removes all edges that are shorter than a
    // threshold `low`. [...] testing before each collapse whether the collapse
    // would produce an edge that is longer than `high`"
    void collapse_short_edges(const double& low, const double& high)
    {
      typedef boost::bimap<
        boost::bimaps::set_of<halfedge_descriptor>,
        boost::bimaps::multiset_of<double, std::less<double> > >  Boost_bimap;
      typedef typename Boost_bimap::value_type                    short_edge;

#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "Collapse short edges (" << low << ", " << high << ")..."
                << std::endl;
#endif
#ifdef CGAL_PMP_REMESHING_VERBOSE_PROGRESS
      std::cout << "Fill bimap...";
      std::cout.flush();
#endif
      double sq_low = low*low;
      double sq_high = high*high;

      Boost_bimap short_edges;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        double sqlen = sqlength(e);
        if( (sqlen < sq_low) && is_collapse_allowed(e) )
          short_edges.insert(short_edge(halfedge(e, mesh_), sqlen));
      }
#ifdef CGAL_PMP_REMESHING_VERBOSE_PROGRESS
      std::cout << "done." << std::endl;
#endif

      unsigned int nb_collapses = 0;
      while (!short_edges.empty())
      {
        //the edge with shortest length
        typename Boost_bimap::right_map::iterator eit = short_edges.right.begin();
        halfedge_descriptor he = eit->second;
        short_edges.right.erase(eit);

#ifdef CGAL_PMP_REMESHING_VERBOSE_PROGRESS
        std::cout << "\r\t(" << short_edges.left.size() << " short edges, ";
        std::cout << nb_collapses << " collapses)";
        std::cout.flush();
#endif

        edge_descriptor e = edge(he, mesh_);
        if (!is_collapse_allowed(e))
          continue; //situation could have changed since it was added to the bimap

        //handle the boundary case :
        //a PATCH_BORDER edge can be collapsed,
        //and an edge incident to PATCH_BORDER can be collapsed,
        //but only if the boundary vertex is kept,
        //so re-insert opposite(he) to collapse it
        if (!is_on_patch(he))
        {
          CGAL_assertion(!protect_constraints_);//is_collapse_allowed returned false
          if (is_on_border(he) || is_on_mesh(he))
          {
            he = opposite(he, mesh_); //he now is PATCH_BORDER
            CGAL_assertion(is_on_patch_border(he));
          }
        }//end if(not on PATCH)

        //let's try to collapse he into vb
        vertex_descriptor va = source(he, mesh_);
        vertex_descriptor vb = target(he, mesh_);

        bool swap_done = false;
        if (is_on_patch_border(va) && !is_on_patch_border(vb))
        {
          he = opposite(he, mesh_);
          va = source(he, mesh_);
          vb = target(he, mesh_);
          swap_done = true;
          CGAL_assertion(is_on_patch_border(vb) && !is_on_patch_border(va));
        }

        if(!collapse_does_not_invert_face(he))
        {
          if (!swap_done//if swap has already been done, don't re-swap
           && collapse_does_not_invert_face(opposite(he, mesh_)))
          {
            he = opposite(he, mesh_);
            va = source(he, mesh_);
            vb = target(he, mesh_);

            if (is_on_patch_border(va) && !is_on_patch_border(vb))
              continue;//we cannot swap again. It would lead to a face inversion
            else if (is_corner(va) && !is_corner(vb))
              continue;//idem
          }
          else
            continue;//both directions invert a face
        }
        CGAL_assertion(collapse_does_not_invert_face(he));
        CGAL_assertion(is_collapse_allowed(e));

        if (degree(va, mesh_) < 3
          || degree(vb, mesh_) < 3
          || !CGAL::Euler::does_satisfy_link_condition(e, mesh_))//necessary to collapse
          continue;

        //check that collapse would not create an edge with length > high
        //iterate on vertices va_i of the one-ring of va
        bool collapse_ok = true;
        BOOST_FOREACH(halfedge_descriptor ha, halfedges_around_target(va, mesh_))
        {
          vertex_descriptor va_i = source(ha, mesh_);
          if (sqlength(vb, va_i) > sq_high)
          {
            collapse_ok = false;
            break;
          }
        }
        //before collapsing va into vb, check that it does not break a corner
        //if it is allowed, perform the collapse
        if (collapse_ok && !is_constrained(va) && !is_corner(va))
        {
          //"collapse va into vb along e"
          // remove edges incident to va and vb, because their lengths will change
          BOOST_FOREACH(halfedge_descriptor ha, halfedges_around_target(va, mesh_))
          {
            short_edges.left.erase(ha);
            short_edges.left.erase(opposite(ha, mesh_));
          }
          BOOST_FOREACH(halfedge_descriptor hb, halfedges_around_target(vb, mesh_))
          {
            short_edges.left.erase(hb);
            short_edges.left.erase(opposite(hb, mesh_));
          }

          //before collapse
          bool mesh_border_case     = is_on_border(he);
          bool mesh_border_case_opp = is_on_border(opposite(he, mesh_));
          halfedge_descriptor ep_p  = prev(opposite(he, mesh_), mesh_);
          halfedge_descriptor epo_p = opposite(ep_p, mesh_);
          halfedge_descriptor en    = next(he, mesh_);
          halfedge_descriptor en_p  = next(opposite(he, mesh_), mesh_);
          Halfedge_status s_ep_p    = status(ep_p);
          Halfedge_status s_epo_p   = status(epo_p);
          Halfedge_status s_ep      = status(prev(he, mesh_));
          Halfedge_status s_epo     = status(opposite(prev(he, mesh_), mesh_));

          // merge halfedge_status to keep the more important on both sides
          //do it before collapse is performed to be sure everything is valid
          if (!mesh_border_case)
            merge_status(en, s_epo, s_ep);
          if (!mesh_border_case_opp)
            merge_status(en_p, s_epo_p, s_ep_p);

          halfedge_and_opp_removed(he);
          if (!mesh_border_case)
            halfedge_and_opp_removed(prev(he, mesh_));
          if (!mesh_border_case_opp)
            halfedge_and_opp_removed(prev(opposite(he, mesh_), mesh_));

          //constrained case
          bool constrained_case = is_constrained(va) || is_constrained(vb);
          if (constrained_case)
          {
            CGAL_assertion(is_constrained(va) ^ is_constrained(vb));//XOR
            set_constrained(va, false);
            set_constrained(vb, false);
          }

          //perform collapse
          Point target_point = get(vpmap_, vb);
          vertex_descriptor vkept = CGAL::Euler::collapse_edge(e, mesh_);
          put(vpmap_, vkept, target_point);
          ++nb_collapses;

          //fix constrained case
          if (constrained_case)//we have made sure that collapse goes to constrained vertex
            set_constrained(vkept, true);

          fix_degenerate_faces(vkept, short_edges, sq_low);

#ifdef CGAL_PMP_REMESHING_DEBUG
          debug_status_map();
          CGAL_assertion(!incident_to_degenerate(halfedge(vkept, mesh_)));
#endif

          //insert new/remaining short edges
          BOOST_FOREACH(halfedge_descriptor ht, halfedges_around_target(vkept, mesh_))
          {
            double sqlen = sqlength(ht);
            if( (sqlen < sq_low) && is_collapse_allowed(edge(ht, mesh_)) )
              short_edges.insert(short_edge(ht, sqlen));
          }
        }//end if(collapse_ok)
      }

#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << " done (" << nb_collapses << " collapses)." << std::endl;
#endif

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("2-edge_collapse.off");
#endif

#ifdef CGAL_PMP_REMESHING_DEBUG
      CGAL_expensive_assertion(is_triangle_mesh(mesh_));
      debug_status_map();
      debug_self_intersections();
      CGAL_assertion(0 == PMP::remove_degenerate_faces(mesh_,
        PMP::parameters::vertex_point_map(vpmap_).geom_traits(GeomTraits())));
#endif
    }

    // PMP book :
    // "equalizes the vertex valences by flipping edges.
    // The target valence is 6 and 4 for interior and boundary vertices, resp.
    // The algo. tentatively flips each edge `e` and checks whether the deviation
    // to the target valences decreases. If not, the edge is flipped back"
    void equalize_valences()
    {
#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "Equalize valences..." << std::endl;
#endif
      unsigned int nb_flips = 0;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        //only the patch edges are allowed to be flipped
        if (!is_flip_allowed(e))
          continue;

        halfedge_descriptor he = halfedge(e, mesh_);
        vertex_descriptor va = source(he, mesh_);
        vertex_descriptor vb = target(he, mesh_);
        vertex_descriptor vc = target(next(he, mesh_), mesh_);
        vertex_descriptor vd = target(next(opposite(he, mesh_), mesh_), mesh_);

        int vva = valence(va), tvva = target_valence(va);
        int vvb = valence(vb), tvvb = target_valence(vb);
        int vvc = valence(vc), tvvc = target_valence(vc);
        int vvd = valence(vd), tvvd = target_valence(vd);
        int deviation_pre = CGAL::abs(vva - tvva)
                          + CGAL::abs(vvb - tvvb)
                          + CGAL::abs(vvc - tvvc)
                          + CGAL::abs(vvd - tvvd);

        CGAL_assertion_code(Halfedge_status s1 = status(he));
        CGAL_assertion_code(Halfedge_status s1o = status(opposite(he, mesh_)));

        Patch_id pid = get_patch_id(face(he, mesh_));

        CGAL::Euler::flip_edge(he, mesh_);
        vva -= 1;
        vvb -= 1;
        vvc += 1;
        vvd += 1;
        ++nb_flips;

#ifdef CGAL_PMP_REMESHING_VERBOSE_PROGRESS
        std::cout << "\r\t(" << nb_flips << " flips)";
        std::cout.flush();
#endif
        CGAL_assertion_code(Halfedge_status s2 = status(he));
        CGAL_assertion_code(Halfedge_status s2o = status(opposite(he, mesh_)));
        CGAL_assertion(s1 == s2   && s1 == PATCH);
        CGAL_assertion(s1o == s2o && s1o == PATCH);
        CGAL_assertion(!is_border(he, mesh_));

        CGAL_assertion(
             (vc == target(he, mesh_) && vd == source(he, mesh_))
          || (vd == target(he, mesh_) && vc == source(he, mesh_)));

        int deviation_post = CGAL::abs(vva - tvva)
                           + CGAL::abs(vvb - tvvb)
                           + CGAL::abs(vvc - tvvc)
                           + CGAL::abs(vvd - tvvd);

        //check that mesh does not become non-triangle,
        //nor has inverted faces
        if (deviation_pre <= deviation_post
          || !check_normals(he)
          || incident_to_degenerate(he)
          || incident_to_degenerate(opposite(he, mesh_))
          || !is_on_triangle(he)
          || !is_on_triangle(opposite(he, mesh_))
          || !check_normals(target(he, mesh_))
          || !check_normals(source(he, mesh_)))
        {
          CGAL::Euler::flip_edge(he, mesh_);
          --nb_flips;

          CGAL_assertion_code(Halfedge_status s3 = status(he));
          CGAL_assertion(s1 == s3);
          CGAL_assertion(!is_border(he, mesh_));
          CGAL_assertion(
               (va == source(he, mesh_) && vb == target(he, mesh_))
            || (vb == source(he, mesh_) && va == target(he, mesh_)));
        }

        set_patch_id(face(he, mesh_), pid);
        set_patch_id(face(opposite(he, mesh_), mesh_), pid);
      }

#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "\r\tdone ("<< nb_flips << " flips)" << std::endl;
#endif

#ifdef CGAL_PMP_REMESHING_DEBUG
      debug_status_map();
      CGAL_assertion(0 == PMP::remove_degenerate_faces(mesh_
                            , PMP::parameters::vertex_point_map(vpmap_)
                            .geom_traits(GeomTraits())));
      debug_self_intersections();
#endif

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("3-edge_flips.off");
#endif
    }

    // PMP book :
    // "applies an iterative smoothing filter to the mesh.
    // The vertex movement has to be constrained to the vertex tangent plane [...]
    // smoothing algorithm with uniform Laplacian weights"
    void tangential_relaxation(const bool relax_constraints/*1d smoothing*/
                             , const unsigned int nb_iterations)
    {
#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "Tangential relaxation (" << nb_iterations << " iter.)...";
      std::cout << std::endl;
#endif
      for (unsigned int nit = 0; nit < nb_iterations; ++nit)
      {
#ifdef CGAL_PMP_REMESHING_VERBOSE_PROGRESS
        std::cout << "\r\t(iteration " << (nit + 1) << " / ";
        std::cout << nb_iterations << ") ";
        std::cout.flush();
#endif
      //todo : use boost::vector_property_map to improve computing time
      typedef std::map<vertex_descriptor, Vector_3> VNormalsMap;
      VNormalsMap vnormals;
      boost::associative_property_map<VNormalsMap> propmap_normals(vnormals);
      std::map<vertex_descriptor, Point> barycenters;

      // at each vertex, compute vertex normal
      // at each vertex, compute barycenter of neighbors
      BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
      {
        if (is_constrained(v) || is_isolated(v))
          continue;

        else if (is_on_patch(v))
        {
        Vector_3 vn = PMP::compute_vertex_normal(v, mesh_
                            , PMP::parameters::vertex_point_map(vpmap_)
                            .geom_traits(GeomTraits()));
        put(propmap_normals, v, vn);

          Vector_3 move = CGAL::NULL_VECTOR;
          unsigned int star_size = 0;
          BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, mesh_))
          {
            move = move + Vector_3(get(vpmap_, v), get(vpmap_, source(h, mesh_)));
            ++star_size;
          }
          CGAL_assertion(star_size > 0); //isolated vertices have already been discarded
          move = (1. / (double)star_size) * move;

          barycenters[v] = get(vpmap_, v) + move;
        }
        else if (relax_constraints
              && !protect_constraints_
              && is_on_patch_border(v)
              && !is_corner(v))
        {
          put(propmap_normals, v, CGAL::NULL_VECTOR);

          std::vector<halfedge_descriptor> border_halfedges;
          BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, mesh_))
          {
            if (is_on_patch_border(h) || is_on_patch_border(opposite(h, mesh_)))
              border_halfedges.push_back(h);
          }
          if (border_halfedges.size() == 2)//others are corner cases
          {
            vertex_descriptor ph0 = source(border_halfedges[0], mesh_);
            vertex_descriptor ph1 = source(border_halfedges[1], mesh_);
            double dot = to_double(Vector_3(get(vpmap_, v), get(vpmap_, ph0))
                                   * Vector_3(get(vpmap_, v), get(vpmap_, ph1)));
            //check squared cosine is < 0.25 (~120 degrees)
            if (0.25 < dot / (sqlength(border_halfedges[0]) * sqlength(border_halfedges[0])))
              barycenters[v] = CGAL::midpoint(midpoint(border_halfedges[0]),
                                              midpoint(border_halfedges[1]));
          }
        }
      }

      // compute moves
      typedef typename std::map<vertex_descriptor, Point>::value_type VP_pair;
      std::map<vertex_descriptor, Point> new_locations;
      BOOST_FOREACH(const VP_pair& vp, barycenters)
      {
        vertex_descriptor v = vp.first;
        Point pv = get(vpmap_, v);
        Vector_3 nv = boost::get(propmap_normals, v);
        Point qv = vp.second; //barycenter at v

        new_locations[v] = qv + (nv * Vector_3(qv, pv)) * nv;
      }

      // perform moves
      BOOST_FOREACH(const VP_pair& vp, new_locations)
      {
        const Point initial_pos = get(vpmap_, vp.first);
        const Vector_3 move(initial_pos, vp.second);

        put(vpmap_, vp.first, vp.second);

        //check that no inversion happened
        double frac = 1.;
        while (frac > 0.03 //5 attempts maximum
           && !check_normals(vp.first)) //if a face has been inverted
        {
          frac = 0.5 * frac;
          put(vpmap_, vp.first, initial_pos + frac * move);//shorten the move by 2
        }
        if (frac <= 0.02)
          put(vpmap_, vp.first, initial_pos);//cancel move
      }

      CGAL_assertion(is_valid(mesh_));
      CGAL_assertion(is_triangle_mesh(mesh_));
      }//end for loop (nit == nb_iterations)

#ifdef CGAL_PMP_REMESHING_DEBUG
      debug_self_intersections();
#endif
#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "done." << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("4-relaxation.off");
#endif
    }

    // PMP book :
    // "maps the vertices back to the surface"
    void project_to_surface()
    {
      //todo : handle the case of boundary vertices
#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "Project to surface...";
      std::cout.flush();
#endif

      BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
      {
        if (is_constrained(v) || is_isolated(v) || !is_on_patch(v))
          continue;
        //note if v is constrained, it has not moved

        Point proj = trees[patch_id_to_index_map[get_patch_id(face(halfedge(v, mesh_), mesh_))]]->closest_point(get(vpmap_, v));
        put(vpmap_, v, proj);
      }
      CGAL_assertion(is_valid(mesh_));
      CGAL_assertion(is_triangle_mesh(mesh_));
#ifdef CGAL_PMP_REMESHING_DEBUG
      debug_self_intersections();
#endif
#ifdef CGAL_PMP_REMESHING_VERBOSE
      std::cout << "done." << std::endl;
#endif

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("5-project.off");
#endif
    }

private:
  Patch_id get_patch_id(const face_descriptor& f) const
  {
    if (f == boost::graph_traits<PM>::null_face())
      return Patch_id(-1);
    return get(patch_ids_map_, f);
  }

  void set_patch_id(const face_descriptor& f, const Patch_id& i)
  {
    put(patch_ids_map_, f, i);
  }

  struct Patch_id_property_map
  {
    typedef boost::readable_property_map_tag     category;
    typedef Patch_id                             value_type;
    typedef Patch_id&                            reference;
    typedef typename Triangle_list::const_iterator key_type;

    const Self* remesher_ptr_;

    Patch_id_property_map()
      : remesher_ptr_(NULL) {}
    Patch_id_property_map(const Self& remesher)
      : remesher_ptr_(&remesher) {}

    friend Patch_id get(const Patch_id_property_map& m, key_type tr_it)
    {
      //tr_it is an iterator from triangles_
      std::size_t id_in_vec = std::distance(
        m.remesher_ptr_->input_triangles().begin(), tr_it);

      CGAL_assertion(id_in_vec < m.remesher_ptr_->input_patch_ids().size());
      CGAL_assertion(*tr_it == m.remesher_ptr_->input_triangles()[id_in_vec]);

      return m.remesher_ptr_->input_patch_ids()[id_in_vec];
    }
  };

  private:
    Triangle_3 triangle(face_descriptor f) const
    {
      halfedge_descriptor h = halfedge(f, mesh_);
      vertex_descriptor v1  = target(h, mesh_);
      vertex_descriptor v2  = target(next(h, mesh_), mesh_);
      vertex_descriptor v3  = target(next(next(h, mesh_), mesh_), mesh_);
      return Triangle_3(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
    }

    double sqlength(const vertex_descriptor& v1,
                    const vertex_descriptor& v2) const
    {
      return to_double(CGAL::squared_distance(get(vpmap_, v1), get(vpmap_, v2)));
    }

    double sqlength(const halfedge_descriptor& h) const
    {
      vertex_descriptor v1 = target(h, mesh_);
      vertex_descriptor v2 = source(h, mesh_);
      return sqlength(v1, v2);
    }

    double sqlength(const edge_descriptor& e) const
    {
      return sqlength(halfedge(e, mesh_));
    }

    Point midpoint(const halfedge_descriptor& he) const
    {
      Point p1 = get(vpmap_, target(he, mesh_));
      Point p2 = get(vpmap_, source(he, mesh_));
      return CGAL::midpoint(p1, p2);
    }

    void dump(const char* filename) const
    {
      std::ofstream out(filename);
      out << mesh_;
      out.close();
    }

    int valence(const vertex_descriptor& v) const
    {
      return static_cast<int>(degree(v, mesh_));
    }

    int target_valence(const vertex_descriptor& v) const
    {
      return (has_border_ && is_border(v, mesh_)) ? 4 : 6;
    }

    bool is_on_triangle(const halfedge_descriptor& h) const
    {
      return h == next(next(next(h, mesh_), mesh_), mesh_);
    }

    bool is_longest_on_faces(const edge_descriptor& e) const
    {
      halfedge_descriptor h = halfedge(e, mesh_);
      halfedge_descriptor hopp = opposite(h, mesh_);

      //check whether h is the longest edge in its associated face
      //overwise refinement will go for an endless loop
      double sqh = sqlength(h);
      return sqh >= sqlength(next(h, mesh_))
          && sqh >= sqlength(next(next(h, mesh_), mesh_))
          //do the same for hopp
          && sqh >= sqlength(next(hopp, mesh_))
          && sqh >= sqlength(next(next(hopp, mesh_), mesh_));
    }

    bool is_constrained(const edge_descriptor& e) const
    {
      return is_on_border(e) || is_on_patch_border(e);
    }

    bool is_split_allowed(const edge_descriptor& e) const
    {
      halfedge_descriptor h = halfedge(e, mesh_);
      halfedge_descriptor hopp = opposite(h, mesh_);

      if (protect_constraints_ && is_constrained(e))
        return false;
      else //allow splitting constraints
      {
        if (is_on_mesh(h) && is_on_mesh(hopp))
          return false;
        else if (is_on_mesh(h) && is_on_border(hopp))
          return false;
        else if (is_on_mesh(hopp) && is_on_border(h))
          return false;
        else
          return true;
      }
    }

    bool is_collapse_allowed(const edge_descriptor& e) const
    {
      halfedge_descriptor he = halfedge(e, mesh_);
      halfedge_descriptor hopp = opposite(he, mesh_);

      if (is_on_mesh(he) && is_on_mesh(hopp))
        return false;

      if (protect_constraints_ && is_constrained(e))
        return false;
      if (is_on_patch(he)) //hopp is also on patch
      {
        CGAL_assertion(is_on_patch(hopp));
        if (is_on_patch_border(target(he, mesh_)) && is_on_patch_border(source(he, mesh_)))
          return false;//collapse would induce pinching the selection
        else
          return (is_collapse_allowed_on_patch(he)
               && is_collapse_allowed_on_patch(hopp));
      }
      else if (is_on_patch_border(he))
        return is_collapse_allowed_on_patch_border(he);
      else if (is_on_patch_border(hopp))
        return is_collapse_allowed_on_patch_border(hopp);
      return false;
    }

    bool is_collapse_allowed_on_patch(const halfedge_descriptor& he) const
    {
      halfedge_descriptor hopp = opposite(he, mesh_);

      if (is_on_patch_border(next(he, mesh_)) && is_on_patch_border(prev(he, mesh_)))
        return false;//too many cases to be handled
      if (is_on_patch_border(next(hopp, mesh_)) && is_on_patch_border(prev(hopp, mesh_)))
        return false;//too many cases to be handled
      else if (is_on_patch_border(next(he, mesh_)))
      {
        //avoid generation of degenerate faces, and self-intersections
        if (source(he, mesh_) ==
          target(next(next_on_patch_border(next(he, mesh_)), mesh_), mesh_))
          return false;
      }
      else if (is_on_patch_border(prev(hopp, mesh_)))
      {
        //avoid generation of degenerate faces, and self-intersections
        if (target(hopp, mesh_) ==
          source(prev(prev_on_patch_border(prev(hopp, mesh_)), mesh_), mesh_))
          return false;
      }
      return true;
    }

    bool is_collapse_allowed_on_patch_border(const halfedge_descriptor& h) const
    {
      CGAL_precondition(is_on_patch_border(h));
      halfedge_descriptor hopp = opposite(h, mesh_);

      if (is_on_patch_border(next(h, mesh_)) && is_on_patch_border(prev(h, mesh_)))
        return false;

      if (is_on_patch_border(hopp))
      {
        if (is_on_patch_border(next(hopp, mesh_)) && is_on_patch_border(prev(hopp, mesh_)))
          return false;
        else if (next_on_patch_border(h) == hopp)
          return false; //isolated patch border
        else
          return true;
      }
      CGAL_assertion(is_on_mesh(hopp) || is_on_border(hopp));
      return true;//we already checked we're not pinching a hole in the patch
    }

    bool is_flip_allowed(const edge_descriptor& e) const
    {
      return is_flip_allowed(halfedge(e, mesh_))
          && is_flip_allowed(opposite(halfedge(e, mesh_), mesh_));
    }

    bool is_flip_allowed(const halfedge_descriptor& h) const
    {
      if (!is_on_patch(h))
        return false;
      if (!is_on_patch_border(target(h, mesh_)))
        return true;
      if ( is_on_patch_border(next(h, mesh_))
        && is_on_patch_border(prev(opposite(h, mesh_), mesh_)))
        return false;
      return true;
    }

    halfedge_descriptor next_on_patch_border(const halfedge_descriptor& h) const
    {
      CGAL_precondition(is_on_patch_border(h));
      CGAL_assertion_code(const Patch_id& pid = get_patch_id(face(h, mesh_)));

      halfedge_descriptor end = opposite(h, mesh_);
      halfedge_descriptor nxt = next(h, mesh_);
      do
      {
        if (is_on_patch_border(nxt))
        { 
          CGAL_assertion(get_patch_id(face(nxt, mesh_)) == pid);
          return nxt;
        }
        nxt = next(opposite(nxt, mesh_), mesh_);
      }
      while (end != nxt);

      CGAL_assertion(get_patch_id(face(nxt, mesh_)) == pid);
      CGAL_assertion(is_on_patch_border(end));
      return end;
    }

    halfedge_descriptor prev_on_patch_border(const halfedge_descriptor& h) const
    {
      CGAL_precondition(is_on_patch_border(h));
      CGAL_assertion_code(const Patch_id& pid = get_patch_id(face(h, mesh_)));

      halfedge_descriptor end = opposite(h, mesh_);
      halfedge_descriptor prv = prev(h, mesh_);
      do
      {
        if (is_on_patch_border(prv))
        {
          CGAL_assertion(get_patch_id(face(prv, mesh_)) == pid);
          return prv;
        }
        prv = prev(opposite(prv, mesh_), mesh_);
      }
      while (end != prv);

      CGAL_assertion(is_on_patch_border(end));
      CGAL_assertion(get_patch_id(face(prv, mesh_)) == pid);
      return end;
    }

    bool collapse_does_not_invert_face(const halfedge_descriptor& h) const
    {
      vertex_descriptor vs = source(h, mesh_);
      vertex_descriptor vt = target(h, mesh_);
      
      //backup source point
      Point ps = get(vpmap_, vs);
      //move source at target
      put(vpmap_, vs, get(vpmap_, vt));

      //collect normals to faces around vs AND vt
      //vertices are at the same location, but connectivity is still be same,
      //with plenty of degenerate triangles (which are common to both stars)
      bool res = check_normals(
                   boost::range::join(
                     halfedges_around_target(h, mesh_),
                     halfedges_around_target(opposite(h, mesh_), mesh_)));

      //restore position
      put(vpmap_, vs, ps);
      return res;
    }

    bool is_constrained(const vertex_descriptor& v) const
    {
      return get(vcmap_, v);
    }
    void set_constrained(const vertex_descriptor& v, const bool b)
    {
      put(vcmap_, v, b);
    }
    bool is_isolated(const vertex_descriptor& v) const
    {
      return halfedges_around_target(v, mesh_).empty();
    }

    bool is_corner(const vertex_descriptor& v) const
    {
      if(! has_border_){
        return false;
      }
      unsigned int nb_incident_features = 0;
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, mesh_))
      {
        if (is_on_border(h) || is_on_patch_border(h))
          ++nb_incident_features;
        if (nb_incident_features > 2)
          return true;
      }
      return (nb_incident_features == 1);
    }

    Vector_3 compute_normal(const face_descriptor& f) const
    {
      if (f == boost::graph_traits<PM>::null_face())
        return CGAL::NULL_VECTOR;

      halfedge_descriptor hd = halfedge(f, mesh_);
      typename boost::property_traits<VertexPointMap>::reference
        p = get(vpmap_, target(hd, mesh_));
      hd = next(hd,mesh_);
      typename boost::property_traits<VertexPointMap>::reference
        q = get(vpmap_, target(hd, mesh_));
      hd = next(hd,mesh_);
      typename boost::property_traits<VertexPointMap>::reference
        r =get(vpmap_, target(hd, mesh_));
      
      if (GeomTraits().collinear_3_object()(p,q,r))
        return CGAL::NULL_VECTOR;
      else
        return PMP::compute_face_normal(f, mesh_, parameters::vertex_point_map(vpmap_));
    }

    template <typename FaceRange>
    void constrain_patch_corners(const FaceRange& face_range)
    {
      boost::container::flat_set<vertex_descriptor> visited;

      BOOST_FOREACH(face_descriptor f, face_range)
      {
        BOOST_FOREACH(halfedge_descriptor h,
                      halfedges_around_face(halfedge(f, mesh_), mesh_))
        {
          vertex_descriptor vt = target(h, mesh_);
          //treat target(h, mesh_)
          if (visited.find(vt) != visited.end())
            continue;//already treated

          if (status(h) == PATCH)//h not on patch boundary
            continue;            //so neither is target(h, mesh_)

          //count incident MESH_BORDER edges
          unsigned int nb_incident_borders = 0;
          BOOST_FOREACH(halfedge_descriptor hv,
                        halfedges_around_target(h, mesh_))
          {
            CGAL_assertion(vt == target(hv, mesh_));
            if ( (status(hv) == PATCH_BORDER && status(opposite(hv, mesh_)) == MESH_BORDER)
              || (status(hv) == MESH_BORDER && status(opposite(hv, mesh_)) == PATCH_BORDER))
            nb_incident_borders++;
          }

          if (nb_incident_borders == 1) //this is a special corner
            set_constrained(vt, true);

          visited.insert(vt);
        }
      }
    }

    template<typename FaceRange>
    void tag_halfedges_status(const FaceRange& face_range)
    {
      //tag MESH,        //h and hopp belong to the mesh, not the patch
      //tag MESH_BORDER  //h belongs to the mesh, face(hopp, pmesh) == null_face()
      BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh_))
      {
        //being part of the border of the mesh is predominant
        if (is_border(h, mesh_)){
          set_status(h, MESH_BORDER); //erase previous value if exists
          has_border_ = true;
        } else {
          set_status(h, MESH);
        }
      }

      //tag PATCH,       //h and hopp belong to the patch to be remeshed
      BOOST_FOREACH(face_descriptor f, face_range)
      {
        BOOST_FOREACH(halfedge_descriptor h,
                      halfedges_around_face(halfedge(f, mesh_), mesh_))
        {
          set_status(h, PATCH);
        }
      }

      internal::Border_constraint_pmap<PM, FaceRange, FaceIndexMap>
        border_map(mesh_, face_range, fimap_);
      //override the border of PATCH
      //tag PATCH_BORDER,//h belongs to the patch, hopp doesn't
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        if (get(ecmap_, e)
          || get(border_map, e)
          || get_patch_id(face(halfedge(e, mesh_), mesh_))
              != get_patch_id(face(opposite(halfedge(e, mesh_), mesh_), mesh_)))
        {
          //deal with h and hopp for borders that are sharp edges to be preserved
          halfedge_descriptor h = halfedge(e, mesh_);
          if (status(h) == PATCH){
            set_status(h, PATCH_BORDER);
            has_border_ = true;
          }

          halfedge_descriptor hopp = opposite(h, mesh_);
          if (status(hopp) == PATCH){
            set_status(hopp, PATCH_BORDER);
            has_border_ = true;
          }
        }
      }
    }


    // for a Surface_mesh::Property_map we make MESH the default value
    Halfedge_status status(const halfedge_descriptor& h) const
    {
      return get(halfedge_status_pmap_,h);
    }

    void set_status(const halfedge_descriptor& h,
                    const Halfedge_status& s)
    {
      put(halfedge_status_pmap_,h,s);
    }

    void merge_status(const halfedge_descriptor& en,
                      const Halfedge_status& s_epo,
                      const Halfedge_status& s_ep)
    {
      //get missing data
      halfedge_descriptor eno = opposite(en, mesh_);
      Halfedge_status s_eno = status(eno);

      if (s_epo == MESH_BORDER && s_eno == MESH_BORDER)
        return;

      if(s_epo == MESH_BORDER
        || s_ep == MESH_BORDER
        || s_epo == PATCH_BORDER
        || s_ep == PATCH_BORDER)
      {
        set_status(en, s_epo);
        set_status(eno, s_ep);
      }
      // else keep current status for en and eno
    }

    template<typename Bimap>
    void fix_degenerate_faces(const vertex_descriptor& v,
                              Bimap& short_edges,
                              const double& sq_low)
    {
      CGAL_assertion_code(std::size_t nb_done = 0);
      boost::unordered_set<halfedge_descriptor> degenerate_faces;
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_target(halfedge(v, mesh_), mesh_))
      {
        if (is_border(h, mesh_))
          continue;
        if (is_degenerate_triangle_face(h, mesh_, vpmap_, GeomTraits()))
          degenerate_faces.insert(h);
      }
      while(!degenerate_faces.empty())
      {
        halfedge_descriptor h = *(degenerate_faces.begin());
        degenerate_faces.erase(degenerate_faces.begin());

        if (!is_degenerate_triangle_face(h, mesh_, vpmap_, GeomTraits()))
          //this can happen when flipping h has consequences further in the mesh
          continue;

        //check that opposite is not also degenerate
        if (degenerate_faces.find(opposite(h, mesh_)) != degenerate_faces.end())
          degenerate_faces.erase(opposite(h, mesh_));

        if(is_border(h, mesh_))
          continue;

        BOOST_FOREACH(halfedge_descriptor hf,
                      halfedges_around_face(h, mesh_))
        {
          vertex_descriptor vc = target(hf, mesh_);
          vertex_descriptor va = target(next(hf, mesh_), mesh_);
          vertex_descriptor vb = target(next(next(hf, mesh_), mesh_), mesh_);
          Vector_3 ab(get(vpmap_,va), get(vpmap_,vb));
          Vector_3 ac(get(vpmap_,va), get(vpmap_,vc));
          if (ab * ac < 0)
          {
            halfedge_descriptor hfo = opposite(hf, mesh_);
            halfedge_descriptor h_ab = prev(hf, mesh_);
            halfedge_descriptor h_ca = next(hf, mesh_);

            short_edges.left.erase(hf);
            short_edges.left.erase(hfo);

            CGAL::Euler::flip_edge(hf, mesh_);
            CGAL_assertion_code(++nb_done);

            //update status
            set_status(h_ab, merge_status(h_ab, hf, hfo));
            set_status(h_ca, merge_status(h_ca, hf, hfo));
            if (is_on_patch(h_ca) || is_on_patch_border(h_ca))
            {
              set_status(hf, PATCH);
              set_status(hfo, PATCH);
            }
#ifdef CGAL_PMP_REMESHING_DEBUG
            debug_status_map();
#endif

            //insert new edges in 'short_edges'
            if (is_collapse_allowed(edge(hf, mesh_)))
            {
              double sqlen = sqlength(hf);
              if (sqlen < sq_low)
                short_edges.insert(typename Bimap::value_type(hf, sqlen));
            }

            if (!is_border(hf, mesh_)
              && is_degenerate_triangle_face(hf, mesh_, vpmap_, GeomTraits()))
              degenerate_faces.insert(hf);
            if (!is_border(hfo, mesh_)
              && is_degenerate_triangle_face(hfo, mesh_, vpmap_, GeomTraits()))
              degenerate_faces.insert(hfo);

            break;
          }
        }
      }
#ifdef CGAL_PMP_REMESHING_DEBUG
      debug_status_map();
#endif
    }

    bool incident_to_degenerate(const halfedge_descriptor& he)
    {
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_target(he, mesh_))
      {
        if (is_border(h, mesh_))
          continue;
        if (is_degenerate_triangle_face(h, mesh_, vpmap_, GeomTraits()))
          return true;
      }
      return false;
    }

    Halfedge_status merge_status(const halfedge_descriptor& h1,
      const halfedge_descriptor& h2,
      const halfedge_descriptor& h3)
    {
      Halfedge_status s1 = status(h1);
      if (s1 == MESH_BORDER) return s1;
      Halfedge_status s2 = status(h2);
      if (s2 == MESH_BORDER) return s2;
      Halfedge_status s3 = status(h3);
      if (s3 == MESH_BORDER) return s3;
      else if (s1 == PATCH_BORDER) return s1;
      else if (s2 == PATCH_BORDER) return s2;
      else if (s3 == PATCH_BORDER) return s3;

      CGAL_assertion(s1 == s2 && s1 == s3);
      return s1;
    }

    bool is_on_patch(const halfedge_descriptor& h) const
    {
      bool res =(status(h) == PATCH);
      CGAL_assertion(res == (status(opposite(h, mesh_)) == PATCH));
      return res;
    }

    bool is_on_patch(const face_descriptor& f) const
    {
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_face(halfedge(f, mesh_), mesh_))
      {
        if (is_on_patch(h) || is_on_patch_border(h))
          return true;
      }
      return false;
    }

    bool is_on_patch(const vertex_descriptor& v) const
    {
      if(! has_border_){
        return true;
      }
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_target(v, mesh_))
      {
        if (!is_on_patch(h))
          return false;
      }
      return true;
    }

public:
    bool is_on_patch_border(const halfedge_descriptor& h) const
    {
      bool res = (status(h) == PATCH_BORDER);
      if (res)
      {
        CGAL_assertion_code(Halfedge_status hs = status(opposite(h, mesh_)));
        CGAL_assertion(hs == MESH_BORDER
                    || hs == MESH
                    || hs == PATCH_BORDER);//when 2 incident patches are remeshed
      }
      return res;
    }
    bool is_on_patch_border(const edge_descriptor& e) const
    {
      return is_on_patch_border(halfedge(e,mesh_))
          || is_on_patch_border(opposite(halfedge(e, mesh_), mesh_));
    }
    bool is_on_patch_border(const vertex_descriptor& v) const
    {
      if(! has_border_){
        return false;
      }
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, mesh_))
      {
        if (is_on_patch_border(h) || is_on_patch_border(opposite(h, mesh_)))
          return true;
      }
      return false;
    }

    bool is_on_border(const halfedge_descriptor& h) const
    {
      bool res = (status(h) == MESH_BORDER);
      CGAL_assertion(res == is_border(h, mesh_));
      CGAL_assertion(res == is_border(next(h, mesh_), mesh_));
      return res;
    }

    bool is_on_border(const edge_descriptor& e) const
    {
      return is_on_border(halfedge(e, mesh_))
          || is_on_border(opposite(halfedge(e, mesh_), mesh_));
    }

    bool is_on_mesh(const halfedge_descriptor& h) const
    {
      return status(h) == MESH;
    }

private:
    void halfedge_added(const halfedge_descriptor& h,
                        const Halfedge_status& s)
    {
        set_status(h, s);
    }

    void halfedge_and_opp_removed(const halfedge_descriptor& h)
    {
      set_status(h,MESH);
      set_status(opposite(h, mesh_),MESH);
    }

    std::size_t nb_valid_halfedges() const
    {
      return static_cast<std::size_t>(
        std::distance(halfedges(mesh_).first, halfedges(mesh_).second));
    }

    void debug_status_map() const
    {
      unsigned int nb_border = 0;
      unsigned int nb_mesh = 0;
      unsigned int nb_patch = 0;
      unsigned int nb_patch_border = 0;

      BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh_))
      {
        if(is_on_patch(h))              nb_patch++;
        else if(is_on_patch_border(h))  nb_patch_border++;
        else if(is_on_mesh(h))          nb_mesh++;
        else if(is_on_border(h))        nb_border++;
        else CGAL_assertion(false);
      }
    }

#ifdef CGAL_PMP_REMESHING_DEBUG
    void debug_self_intersections() const
    {
      std::cout << "Test self intersections...";
      std::vector<std::pair<face_descriptor, face_descriptor> > facets;
      PMP::self_intersections(
        mesh_,
        std::back_inserter(facets),
        PMP::parameters::vertex_point_map(vpmap_));
      //CGAL_assertion(facets.empty());
      std::cout << "done ("<< facets.size() <<" facets)." << std::endl;
    }

    void debug_self_intersections(const vertex_descriptor& v) const
    {
      std::cout << "Test self intersections...";
      std::vector<std::pair<face_descriptor, face_descriptor> > facets;
      PMP::self_intersections(
        faces_around_target(halfedge(v, mesh_), mesh_),
        mesh_,
        std::back_inserter(facets),
        PMP::parameters::vertex_point_map(vpmap_));
      //CGAL_assertion(facets.empty());
      std::cout << "done ("<< facets.size() <<" facets)." << std::endl;
    }
#endif

    //check whether the normals to faces incident to v
    //have all their 2 by 2 dot products > 0
    bool check_normals(const vertex_descriptor& v) const
    {
      return check_normals(halfedges_around_target(halfedge(v, mesh_), mesh_));
    }

    template <typename HalfedgeRange>
    bool check_normals(const HalfedgeRange& hedges) const
    {
      typedef std::multimap<Patch_id, Vector_3>   Normals_multimap;
      typedef typename Normals_multimap::iterator Normals_iterator;

      Normals_multimap normals_per_patch;
      BOOST_FOREACH(halfedge_descriptor hd, hedges)
      {
        Vector_3 n = compute_normal(face(hd, mesh_));
        if (n == CGAL::NULL_VECTOR) //for degenerate faces
          continue;
        Patch_id pid = get_patch_id(face(hd, mesh_));

        normals_per_patch.insert(std::make_pair(pid, n));
      }

      //on each surface patch,
      //check all normals have same orientation
      for (Normals_iterator it = normals_per_patch.begin();
        it != normals_per_patch.end();/*done inside loop*/)
      {
        std::vector<Vector_3> normals;
        std::pair<Normals_iterator, Normals_iterator> n_range
          = normals_per_patch.equal_range((*it).first);
        for (Normals_iterator iit = n_range.first; iit != n_range.second; ++iit)
          normals.push_back((*iit).second);

        if (!check_orientation(normals))
          return false;

        it = n_range.second;
      }
      return true;
    }

    bool check_normals(const halfedge_descriptor& h) const
    {
      if (!is_on_patch(h))
        return true;//nothing to say
      Vector_3 n = compute_normal(face(h, mesh_));
      Vector_3 no = compute_normal(face(opposite(h, mesh_), mesh_));
      return n * no > 0.;
    }

    bool check_orientation(const std::vector<Vector_3>& normals) const
    {
      if (normals.size() < 2)
        return true;
      for (std::size_t i = 1; i < normals.size(); ++i)/*start at 1 on purpose*/
      {
        double dot = to_double(normals[i - 1] * normals[i]);
        if (dot <= 0.)
          return false;
      }
      return true;
    }

  public:
    const Triangle_list& input_triangles() const {
      return input_triangles_;
    }

    const Patch_id_list& input_patch_ids() const {
      return input_patch_ids_;
    }

    void update_constraints_property_map()
    {
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        if (is_on_patch_border(halfedge(e, mesh_))
          || is_on_patch_border(opposite(halfedge(e, mesh_), mesh_)))
          put(ecmap_, e, true);
        else
          put(ecmap_, e, false);
      }
    }

  private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    bool build_tree_;
    bool has_border_;
    std::vector<AABB_tree*> trees;
    Patch_id_to_index_map patch_id_to_index_map;
    Triangle_list input_triangles_;
    Patch_id_list input_patch_ids_;
    Halfedge_status_pmap halfedge_status_pmap_;
    bool protect_constraints_;
    FacePatchMap patch_ids_map_;
    EdgeIsConstrainedMap ecmap_;
    VertexIsConstrainedMap vcmap_;
    FaceIndexMap fimap_;

  };//end class Incremental_remesher
}//end namespace internal
}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
