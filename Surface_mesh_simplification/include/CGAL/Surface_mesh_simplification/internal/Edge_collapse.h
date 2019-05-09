// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Modifiable_priority_queue.h>

#include <boost/scoped_array.hpp>

namespace CGAL {
namespace Surface_mesh_simplification {

// Implementation of the vertex-pair collapse triangulated surface mesh simplification algorithm
template<class TM_,
         class ShouldStop_,
         class VertexIndexMap_,
         class VertexPointMap_,
         class EdgeIndexMap_,
         class EdgeIsConstrainedMap_,
         class GetCost_,
         class GetPlacement_,
         class VisitorT_>
class EdgeCollapse
{
public:
  typedef TM_                                                             TM;
  typedef ShouldStop_                                                     ShouldStop;
  typedef VertexIndexMap_                                                 VertexIndexMap;
  typedef VertexPointMap_                                                 VertexPointMap;
  typedef EdgeIndexMap_                                                   EdgeIndexMap;
  typedef EdgeIsConstrainedMap_                                           EdgeIsConstrainedMap;
  typedef GetCost_                                                        GetCost;
  typedef GetPlacement_                                                   GetPlacement;
  typedef VisitorT_                                                       VisitorT;

  typedef EdgeCollapse                                                    Self;

  typedef Edge_profile<TM, VertexPointMap>                                Profile;

  typedef boost::graph_traits<TM>                                         GraphTraits;
  typedef boost::graph_traits<const TM>                                   ConstGraphTraits;

  typedef typename GraphTraits::vertex_descriptor                         vertex_descriptor;
  typedef typename GraphTraits::vertex_iterator                           vertex_iterator;
  typedef typename GraphTraits::halfedge_descriptor                       halfedge_descriptor;
  typedef typename GraphTraits::halfedge_iterator                         halfedge_iterator;
  typedef CGAL::Halfedge_around_source_iterator<TM>                       out_edge_iterator;
  typedef CGAL::Halfedge_around_target_iterator<TM>                       in_edge_iterator;
  typedef typename GraphTraits::traversal_category                        traversal_category;
  typedef typename GraphTraits::edges_size_type                           size_type;

  typedef typename GraphTraits::edge_iterator                             edge_iterator;
  typedef VertexPointMap                                                  Vertex_point_pmap;
  typedef typename boost::property_traits<Vertex_point_pmap>::value_type  Point;
  typedef typename Kernel_traits<Point>::Kernel                           Traits;

  typedef typename Traits::Equal_3                                        Equal_3;

  typedef typename Traits::Vector_3                                       Vector;
  typedef typename Traits::FT                                             FT;

  typedef boost::optional<FT>                                             Cost_type;
  typedef boost::optional<Point>                                          Placement_type;

  struct Compare_id
  {
    Compare_id() : mAlgorithm(0) {}
    Compare_id(const Self* aAlgorithm) : mAlgorithm(aAlgorithm) {}

    bool operator()(const halfedge_descriptor& a, const halfedge_descriptor& b) const {
      return mAlgorithm->get_halfedge_id(a) < mAlgorithm->get_halfedge_id(b);
    }

    Self const* mAlgorithm;
  };

  struct Compare_cost
  {
    Compare_cost() : mAlgorithm(0) {}
    Compare_cost(const Self* aAlgorithm) : mAlgorithm(aAlgorithm) {}

    bool operator()(const halfedge_descriptor& a, const halfedge_descriptor& b) const
    {
      // NOTE: A cost is a boost::optional<> value.
      // Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
      // In consequence, edges with undefined costs will be promoted to the top of the priority queue and poped out first.
      return mAlgorithm->get_data(a).cost() < mAlgorithm->get_data(b).cost();
    }

    Self const* mAlgorithm;
  };

  struct edge_id
    : boost::put_get_helper<size_type, edge_id>
  {
    typedef boost::readable_property_map_tag category;
    typedef size_type                        value_type;
    typedef size_type                        reference;
    typedef halfedge_descriptor              key_type;

    edge_id() : mAlgorithm(0) {}
    edge_id(Self const* aAlgorithm) : mAlgorithm(aAlgorithm) {}

    size_type operator[](const halfedge_descriptor e) const { return mAlgorithm->get_edge_id(e); }

    Self const* mAlgorithm;
  };

  typedef Modifiable_priority_queue<halfedge_descriptor, Compare_cost, edge_id>     PQ;
  typedef typename PQ::handle                                                       pq_handle;

  // An Edge_data is associated with EVERY _ edge in the mesh (collapsable or not).
  // It relates the edge with the PQ-handle needed to update the priority queue
  // It also relates the edge with a policy-based cache
  class Edge_data
  {
  public :
    Edge_data() : mPQHandle() {}

    const Cost_type& cost() const { return mCost; }
    Cost_type& cost() { return mCost; }

    pq_handle PQ_handle() const { return mPQHandle;}
    bool is_in_PQ() const { return mPQHandle != PQ::null_handle(); }
    void set_PQ_handle(pq_handle h) { mPQHandle = h; }
    void reset_PQ_handle() { mPQHandle = PQ::null_handle(); }

  private:
    Cost_type mCost;
    pq_handle mPQHandle;
  };

  typedef Edge_data*                                                                Edge_data_ptr;
  typedef boost::scoped_array<Edge_data>                                            Edge_data_array;

public:
  EdgeCollapse(TM& aSurface,
               const ShouldStop& aShouldStop,
               const VertexIndexMap& aVertex_index_map,
               const VertexPointMap& aVertex_point_map,
               const EdgeIndexMap& aEdge_index_map,
               const EdgeIsConstrainedMap& aEdge_is_constrained_map,
               const GetCost& aGetCost,
               const GetPlacement& aGetPlacement,
               VisitorT aVisitor);

  int run();

private:
  void collect();
  void loop();
  bool is_collapse_topologically_valid(const Profile& aProfile);
  bool is_tetrahedron(const halfedge_descriptor h1);
  bool is_open_triangle(const halfedge_descriptor h1);
  bool is_collapse_geometrically_valid(const Profile& aProfile, Placement_type aPlacement);
  void collapse(const Profile& aProfile, Placement_type aPlacement);
  void update_neighbors(const vertex_descriptor aKeptV);

  Profile create_profile(const halfedge_descriptor aEdge) {
    return Profile(aEdge, mSurface, Vertex_index_map, Vertex_point_map, Edge_index_map, m_has_border);
  }

  size_type get_halfedge_id(const halfedge_descriptor aEdge) const { return Edge_index_map[aEdge]; }
  size_type get_edge_id(const halfedge_descriptor aEdge) const { return get_halfedge_id(aEdge) / 2; }

  bool is_primary_edge(const halfedge_descriptor aEdge) const { return (get_halfedge_id(aEdge) % 2) == 0; }
  halfedge_descriptor primary_edge(const halfedge_descriptor aEdge) {
    return is_primary_edge(aEdge) ? aEdge : opposite(aEdge, mSurface);
  }

  bool is_border(const vertex_descriptor aV) const;
  bool is_border(const halfedge_descriptor aEdge) const {
    return face(aEdge, mSurface) == boost::graph_traits<TM>::null_face();
  }

  bool is_constrained(const vertex_descriptor aVertex) const;
  bool is_constrained(const halfedge_descriptor aEdge) const {
    return get(Edge_is_constrained_map, edge(aEdge, mSurface));
  }

  bool is_border_or_constrained(const vertex_descriptor aV) const;
  bool is_border_or_constrained(const halfedge_descriptor aEdge) const { return is_border(aEdge) || is_constrained(aEdge); }

  bool is_edge_a_border(const halfedge_descriptor aEdge) const
  {
    return is_border(aEdge) || is_border(opposite(aEdge, mSurface));
  }

  bool are_shared_triangles_valid(const Point& p0, const Point& p1, const Point& p2, const Point& p3) const;

  halfedge_descriptor find_connection(const vertex_descriptor v0, const vertex_descriptor v1) const;
  vertex_descriptor find_exterior_link_triangle_3rd_vertex(const halfedge_descriptor e, const vertex_descriptor v0, const vertex_descriptor v1) const;

  Edge_data& get_data(const halfedge_descriptor aEdge) const
  {
    CGAL_assertion(is_primary_edge(aEdge));
    return mEdgeDataArray[get_edge_id(aEdge)];
  }

  typename boost::property_traits<VertexPointMap>::reference
  get_point(const vertex_descriptor aV) const { return get(Vertex_point_map, aV); }

  boost::tuple<vertex_descriptor,vertex_descriptor> get_vertices(const halfedge_descriptor aEdge) const
  {
    vertex_descriptor p, q;
    p = source(aEdge, mSurface);
    q = target(aEdge, mSurface);
    return boost::make_tuple(p,q);
  }

  boost::tuple<vertex_descriptor,vertex_descriptor> get_vertices(const halfedge_descriptor aEdge)
  {
    vertex_descriptor p,q;
    p = source(aEdge, mSurface);
    q = target(aEdge, mSurface);
    return boost::make_tuple(p,q);
  }

  std::string vertex_to_string(const vertex_descriptor v) const
  {
    const Point& p = get_point(v);
    return boost::str(boost::format("[V%1%:%2%]") % get(Vertex_index_map,v) % xyz_to_string(p));
  }

  std::string edge_to_string(const halfedge_descriptor aEdge) const
  {
    vertex_descriptor p, q;
    boost::tie(p,q) = get_vertices(aEdge);
    return boost::str(boost::format("{E%1% %2%->%3%}%4%") % get(Edge_index_map, aEdge) % vertex_to_string(p) % vertex_to_string(q) % (is_border(aEdge) ? " (BORDER)" : (is_border(opposite(aEdge, mSurface)) ? " (~BORDER)": "")));
  }

  Cost_type get_cost(const Profile& aProfile) const {
    return Get_cost(aProfile, get_placement(aProfile));
  }

  Placement_type get_placement(const Profile& aProfile) const {
    return Get_placement(aProfile);
  }

  void insert_in_PQ(const halfedge_descriptor aEdge, Edge_data& aData)
  {
    CGAL_SURF_SIMPL_TEST_assertion(is_primary_edge(aEdge));
    CGAL_SURF_SIMPL_TEST_assertion(!aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(!mPQ->contains(aEdge));

    aData.set_PQ_handle(mPQ->push(aEdge));

    CGAL_SURF_SIMPL_TEST_assertion(aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(mPQ->contains(aEdge));
  }

  void update_in_PQ(const halfedge_descriptor aEdge, Edge_data& aData)
  {
    CGAL_SURF_SIMPL_TEST_assertion(is_primary_edge(aEdge));
    CGAL_SURF_SIMPL_TEST_assertion(aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(mPQ->contains(aEdge));

    aData.set_PQ_handle(mPQ->update(aEdge,aData.PQ_handle()));

    CGAL_SURF_SIMPL_TEST_assertion(aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(mPQ->contains(aEdge));
  }

  void remove_from_PQ(const halfedge_descriptor aEdge, Edge_data& aData)
  {
    CGAL_SURF_SIMPL_TEST_assertion(is_primary_edge(aEdge));
    CGAL_SURF_SIMPL_TEST_assertion(aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(mPQ->contains(aEdge));

    aData.set_PQ_handle(mPQ->erase(aEdge,aData.PQ_handle()));

    CGAL_SURF_SIMPL_TEST_assertion(!aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(!mPQ->contains(aEdge));
  }

  boost::optional<halfedge_descriptor> pop_from_PQ()
  {
    boost::optional<halfedge_descriptor> rEdge = mPQ->extract_top();
    if(rEdge)
    {
      CGAL_SURF_SIMPL_TEST_assertion(is_primary_edge(*rEdge));
      CGAL_SURF_SIMPL_TEST_assertion(get_data(*rEdge).is_in_PQ());

      get_data(*rEdge).reset_PQ_handle();

      CGAL_SURF_SIMPL_TEST_assertion(!get_data(*rEdge).is_in_PQ());
      CGAL_SURF_SIMPL_TEST_assertion(!mPQ->contains(*rEdge));
    }
    return rEdge;
  }

  /// Functions to ensure the backward compatibility before addition of the constrained edge map
  template<class AEdgeIsConstrainedMap>
  vertex_descriptor
  halfedge_collapse_bk_compatibility(const halfedge_descriptor pq,
                                     AEdgeIsConstrainedMap aEdge_is_constrained_map)
  {
    return CGAL::Euler::collapse_edge(edge(pq, mSurface), mSurface, aEdge_is_constrained_map);
  }

  template<class TM>
  vertex_descriptor halfedge_collapse_bk_compatibility(const halfedge_descriptor pq,
                                                       No_constrained_edge_map<TM>)
  {
    vertex_descriptor vd = CGAL::Euler::collapse_edge(edge(pq, mSurface), mSurface);
    return vd;
  }

  /// We wrap this test to avoid penalizing runtime when no constraints are present
  template<class AEdgeIsConstrainedMap>
  bool is_edge_adjacent_to_a_constrained_edge(const Profile& aProfile, AEdgeIsConstrainedMap)
  {
    return is_constrained(aProfile.v0()) && is_constrained(aProfile.v1());
  }

  template<class TM>
  bool is_edge_adjacent_to_a_constrained_edge(const halfedge_descriptor, No_constrained_edge_map<TM>)
  {
    return false;
  }

private:
  TM& mSurface;
  const ShouldStop& Should_stop;
  const VertexIndexMap& Vertex_index_map;
  const VertexPointMap& Vertex_point_map;
  const EdgeIndexMap& Edge_index_map;
  const EdgeIsConstrainedMap& Edge_is_constrained_map;
  const GetCost& Get_cost;
  const GetPlacement& Get_placement;
  VisitorT Visitor;
  bool m_has_border;

private:
  Edge_data_array mEdgeDataArray;

  boost::scoped_ptr<PQ> mPQ;

  size_type mInitialEdgeCount;
  size_type mCurrentEdgeCount;

  FT mcMaxDihedralAngleCos2;

  CGAL_SMS_DEBUG_CODE(unsigned mStep;)
};

template<class M, class SP, class VIM, class VPM,class EIM, class ECTM, class CF, class PF, class V>
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
EdgeCollapse(TM& aSurface,
             const ShouldStop& aShould_stop,
             const VertexIndexMap& aVertex_index_map,
             const VertexPointMap& aVertex_point_map,
             const EdgeIndexMap& aEdge_index_map,
             const EdgeIsConstrainedMap& aEdge_is_constrained_map,
             const GetCost& aGet_cost,
             const GetPlacement& aGet_placement,
             VisitorT aVisitor)
  :
    mSurface(aSurface),
    Should_stop(aShould_stop),
    Vertex_index_map(aVertex_index_map),
    Vertex_point_map(aVertex_point_map),
    Edge_index_map(aEdge_index_map),
    Edge_is_constrained_map(aEdge_is_constrained_map),
    Get_cost(aGet_cost),
    Get_placement(aGet_placement),
    Visitor(aVisitor),
    m_has_border(false)
{
  const FT cMaxDihedralAngleCos = std::cos(1.0 * CGAL_PI / 180.0);

  mcMaxDihedralAngleCos2 = cMaxDihedralAngleCos * cMaxDihedralAngleCos;

  for(halfedge_descriptor ed : halfedges(mSurface))
  {
    if(is_border(ed))
    {
      m_has_border = true;
      break;
    }
  }

  CGAL_SMS_TRACE(0, "EdgeCollapse of TM with " << (num_edges(aSurface)/2) << " edges");

  CGAL_SMS_DEBUG_CODE(mStep = 0;)

#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
  for(vertex_descriptor vd : vertices(mSurface))
    CGAL_SMS_TRACE(1, vertex_to_string(vd));

  for(halfedge_descriptor ed : halfedges(mSurface))
    CGAL_SMS_TRACE(1, edge_to_string(ed));
#endif
}

template<class M,class SP, class VIM, class VPM,class EIM,class ECTM, class CF,class PF,class V>
int
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
run()
{
  CGAL_SURF_SIMPL_TEST_assertion(mSurface.is_valid() && mSurface.is_pure_triangle());

  Visitor.OnStarted(mSurface);

  // First collect all candidate edges in a PQ
  collect();

  // Then proceed to collapse each edge in turn
  loop();

  CGAL_SMS_TRACE(0, "Finished: " << (mInitialEdgeCount - mCurrentEdgeCount) << " edges removed.");

  int r = int(mInitialEdgeCount - mCurrentEdgeCount);

  Visitor.OnFinished(mSurface);

  return r;
}

template<class M,class SP, class VIM, class VPM,class EIM,class ECTM, class CF,class PF,class V>
void
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
collect()
{
  CGAL_SMS_TRACE(0, "collecting edges...");

  // loop over all the _undirected_ edges in the surface putting them in the PQ

  Equal_3 equal_points = Traits().equal_3_object();

  size_type lSize = num_edges(mSurface);

  mInitialEdgeCount = mCurrentEdgeCount = static_cast<size_type>(
                                            std::distance(boost::begin(edges(mSurface)),
                                                          boost::end(edges(mSurface))));;

  mEdgeDataArray.reset(new Edge_data[lSize]);

  mPQ.reset(new PQ(lSize, Compare_cost(this), edge_id(this)));

  std::size_t id = 0;
  CGAL_SURF_SIMPL_TEST_assertion_code(size_type lInserted = 0);
  CGAL_SURF_SIMPL_TEST_assertion_code(size_type lNotInserted = 0);

  std::set<halfedge_descriptor> zero_length_edges;

  edge_iterator eb, ee;
  for(boost::tie(eb,ee) = edges(mSurface); eb!=ee; ++eb, id+=2)
  {
    halfedge_descriptor lEdge = halfedge(*eb, mSurface);

    if(is_constrained(lEdge))
      continue; // no not insert constrainted edges

    // AF CGAL_assertion(get_halfedge_id(lEdge) == id);
    // AF CGAL_assertion(get_halfedge_id(opposite(lEdge, mSurface)) == id+1);

    const Profile& lProfile = create_profile(lEdge);
    if(!equal_points(lProfile.p0(), lProfile.p1()))
    {
      Edge_data& lData = get_data(lEdge);

      lData.cost() = get_cost(lProfile);
      insert_in_PQ(lEdge, lData);

      Visitor.OnCollected(lProfile, lData.cost());

      CGAL_SURF_SIMPL_TEST_assertion_code(++lInserted);
    }
    else
    {
      zero_length_edges.insert(primary_edge(lEdge));
      CGAL_SURF_SIMPL_TEST_assertion_code(++lNotInserted);
    }

    CGAL_SMS_TRACE(2,edge_to_string(lEdge));
  }

  CGAL_SURF_SIMPL_TEST_assertion(lInserted + lNotInserted == mInitialEdgeCount);

  for(halfedge_descriptor hd : zero_length_edges)
  {
    const Profile& lProfile = create_profile(hd);

    if(!is_collapse_topologically_valid(lProfile))
      continue;

    // edges of length 0 removed no longer need to be treated
    if(lProfile.left_face_exists())
    {
      halfedge_descriptor lEdge_to_remove = is_constrained(lProfile.vL_v0()) ?
                                              primary_edge(lProfile.v1_vL()) :
                                              primary_edge(lProfile.vL_v0());
      zero_length_edges.erase(lEdge_to_remove);
      Edge_data& lData = get_data(lEdge_to_remove);

      if(lData.is_in_PQ())
      {
        CGAL_SMS_TRACE(2,"Removing E" << get(Edge_index_map,lEdge_to_remove) << " from PQ");
        remove_from_PQ(lEdge_to_remove, lData);
      }

      --mCurrentEdgeCount;
    }

    if(lProfile.right_face_exists())
    {
      halfedge_descriptor lEdge_to_remove = is_constrained(lProfile.vR_v1()) ?
                                              primary_edge(lProfile.v0_vR()) :
                                              primary_edge(lProfile.vR_v1());
      zero_length_edges.erase(lEdge_to_remove);
      Edge_data& lData = get_data(lEdge_to_remove);

      if(lData.is_in_PQ())
      {
        CGAL_SMS_TRACE(2,"Removing E" << get(Edge_index_map,lEdge_to_remove) << " from PQ");
        remove_from_PQ(lEdge_to_remove, lData);
      }

      --mCurrentEdgeCount;
    }

    --mCurrentEdgeCount;

    //the placement is trivial, it's always the point itself
    Placement_type lPlacement = lProfile.p0();
    vertex_descriptor rResult = halfedge_collapse_bk_compatibility(lProfile.v0_v1(), Edge_is_constrained_map);
    put(Vertex_point_map, rResult, *lPlacement);
    Visitor.OnCollapsed(lProfile, rResult);
  }

  CGAL_SMS_TRACE(0, "Initial edge count: " << mInitialEdgeCount);
}

template<class M,class SP, class VIM, class VPM,class EIM,class ECTM, class CF,class PF,class V>
void
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
loop()
{
  CGAL_SMS_TRACE(0, "Collapsing edges...");

  CGAL_SURF_SIMPL_TEST_assertion_code(size_type lloop_watchdog = 0);
  CGAL_SURF_SIMPL_TEST_assertion_code(size_type lNonCollapsableCount = 0);

  // Pops and processes each edge from the PQ

  boost::optional<halfedge_descriptor> lEdge;
#ifdef CGAL_SURF_SIMPL_INTERMEDIATE_STEPS_PRINTING
  int i_rm=0;
#endif

  while((lEdge = pop_from_PQ()))
  {
    CGAL_SURF_SIMPL_TEST_assertion(lloop_watchdog++ < mInitialEdgeCount);

    CGAL_SMS_TRACE(1, "Popped " << edge_to_string(*lEdge));
    CGAL_assertion(!is_constrained(*lEdge));

    const Profile& lProfile = create_profile(*lEdge);
    Cost_type lCost = get_data(*lEdge).cost();

    Visitor.OnSelected(lProfile, lCost, mInitialEdgeCount, mCurrentEdgeCount);

    if(lCost)
    {
      if(Should_stop(*lCost, lProfile, mInitialEdgeCount, mCurrentEdgeCount))
      {
        Visitor.OnStopConditionReached(lProfile);

        CGAL_SMS_TRACE(0, "Stop condition reached with InitialEdgeCount=" << mInitialEdgeCount
                            << " CurrentEdgeCount=" << mCurrentEdgeCount
                            << " Current Edge: " << edge_to_string(*lEdge));
        break;
      }

      if(is_collapse_topologically_valid(lProfile))
      {
        // The external function Get_new_vertex_point() is allowed to return an absent point if there is no way to place the vertex
        // satisfying its constraints. In that case the remaining vertex is simply left unmoved.
        Placement_type lPlacement = get_placement(lProfile);

        if(is_collapse_geometrically_valid(lProfile,lPlacement))
        {
#ifdef CGAL_SURF_SIMPL_INTERMEDIATE_STEPS_PRINTING
          std::cout << "step " << i_rm << " " << source(*lEdge, mSurface)->point() << " " << target(*lEdge, mSurface)->point() << "\n";
#endif
          collapse(lProfile, lPlacement);

#ifdef CGAL_SURF_SIMPL_INTERMEDIATE_STEPS_PRINTING
          std::stringstream sstr;
          sstr << "debug/P-";
          if(i_rm<10) sstr << "0";
          if(i_rm<100) sstr << "0";
          sstr <<  i_rm  << ".off";
          std::ofstream out(sstr.str().c_str());
          out << mSurface;
          ++i_rm;
#endif
        }
      }
      else
      {
        CGAL_SURF_SIMPL_TEST_assertion_code(lNonCollapsableCount++);

        Visitor.OnNonCollapsable(lProfile);

        CGAL_SMS_TRACE(1, edge_to_string(*lEdge) << " NOT Collapsable" );
      }
    }
    else
    {
      CGAL_SMS_TRACE(1,edge_to_string(*lEdge) << " uncomputable cost." );
    }
  }
}

template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
bool
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
is_border(const vertex_descriptor aV) const
{
  bool rR = false;

  for(halfedge_descriptor lEdge : halfedges_around_target(aV, mSurface))
  {
    if(is_edge_a_border(lEdge))
    {
      rR = true;
      break;
    }
  }

  return rR;
}

template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
bool
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
is_border_or_constrained(const vertex_descriptor aV) const
{
  for(halfedge_descriptor lEdge : halfedges_around_target(aV, mSurface))
  {
    if(is_edge_a_border(lEdge) || is_constrained(lEdge))
      return true;
  }

  return false;
}

template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
bool
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
is_constrained(const vertex_descriptor aV) const
{
  for(halfedge_descriptor lEdge : halfedges_around_target(aV, mSurface))
    if(is_constrained(lEdge))
      return true;

  return false;
}

// Some edges are NOT collapsable: doing so would break the topological consistency of the mesh.
// This function returns true if a edge 'p->q' can be collapsed.
//
// An edge p->q can be collapsed iff it satisfies the "link condition"
// (as described in the "Mesh Optimization" article of Hoppe et al (1993))
//
// The link condition is as follows: for every vertex 'k' adjacent to both 'p and 'q',
// "p,k,q" is a facet of the mesh.
//
template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
bool
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
is_collapse_topologically_valid(const Profile& aProfile)
{
  bool rR = true;

  CGAL_SMS_TRACE(3,"Testing topological collapsabilty of p_q=V" << get(Vertex_index_map,aProfile.v0()) << "(%" << degree(aProfile.v0(), mSurface) << ")"
                 << "->V" << get(Vertex_index_map,aProfile.v1()) << "(%" << degree(aProfile.v1(), mSurface) << ")");

  CGAL_SMS_TRACE(4, "is p_q border:" << aProfile.is_v0_v1_a_border());
  CGAL_SMS_TRACE(4, "is q_q border:" << aProfile.is_v1_v0_a_border());

  out_edge_iterator eb1, ee1;
  out_edge_iterator eb2, ee2;

  CGAL_SMS_TRACE(4, "  t=V" << (aProfile.left_face_exists() ? get(Vertex_index_map,aProfile.vL()) : -1)
                      << "(%" << (aProfile.left_face_exists() ? degree(aProfile.vL(), mSurface) : 0) << ")");
  CGAL_SMS_TRACE(4, "  b=V" << (aProfile.right_face_exists() ? get(Vertex_index_map,aProfile.vR()) : -1)
                      << "(%" << (aProfile.right_face_exists() ? degree(aProfile.vR(), mSurface) :0) << ")");

  // Simple tests handling the case of non-manifold situations at a vertex or edge (pinching)
  // (even if we advertise one should not use a surface mesh with such features)
  if(aProfile.left_face_exists())
  {
    if(CGAL::is_border(opposite(aProfile.v1_vL(), mSurface), mSurface) &&
       CGAL::is_border(opposite(aProfile.vL_v0(), mSurface), mSurface))
      return false;

    if(aProfile.right_face_exists() &&
       CGAL::is_border(opposite(aProfile.vR_v1(), mSurface), mSurface) &&
       CGAL::is_border(opposite(aProfile.v0_vR(), mSurface), mSurface))
      return false;
  }
  else
  {
    if(aProfile.right_face_exists())
    {
      if(CGAL::is_border(opposite(aProfile.vR_v1(), mSurface), mSurface) &&
         CGAL::is_border(opposite(aProfile.v0_vR(), mSurface), mSurface))
        return false;
    }
    else
      return false;
  }

  // The following loop checks the link condition for v0_v1.
  // Specifically, that for every vertex 'k' adjacent to both 'p and 'q', 'pkq' is a face of the mesh.
  //
  for(boost::tie(eb1,ee1) = halfedges_around_source(aProfile.v0(), mSurface); rR && eb1 != ee1; ++eb1)
  {
    halfedge_descriptor v0_k = *eb1;

    if(v0_k != aProfile.v0_v1())
    {
      vertex_descriptor k = target(v0_k, mSurface);

      for(boost::tie(eb2,ee2) = halfedges_around_source(k, mSurface); rR && eb2 != ee2; ++eb2)
      {
        halfedge_descriptor k_v1 = *eb2;

        if(target(k_v1, mSurface) == aProfile.v1())
        {
          // At this point we know p-q-k are connected and we need to determine if this triangle is a face of the mesh.
          //
          // Since the mesh is known to be triangular there are at most two faces sharing the edge p-q.
          //
          // If p->q is NOT a border edge, the top face is p->q->t where t is target(next(p->q))
          // If q->p is NOT a border edge, the bottom face is q->p->b where b is target(next(q->p))
          //
          // If k is either t or b then p-q-k *might* be a face of the mesh. It won't be if k==t but p->q is border
          // or k==b but q->b is a border (because in that case even though there exists triangles p->q->t (or q->p->b)
          // they are holes, not faces)
          //

          bool lIsFace = (aProfile.vL() == k && aProfile.left_face_exists()) ||
                         (aProfile.vR() == k && aProfile.right_face_exists());

          CGAL_SURF_SIMPL_TEST_assertion_code(
            if(lIsFace)
            {
              // Is k_v1 the halfedge bounding the face 'k-v1-v0'?
              if(!is_border(k_v1) && target(next(k_v1, mSurface), mSurface) == aProfile.v0())
              {
                CGAL_SURF_SIMPL_TEST_assertion(!is_border(k_v1));
                CGAL_SURF_SIMPL_TEST_assertion(target(k_v1, mSurface) == aProfile.v1());
                CGAL_SURF_SIMPL_TEST_assertion(target(next(k_v1, mSurface), mSurface) == aProfile.v0());
                CGAL_SURF_SIMPL_TEST_assertion(target(next(next(k_v1, mSurface), mSurface), mSurface) == k);
              }
              else // or is it the opposite?
              {
                halfedge_descriptor v1_k = opposite(k_v1, mSurface);
                CGAL_SURF_SIMPL_TEST_assertion(!is_border(v1_k));
                CGAL_SURF_SIMPL_TEST_assertion(target(v1_k, mSurface) == k);
                CGAL_SURF_SIMPL_TEST_assertion(target(next(v1_k, mSurface), mSurface) == aProfile.v0());
                CGAL_SURF_SIMPL_TEST_assertion(target(next(next(v1_k, mSurface), mSurface), mSurface) == aProfile.v1());
              }
            }
          );

          if(!lIsFace)
          {
            CGAL_SMS_TRACE(3,"  k=V" << get(Vertex_index_map,k) << " IS NOT in a face with p-q. NON-COLLAPSABLE edge.");
            rR = false;
            break;
          }
          else
          {
            CGAL_SMS_TRACE(4,"  k=V" << get(Vertex_index_map,k) << " is in a face with p-q");
          }
        }
      }
    }
  }

  if(rR)
  {
    /// ensure two constrained edges cannot get merged
    if(is_edge_adjacent_to_a_constrained_edge(aProfile, Edge_is_constrained_map))
      return false;

    if(aProfile.is_v0_v1_a_border())
    {
      if(is_open_triangle(aProfile.v0_v1()))
      {
        rR = false;
        CGAL_SMS_TRACE(3,"  p-q belongs to an open triangle. NON-COLLAPSABLE edge.");
      }
    }
    else if(aProfile.is_v1_v0_a_border())
    {
      if(is_open_triangle(aProfile.v1_v0()))
      {
        rR = false;
        CGAL_SMS_TRACE(3,"  p-q belongs to an open triangle. NON-COLLAPSABLE edge.");
      }
    }
    else
    {
      if(is_border(aProfile.v0()) && is_border(aProfile.v1()))
      {
        rR = false;
        CGAL_SMS_TRACE(3,"  both p and q are boundary vertices but p-q is not. NON-COLLAPSABLE edge.");
      }
      else
      {
        bool lTetra = is_tetrahedron(aProfile.v0_v1());

        //CGAL_SURF_SIMPL_TEST_assertion(lTetra == mSurface.is_tetrahedron(aProfile.v0_v1()));

        if(lTetra)
        {
          rR = false;
          CGAL_SMS_TRACE(3,"  p-q belongs to a tetrahedron. NON-COLLAPSABLE edge.");
        }

        if(next(aProfile.v0_v1(), mSurface) == opposite(prev(aProfile.v1_v0(), mSurface), mSurface) &&
           prev(aProfile.v0_v1(), mSurface) == opposite(next(aProfile.v1_v0(), mSurface), mSurface))
        {
          CGAL_SMS_TRACE(3,"  degenerate volume.");
          return false;
        }
      }
    }
  }

  return rR;
}

template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
bool
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
is_tetrahedron(const halfedge_descriptor h1)
{
  return CGAL::is_tetrahedron(h1, mSurface);
}

template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
bool
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
is_open_triangle(const halfedge_descriptor h1)
{
  bool rR = false;

  halfedge_descriptor h2 = next(h1, mSurface);
  halfedge_descriptor h3 = next(h2, mSurface);

  // First check if it is a triangle
  if(next(h3, mSurface) == h1)
  {
    // Now check if it is open
    CGAL_SMS_TRACE(4,"  p-q is a border edge... checking E" << get(Edge_index_map, h2) << " and E" << get(Edge_index_map,h3));

    rR = is_border(h2) && is_border(h3);

    CGAL_SURF_SIMPL_TEST_assertion(rR == (is_border(h1) &&
                                          is_border(next(h1, mSurface)) &&
                                          is_border(next(next(h1, mSurface), mSurface))));
  }

  return rR;
}

// Given triangles 'p0,p1,p2' and 'p0,p2,p3', both shared along edge 'v0-v2',
// determine if they are geometrically valid: that is, the ratio of their
// respective areas is no greater than a max value and the internal
// dihedral angle formed by their supporting planes is no greater than
// a given threshold
template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
bool
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
are_shared_triangles_valid(const Point& p0, const Point& p1, const Point& p2, const Point& p3) const
{
  bool rR = false;

  Vector e01 = Traits().construct_vector_3_object()(p0, p1);
  Vector e02 = Traits().construct_vector_3_object()(p0, p2);
  Vector e03 = Traits().construct_vector_3_object()(p0, p3);

  Vector n012 = Traits().construct_cross_product_vector_3_object()(e01, e02);
  Vector n023 = Traits().construct_cross_product_vector_3_object()(e02, e03);

  FT l012 = Traits().compute_scalar_product_3_object()(n012, n012);
  FT l023 = Traits().compute_scalar_product_3_object()(n023, n023);

  FT larger  = (std::max)(l012, l023);
  FT smaller = (std::min)(l012, l023);

  const double cMaxAreaRatio = 1e8;

  CGAL_SMS_TRACE(4,"    Testing validity of shared triangles:"
                 << "\n      p0=" << xyz_to_string(p0) << "\n      p1=" << xyz_to_string(p1) << "\n      p2=" << xyz_to_string(p2) << "\n      p3=" << xyz_to_string(p3)
                 << "\n      e01=" << xyz_to_string(e01) << "\n      e02=" << xyz_to_string(e02) << "\n      e03=" << xyz_to_string(e03)
                 << "\n      n012=" << xyz_to_string(n012) << "\n      n023=" << xyz_to_string(n023)
                 << "\n      l012=" << n_to_string(l012) << "\n      l023=" << n_to_string(l023));

  if(larger < cMaxAreaRatio * smaller)
  {
    FT l0123 = Traits().compute_scalar_product_3_object()(n012, n023);
    CGAL_SMS_TRACE(4,"\n      l0123=" << n_to_string(l0123));

    if(CGAL_NTS is_positive(l0123))
    {
      rR = true;
    }
    else
    {
      CGAL_SMS_TRACE(4,"\n      lhs: " << n_to_string((l0123 * l0123) / (l012 * l023)) << " <= rhs: " << mcMaxDihedralAngleCos2);

      if((l0123 * l0123) <= mcMaxDihedralAngleCos2 * (l012 * l023))
      {
        rR = true;
      }
    }
  }

  return rR;
}

// Returns the directed halfedge connecting v0 to v1, if exists.
template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
typename EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::halfedge_descriptor
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
find_connection(const vertex_descriptor v0,
                const vertex_descriptor v1) const
{
  for(halfedge_descriptor out : halfedges_around_source(v0, mSurface))
  {
    if(target(out, mSurface) == v1)
      return out;
  }

  return GraphTraits::null_halfedge();
}

// Given the edge 'e' around the link for the collapsinge edge "v0-v1", finds the vertex that makes a triangle adjacent to 'e' but exterior to the link (i.e not containing v0 nor v1)
// If 'e' is a null handle OR 'e' is a border edge, there is no such triangle and a null handle is returned.
template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
typename EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::vertex_descriptor
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
find_exterior_link_triangle_3rd_vertex(const halfedge_descriptor e,
                                       const vertex_descriptor v0,
                                       const vertex_descriptor v1) const
{
  vertex_descriptor r;

  if(handle_assigned(e))
  {
    vertex_descriptor ra = target(next(e, mSurface), mSurface);
    vertex_descriptor rb = source(prev(e, mSurface), mSurface);

    if(ra == rb && ra != v0 && ra != v1)
    {
      r = ra;
    }
    else
    {
      ra = target(next(opposite(e, mSurface), mSurface), mSurface);
      rb = source(prev(opposite(e, mSurface), mSurface), mSurface);

      if(ra == rb && ra != v0 && ra != v1)
      {
        r = ra;
      }
    }
  }

  return r;
}

// A collapse is geometrically valid if, in the resulting local mesh no two adjacent triangles form an internal dihedral angle
// greater than a fixed threshold (i.e. triangles do not "fold" into each other)
//
template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
bool
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
is_collapse_geometrically_valid(const Profile& aProfile, Placement_type k0)
{
  bool rR = true;

  CGAL_SMS_TRACE(3,"Testing geometrical collapsabilty of v0-v1=E" << get(Edge_index_map, aProfile.v0_v1()));
  if(k0)
  {
    // Use the current link to extract all local triangles incident to 'vx' in the collapsed mesh (which at this point doesn't exist yet)
    //
    typedef typename Profile::vertex_descriptor_vector::const_iterator link_iterator;
    link_iterator linkb = aProfile.link().begin();
    link_iterator linke = aProfile.link().end();
    link_iterator linkl = std::prev(linke);

    for(link_iterator l=linkb; l!=linke && rR; ++l)
    {
      link_iterator pv = (l == linkb ? linkl : std::prev(l));
      link_iterator nx = (l == linkl ? linkb : std::next(l));

      // k0,k1 and k3 are three consecutive vertices along the link.
      vertex_descriptor k1 = *pv;
      vertex_descriptor k2 = *l;
      vertex_descriptor k3 = *nx;

      CGAL_SMS_TRACE(4, "  Screening link vertices k1=V" << get(Vertex_index_map, k1)
                          << " k2=V" << get(Vertex_index_map, k2)
                          << " k3=V" << get(Vertex_index_map, k3));

      halfedge_descriptor e12 = find_connection(k1,k2);
      halfedge_descriptor e23 = k3 != k1 ? find_connection(k2,k3) : GraphTraits::null_halfedge();

      // If 'k1-k2-k3' are connected there will be two adjacent triangles 'k0,k1,k2' and 'k0,k2,k3' after the collapse.
      if(handle_assigned(e12) && handle_assigned(e23))
      {
        CGAL_SMS_TRACE(4,"    Link triangles shared");

        if(!are_shared_triangles_valid(*k0, get_point(k1), get_point(k2), get_point(k3)))
        {
          CGAL_SMS_TRACE(3, "    Triangles VX-V" << get(Vertex_index_map, k1)
                              << "-V" << get(Vertex_index_map, k2)
                              << " and VX-V" << get(Vertex_index_map, k3)
                              << " are not geometrically valid. Collapse rejected");
          rR = false;
        }
      }

      if(rR)
      {
        // Also check the triangles 'k0,k1,k2' and it's adjacent along e12: 'k4,k2,k1', if exist
        vertex_descriptor k4 = find_exterior_link_triangle_3rd_vertex(e12,aProfile.v0(), aProfile.v1());

        // There is indeed a triangle shared along e12
        if(handle_assigned(k4))
        {
          CGAL_SMS_TRACE(4, "    Found exterior link triangle shared along E" << get(Edge_index_map, e12)
                              << " with third vertex: V" << get(Vertex_index_map, k4));

          if(!are_shared_triangles_valid(get_point(k1), get_point(k4), get_point(k2), *k0))
          {
            CGAL_SMS_TRACE(3, "    Triangles V" << get(Vertex_index_map, k1)
                                << "-V" << get(Vertex_index_map, k4) << " and V"
                                << get(Vertex_index_map, k2) << "-VX are not geometrically valid. Collapse rejected");
            rR = false;
          }
        }
      }

      if(rR)
      {
        // And finally, check the triangles 'k0,k2,k3' and it's adjacent e23: 'k5,k3,k2' if exist
        vertex_descriptor k5 = find_exterior_link_triangle_3rd_vertex(e23, aProfile.v0(), aProfile.v1());

        // There is indeed a triangle shared along e12
        if(handle_assigned(k5))
        {
          CGAL_SMS_TRACE(4, "    Found exterior link triangle shared along E" << get(Edge_index_map, e23)
                             << " with third vertex: V" << get(Vertex_index_map, k5));

          if(!are_shared_triangles_valid(get_point(k2), get_point(k5), get_point(k3), *k0))
          {
            CGAL_SMS_TRACE(3, "    Triangles V" << get(Vertex_index_map, k2)
                                << "-V" << get(Vertex_index_map, k5)
                                << " and V" << get(Vertex_index_map, k3) << "-VX are not geometrically valid. Collapse rejected");
            rR = false;
          }
        }
      }
    }
  }

  return rR;
}

template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
void
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
collapse(const Profile& aProfile,
         Placement_type aPlacement)
{
  CGAL_SMS_TRACE(1,"S" << mStep << ". Collapsing " << edge_to_string(aProfile.v0_v1()));

  vertex_descriptor rResult;

  Visitor.OnCollapsing(aProfile,aPlacement);

  --mCurrentEdgeCount;

  CGAL_SURF_SIMPL_TEST_assertion_code(
    size_type lResultingVertexCount = mSurface.size_of_vertices();
    size_type lResultingEdgeCount = mSurface.size_of_halfedges() / 2;
  );

  // If the top/bottom facets exists, they are removed and the edges v0vt and Q-B along with them.
  // In that case their corresponding pairs must be pop off the queue
  if(aProfile.left_face_exists())
  {
    halfedge_descriptor lV0VL = primary_edge(aProfile.vL_v0());
    if(is_constrained(lV0VL)) // make sure a constrained edge will not disappear
      lV0VL=primary_edge(aProfile.v1_vL());

    CGAL_SMS_TRACE(3, "V0VL E" << get(Edge_index_map,lV0VL)
                       << "(V" <<  get(Vertex_index_map, source(lV0VL, mSurface))
                       << "->V" << get(Vertex_index_map,target(lV0VL, mSurface)) << ")");

    Edge_data& lData = get_data(lV0VL);
    if(lData.is_in_PQ())
    {
      CGAL_SMS_TRACE(2,"Removing E" << get(Edge_index_map, lV0VL) << " from PQ");
      remove_from_PQ(lV0VL, lData);
    }

    --mCurrentEdgeCount;
    CGAL_SURF_SIMPL_TEST_assertion_code(--lResultingEdgeCount);
  }

  if(aProfile.right_face_exists())
  {
    halfedge_descriptor lVRV1 = primary_edge(aProfile.vR_v1());
    if(is_constrained(lVRV1)) //make sure a constrained edge will not disappear
      lVRV1=primary_edge(aProfile.v0_vR());

    CGAL_SMS_TRACE(3, "V1VRE" << get(Edge_index_map, lVRV1)
                        << "(V" <<  get(Vertex_index_map, source(lVRV1, mSurface))
                        << "->V" << get(Vertex_index_map, target(lVRV1, mSurface)) << ")");

    Edge_data& lData = get_data(lVRV1);
    if(lData.is_in_PQ())
    {
      CGAL_SMS_TRACE(2,"Removing E" << get(Edge_index_map, lVRV1) << " from PQ");
      remove_from_PQ(lVRV1, lData);
    }

    --mCurrentEdgeCount;
    CGAL_SURF_SIMPL_TEST_assertion_code(--lResultingEdgeCount);
  }

  CGAL_SMS_TRACE(1, "Removing:\n  v0v1: E" << get(Edge_index_map, aProfile.v0_v1())
                      << "(V" << get(Vertex_index_map, aProfile.v0())
                      << "->V" << get(Vertex_index_map, aProfile.v1()) << ")");

  // Perform the actuall collapse.
  // This is an external function.
  // It's REQUIRED to remove ONLY 1 vertex (P or Q) and edges PQ,PT and QB
  // (PT and QB are removed if they are not null).
  // All other edges must be kept.
  // All directed edges incident to vertex removed are relink to the vertex kept.
  rResult = halfedge_collapse_bk_compatibility(aProfile.v0_v1(), Edge_is_constrained_map);

  CGAL_SURF_SIMPL_TEST_assertion_code(--lResultingEdgeCount);
  CGAL_SURF_SIMPL_TEST_assertion_code(--lResultingVertexCount);

  CGAL_SURF_SIMPL_TEST_assertion(lResultingEdgeCount * 2 == mSurface.size_of_halfedges());
  CGAL_SURF_SIMPL_TEST_assertion(lResultingVertexCount == mSurface.size_of_vertices());
  CGAL_SURF_SIMPL_TEST_assertion(mSurface.is_valid() && mSurface.is_pure_triangle());

  CGAL_SMS_TRACE(1,"V" << get(Vertex_index_map,rResult) << " kept.");

#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
  for(halfedge_descriptor hd : halfedges_around_source(rResult, mSurface))
    CGAL_SMS_TRACE(2, edge_to_string(hd));
#endif

  if(aPlacement)
  {
    CGAL_SMS_TRACE(1,"New vertex point: " << xyz_to_string(*aPlacement));
    put(Vertex_point_map,rResult,*aPlacement);
  }

  Visitor.OnCollapsed(aProfile,rResult);

  update_neighbors(rResult);

  CGAL_SMS_DEBUG_CODE(++mStep;)
}

template<class M, class SP, class VIM, class VPM, class EIM, class ECTM, class CF, class PF, class V>
void
EdgeCollapse<M,SP,VIM,VPM,EIM,ECTM,CF,PF,V>::
update_neighbors(const vertex_descriptor aKeptV)
{
  CGAL_SMS_TRACE(3,"Updating cost of neighboring edges...");

  // (A) collect all edges to update their cost: all those around each vertex adjacent to the vertex kept
  typedef std::set<halfedge_descriptor, Compare_id>                       edges;

  edges lToUpdate(Compare_id(this));
  edges lToInsert(Compare_id(this));

  // (A.1) loop around all vertices adjacent to the vertex kept
  for(halfedge_descriptor lEdge1 : halfedges_around_target(aKeptV, mSurface))
  {
    vertex_descriptor lAdj_k = source(lEdge1, mSurface);

    // (A.2) loop around all edges incident on each adjacent vertex
    for(halfedge_descriptor lEdge2 : halfedges_around_target(lAdj_k, mSurface))
    {
      lEdge2 = primary_edge(lEdge2);

      Edge_data& lData2 = get_data(lEdge2);
      CGAL_SMS_TRACE(4,"Inedge around V" << get(Vertex_index_map, lAdj_k) << edge_to_string(lEdge2));

      // Only edges still in the PQ needs to be updated, the other needs to be re-inserted
      if(lData2.is_in_PQ())
        lToUpdate.insert(lEdge2);
      else
        lToInsert.insert(lEdge2);
    }
  }

  // (B) Proceed to update the costs.
  for(halfedge_descriptor lEdge : lToUpdate)
  {
    Edge_data& lData = get_data(lEdge);
    const Profile& lProfile = create_profile(lEdge);
    lData.cost() = get_cost(lProfile);

    CGAL_SMS_TRACE(3, edge_to_string(lEdge) << " updated in the PQ");

    update_in_PQ(lEdge, lData);
  }

  // (C) Insert ignored edges
  //
  // I think that this should be done for edges eliminated because of the geometric criteria
  // and not the topological one.However maintaining such a set might be more expensive
  // and hard to be safe ...
  for(halfedge_descriptor lEdge : lToInsert)
  {
    if(is_constrained(lEdge))
      continue; //do not insert constrained edges

    Edge_data& lData = get_data(lEdge);
    const Profile& lProfile = create_profile(lEdge);
    lData.cost() = get_cost(lProfile);

    CGAL_SMS_TRACE(3, edge_to_string(lEdge) << " re-inserted in the PQ");
    insert_in_PQ(lEdge, lData);
  }
}

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_H //
