// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
namespace internal {

BOOST_MPL_HAS_XXX_TRAIT_DEF(Update_tag)

template <typename Cost_oracle,
          bool has_Update_tag = has_Update_tag<Cost_oracle>::value>
struct Oracles_require_updates :
  public CGAL::Boolean_tag<Cost_oracle::Update_tag::value>
  // when Mesh_domain has the nested type Has_features
{ };

template <typename Cost_oracle>
struct Oracles_require_updates<Cost_oracle, false> : public CGAL::Tag_false { };

template <typename EdgeCollapse,
          bool Tag = internal::Oracles_require_updates<typename EdgeCollapse::Get_cost>::value>
struct Oracles_initializer
{
  Oracles_initializer(const EdgeCollapse& e) : e(e) { }
  void operator()() const { }
  const EdgeCollapse& e;
};

template <typename EdgeCollapse>
struct Oracles_initializer<EdgeCollapse, true>
{
  Oracles_initializer(const EdgeCollapse& e) : e(e) { }
  void operator()() const { e.get_cost().initialize(e.mesh(), e.vpm(), e.geom_traits()); }
  const EdgeCollapse& e;
};

template <typename EdgeCollapse,
          bool Tag = internal::Oracles_require_updates<typename EdgeCollapse::Get_cost>::value>
struct After_collapse_oracles_updater
{
  After_collapse_oracles_updater(const EdgeCollapse& e) : e(e) { }
  void operator()(const typename EdgeCollapse::Profile& /*profile*/,
                  const typename boost::graph_traits<typename EdgeCollapse::Triangle_mesh>::vertex_descriptor /*v_kept*/) const { };
  const EdgeCollapse& e;
};

template <typename EdgeCollapse>
struct After_collapse_oracles_updater<EdgeCollapse, true>
{
  After_collapse_oracles_updater(const EdgeCollapse& e) : e(e) { }
  void operator()(const typename EdgeCollapse::Profile& profile,
                  const typename boost::graph_traits<typename EdgeCollapse::Triangle_mesh>::vertex_descriptor v_kept) const {
    e.get_cost().update_after_collapse(profile, v_kept);
  }
  const EdgeCollapse& e;
};

} // namespace internal

// Implementation of the vertex-pair collapse triangulated surface mesh simplification algorithm
template<class TM_,
         class GeomTraits_,
         class ShouldStop_,
         class VertexIndexMap_,
         class VertexPointMap_,
         class HalfedgeIndexMap_,
         class EdgeIsConstrainedMap_,
         class GetCost_,
         class GetPlacement_,
         class VisitorT_>
class EdgeCollapse
{
  typedef EdgeCollapse                                                    Self;

public:
  typedef TM_                                                             Triangle_mesh;
  typedef GeomTraits_                                                     Geom_traits;
  typedef ShouldStop_                                                     Should_stop;
  typedef VertexIndexMap_                                                 Vertex_index_map;
  typedef VertexPointMap_                                                 Vertex_point_map;
  typedef HalfedgeIndexMap_                                               Halfedge_index_map;
  typedef EdgeIsConstrainedMap_                                           Edge_is_constrained_map;
  typedef GetCost_                                                        Get_cost;
  typedef GetPlacement_                                                   Get_placement;
  typedef VisitorT_                                                       Visitor;

  typedef Edge_profile<Triangle_mesh, Vertex_point_map, Geom_traits>      Profile;

  typedef boost::graph_traits<Triangle_mesh>                              Graph_traits;
  typedef typename Graph_traits::vertex_descriptor                        vertex_descriptor;
  typedef typename Graph_traits::vertex_iterator                          vertex_iterator;
  typedef typename Graph_traits::halfedge_descriptor                      halfedge_descriptor;
  typedef typename Graph_traits::halfedge_iterator                        halfedge_iterator;
  typedef typename Graph_traits::edge_descriptor                          edge_descriptor;
  typedef typename Graph_traits::edge_iterator                            edge_iterator;
  typedef CGAL::Halfedge_around_source_iterator<Triangle_mesh>            out_edge_iterator;
  typedef typename Graph_traits::edges_size_type                          size_type;

  typedef typename boost::property_traits<Vertex_point_map>::value_type   Point;

  typedef typename Geom_traits::FT                                        FT;
  typedef typename Geom_traits::Vector_3                                  Vector;
  typedef typename Geom_traits::Equal_3                                   Equal_3;

  typedef boost::optional<FT>                                             Cost_type;
  typedef boost::optional<Point>                                          Placement_type;

  struct Compare_id
  {
    Compare_id() : m_algorithm(0) {}
    Compare_id(const Self* algorithm) : m_algorithm(algorithm) {}

    bool operator()(const halfedge_descriptor a, const halfedge_descriptor b) const {
      return m_algorithm->get_halfedge_id(a) < m_algorithm->get_halfedge_id(b);
    }

    const Self* m_algorithm;
  };

  struct Compare_cost
  {
    Compare_cost() : m_algorithm(0) {}
    Compare_cost(const Self* algorithm) : m_algorithm(algorithm) {}

    bool operator()(const halfedge_descriptor a, const halfedge_descriptor b) const
    {
      // NOTE: A cost is a boost::optional<> value.
      // Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
      // In consequence, edges with undefined costs will be promoted to the top of the priority queue and poped out first.
      return m_algorithm->get_data(a).cost() < m_algorithm->get_data(b).cost();
    }

    const Self* m_algorithm;
  };

  struct edge_id
    : boost::put_get_helper<size_type, edge_id>
  {
    typedef boost::readable_property_map_tag category;
    typedef size_type                        value_type;
    typedef size_type                        reference;
    typedef halfedge_descriptor              key_type;

    edge_id() : m_algorithm(0) {}
    edge_id(const Self* algorithm) : m_algorithm(algorithm) {}

    size_type operator[](const halfedge_descriptor e) const { return m_algorithm->get_edge_id(e); }

    const Self* m_algorithm;
  };

  typedef Modifiable_priority_queue<halfedge_descriptor, Compare_cost, edge_id>     PQ;
  typedef typename PQ::handle                                                       PQ_handle;

  // An Edge_data is associated with EVERY _ edge in the mesh (collapsable or not).
  // It relates the edge with the PQ-handle needed to update the priority queue
  // It also relates the edge with a policy-based cache
  class Edge_data
  {
  public :
    Edge_data() : m_PQ_h() {}

    const Cost_type& cost() const { return m_cost; }
    Cost_type& cost() { return m_cost; }

    PQ_handle queue_handle() const { return m_PQ_h;}
    bool is_in_PQ() const { return m_PQ_h != PQ::null_handle(); }
    void set_PQ_handle(PQ_handle h) { m_PQ_h = h; }
    void reset_queue_handle() { m_PQ_h = PQ::null_handle(); }

  private:
    Cost_type m_cost;
    PQ_handle m_PQ_h;
  };

  typedef Edge_data*                                                                Edge_data_ptr;
  typedef boost::scoped_array<Edge_data>                                            Edge_data_array;

public:
  EdgeCollapse(Triangle_mesh& tmesh,
               const Geom_traits& traits,
               const Should_stop& should_stop,
               const Vertex_index_map& vim,
               const Vertex_point_map& vpm,
               const Halfedge_index_map& him,
               const Edge_is_constrained_map& ecm,
               const Get_cost& aGetCost,
               const Get_placement& aGetPlacement,
               Visitor visitor);

  int run();

public:
  const Triangle_mesh& mesh() const { return m_tm; }
  const Geom_traits& geom_traits() const { return m_traits; }
  const Get_cost& get_cost() const { return m_get_cost; }
  const Vertex_point_map& vpm() const { return m_vpm; }

private:
  void collect();
  void loop();

  bool is_collapse_topologically_valid(const Profile& profile);
  bool is_tetrahedron(const halfedge_descriptor h);
  bool is_open_triangle(const halfedge_descriptor h1);
  bool is_collapse_geometrically_valid(const Profile& profile, Placement_type placement);
  void collapse(const Profile& profile, Placement_type placement);
  void update_neighbors(const vertex_descriptor v_kept);

  Profile create_profile(const halfedge_descriptor h) {
    return Profile(h, m_tm, m_traits, m_vim, m_vpm, m_him, m_has_border);
  }

  size_type get_halfedge_id(const halfedge_descriptor h) const { return get(m_him, h); }
  size_type get_edge_id(const halfedge_descriptor h) const { return get_halfedge_id(h) / 2; }

  bool is_primary_edge(const halfedge_descriptor h) const { return (get_halfedge_id(h) % 2) == 0; }
  halfedge_descriptor primary_edge(const halfedge_descriptor h) {
    return is_primary_edge(h) ? h : opposite(h, m_tm);
  }

  bool is_constrained(const vertex_descriptor v) const;
  bool is_constrained(const halfedge_descriptor h) const {
    return get(m_ecm, edge(h, m_tm));
  }

  bool is_border_or_constrained(const vertex_descriptor v) const;
  bool is_border_or_constrained(const halfedge_descriptor h) const { return is_border(h) || is_constrained(h); }

  bool are_shared_triangles_valid(const Point& p0, const Point& p1, const Point& p2, const Point& p3) const;

  halfedge_descriptor find_connection(const vertex_descriptor v0, const vertex_descriptor v1) const;
  vertex_descriptor find_exterior_link_triangle_3rd_vertex(const halfedge_descriptor e, const vertex_descriptor v0, const vertex_descriptor v1) const;

  Edge_data& get_data(const halfedge_descriptor h) const
  {
    CGAL_assertion(is_primary_edge(h));
    return m_edge_data[get_edge_id(h)];
  }

  typename boost::property_traits<Vertex_point_map>::reference
  get_point(const vertex_descriptor v) const { return get(m_vpm, v); }

  boost::tuple<vertex_descriptor, vertex_descriptor> get_vertices(const halfedge_descriptor h) const
  {
    vertex_descriptor p, q;
    p = source(h, m_tm);
    q = target(h, m_tm);
    return boost::make_tuple(p, q);
  }

  std::string vertex_to_string(const vertex_descriptor v) const
  {
    const Point& p = get_point(v);
    return boost::str(boost::format("[V%1%:%2%]") % get(m_vim,v) % xyz_to_string(p));
  }

  std::string edge_to_string(const halfedge_descriptor h) const
  {
    vertex_descriptor p, q;
    boost::tie(p,q) = get_vertices(h);
    return boost::str(boost::format("{E%1% %2%->%3%}%4%") % get_edge_id(h) % vertex_to_string(p) % vertex_to_string(q) % (is_border(h, m_tm) ? " (BORDER)" : (is_border(opposite(h, m_tm), m_tm) ? " (~BORDER)": "")));
  }

  Cost_type cost(const Profile& profile) const {
    return m_get_cost(profile, get_placement(profile));
  }

  Placement_type get_placement(const Profile& profile) const {
    return m_get_placement(profile);
  }

  void insert_in_PQ(const halfedge_descriptor h, Edge_data& data)
  {
    CGAL_assertion(is_primary_edge(h));
    CGAL_expensive_assertion(!data.is_in_PQ());
    CGAL_expensive_assertion(!mPQ->contains(h));

    data.set_PQ_handle(mPQ->push(h));

    CGAL_expensive_assertion(data.is_in_PQ());
    CGAL_expensive_assertion(mPQ->contains(h));
  }

  void update_in_PQ(const halfedge_descriptor h, Edge_data& data)
  {
    CGAL_assertion(is_primary_edge(h));
    CGAL_expensive_assertion(data.is_in_PQ());
    CGAL_expensive_assertion(mPQ->contains(h));

    data.set_PQ_handle(mPQ->update(h, data.queue_handle()));

    CGAL_assertion(data.is_in_PQ());
    CGAL_expensive_assertion(mPQ->contains(h));
  }

  void remove_from_PQ(const halfedge_descriptor h, Edge_data& data)
  {
    CGAL_assertion(is_primary_edge(h));
    CGAL_expensive_assertion(data.is_in_PQ());
    CGAL_expensive_assertion(mPQ->contains(h));

    data.set_PQ_handle(mPQ->erase(h, data.queue_handle()));

    CGAL_expensive_assertion(!data.is_in_PQ());
    CGAL_expensive_assertion(!mPQ->contains(h));
  }

  boost::optional<halfedge_descriptor> pop_from_PQ()
  {
    boost::optional<halfedge_descriptor> opt_h = mPQ->extract_top();
    if(opt_h)
    {
      CGAL_assertion(is_primary_edge(*opt_h));
      CGAL_expensive_assertion(get_data(*opt_h).is_in_PQ());

      get_data(*opt_h).reset_queue_handle();

      CGAL_expensive_assertion(!get_data(*opt_h).is_in_PQ());
      CGAL_expensive_assertion(!mPQ->contains(*opt_h));
    }
    return opt_h;
  }

  /// Functions to ensure the backward compatibility before addition of the constrained edge map
  template<class IsConstrainedMap_>
  vertex_descriptor halfedge_collapse_bk_compatibility(const halfedge_descriptor h,
                                                       const IsConstrainedMap_& ecm)
  {
    return CGAL::Euler::collapse_edge(edge(h, m_tm), m_tm, ecm);
  }
  template<class TM>
  vertex_descriptor halfedge_collapse_bk_compatibility(const halfedge_descriptor h,
                                                       const No_constrained_edge_map<TM>&)
  {
    return CGAL::Euler::collapse_edge(edge(h, m_tm), m_tm);
  }

  /// We wrap this test to avoid penalizing runtime when no constraints are present
  template<class ECM_>
  bool is_edge_adjacent_to_a_constrained_edge(const Profile& profile, const ECM_&)
  {
    return is_constrained(profile.v0()) && is_constrained(profile.v1());
  }
  template<class TM>
  bool is_edge_adjacent_to_a_constrained_edge(const halfedge_descriptor, const No_constrained_edge_map<TM>&)
  {
    return false;
  }

private:
  Triangle_mesh& m_tm;
  const Geom_traits& m_traits;
  const Should_stop& m_should_stop;
  const Vertex_index_map& m_vim;
  const Vertex_point_map& m_vpm;
  const Halfedge_index_map& m_him;
  const Edge_is_constrained_map& m_ecm;
  const Get_cost& m_get_cost;
  const Get_placement& m_get_placement;
  Visitor m_visitor;
  bool m_has_border;

private:
  Edge_data_array m_edge_data;

  boost::scoped_ptr<PQ> mPQ;

  size_type m_initial_edge_count;
  size_type m_current_edge_count;

  FT m_max_dihedral_angle_squared_cos;

  CGAL_SMS_DEBUG_CODE(unsigned m_step;)
};

template<class TM, class GT, class SP, class VIM, class VPM,class HIM, class ECM, class CF, class PF, class V>
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
EdgeCollapse(Triangle_mesh& tmesh,
             const Geom_traits& traits,
             const Should_stop& should_stop,
             const Vertex_index_map& vim,
             const Vertex_point_map& vpm,
             const Halfedge_index_map& him,
             const Edge_is_constrained_map& ecm,
             const Get_cost& get_cost,
             const Get_placement& get_placement,
             Visitor visitor)
  :
    m_tm(tmesh),
    m_traits(traits),
    m_should_stop(should_stop),
    m_vim(vim),
    m_vpm(vpm),
    m_him(him),
    m_ecm(ecm),
    m_get_cost(get_cost),
    m_get_placement(get_placement),
    m_visitor(visitor),
    m_has_border(!is_closed(tmesh))
{
  m_max_dihedral_angle_squared_cos = CGAL::square(std::cos(1.0 * CGAL_PI / 180.0));

  CGAL_SMS_TRACE(0, "EdgeCollapse of TM with " << (num_edges(tmesh)/2) << " edges");

  CGAL_SMS_DEBUG_CODE(m_step = 0;)

#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
  for(vertex_descriptor vd : vertices(m_tm))
    CGAL_SMS_TRACE(1, vertex_to_string(vd));

  for(halfedge_descriptor ed : halfedges(m_tm))
    CGAL_SMS_TRACE(1, edge_to_string(ed));
#endif
}

template<class TM, class GT, class SP, class VIM, class VPM,class HIM, class ECM, class CF, class PF, class V>
int
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
run()
{
  CGAL_precondition(is_valid_polygon_mesh(m_tm) && CGAL::is_triangle_mesh(m_tm));

  m_visitor.OnStarted(m_tm);

  // this is similar to the visitor, but for the cost/stop/placement oracles
  internal::Oracles_initializer<Self>(*this)();

  // First collect all candidate edges in a PQ
  collect();

  // Then proceed to collapse each edge in turn
  loop();

  CGAL_SMS_TRACE(0, "Finished: " << (m_initial_edge_count - m_current_edge_count) << " edges removed.");

  int r = int(m_initial_edge_count - m_current_edge_count);

  m_visitor.OnFinished(m_tm);

  return r;
}

template<class TM, class GT, class SP, class VIM, class VPM,class HIM, class ECM, class CF, class PF, class V>
void
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
collect()
{
  CGAL_SMS_TRACE(0, "collecting edges...");

  // loop over all the _undirected_ edges in the surface putting them in the PQ

  const size_type ne = num_edges(m_tm); // if the mesh has garbage, you might have "ne > edges(tm).size()"
  m_initial_edge_count = m_current_edge_count = size_type(edges(m_tm).size());

  m_edge_data.reset(new Edge_data[ne]);
  mPQ.reset(new PQ(ne, Compare_cost(this), edge_id(this)));

  CGAL_assertion_code(size_type num_inserted = 0);
  CGAL_assertion_code(size_type num_not_inserted = 0);

  std::set<halfedge_descriptor> zero_length_edges;

  for(edge_descriptor e : edges(m_tm))
  {
    const halfedge_descriptor h = halfedge(e, m_tm);

    if(is_constrained(h))
    {
      CGAL_assertion_code(++num_not_inserted);
      continue; // no not insert constrainted edges
    }

    const Profile profile = create_profile(h);
    if(!m_traits.equal_3_object()(profile.p0(), profile.p1()))
    {
      Edge_data& data = get_data(h);

      data.cost() = cost(profile);
      insert_in_PQ(h, data);

      m_visitor.OnCollected(profile, data.cost());

      CGAL_assertion_code(++num_inserted);
    }
    else
    {
      zero_length_edges.insert(primary_edge(h));
      CGAL_assertion_code(++num_not_inserted);
    }

    CGAL_SMS_TRACE(2, edge_to_string(h));
  }

  CGAL_assertion(num_inserted + num_not_inserted == m_initial_edge_count);

  for(halfedge_descriptor hd : zero_length_edges)
  {
    const Profile profile = create_profile(hd);

    if(!is_collapse_topologically_valid(profile))
      continue;

    // edges of length 0 removed no longer need to be treated
    if(profile.left_face_exists())
    {
      halfedge_descriptor h_to_remove = is_constrained(profile.vL_v0()) ?
                                          primary_edge(profile.v1_vL()) :
                                          primary_edge(profile.vL_v0());
      zero_length_edges.erase(h_to_remove);
      Edge_data& data = get_data(h_to_remove);

      if(data.is_in_PQ())
      {
        CGAL_SMS_TRACE(2, "Removing E" << get_edge_id(h_to_remove) << " from PQ");
        remove_from_PQ(h_to_remove, data);
      }

      --m_current_edge_count;
    }

    if(profile.right_face_exists())
    {
      halfedge_descriptor h_to_remove = is_constrained(profile.vR_v1()) ?
                                          primary_edge(profile.v0_vR()) :
                                          primary_edge(profile.vR_v1());
      zero_length_edges.erase(h_to_remove);
      Edge_data& data = get_data(h_to_remove);

      if(data.is_in_PQ())
      {
        CGAL_SMS_TRACE(2, "Removing E" << get_edge_id(h_to_remove) << " from PQ");
        remove_from_PQ(h_to_remove, data);
      }

      --m_current_edge_count;
    }

    --m_current_edge_count;

    //the placement is trivial, it's always the point itself
    Placement_type placement = profile.p0();
    vertex_descriptor v = halfedge_collapse_bk_compatibility(profile.v0_v1(), m_ecm);
    put(m_vpm, v, *placement);

    m_visitor.OnCollapsed(profile, v);
    internal::After_collapse_oracles_updater<Self>(*this)(profile, v);
  }

  CGAL_SMS_TRACE(0, "Initial edge count: " << m_initial_edge_count);
}

template<class TM, class GT, class SP, class VIM, class VPM,class HIM, class ECM, class CF, class PF, class V>
void
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
loop()
{
  CGAL_SMS_TRACE(0, "Collapsing edges...");

  CGAL_assertion_code(size_type non_collapsable_count = 0);

  // Pops and processes each edge from the PQ

  boost::optional<halfedge_descriptor> opt_h;

#ifdef CGAL_SURF_SIMPL_INTERMEDIATE_STEPS_PRINTING
  int i_rm = 0;
#endif

  while((opt_h = pop_from_PQ()))
  {
    CGAL_SMS_TRACE(1, "Popped " << edge_to_string(*opt_h));
    CGAL_assertion(!is_constrained(*opt_h));

    const Profile profile = create_profile(*opt_h);
    Cost_type cost = get_data(*opt_h).cost();

    m_visitor.OnSelected(profile, cost, m_initial_edge_count, m_current_edge_count);

    if(cost)
    {
      if(m_should_stop(*cost, profile, m_initial_edge_count, m_current_edge_count))
      {
        m_visitor.OnStopConditionReached(profile);

        CGAL_SMS_TRACE(0, "Stop condition reached with initial edge count=" << m_initial_edge_count
                            << " current edge count=" << m_current_edge_count
                            << " current edge: " << edge_to_string(*opt_h));
        break;
      }

      if(is_collapse_topologically_valid(profile))
      {
        Placement_type placement = get_placement(profile);

        if(is_collapse_geometrically_valid(profile, placement))
        {
#ifdef CGAL_SURF_SIMPL_INTERMEDIATE_STEPS_PRINTING
          std::cout << "step " << i_rm << " " << get(m_vpm, source(*h, m_tm))
                                       << " " << get(m_vpm, target(*h, m_tm)) << "\n";
#endif
          collapse(profile, placement);

#ifdef CGAL_SURF_SIMPL_INTERMEDIATE_STEPS_PRINTING
          std::stringstream sstr;
          sstr << "debug/P-";
          if(i_rm<10) sstr << "0";
          if(i_rm<100) sstr << "0";
          sstr << i_rm << ".off";
          std::ofstream out(sstr.str().c_str());
          out << m_tm;
          ++i_rm;
#endif
        }
      }
      else
      {
        CGAL_assertion_code(++non_collapsable_count);

        m_visitor.OnNonCollapsable(profile);

        CGAL_SMS_TRACE(1, edge_to_string(*opt_h) << " NOT Collapsable" );
      }
    }
    else
    {
      CGAL_SMS_TRACE(1, edge_to_string(*opt_h) << " uncomputable cost." );
    }
  }
}

template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
bool
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
is_border_or_constrained(const vertex_descriptor v) const
{
  for(halfedge_descriptor h : halfedges_around_target(v, m_tm))
  {
    if(is_border_edge(h) || is_constrained(h))
      return true;
  }

  return false;
}

template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
bool
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
is_constrained(const vertex_descriptor v) const
{
  for(halfedge_descriptor h : halfedges_around_target(v, m_tm))
    if(is_constrained(h))
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
template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
bool
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
is_collapse_topologically_valid(const Profile& profile)
{
  bool res = true;

  CGAL_SMS_TRACE(3,"Testing topological collapsabilty of p_q=V" << get(m_vim,profile.v0()) << "(%" << degree(profile.v0(), m_tm) << ")"
                 << "->V" << get(m_vim, profile.v1()) << "(%" << degree(profile.v1(), m_tm) << ")");

  CGAL_SMS_TRACE(4, "is p_q border:" << profile.is_v0_v1_a_border());
  CGAL_SMS_TRACE(4, "is q_q border:" << profile.is_v1_v0_a_border());

  out_edge_iterator eb1, ee1;
  out_edge_iterator eb2, ee2;

  CGAL_SMS_TRACE(4, "  t=V" << (profile.left_face_exists() ? get(m_vim,profile.vL()) : -1)
                            << "(%" << (profile.left_face_exists() ? degree(profile.vL(), m_tm) : 0) << ")");
  CGAL_SMS_TRACE(4, "  b=V" << (profile.right_face_exists() ? get(m_vim,profile.vR()) : -1)
                            << "(%" << (profile.right_face_exists() ? degree(profile.vR(), m_tm) :0) << ")");

  // Simple tests handling the case of non-manifold situations at a vertex or edge (pinching)
  // (even if we advertise one should not use a surface mesh with such features)
  if(profile.left_face_exists())
  {
    if(CGAL::is_border(opposite(profile.v1_vL(), m_tm), m_tm) &&
       CGAL::is_border(opposite(profile.vL_v0(), m_tm), m_tm))
      return false;

    if(profile.right_face_exists() &&
       CGAL::is_border(opposite(profile.vR_v1(), m_tm), m_tm) &&
       CGAL::is_border(opposite(profile.v0_vR(), m_tm), m_tm))
      return false;
  }
  else
  {
    if(profile.right_face_exists())
    {
      if(CGAL::is_border(opposite(profile.vR_v1(), m_tm), m_tm) &&
         CGAL::is_border(opposite(profile.v0_vR(), m_tm), m_tm))
        return false;
    }
    else
      return false;
  }

  // The following loop checks the link condition for v0_v1.
  // Specifically, that for every vertex 'k' adjacent to both 'p and 'q', 'pkq' is a face of the mesh.
  //
  for(boost::tie(eb1,ee1) = halfedges_around_source(profile.v0(), m_tm); res && eb1 != ee1; ++eb1)
  {
    halfedge_descriptor v0_k = *eb1;

    if(v0_k != profile.v0_v1())
    {
      vertex_descriptor k = target(v0_k, m_tm);

      for(boost::tie(eb2,ee2) = halfedges_around_source(k, m_tm); res && eb2 != ee2; ++eb2)
      {
        halfedge_descriptor k_v1 = *eb2;

        if(target(k_v1, m_tm) == profile.v1())
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

          bool is_face = (profile.vL() == k && profile.left_face_exists()) ||
                         (profile.vR() == k && profile.right_face_exists());

          CGAL_assertion_code(
            if(is_face)
            {
              // Is k_v1 the halfedge bounding the face 'k-v1-v0'?
              if(!is_border(k_v1, m_tm) && target(next(k_v1, m_tm), m_tm) == profile.v0())
              {
                CGAL_assertion(target(k_v1, m_tm) == profile.v1());
                CGAL_assertion(target(next(k_v1, m_tm), m_tm) == profile.v0());
                CGAL_assertion(target(next(next(k_v1, m_tm), m_tm), m_tm) == k);
              }
              else // or is it the opposite?
              {
                halfedge_descriptor v1_k = opposite(k_v1, m_tm);
                CGAL_assertion(!is_border(v1_k, m_tm));
                CGAL_assertion(target(v1_k, m_tm) == k);
                CGAL_assertion(target(next(v1_k, m_tm), m_tm) == profile.v0());
                CGAL_assertion(target(next(next(v1_k, m_tm), m_tm), m_tm) == profile.v1());
              }
            }
          );

          if(!is_face)
          {
            CGAL_SMS_TRACE(3,"  k=V" << get(m_vim,k) << " IS NOT in a face with p-q. NON-COLLAPSABLE edge.");
            res = false;
            break;
          }
          else
          {
            CGAL_SMS_TRACE(4,"  k=V" << get(m_vim,k) << " is in a face with p-q");
          }
        }
      }
    }
  }

  if(res)
  {
    /// ensure two constrained edges cannot get merged
    if(is_edge_adjacent_to_a_constrained_edge(profile, m_ecm))
      return false;

    if(profile.is_v0_v1_a_border())
    {
      if(is_open_triangle(profile.v0_v1()))
      {
        res = false;
        CGAL_SMS_TRACE(3,"  p-q belongs to an open triangle. NON-COLLAPSABLE edge.");
      }
    }
    else if(profile.is_v1_v0_a_border())
    {
      if(is_open_triangle(profile.v1_v0()))
      {
        res = false;
        CGAL_SMS_TRACE(3,"  p-q belongs to an open triangle. NON-COLLAPSABLE edge.");
      }
    }
    else
    {
      if(is_border(profile.v0(), m_tm) && is_border(profile.v1(), m_tm))
      {
        res = false;
        CGAL_SMS_TRACE(3,"  both p and q are boundary vertices but p-q is not. NON-COLLAPSABLE edge.");
      }
      else
      {
        bool tetra = is_tetrahedron(profile.v0_v1());

        //CGAL_assertion(tetra == m_tm.is_tetrahedron(profile.v0_v1()));

        if(tetra)
        {
          res = false;
          CGAL_SMS_TRACE(3,"  p-q belongs to a tetrahedron. NON-COLLAPSABLE edge.");
        }

        if(next(profile.v0_v1(), m_tm) == opposite(prev(profile.v1_v0(), m_tm), m_tm) &&
           prev(profile.v0_v1(), m_tm) == opposite(next(profile.v1_v0(), m_tm), m_tm))
        {
          CGAL_SMS_TRACE(3,"  degenerate volume.");
          return false;
        }
      }
    }
  }

  return res;
}

template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
bool
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
is_tetrahedron(const halfedge_descriptor h)
{
  return CGAL::is_tetrahedron(h, m_tm);
}

template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
bool
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
is_open_triangle(const halfedge_descriptor h1)
{
  bool res = false;

  halfedge_descriptor h2 = next(h1, m_tm);
  halfedge_descriptor h3 = next(h2, m_tm);

  // First check if it is a triangle
  if(next(h3, m_tm) == h1)
  {
    // Now check if it is open
    CGAL_SMS_TRACE(4,"  p-q is a border edge... checking E" << get_edge_id(h2) << " and E" << get_edge_id(h3));

    res = is_border(h2, m_tm) && is_border(h3, m_tm);

    CGAL_assertion(res == (is_border(h1, m_tm) &&
                           is_border(next(h1, m_tm), m_tm) &&
                           is_border(next(next(h1, m_tm), m_tm), m_tm)));
  }

  return res;
}

// Given triangles 'p0,p1,p2' and 'p0,p2,p3', both shared along edge 'v0-v2',
// determine if they are geometrically valid: that is, the ratio of their
// respective areas is no greater than a max value and the internal
// dihedral angle formed by their supporting planes is no greater than
// a given threshold
template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
bool
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
are_shared_triangles_valid(const Point& p0, const Point& p1, const Point& p2, const Point& p3) const
{
  bool res = false;

  Vector e01 = m_traits.construct_vector_3_object()(p0, p1);
  Vector e02 = m_traits.construct_vector_3_object()(p0, p2);
  Vector e03 = m_traits.construct_vector_3_object()(p0, p3);

  Vector n012 = m_traits.construct_cross_product_vector_3_object()(e01, e02);
  Vector n023 = m_traits.construct_cross_product_vector_3_object()(e02, e03);

  FT l012 = m_traits.compute_scalar_product_3_object()(n012, n012);
  FT l023 = m_traits.compute_scalar_product_3_object()(n023, n023);

  FT larger = (std::max)(l012, l023);
  FT smaller = (std::min)(l012, l023);

  const FT max_area_ratio = 1e8;

  CGAL_SMS_TRACE(4,"    Testing validity of shared triangles:"
                 << "\n      p0=" << xyz_to_string(p0) << "\n      p1=" << xyz_to_string(p1) << "\n      p2=" << xyz_to_string(p2) << "\n      p3=" << xyz_to_string(p3)
                 << "\n      e01=" << xyz_to_string(e01) << "\n      e02=" << xyz_to_string(e02) << "\n      e03=" << xyz_to_string(e03)
                 << "\n      n012=" << xyz_to_string(n012) << "\n      n023=" << xyz_to_string(n023)
                 << "\n      l012=" << n_to_string(l012) << "\n      l023=" << n_to_string(l023));

  if(larger < max_area_ratio * smaller)
  {
    FT l0123 = m_traits.compute_scalar_product_3_object()(n012, n023);
    CGAL_SMS_TRACE(4,"\n      l0123=" << n_to_string(l0123));

    if(CGAL_NTS is_positive(l0123))
    {
      res = true;
    }
    else
    {
      CGAL_SMS_TRACE(4,"\n      lhs: " << n_to_string((l0123 * l0123) / (l012 * l023)) << " <= rhs: " << m_max_dihedral_angle_squared_cos);

      if((l0123 * l0123) <= m_max_dihedral_angle_squared_cos * (l012 * l023))
      {
        res = true;
      }
    }
  }

  return res;
}

// Returns the directed halfedge connecting v0 to v1, if exists.
template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
typename EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::halfedge_descriptor
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
find_connection(const vertex_descriptor v0,
                const vertex_descriptor v1) const
{
  for(halfedge_descriptor out : halfedges_around_source(v0, m_tm))
  {
    if(target(out, m_tm) == v1)
      return out;
  }

  return Graph_traits::null_halfedge();
}

// Given the edge 'e' around the link for the collapsinge edge "v0-v1", finds the vertex that makes a triangle adjacent to 'e' but exterior to the link (i.e not containing v0 nor v1)
// If 'e' is a null handle OR 'e' is a border edge, there is no such triangle and a null handle is returned.
template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
typename EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::vertex_descriptor
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
find_exterior_link_triangle_3rd_vertex(const halfedge_descriptor e,
                                       const vertex_descriptor v0,
                                       const vertex_descriptor v1) const
{
  vertex_descriptor r;

  if(handle_assigned(e))
  {
    vertex_descriptor ra = target(next(e, m_tm), m_tm);
    vertex_descriptor rb = source(prev(e, m_tm), m_tm);

    if(ra == rb && ra != v0 && ra != v1)
    {
      r = ra;
    }
    else
    {
      ra = target(next(opposite(e, m_tm), m_tm), m_tm);
      rb = source(prev(opposite(e, m_tm), m_tm), m_tm);

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
template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
bool
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
is_collapse_geometrically_valid(const Profile& profile, Placement_type k0)
{
  bool res = false;

  CGAL_SMS_TRACE(3,"Testing geometrical collapsabilty of v0-v1=E" << get_edge_id(profile.v0_v1()));
  if(k0)
  {
    res = true;

    // Use the current link to extract all local triangles incident to 'vx' in the collapsed mesh
    // (which at this point doesn't exist yet)
    typedef typename std::vector<vertex_descriptor>::const_iterator link_iterator;

    link_iterator linkb = profile.link().begin();
    link_iterator linke = profile.link().end();
    link_iterator linkl = std::prev(linke);

    for(link_iterator l=linkb; l!=linke && res; ++l)
    {
      link_iterator pv = (l == linkb ? linkl : std::prev(l));
      link_iterator nx = (l == linkl ? linkb : std::next(l));

      // k0,k1 and k3 are three consecutive vertices along the link.
      vertex_descriptor k1 = *pv;
      vertex_descriptor k2 = *l;
      vertex_descriptor k3 = *nx;

      CGAL_SMS_TRACE(4, "  Screening link vertices k1=V" << get(m_vim, k1)
                          << " k2=V" << get(m_vim, k2)
                          << " k3=V" << get(m_vim, k3));

      halfedge_descriptor e12 = find_connection(k1,k2);
      halfedge_descriptor e23 = (k3 != k1) ? find_connection(k2,k3) : Graph_traits::null_halfedge();

      // If 'k1-k2-k3' are connected there will be two adjacent triangles 'k0,k1,k2' and 'k0,k2,k3' after the collapse.
      if(handle_assigned(e12) && handle_assigned(e23))
      {
        CGAL_SMS_TRACE(4,"    Link triangles shared");

        if(!are_shared_triangles_valid(*k0, get_point(k1), get_point(k2), get_point(k3)))
        {
          CGAL_SMS_TRACE(3, "    Triangles VX-V" << get(m_vim, k1)
                              << "-V" << get(m_vim, k2)
                              << " and VX-V" << get(m_vim, k3)
                              << " are not geometrically valid. Collapse rejected");
          res = false;
        }
      }

      if(res)
      {
        // Also check the triangles 'k0,k1,k2' and it's adjacent along e12: 'k4,k2,k1', if exist
        vertex_descriptor k4 = find_exterior_link_triangle_3rd_vertex(e12, profile.v0(), profile.v1());

        // There is indeed a triangle shared along e12
        if(handle_assigned(k4))
        {
          CGAL_SMS_TRACE(4, "    Found exterior link triangle shared along E" << get_edge_id(e12)
                              << " with third vertex: V" << get(m_vim, k4));

          if(!are_shared_triangles_valid(get_point(k1), get_point(k4), get_point(k2), *k0))
          {
            CGAL_SMS_TRACE(3, "    Triangles V" << get(m_vim, k1)
                                << "-V" << get(m_vim, k4) << " and V"
                                << get(m_vim, k2) << "-VX are not geometrically valid. Collapse rejected");
            res = false;
          }
        }
      }

      if(res)
      {
        // And finally, check the triangles 'k0,k2,k3' and it's adjacent e23: 'k5,k3,k2' if exist
        vertex_descriptor k5 = find_exterior_link_triangle_3rd_vertex(e23, profile.v0(), profile.v1());

        // There is indeed a triangle shared along e12
        if(handle_assigned(k5))
        {
          CGAL_SMS_TRACE(4, "    Found exterior link triangle shared along E" << get_edge_id(e23)
                             << " with third vertex: V" << get(m_vim, k5));

          if(!are_shared_triangles_valid(get_point(k2), get_point(k5), get_point(k3), *k0))
          {
            CGAL_SMS_TRACE(3, "    Triangles V" << get(m_vim, k2)
                                << "-V" << get(m_vim, k5)
                                << " and V" << get(m_vim, k3) << "-VX are not geometrically valid. Collapse rejected");
            res = false;
          }
        }
      }
    }
  }

  return res;
}

template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
void
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
collapse(const Profile& profile,
         Placement_type placement)
{
  CGAL_SMS_TRACE(1, "S" << m_step << ". Collapsing " << edge_to_string(profile.v0_v1()));

  vertex_descriptor v_res;

  m_visitor.OnCollapsing(profile, placement);

  --m_current_edge_count;

  CGAL_assertion_code(
    size_type resulting_vertex_count = size_type(vertices(m_tm).size());
    size_type result_edge_count = size_type(edges(m_tm).size());
  );

  // If the top/bottom facets exists, they are removed and the edges v0vt and Q-B along with them.
  // In that case their corresponding pairs must be pop off the queue
  if(profile.left_face_exists())
  {
    halfedge_descriptor h_v0vL = primary_edge(profile.vL_v0());
    if(is_constrained(h_v0vL)) // make sure a constrained edge will not disappear
      h_v0vL = primary_edge(profile.v1_vL());

    CGAL_SMS_TRACE(3, "V0VL E" << get_edge_id(h_v0vL)
                       << "(V" << get(m_vim, source(h_v0vL, m_tm))
                       << "->V" << get(m_vim,target(h_v0vL, m_tm)) << ")");

    Edge_data& data = get_data(h_v0vL);
    if(data.is_in_PQ())
    {
      CGAL_SMS_TRACE(2, "Removing E" << get_edge_id(h_v0vL) << " from PQ");
      remove_from_PQ(h_v0vL, data);
    }

    --m_current_edge_count;
    CGAL_assertion_code(--result_edge_count);
  }

  if(profile.right_face_exists())
  {
    halfedge_descriptor lVRV1 = primary_edge(profile.vR_v1());
    if(is_constrained(lVRV1)) // make sure a constrained edge will not disappear
      lVRV1 = primary_edge(profile.v0_vR());

    CGAL_SMS_TRACE(3, "V1VRE" << get_edge_id(lVRV1)
                        << "(V" << get(m_vim, source(lVRV1, m_tm))
                        << "->V" << get(m_vim, target(lVRV1, m_tm)) << ")");

    Edge_data& data = get_data(lVRV1);
    if(data.is_in_PQ())
    {
      CGAL_SMS_TRACE(2, "Removing E" << get_edge_id(lVRV1) << " from PQ");
      remove_from_PQ(lVRV1, data);
    }

    --m_current_edge_count;
    CGAL_assertion_code(--result_edge_count);
  }

  CGAL_SMS_TRACE(1, "Removing:\n  v0v1: E" << get_edge_id(profile.v0_v1())
                      << "(V" << get(m_vim, profile.v0())
                      << "->V" << get(m_vim, profile.v1()) << ")");

  // Perform the actuall collapse.
  // This is an external function.
  // It's REQUIRED to remove ONLY 1 vertex (P or Q) and edges PQ, PT and QB
  // (PT and QB are removed if they are not null).
  // All other edges must be kept.
  // All directed edges incident to vertex removed are relink to the vertex kept.
  v_res = halfedge_collapse_bk_compatibility(profile.v0_v1(), m_ecm);

  CGAL_assertion_code(--result_edge_count);
  CGAL_assertion_code(--resulting_vertex_count);

  CGAL_assertion(result_edge_count == edges(m_tm).size());
  CGAL_assertion(resulting_vertex_count == vertices(m_tm).size());
  CGAL_expensive_assertion(is_valid_polygon_mesh(m_tm) && CGAL::is_triangle_mesh(m_tm));

  CGAL_SMS_TRACE(1, "V" << get(m_vim, v_res) << " kept.");

#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
  for(halfedge_descriptor hd : halfedges_around_source(v_res, m_tm))
    CGAL_SMS_TRACE(2, edge_to_string(hd));
#endif

  if(placement)
  {
    CGAL_SMS_TRACE(1, "New vertex point: " << xyz_to_string(*placement));
    put(m_vpm, v_res, *placement);
  }

  m_visitor.OnCollapsed(profile, v_res);
  internal::After_collapse_oracles_updater<Self>(*this)(profile, v_res);

  update_neighbors(v_res);

  CGAL_SMS_DEBUG_CODE(++m_step;)
}

template<class TM, class GT, class SP, class VIM, class VPM, class HIM, class ECM, class CF, class PF, class V>
void
EdgeCollapse<TM,GT,SP,VIM,VPM,HIM,ECM,CF,PF,V>::
update_neighbors(const vertex_descriptor v_kept)
{
  CGAL_SMS_TRACE(3,"Updating cost of neighboring edges...");

  // (A) collect all edges to update their cost: all those around each vertex adjacent to the vertex kept
  typedef std::set<halfedge_descriptor, Compare_id>                       Edge_set;

  Edge_set edges_to_update(Compare_id(this));
  Edge_set edges_to_insert(Compare_id(this));

  // (A.1) loop around all vertices adjacent to the vertex kept
  for(halfedge_descriptor h : halfedges_around_target(v_kept, m_tm))
  {
    vertex_descriptor v_adj = source(h, m_tm);

    // (A.2) loop around all edges incident on each adjacent vertex
    for(halfedge_descriptor h2 : halfedges_around_target(v_adj, m_tm))
    {
      h2 = primary_edge(h2);

      Edge_data& data2 = get_data(h2);
      CGAL_SMS_TRACE(4,"Inedge around V" << get(m_vim, v_adj) << edge_to_string(h2));

      // Only edges still in the PQ needs to be updated, the other needs to be re-inserted
      if(data2.is_in_PQ())
        edges_to_update.insert(h2);
      else
        edges_to_insert.insert(h2);
    }
  }

  // (B) Proceed to update the costs.
  for(halfedge_descriptor h : edges_to_update)
  {
    Edge_data& data = get_data(h);
    const Profile& profile = create_profile(h);
    data.cost() = cost(profile);

    CGAL_SMS_TRACE(3, edge_to_string(h) << " updated in the PQ");

    update_in_PQ(h, data);
  }

  // (C) Insert ignored edges
  //
  // I think that this should be done for edges eliminated because of the geometric criteria
  // and not the topological one.However maintaining such a set might be more expensive
  // and hard to be safe ...
  for(halfedge_descriptor h : edges_to_insert)
  {
    if(is_constrained(h))
      continue; //do not insert constrained edges

    Edge_data& data = get_data(h);
    const Profile& profile = create_profile(h);
    data.cost() = cost(profile);

    CGAL_SMS_TRACE(3, edge_to_string(h) << " re-inserted in the PQ");
    insert_in_PQ(h, data);
  }
}

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_H //
