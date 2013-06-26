#ifndef MEAN_CURVATURE_SKELETON_H
#define MEAN_CURVATURE_SKELETON_H

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

// Compute cotangent Laplacian
#include <CGAL/internal/Mean_curvature_skeleton/Weights.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/internal/Mean_curvature_skeleton/Edge_minimum_length_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>

// Map used to mark edges as fixed
#include <CGAL/Unique_hash_map.h>

namespace SMS = CGAL::Surface_mesh_simplification;

namespace CGAL {
namespace internal {

template<class Polyhedron, class edge_descriptor, class Point>
edge_descriptor mesh_split(Polyhedron *polyhedron, edge_descriptor ei, Point pn)
{
  edge_descriptor en = polyhedron->split_edge(ei);
  en->vertex()->point() = pn;
  polyhedron->split_facet(en, ei->next());

  edge_descriptor ej = en->opposite();
  if (!(ej->is_border()))
  {
    polyhedron->split_facet(en, ej->next());
  }

  return en;
}

}
}

namespace CGAL {

template <class Polyhedron, class SparseLinearAlgebraTraits_d,
          class PolyhedronVertexIndexMap, class PolyhedronEdgeIndexMap>
class Mean_curvature_skeleton
{
// Public types
public:

  // Geometric types
  typedef typename Polyhedron::Traits         Kernel;
  typedef typename Kernel::Vector_3           Vector;
  typedef typename Kernel::Point_3            Point;

  // Repeat Polyhedron types
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor	         vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator            vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator              edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator           in_edge_iterator;
  typedef typename internal::Cotangent_weight<Polyhedron>                      Weight_calculator;

  // Data members.
private:

  Polyhedron* polyhedron;
  PolyhedronVertexIndexMap vertex_id_pmap;
  PolyhedronEdgeIndexMap edge_id_pmap;

  Weight_calculator weight_calculator;
  std::vector<double> edge_weight;
  SparseLinearAlgebraTraits_d m_solver;

  double omega_L;
  double omega_H;
  double edgelength_TH;
  double TH_ALPHA;

  int vertex_id_count;

  std::map<size_t, bool> is_vertex_fixed_map;
  std::vector<double> halfedge_angle;

  //
  // BGL property map which indicates whether an edge is border OR is marked as non-removable
  //
  class Constrains_map : public boost::put_get_helper<bool, Constrains_map>
  {
  public:

    typedef boost::readable_property_map_tag                                category;
    typedef bool                                                            value_type;
    typedef bool                                                            reference;
    typedef typename boost::graph_traits<Polyhedron const>::edge_descriptor key_type;

    Constrains_map() : mConstrains(false) {}

    reference operator[](key_type const& e) const
    {
      return e->is_border() || is_constrained(e);
    }

    void set_is_constrained (key_type const& e, bool is)
    {
      mConstrains[e] = is;
    }

    bool is_constrained(key_type const& e) const
    {
      return mConstrains.is_defined(e) ? mConstrains[e] : false;
    }

  private:

    CGAL::Unique_hash_map<key_type, bool> mConstrains ;

  };

  // Public methods
public:

  // The constructor gets the polyhedron that we will model
  Mean_curvature_skeleton(Polyhedron* P,
                          PolyhedronVertexIndexMap Vertex_index_map,
                          PolyhedronEdgeIndexMap Edge_index_map,
                          double omega_L, double omega_H, double edgelength_TH,
                          Weight_calculator weight_calculator = Weight_calculator()
                          )
    :polyhedron(P), vertex_id_pmap(Vertex_index_map), edge_id_pmap(Edge_index_map),
      omega_L(omega_L), omega_H(omega_H), edgelength_TH(edgelength_TH), TH_ALPHA(110),
      weight_calculator(weight_calculator)
  {
    TH_ALPHA *= (M_PI / 180.0);

    // initialize index maps
    vertex_iterator vb, ve;
    vertex_id_count = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, vertex_id_count++);
    }

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    is_vertex_fixed_map.clear();
  }

  // Release resources
  ~Mean_curvature_skeleton(void)
  {
  }

  // compute cotangent weights of all edges
  void compute_edge_weight()
  {
    edge_weight.reserve(boost::num_edges(*polyhedron));
    edge_iterator eb, ee;
    for(boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_weight.push_back(this->weight_calculator(*eb, *polyhedron));
    }
  }

  void assemble_LHS(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    int nver = boost::num_vertices(*polyhedron);

    // initialize the Laplacian matrix
    for (int i = 0; i < nver; i++)
    {
      A.set_coef(i, i, 0.0, true);
      A.set_coef(i + nver, i, omega_H, true);
    }

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      int i = boost::get(vertex_id_pmap, *vb);
      double diagonal = 0;
      in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = boost::in_edges(*vb, *polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, *polyhedron);
        double wij = edge_weight[boost::get(edge_id_pmap, *e)] * 2.0;
        int j = boost::get(vertex_id_pmap, vj);
        A.set_coef(i, j, wij * omega_L, true);
        diagonal += -wij;
      }
      A.set_coef(i, i, diagonal);
    }
  }

  void assemble_RHS(typename SparseLinearAlgebraTraits_d::Vector& Bx,
                    typename SparseLinearAlgebraTraits_d::Vector& By,
                    typename SparseLinearAlgebraTraits_d::Vector& Bz)
  {
    // assemble right columns of linear system
    int nver = boost::num_vertices(*polyhedron);
    vertex_iterator vb, ve;
    for (int i = 0; i < nver; i++)
    {
      Bx[i] = 0;
      By[i] = 0;
      Bz[i] = 0;
    }
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      vertex_descriptor vi = *vb;
      int i = boost::get(vertex_id_pmap, vi);
      Bx[i + nver] = vi->point().x() * omega_H;
      By[i + nver] = vi->point().y() * omega_H;
      Bz[i + nver] = vi->point().z() * omega_H;
    }
  }

  void contract_geometry()
  {
    compute_edge_weight();

    // Assemble linear system At * A * X = At * B
    int nver = boost::num_vertices(*polyhedron);
    typename SparseLinearAlgebraTraits_d::Matrix A(nver * 2, nver);
    assemble_LHS(A);

    typename SparseLinearAlgebraTraits_d::Vector X(nver), Bx(nver * 2);
    typename SparseLinearAlgebraTraits_d::Vector Y(nver), By(nver * 2);
    typename SparseLinearAlgebraTraits_d::Vector Z(nver), Bz(nver * 2);
    assemble_RHS(Bx, By, Bz);

    // solve "At * A * X = At * B".
    double D;
    m_solver.pre_factor_non_symmetric(A, D);
    m_solver.linear_solver_non_symmetric(A, Bx, X);
    m_solver.linear_solver_non_symmetric(A, By, Y);
    m_solver.linear_solver_non_symmetric(A, Bz, Z);

    // copy to mesh
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      vertex_descriptor vi = *vb;
      int i = boost::get(vertex_id_pmap, vi);
      Point p(X[i], Y[i], Z[i]);
      vi->point() = p;
    }
  }

  int collapse_short_edges(double edgelength_TH)
  {
    Constrains_map constrains_map;

    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, *polyhedron);
      vertex_descriptor vj = boost::target(*eb, *polyhedron);
      size_t vi_idx = boost::get(vertex_id_pmap, vi);
      size_t vj_idx = boost::get(vertex_id_pmap, vj);

      if (is_vertex_fixed_map.find(vi_idx) == is_vertex_fixed_map.end()
       || is_vertex_fixed_map.find(vj_idx) == is_vertex_fixed_map.end())
      {
        continue;
      }
      // if both vertices are fixed, the edge is fixed
      if (is_vertex_fixed_map[vi_idx] && is_vertex_fixed_map[vj_idx])
      {
        constrains_map.set_is_constrained(*eb, true);
      }
    }

    // This is a stop predicate (defines when the algorithm terminates).
    // The simplification stops when the length of all edges is greater than the minimum threshold.
    CGAL::internal::Minimum_length_predicate<Polyhedron> stop(edgelength_TH);

    int r = SMS::edge_collapse
                (*polyhedron
                ,stop
                ,CGAL::get_cost     (SMS::Edge_length_cost  <Polyhedron>())
                      .get_placement(SMS::Midpoint_placement<Polyhedron>())
                      .edge_is_border_map(constrains_map)
                );
    return r;
  }

  void compute_incident_angle()
  {
    halfedge_angle.clear();
    int ne = boost::num_edges(*polyhedron);
    halfedge_angle.resize(ne, 0);

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      int e_id = boost::get(edge_id_pmap, *eb);
      edge_descriptor ed = *eb;

      if (ed->is_border())
      {
        halfedge_angle[e_id] = -1;
      }
      else
      {
        vertex_descriptor vi = boost::source(ed, *polyhedron);
        vertex_descriptor vj = boost::target(ed, *polyhedron);
        edge_descriptor ed_next = ed->next();
        vertex_descriptor vk = boost::target(ed_next, *polyhedron);
        Point pi = vi->point();
        Point pj = vj->point();
        Point pk = vk->point();

        double dis2_ij = squared_distance(pi, pj);
        double dis2_ik = squared_distance(pi, pk);
        double dis2_jk = squared_distance(pj, pk);
        double dis_ij = sqrt(dis2_ij);
        double dis_ik = sqrt(dis2_ik);
        double dis_jk = sqrt(dis2_jk);

        /// A degenerate triangle will never undergo a split (but rather a collapse...)
        if (dis_ij < edgelength_TH || dis_ik < edgelength_TH || dis_jk < edgelength_TH)
        {
          halfedge_angle[e_id] = -1;
        }
        else
        {
          halfedge_angle[e_id] = acos((dis2_ik + dis2_jk - dis2_ij) / (2.0 * dis_ik * dis_jk));
        }
      }
    }
  }

  Point project_vertex(const Point& ps, const Point& pt, const Point& pk)
  {
    CGAL::internal::Vector vec_st = CGAL::internal::Vector(ps, pt);
    CGAL::internal::Vector vec_sk = CGAL::internal::Vector(ps, pk);

    vec_st.normalize();
    double len = vec_st.length();
    double t = vec_st.dot(vec_sk);
    Point st = Point(vec_st[0] * t, vec_st[1] * t, vec_st[2] * t);
    Point pn = Point(ps[0] + st[0], ps[1] + st[1], ps[2] + st[2]);
    return pn;
  }

  bool split_flat_triangle()
  {
    compute_incident_angle();

    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ei = *eb;
      edge_descriptor ej = ei->opposite();
      int ei_id = boost::get(edge_id_pmap, ei);
      int ej_id = boost::get(edge_id_pmap, ej);

      vertex_descriptor vs = boost::source(ei, *polyhedron);
      vertex_descriptor vt = boost::target(ei, *polyhedron);
      size_t vs_id = boost::get(vertex_id_pmap, vs);
      size_t vt_id = boost::get(vertex_id_pmap, vt);
      Point ps = vs->point();
      Point pt = vt->point();

      if (is_vertex_fixed_map.find(vs_id) != is_vertex_fixed_map.end()
       && is_vertex_fixed_map.find(vt_id) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[vs_id] && is_vertex_fixed_map[vt_id])
        {
          continue;
        }
      }

      // for border edge, the angle is -1
      double angle_i = halfedge_angle[ei_id];
      double angle_j = halfedge_angle[ej_id];
      if (angle_i < TH_ALPHA || angle_j < TH_ALPHA)
      {
        continue;
      }

      edge_descriptor ek;
      if (angle_i > angle_j)
      {
        ek = ei->next();
      }
      else
      {
        ek = ej->next();
      }
      vertex_descriptor vk = boost::target(ek, *polyhedron);
      Point pk = vk->point();
      Point pn = project_vertex(ps, pt, pk);
      edge_descriptor en = CGAL::internal::mesh_split(polyhedron, ei, pn);
      // set id for new vertex
      boost::put(vertex_id_pmap, en->vertex(), vertex_id_count++);
      return true;
    }
    return false;
  }

  int iteratively_split_triangles()
  {
    int num_splits = 0;
    while(true)
    {
      if (split_flat_triangle())
      {
        num_splits++;
      }
      else
      {
        break;
      }
    }
    return num_splits;
  }

  void updateTopology()
  {
    int num_collapses = collapse_short_edges(edgelength_TH);
    std::cout << "collapse " << num_collapses << " edges.\n";

    int num_splits = iteratively_split_triangles();
    std::cout << "split " << num_splits << " edges.\n";
  }
};

} //namespace CGAL

#endif // MEAN_CURVATURE_SKELETON_H
