#ifndef MEAN_CURVATURE_SKELETON_H
#define MEAN_CURVATURE_SKELETON_H

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/copy.hpp>

// Compute cotangent Laplacian
#include <CGAL/internal/Mean_curvature_skeleton/Weights.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/internal/Mean_curvature_skeleton/Edge_minimum_length_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// Skip the geometric test
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Geometric_test_skipper.h>

// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>

// Map used to mark edges as fixed
#include <CGAL/Unique_hash_map.h>

// Curve skeleton data structure
#include <CGAL/Curve_skeleton.h>

#include <queue>

namespace SMS = CGAL::Surface_mesh_simplification;

namespace CGAL {
namespace internal {

template<class Polyhedron, class edge_descriptor, class Point>
edge_descriptor mesh_split(Polyhedron *polyhedron, edge_descriptor ei, Point pn)
{
  edge_descriptor en = polyhedron->split_edge(ei);
  en->vertex()->point() = pn;
  polyhedron->split_facet(en, ei->next());

  en->id() = -1;
  en->opposite()->id() = -1;
  ei->id() = -1;
  ei->opposite()->id() = -1;
  en->next()->id() = -1;
  en->next()->opposite()->id() = -1;
  en->next()->next()->id() = -1;
  ei->next()->id() = -1;
  edge_descriptor ej = en->opposite();
  if (!(ej->is_border()))
  {
    polyhedron->split_facet(ei->opposite(), ej->next());
    ej->next()->id() = -1;
    edge_descriptor ei_op_next = ei->opposite()->next();
    ei_op_next->id() = -1;
    ei_op_next->opposite()->id() = -1;
    ei_op_next->next()->id() = -1;
  }

  return en;
}

} //namespace internal
} //namespace CGAL

namespace CGAL {

template <class Polyhedron, class SparseLinearAlgebraTraits_d,
          class PolyhedronVertexIndexMap, class PolyhedronEdgeIndexMap,
          class Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> >
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
  typedef typename Polyhedron::Vertex_handle                                   Vertex_handle;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator              edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator           in_edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::out_edge_iterator		       out_edge_iterator;
  typedef typename Polyhedron::Facet_iterator                                  Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator                Halfedge_facet_circulator;
  typedef typename internal::Cotangent_weight<Polyhedron,
  internal::Cotangent_value_minimum_zero<Polyhedron,
  internal::Cotangent_value_Meyer_secure<Polyhedron> > >                       Weight_calculator;

  // Skeleton types
  typedef Curve_skeleton<Polyhedron, Graph,
  PolyhedronVertexIndexMap, PolyhedronEdgeIndexMap>                            Skeleton;

  // Mesh simplification types
  typedef SMS::Edge_profile<Polyhedron>                                        Profile;

// Data members
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
  double zero_TH;
  double area_TH;

  int vertex_id_count;
  int max_id;

  std::map<size_t, bool> is_vertex_fixed_map;
  std::map<int, int> new_id;
  std::vector<double> halfedge_angle;

  Graph g;
  std::vector<Point> points;
  // record the correspondence between final surface and original surface points
  std::map<int, std::vector<int> > correspondence;
  // record the correspondence between skeletal points and original surface points
  std::vector<std::vector<int> > skeleton_to_surface;

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

    CGAL::Unique_hash_map<key_type, bool> mConstrains;
  };

  struct Track_vertex_visitor : SMS::Edge_collapse_visitor_base<Polyhedron>
  {
    Track_vertex_visitor(std::map<int, std::vector<int> >* corr, int max_id) :
      corr(corr), max_id(max_id){}

    // Called AFTER each edge has been collapsed
    void OnCollapsed(Profile const& edge, Vertex_handle v)
    {
//      std::cerr << "before onCollapse\n";
      Vertex_handle v0 = edge.v0();
      Vertex_handle v1 = edge.v1();
      int id0 = v0->id();
      int id1 = v1->id();
      int vid = v->id();
      int from, to;
      if (id0 == vid)
      {
        from = id1;
        to = id0;
      }
      else if (id1 == vid)
      {
        from = id0;
        to = id1;
      }
      else
      {
        std::cerr << "very wrong\n";
      }
      if ((*corr).find(to) == (*corr).end())
      {
        (*corr)[to] = std::vector<int>();
      }
      // only track vertex in original mesh
      if (from < max_id)
      {
        (*corr)[to].push_back(from);
      }
      std::map<int, std::vector<int> >::iterator iter = (*corr).find(from);
      if (iter != (*corr).end())
      {
        for (size_t i = 0; i < (iter->second).size(); i++)
        {
          (*corr)[to].push_back((iter->second)[i]);
        }
        (iter->second).clear();
        (*corr).erase(iter);
      }
//      std::cerr << "after onCollapse\n";
    }

    std::map<int, std::vector<int> >* corr;
    int max_id;
  };

// Public methods
public:

  // The constructor gets the polyhedron that we will model
  Mean_curvature_skeleton(Polyhedron* P,
                          PolyhedronVertexIndexMap Vertex_index_map,
                          PolyhedronEdgeIndexMap Edge_index_map,
                          double omega_L, double omega_H,
                          double edgelength_TH, double zero_TH, double area_TH = 1e-5,
                          Weight_calculator weight_calculator = Weight_calculator()
                          )
    :polyhedron(P), vertex_id_pmap(Vertex_index_map), edge_id_pmap(Edge_index_map),
      omega_L(omega_L), omega_H(omega_H), edgelength_TH(edgelength_TH), TH_ALPHA(110),
      weight_calculator(weight_calculator), zero_TH(zero_TH), area_TH(area_TH)
  {
    TH_ALPHA *= (M_PI / 180.0);

    // initialize index maps
    vertex_iterator vb, ve;
    vertex_id_count = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, vertex_id_count++);
    }
    max_id = vertex_id_count;

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    is_vertex_fixed_map.clear();
    correspondence.clear();
  }

  // Release resources
  ~Mean_curvature_skeleton(void)
  {
  }

  void set_omega_L(double value)
  {
    omega_L = value;
  }

  void set_omega_H(double value)
  {
    omega_H = value;
  }

  void set_edgelength_TH(double value)
  {
    edgelength_TH = value;
  }

  void set_TH_ALPHA(double value)
  {
    TH_ALPHA = value;
  }

  void set_zero_TH(double value)
  {
    zero_TH = value;
  }

  Polyhedron* get_polyhedron()
  {
    return polyhedron;
  }

  void get_fixed_points(std::vector<Point>& fixed_points)
  {
    fixed_points.clear();
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[id])
        {
          vertex_descriptor vd = *vb;
          fixed_points.push_back(vd->point());
        }
      }
    }
  }

  void get_non_fixed_points(std::vector<Point>& non_fixed_points)
  {
    non_fixed_points.clear();
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) == is_vertex_fixed_map.end())
      {
          vertex_descriptor vd = *vb;
          non_fixed_points.push_back(vd->point());
      }
    }
  }

  double get_triangle_area(vertex_descriptor v1,
                           vertex_descriptor v2,
                           vertex_descriptor v3)
  {
    Point p1 = v1->point();
    Point p2 = v2->point();
    Point p3 = v3->point();
    Vector v12(p1, p2);
    Vector v13(p1, p3);
    return sqrt(cross_product(v12, v13).squared_length()) * 0.5;
  }

  double get_surface_area()
  {
    double total_area = 0;
    for (Facet_iterator i = polyhedron->facets_begin(); i != polyhedron->facets_end(); ++i)
    {
      Halfedge_facet_circulator j = i->facet_begin();
      vertex_descriptor v1 = j->vertex();
      ++j;
      vertex_descriptor v2 = j->vertex();
      ++j;
      vertex_descriptor v3 = j->vertex();
      total_area += get_triangle_area(v1, v2, v3);
    }
    return total_area;
  }

  // compute cotangent weights of all edges
  void compute_edge_weight()
  {
    edge_weight.clear();
    edge_weight.reserve(boost::num_edges(*polyhedron));
    edge_iterator eb, ee;
    for(boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_weight.push_back(this->weight_calculator(*eb, *polyhedron));
    }
  }

  void assemble_LHS(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
//    std::cerr << "start LHS\n";
    int nver = boost::num_vertices(*polyhedron);

    vertex_iterator vb, ve;
    // initialize the Laplacian matrix
    int cnt_fix = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (new_id.find(id) == new_id.end())
      {
        std::cerr << "id does not exist!\n";
      }

      int i = new_id[id];
      if (i >= nver)
      {
        std::cerr << "id is too large\n";
      }
      if (i < 0)
      {
        std::cerr << "id is too small\n";
      }
      // if the vertex is fixed
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end()
          && is_vertex_fixed_map[id])
      {
        cnt_fix++;
        A.set_coef(i + nver, i, 1.0 / zero_TH, true);
      }
      else
      {
        A.set_coef(i + nver, i, omega_H, true);
      }
    }

    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      int i = new_id[id];
      double L = omega_L;
      // if the vertex is fixed
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end()
          && is_vertex_fixed_map[id])
      {
        L = 0;
      }
      double diagonal = 0;
      in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = boost::in_edges(*vb, *polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, *polyhedron);
        double wij = edge_weight[boost::get(edge_id_pmap, *e)] * 2.0;
        int jd = boost::get(vertex_id_pmap, vj);
        int j = new_id[jd];
        A.set_coef(i, j, wij * L, true);
        diagonal += -wij;
      }
      A.set_coef(i, i, diagonal, true);
    }
//    std::cerr << "fix " << cnt_fix << " vertices\n";
//    std::cerr << "end LHS\n";
  }

  void assemble_RHS(typename SparseLinearAlgebraTraits_d::Vector& Bx,
                    typename SparseLinearAlgebraTraits_d::Vector& By,
                    typename SparseLinearAlgebraTraits_d::Vector& Bz)
  {
//    std::cerr << "start RHS\n";
    // assemble right columns of linear system
    int nver = boost::num_vertices(*polyhedron);
    vertex_iterator vb, ve;
    for (int i = 0; i < nver; i++)
    {
      Bx[i] = 0;
      By[i] = 0;
      Bz[i] = 0;
    }
    int cnt_fix = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      vertex_descriptor vi = *vb;
      int id = boost::get(vertex_id_pmap, vi);
      int i = new_id[id];

      double omega;
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end()
          && is_vertex_fixed_map[id])
      {
        cnt_fix++;
        omega = 1.0 / zero_TH;
      }
      else
      {
        omega = omega_H;
      }
      Bx[i + nver] = vi->point().x() * omega;
      By[i + nver] = vi->point().y() * omega;
      Bz[i + nver] = vi->point().z() * omega;
    }
//    std::cerr << "end RHS\n";
  }

  void update_vertex_id()
  {
//    std::cerr << "start update id\n";
    new_id.clear();
    int cnt = 0;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      new_id[id] = cnt++;
    }
//    std::cerr << "end update id\n";
  }

  void contract_geometry()
  {
//    std::cerr << "before contract geometry";
    update_vertex_id();

    compute_edge_weight();

    // Assemble linear system At * A * X = At * B
    int nver = boost::num_vertices(*polyhedron);
    typename SparseLinearAlgebraTraits_d::Matrix A(nver * 2, nver);
    assemble_LHS(A);

    typename SparseLinearAlgebraTraits_d::Vector X(nver), Bx(nver * 2);
    typename SparseLinearAlgebraTraits_d::Vector Y(nver), By(nver * 2);
    typename SparseLinearAlgebraTraits_d::Vector Z(nver), Bz(nver * 2);
    assemble_RHS(Bx, By, Bz);

//    std::cerr << "before solve\n";
    // solve "At * A * X = At * B".
    double D;
    m_solver.pre_factor_non_symmetric(A, D);
    m_solver.linear_solver_non_symmetric(A, Bx, X);
    m_solver.linear_solver_non_symmetric(A, By, Y);
    m_solver.linear_solver_non_symmetric(A, Bz, Z);
//    std::cerr << "after solve\n";

    // copy to mesh
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      vertex_descriptor vi = *vb;
      int id = boost::get(vertex_id_pmap, vi);
      int i = new_id[id];
      Point p(X[i], Y[i], Z[i]);
      vi->point() = p;
    }
//    std::cerr << "leave contract geometry\n";
  }

  int collapse_short_edges()
  {
    Constrains_map constrains_map;

    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, *polyhedron);
      vertex_descriptor vj = boost::target(*eb, *polyhedron);
      size_t vi_idx = boost::get(vertex_id_pmap, vi);
      size_t vj_idx = boost::get(vertex_id_pmap, vj);

      if (is_vertex_fixed_map.find(vi_idx) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[vi_idx])
        {
          constrains_map.set_is_constrained(*eb, true);
        }
      }
      if (is_vertex_fixed_map.find(vj_idx) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[vj_idx])
        {
          constrains_map.set_is_constrained(*eb, true);
        }
      }

//      if (is_vertex_fixed_map.find(vi_idx) != is_vertex_fixed_map.end()
//       && is_vertex_fixed_map.find(vj_idx) != is_vertex_fixed_map.end())
//      {
//        if (is_vertex_fixed_map[vi_idx] && is_vertex_fixed_map[vj_idx])
//        {
//          constrains_map.set_is_constrained(*eb, true);
//        }
//      }

    }

    int edge_id = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, edge_id++);
    }

    // This is a stop predicate (defines when the algorithm terminates).
    // The simplification stops when the length of all edges is greater than the minimum threshold.
    CGAL::internal::Minimum_length_predicate<Polyhedron> stop(edgelength_TH);

    // midpoint placement without geometric test
    SMS::Geometric_test_skipper< SMS::Midpoint_placement<Polyhedron> > placement;

    Track_vertex_visitor vis(&correspondence, max_id);

    int r = SMS::edge_collapse
                (*polyhedron
                ,stop
                ,CGAL::get_cost(SMS::Edge_length_cost<Polyhedron>())
                      .get_placement(placement)
                      .visitor(vis)
                      .edge_is_border_map(constrains_map)
                );

    return r;
  }

  int iteratively_collapse_edges()
  {
    int num_collapses = 0;
    while (true)
    {
      int cnt = collapse_short_edges();
      if (cnt == 0)
      {
        break;
      }
      else
      {
//        std::cerr << "collapse " << cnt << "\n";
        num_collapses += cnt;
      }
    }
    return num_collapses;
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

      // for border edge, the angle is -1
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
        if (dis_ij < zero_TH || dis_ik < zero_TH || dis_jk < zero_TH)
        {
          halfedge_angle[e_id] = -1;
        }
        else
        {
          halfedge_angle[e_id] = acos((dis2_ik + dis2_jk - dis2_ij) / (2.0 * dis_ik * dis_jk));
          if (halfedge_angle[e_id] > M_PI)
          {
            std::cerr << "angle too large\n";
          }
          if (halfedge_angle[e_id] < 0)
          {
            std::cerr << "angle too small\n";
          }
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

  int split_flat_triangle()
  {
//    std::cerr << "TH_ALPHA " << TH_ALPHA << "\n";
//    std::cerr << "short " << zero_TH << "\n";
    int ne = boost::num_edges(*polyhedron);
    compute_incident_angle();

    int cnt = 0;
    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ei = *eb;
      edge_descriptor ej = ei->opposite();
      int ei_id = boost::get(edge_id_pmap, ei);
      int ej_id = boost::get(edge_id_pmap, ej);
      if (ei_id < 0 || ei_id >= ne
       || ej_id < 0 || ej_id >= ne)
      {
        continue;
      }

      vertex_descriptor vs = boost::source(ei, *polyhedron);
      vertex_descriptor vt = boost::target(ei, *polyhedron);
      size_t vs_id = boost::get(vertex_id_pmap, vs);
      size_t vt_id = boost::get(vertex_id_pmap, vt);
      Point ps = vs->point();
      Point pt = vt->point();

      if (is_vertex_fixed_map.find(vs_id) != is_vertex_fixed_map.end()
       || is_vertex_fixed_map.find(vt_id) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[vs_id] || is_vertex_fixed_map[vt_id])
        {
          continue;
        }
      }

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
      cnt++;
    }
    return cnt;
  }

  int iteratively_split_triangles()
  {
//    std::cerr << "before split\n";
    int num_splits = 0;
    while (true)
    {
      int cnt = split_flat_triangle();
      if (cnt == 0)
      {
        break;
      }
      else
      {
//        std::cerr << "split " << cnt << "\n";
        num_splits += cnt;
      }
    }
//    std::cerr << "after split\n";
    return num_splits;
  }

  int update_topology()
  {
//    std::cerr << "before collapse edges\n";
    int num_collapses = iteratively_collapse_edges();
    std::cerr << "collapse " << num_collapses << " edges.\n";

    int num_splits = iteratively_split_triangles();
    std::cerr << "split " << num_splits << " edges.\n";

    return num_collapses + num_splits;
  }

  bool is_vertex_degenerate(vertex_descriptor root)
  {
    std::set<edge_descriptor> edge_visited;
    std::map<vertex_descriptor, int> vertex_visited;

    std::map<vertex_descriptor, int> D;
    std::queue<vertex_descriptor> Q;
    Q.push(root);
    D[root] = 0;
    vertex_visited[root] = 0;

    int dist_v;
    double max_distance = 0.0;
    // size of k-ring
    int k = 2;
    while (!Q.empty() && (dist_v = D[Q.front()]) < k)
    {
      vertex_descriptor v = Q.front();
      Q.pop();

      out_edge_iterator e, e_end;
      for(boost::tie(e, e_end) = boost::out_edges(v, *polyhedron); e != e_end; e++)
      {
        edge_descriptor ed = *e;
        if (edge_visited.find(ed) != edge_visited.end())
        {
          continue;
        }

        vertex_descriptor new_v = boost::target(ed, *polyhedron);
        if (vertex_visited.find(new_v) != vertex_visited.end())
        {
          if (vertex_visited[new_v] != dist_v)
          {
//            std::cerr << vertex_visited[new_v] << " " << dist_v << "\n";
            return true;
          }
        }
        edge_visited.insert(ed);
        edge_visited.insert(ed->opposite());
        edge_visited.insert(ed->next());
        edge_visited.insert(ed->next()->opposite());
        vertex_visited[new_v] = dist_v;

        if (D.insert(std::make_pair(new_v, dist_v + 1)).second)
        {
          max_distance = (std::max)((new_v->point() - root->point()).squared_length(), max_distance);
          Q.push(new_v);
        }
      }
    }
    // now Q contains all nonprocessed
//    while (!Q.empty())
//    {
//      vertex_descriptor v = Q.front();
//      Q.pop();

//      out_edge_iterator e, e_end;
//      for (boost::tie(e, e_end) = boost::out_edges(v, *polyhedron); e != e_end; e++)
//      {
//        vertex_descriptor new_v = boost::target(*e, *polyhedron);
//        double distance = (new_v->point() - root->point()).squared_length();
//        if (distance < max_distance)
//        {
//          if (D.insert(std::make_pair(new_v, dist_v + 1)).second)
//          {
//            Q.push(new_v);
//          }
//        }
//      }
//    }
    return false;
  }

  int detect_degeneracies_in_disk()
  {
    int num_fixed = 0;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);

      if (is_vertex_fixed_map.find(idx) == is_vertex_fixed_map.end() || !is_vertex_fixed_map[idx])
      {
        bool willbefixed = is_vertex_degenerate(v);
        if (willbefixed)
        {
//          std::cerr << "detect " << idx << "\n";
          is_vertex_fixed_map[idx] = willbefixed;
          num_fixed++;
        }
      }
    }
    std::cerr << "fixed " << num_fixed << " vertices.\n";
    return num_fixed;
  }

  // TODO: check if the local neighborhood is a disk
  int detect_degeneracies()
  {
    int num_fixed = 0;
    double elength_fixed = edgelength_TH;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);
//      std::cerr << v->point() << "\n";
      if (is_vertex_fixed_map.find(idx) == is_vertex_fixed_map.end() || !is_vertex_fixed_map[idx])
      {
        bool willbefixed = false;
        int bad_counter = 0;

        in_edge_iterator eb, ee;
        for (boost::tie(eb, ee) = boost::in_edges(v, *polyhedron); eb != ee; eb++)
        {
          edge_descriptor edge = *eb;
          vertex_descriptor v0 = boost::source(edge, *polyhedron);
          vertex_descriptor v1 = boost::target(edge, *polyhedron);
          double length = sqrt(squared_distance(v0->point(), v1->point()));
//          std::cerr << length << "\n";
          if (length < elength_fixed)
          {
            if (!is_collapse_ok(edge))
            {
              bad_counter++;
            }
          }
        }
//        std::cerr << "bad " << bad_counter << "\n";
        willbefixed = (bad_counter >= 2);
        if (willbefixed)
        {
//          std::cerr << "detect " << idx << "\n";
          is_vertex_fixed_map[idx] = willbefixed;
          num_fixed++;
        }
      }
    }
    std::cerr << "fixed " << num_fixed << " vertices.\n";
    return num_fixed;
  }

  bool is_collapse_ok(edge_descriptor v0v1)
  {
    edge_descriptor v1v0 = v0v1->opposite();
    vertex_descriptor v0 = boost::target(v1v0, *polyhedron);
    vertex_descriptor v1 = boost::source(v1v0, *polyhedron);

    vertex_descriptor vv, vl, vr;
    edge_descriptor  h1, h2;

    // the edges v1-vl and vl-v0 must not be both boundary edges
    if (!(v0v1->is_border()))
    {
      vl = boost::target(v0v1->next(), *polyhedron);
      h1 = v0v1->next();
      h2 = h1->next();
      if (h1->opposite()->is_border() && h2->opposite()->is_border())
      {
        return false;
      }
    }

    // the edges v0-vr and vr-v1 must not be both boundary edges
    if (!(v1v0->is_border()))
    {
      vr = boost::target(v1v0->next(), *polyhedron);
      h1 = v1v0->next();
      h2 = h1->next();
      if (h1->opposite()->is_border() && h2->opposite()->is_border())
      {
        return false;
      }
    }

    // if vl and vr are equal or both invalid -> fail
    if (vl == vr)
    {
      return false;
    }

    // edge between two boundary vertices should be a boundary edge
    if (is_border(v0) && is_border(v1) &&
        !(v0v1->is_border()) && !(v1v0->is_border()))
    {
      return false;
    }

    // test intersection of the one-rings of v0 and v1
    in_edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::in_edges(v0, *polyhedron); eb != ee; eb++)
    {
      vv = boost::source(*eb, *polyhedron);
      if (vv != v1 && vv != vl && vv != vr)
      {
        if (find_halfedge(vv, v1))
        {
          return false;
        }
      }
    }

    // passed all tests
    return true;
  }

  bool find_halfedge(vertex_descriptor vi, vertex_descriptor vj)
  {
    in_edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::in_edges(vj, *polyhedron); eb != ee; eb++)
    {
      vertex_descriptor vv = boost::source(*eb, *polyhedron);
      if (vv == vi)
      {
        return true;
      }
    }
    return false;
  }

  bool is_border(vertex_descriptor aV)
  {
    bool rR = false;

    in_edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::in_edges(aV, *polyhedron); eb != ee; eb++)
    {
      edge_descriptor lEdge = *eb;
      if (is_undirected_edge_a_border(lEdge))
      {
        rR = true;
        break;
      }
    }

    return rR;
  }

  bool is_undirected_edge_a_border(edge_descriptor aEdge)
  {
    return aEdge->is_border() || aEdge->opposite()->is_border();
  }

  void contract()
  {
    contract_geometry();
    update_topology();
    detect_degeneracies();
//    double area = get_surface_area();
//    std::cout << "area " << area << "\n";
//    detect_degeneracies_in_disk();
  }

  void run_to_converge()
  {
//    double last_area = 0;
    while (true)
    {
      contract_geometry();
      int num_events = update_topology();
      detect_degeneracies();
//      double area = get_surface_area();
//      std::cout << "area " << area << "\n";
//      if (fabs(last_area - area) < area_TH)
//      {
//        break;
//      }
      if (num_events == 0)
      {
        break;
      }
//      last_area = area;
    }
  }

  void convert_to_skeleton()
  {
    Skeleton skeleton(polyhedron);
    std::vector<std::vector<int> > record;
    skeleton.extract_skeleton(g, points, record);

    skeleton_to_surface.resize(record.size());
    for (size_t i = 0; i < record.size(); i++)
    {
      for (size_t j = 0; j < record[i].size(); j++)
      {
        int id = record[i][j];
        if (correspondence.find(id) != correspondence.end())
        {
          skeleton_to_surface[i].insert(skeleton_to_surface[i].end(),
                                        correspondence[id].begin(),
                                        correspondence[id].end());
        }

        if (id < max_id)
        {
          skeleton_to_surface[i].push_back(id);
        }
      }
    }
    int cnt = 0;
    for (size_t i = 0; i < skeleton_to_surface.size(); i++)
    {
      cnt += skeleton_to_surface[i].size();
    }
    std::cout << "tracked " << cnt << " vertices\n";
  }

  void get_skeleton(Graph& g, std::vector<Point>& points)
  {
    g = this->g;
    points = this->points;
  }

  void get_correspondent_vertices(std::vector<std::vector<int> >& corr)
  {
    corr = skeleton_to_surface;
  }
};

} //namespace CGAL

#endif // MEAN_CURVATURE_SKELETON_H
