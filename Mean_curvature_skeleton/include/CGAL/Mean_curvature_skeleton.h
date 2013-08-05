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

// Compute the vertex normal
#include <CGAL/internal/Mean_curvature_skeleton/get_normal.h>

// Low level collapse function
#include <CGAL/Surface_mesh_simplification/halfedge_collapse_Polyhedron_3.h>

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

// For Voronoi diagram
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

// For debugging macro
#include <CGAL/internal/Mean_curvature_skeleton/Debug.h>

// For mesh_split
#include <CGAL/internal/Mean_curvature_skeleton/Utility.h>

#include <queue>

namespace SMS = CGAL::Surface_mesh_simplification;

namespace CGAL {

template <class Polyhedron,
          class SparseLinearAlgebraTraits_d,
          class PolyhedronVertexIndexMap,
          class PolyhedronEdgeIndexMap,
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

  typedef typename Polyhedron::Face_handle                                     Face_handle;
  typedef typename Polyhedron::Facet_iterator                                  Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator                Halfedge_facet_circulator;

  // Cotangent weight calculator
  typedef typename internal::Cotangent_weight<Polyhedron,
  internal::Cotangent_value_minimum_zero<Polyhedron,
  internal::Cotangent_value_Meyer_secure<Polyhedron> > >                       Weight_calculator;

  // Repeat Graph types
  typedef typename boost::graph_traits<Graph>::in_edge_iterator                in_edge_iter;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator               out_edge_iter;
  typedef typename boost::graph_traits<Graph>::edge_iterator                   edge_iter;
  typedef typename boost::graph_traits<Graph>::edge_descriptor                 edge_desc;

  // Skeleton types
  typedef Curve_skeleton<Polyhedron, Graph,
  PolyhedronVertexIndexMap, PolyhedronEdgeIndexMap>                            Skeleton;

  // Mesh simplification types
  typedef SMS::Edge_profile<Polyhedron>                                        Profile;

  // Repeat Triangulation types
  typedef CGAL::Exact_predicates_exact_constructions_kernel                    K;
  typedef K::Vector_3                                                          Exact_vector;
  typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>             Vb;
  typedef CGAL::Triangulation_data_structure_3<Vb>                             Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds>                               Delaunay;
  typedef Delaunay::Point                                                      Exact_point;
  typedef Delaunay::Cell_handle                                                Cell_handle;
  typedef Delaunay::Vertex_handle                                              TriVertex_handle;
  typedef Delaunay::Locate_type                                                Locate_type;
  typedef Delaunay::Finite_vertices_iterator                                   Finite_vertices_iterator;
  typedef Delaunay::Finite_edges_iterator                                      Finite_edges_iterator;
  typedef Delaunay::Finite_facets_iterator                                     Finite_facets_iterator;
  typedef Delaunay::Finite_cells_iterator                                      Finite_cells_iterator;

// Data members
private:

  Polyhedron& polyhedron;
  PolyhedronVertexIndexMap vertex_id_pmap;
  PolyhedronEdgeIndexMap edge_id_pmap;

  double omega_H;
  double omega_P;
  double edgelength_TH;
  double TH_ALPHA;
  double zero_TH;
  double area_TH;
  double original_area;
  int iteration_TH;

  Weight_calculator weight_calculator;
  std::vector<double> edge_weight;
  SparseLinearAlgebraTraits_d m_solver;

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

  bool is_medially_centered;
  // record the corresponding pole of a point
  std::map<int, int> poles;
  // the normal of surface points
  std::vector<Vector> normals;
  // the dual of a cell in Triangulation(a Voronoi point)
  std::vector<Point> cell_dual;

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
    Track_vertex_visitor(){}

    Track_vertex_visitor(std::map<int, std::vector<int> >* corr, int max_id) :
      corr(corr), max_id(max_id), is_medially_centered(false){}

    Track_vertex_visitor(std::map<int, std::vector<int> >* corr,
                         std::map<int, int>* poles,
                         std::vector<Point>* cell_dual,
                         int max_id) :
      corr(corr), max_id(max_id), is_medially_centered(true),
      poles(poles), cell_dual(cell_dual)
      {}

    // Called AFTER each edge has been collapsed
    void OnCollapsed(Profile const& edge, Vertex_handle v)
    {
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
        for (size_t i = 0; i < (iter->second).size(); ++i)
        {
          (*corr)[to].push_back((iter->second)[i]);
        }
        (iter->second).clear();
        (*corr).erase(iter);
      }

      if (is_medially_centered)
      {
        Point pole0 = Point(to_double((*cell_dual)[(*poles)[id0]].x()),
                            to_double((*cell_dual)[(*poles)[id0]].y()),
                            to_double((*cell_dual)[(*poles)[id0]].z()));
        Point pole1 = Point(to_double((*cell_dual)[(*poles)[id1]].x()),
                            to_double((*cell_dual)[(*poles)[id1]].y()),
                            to_double((*cell_dual)[(*poles)[id1]].z()));
        Point p1 = v1->point();
        double dis_to_pole0 = sqrt(squared_distance(pole0, p1));
        double dis_to_pole1 = sqrt(squared_distance(pole1, p1));
        if (dis_to_pole0 < dis_to_pole1)
        {
          (*poles)[id1] = (*poles)[id0];
        }
        std::map<int, int>::iterator pole_iter = (*poles).find(id0);
        (*poles).erase(pole_iter);
      }
    }

    std::map<int, std::vector<int> >* corr;
    int max_id;

    bool is_medially_centered;
    std::map<int, int>* poles;
    std::vector<Point>* cell_dual;
  };

// Public methods
public:

  // The constructor gets the polyhedron that we will model
  Mean_curvature_skeleton(Polyhedron& P,
                          PolyhedronVertexIndexMap Vertex_index_map,
                          PolyhedronEdgeIndexMap Edge_index_map,
                          double omega_H,
                          double edgelength_TH,
                          double area_TH = 1e-5,
                          double zero_TH = 1e-7,
                          Weight_calculator weight_calculator = Weight_calculator()
                          )
    :polyhedron(P),
     vertex_id_pmap(Vertex_index_map),
     edge_id_pmap(Edge_index_map),
     omega_H(omega_H),
     edgelength_TH(edgelength_TH),
     TH_ALPHA(110),
     zero_TH(zero_TH),
     area_TH(area_TH),
     weight_calculator(weight_calculator),
     is_medially_centered(false)
  {
    TH_ALPHA *= (M_PI / 180.0);
    double area = get_surface_area();
    area_TH = 0.0001 * area;
    original_area = get_surface_area();
    iteration_TH = 500;

    // initialize index maps
    vertex_iterator vb, ve;
    vertex_id_count = 0;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, vertex_id_count++);
    }
    max_id = vertex_id_count;

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    is_vertex_fixed_map.clear();
    correspondence.clear();

    if (is_medially_centered)
    {
      compute_voronoi_pole();
    }
  }

  Mean_curvature_skeleton(Polyhedron& P,
                          PolyhedronVertexIndexMap Vertex_index_map,
                          PolyhedronEdgeIndexMap Edge_index_map,
                          double omega_H,
                          double omega_P,
                          double edgelength_TH,
                          bool is_medially_centered,
                          double area_TH = 1e-5,
                          double zero_TH = 1e-7,
                          Weight_calculator weight_calculator = Weight_calculator()
                          )
    :polyhedron(P),
     vertex_id_pmap(Vertex_index_map),
     edge_id_pmap(Edge_index_map),
     omega_H(omega_H),
     omega_P(omega_P),
     edgelength_TH(edgelength_TH),
     TH_ALPHA(110),
     zero_TH(zero_TH),
     area_TH(area_TH),
     weight_calculator(weight_calculator),
     is_medially_centered(is_medially_centered)
  {
    TH_ALPHA *= (M_PI / 180.0);
    double area = get_surface_area();
    area_TH = 0.0001 * area;
    original_area = get_surface_area();
    iteration_TH = 500;

    // initialize index maps
    vertex_iterator vb, ve;
    vertex_id_count = 0;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, vertex_id_count++);
    }
    max_id = vertex_id_count;

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    is_vertex_fixed_map.clear();
    correspondence.clear();

    if (is_medially_centered)
    {
      compute_voronoi_pole();
    }
  }

  // Release resources
  ~Mean_curvature_skeleton(void)
  {
  }

  void set_omega_H(double value)
  {
    omega_H = value;
  }

  void set_omega_P(double value)
  {
    omega_P = value;
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

  void set_medially_centered(bool value)
  {
    is_medially_centered = value;
  }

  Polyhedron& get_polyhedron()
  {
    return polyhedron;
  }

  void get_fixed_points(std::vector<Point>& fixed_points)
  {
    fixed_points.clear();
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        vertex_descriptor vd = *vb;
        fixed_points.push_back(vd->point());
      }
    }
  }

  void get_non_fixed_points(std::vector<Point>& non_fixed_points)
  {
    non_fixed_points.clear();
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
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
    for (Facet_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
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
    edge_weight.reserve(boost::num_edges(polyhedron));
    edge_iterator eb, ee;
    for(boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      edge_weight.push_back(this->weight_calculator(*eb, polyhedron));
    }
  }

  void assemble_LHS(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    MCFSKEL_DEBUG(std::cerr << "start LHS\n";)

    int nver = boost::num_vertices(polyhedron);

    vertex_iterator vb, ve;
    // initialize the Laplacian matrix
    int cnt_fix = 0;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);

      int i = new_id[id];
      // if the vertex is fixed
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        cnt_fix++;
        A.set_coef(i + nver, i, 1.0 / zero_TH, true);
      }
      else
      {
        A.set_coef(i + nver, i, omega_H, true);
        if (is_medially_centered)
        {
          if (id < max_id)
          {
            A.set_coef(i + nver * 2, i, omega_P, true);
          }
        }
      }
    }

    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      int i = new_id[id];
      double L = 1.0;
      // if the vertex is fixed
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        L = 0;
      }
      double diagonal = 0;
      in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = boost::in_edges(*vb, polyhedron); e != e_end; ++e)
      {
        vertex_descriptor vj = boost::source(*e, polyhedron);
        double wij = edge_weight[boost::get(edge_id_pmap, *e)] * 2.0;
        int jd = boost::get(vertex_id_pmap, vj);
        int j = new_id[jd];
        A.set_coef(i, j, wij * L, true);
        diagonal += -wij;
      }
      A.set_coef(i, i, diagonal, true);
    }

    MCFSKEL_DEBUG(std::cerr << "end LHS\n";)
  }

  void assemble_RHS(typename SparseLinearAlgebraTraits_d::Vector& Bx,
                    typename SparseLinearAlgebraTraits_d::Vector& By,
                    typename SparseLinearAlgebraTraits_d::Vector& Bz)
  {
    MCFSKEL_DEBUG(std::cerr << "start RHS\n";)
    // assemble right columns of linear system
    int nver = boost::num_vertices(polyhedron);
    vertex_iterator vb, ve;
    for (int i = 0; i < nver; ++i)
    {
      Bx[i] = 0;
      By[i] = 0;
      Bz[i] = 0;
    }

    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor vi = *vb;
      int id = boost::get(vertex_id_pmap, vi);
      int i = new_id[id];

      double oh, op = 0.0;
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        oh = 1.0 / zero_TH;
      }
      else
      {
        oh = omega_H;
        if (id < max_id)
        {
          op = omega_P;
        }
      }
      Bx[i + nver] = vi->point().x() * oh;
      By[i + nver] = vi->point().y() * oh;
      Bz[i + nver] = vi->point().z() * oh;
      if (is_medially_centered)
      {
        double x = to_double(cell_dual[poles[id]].x());
        double y = to_double(cell_dual[poles[id]].y());
        double z = to_double(cell_dual[poles[id]].z());
        Bx[i + nver * 2] = x * op;
        By[i + nver * 2] = y * op;
        Bz[i + nver * 2] = z * op;
      }
    }
    MCFSKEL_DEBUG(std::cerr << "end RHS\n";)
  }

  void update_vertex_id()
  {
    new_id.clear();
    int cnt = 0;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      new_id[id] = cnt++;
    }
  }

  void contract_geometry()
  {
    MCFSKEL_DEBUG(std::cerr << "before contract geometry";)

    update_vertex_id();

    compute_edge_weight();

    int nver = boost::num_vertices(polyhedron);
    int nrows;
    if (is_medially_centered)
    {
      nrows = nver * 3;
    }
    else
    {
      nrows = nver * 2;
    }
    // Assemble linear system At * A * X = At * B
    typename SparseLinearAlgebraTraits_d::Matrix A(nrows, nver);
    assemble_LHS(A);

    typename SparseLinearAlgebraTraits_d::Vector X(nver), Bx(nrows);
    typename SparseLinearAlgebraTraits_d::Vector Y(nver), By(nrows);
    typename SparseLinearAlgebraTraits_d::Vector Z(nver), Bz(nrows);
    assemble_RHS(Bx, By, Bz);

    MCFSKEL_DEBUG(std::cerr << "before solve\n";)

    // solve "At * A * X = At * B".
    double D;
    m_solver.pre_factor_non_symmetric(A, D);
    m_solver.linear_solver_non_symmetric(A, Bx, X);
    m_solver.linear_solver_non_symmetric(A, By, Y);
    m_solver.linear_solver_non_symmetric(A, Bz, Z);

    MCFSKEL_DEBUG(std::cerr << "after solve\n";)

    // copy to mesh
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor vi = *vb;
      int id = boost::get(vertex_id_pmap, vi);
      int i = new_id[id];
      Point p(X[i], Y[i], Z[i]);
      vi->point() = p;
    }

    MCFSKEL_DEBUG(std::cerr << "leave contract geometry\n";)
  }

  int collapse_short_edges()
  {
    Constrains_map constrains_map;

    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, polyhedron);
      vertex_descriptor vj = boost::target(*eb, polyhedron);
      size_t vi_idx = boost::get(vertex_id_pmap, vi);
      size_t vj_idx = boost::get(vertex_id_pmap, vj);

      if (is_vertex_fixed_map.find(vi_idx) != is_vertex_fixed_map.end()
       && is_vertex_fixed_map.find(vj_idx) != is_vertex_fixed_map.end())
      {
        constrains_map.set_is_constrained(*eb, true);
      }
    }

    int edge_id = 0;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, edge_id++);
    }

    // This is a stop predicate (defines when the algorithm terminates).
    // The simplification stops when the length of all edges is greater than the minimum threshold.
    CGAL::internal::Minimum_length_predicate<Polyhedron> stop(edgelength_TH);

    // midpoint placement without geometric test
    SMS::Geometric_test_skipper< SMS::Midpoint_placement<Polyhedron> > placement;

    Track_vertex_visitor vis;
    if (is_medially_centered)
    {
      vis = Track_vertex_visitor(&correspondence, &poles, &cell_dual, max_id);
    }
    else
    {
      vis = Track_vertex_visitor(&correspondence, max_id);
    }

    int r = SMS::edge_collapse
                (polyhedron
                ,stop
                ,CGAL::get_cost(SMS::Edge_length_cost<Polyhedron>())
                      .get_placement(placement)
                      .visitor(vis)
                      .edge_is_border_map(constrains_map)
                );

    return r;
  }

  void track_correspondence(vertex_descriptor v0, vertex_descriptor v1,
                            vertex_descriptor v)
  {
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

    if (correspondence.find(to) == correspondence.end())
    {
      correspondence[to] = std::vector<int>();
    }
    // only track vertex in original mesh
    if (from < max_id)
    {
      correspondence[to].push_back(from);
    }
    std::map<int, std::vector<int> >::iterator iter = correspondence.find(from);
    if (iter != correspondence.end())
    {
      for (size_t i = 0; i < (iter->second).size(); ++i)
      {
        correspondence[to].push_back((iter->second)[i]);
      }
      (iter->second).clear();
      correspondence.erase(iter);
    }

    if (is_medially_centered)
    {
      Point pole0 = Point(to_double(cell_dual[poles[id0]].x()),
                          to_double(cell_dual[poles[id0]].y()),
                          to_double(cell_dual[poles[id0]].z()));
      Point pole1 = Point(to_double(cell_dual[poles[id1]].x()),
                          to_double(cell_dual[poles[id1]].y()),
                          to_double(cell_dual[poles[id1]].z()));
      Point p1 = v1->point();
      double dis_to_pole0 = sqrt(squared_distance(pole0, p1));
      double dis_to_pole1 = sqrt(squared_distance(pole1, p1));
      if (dis_to_pole0 < dis_to_pole1)
      {
        poles[id1] = poles[id0];
      }
      std::map<int, int>::iterator pole_iter = poles.find(id0);
      poles.erase(pole_iter);
    }
  }

  int collapse_edges(Constrains_map& constrains_map)
  {
    std::vector<edge_descriptor> edges;
    edges.reserve(boost::num_edges(polyhedron));
    edge_iterator eb, ee;

    boost::tie(eb, ee) = boost::edges(polyhedron);
    std::copy(eb, ee, std::back_inserter(edges));

    int cnt = 0;
    for (size_t i = 0; i < edges.size(); ++i)
    {
      edge_descriptor ed = edges[i];
      if (constrains_map.is_constrained(ed))
      {
        continue;
      }

      vertex_descriptor vi = boost::source(ed, polyhedron);
      vertex_descriptor vj = boost::target(ed, polyhedron);
      double edge_length = sqrt(squared_distance(vi->point(), vj->point()));
      if (is_collapse_ok(ed) && edge_length < edgelength_TH)
      {
        Point p = midpoint(
          boost::get(vertex_point, polyhedron, boost::source(ed, polyhedron)),
          boost::get(vertex_point, polyhedron, boost::target(ed, polyhedron)));

        // invalidate the edges that will be collapsed
        // since the mesh is closed, 6 halfedges will be collapsed
        constrains_map.set_is_constrained(ed, true);
        constrains_map.set_is_constrained(ed->opposite(), true);
        constrains_map.set_is_constrained(ed->prev(), true);
        constrains_map.set_is_constrained(ed->prev()->opposite(), true);
        constrains_map.set_is_constrained(ed->opposite()->prev(), true);
        constrains_map.set_is_constrained(ed->opposite()->prev()->opposite(), true);

        vertex_descriptor v = Surface_mesh_simplification::halfedge_collapse(ed, polyhedron);
        boost::put(vertex_point, polyhedron, v, p);

        track_correspondence(vi, vj, v);

        cnt++;
      }
    }

    return cnt;
  }

  void init_constraint_map(Constrains_map& constrains_map)
  {
    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, polyhedron);
      vertex_descriptor vj = boost::target(*eb, polyhedron);
      size_t vi_idx = boost::get(vertex_id_pmap, vi);
      size_t vj_idx = boost::get(vertex_id_pmap, vj);

      if (is_vertex_fixed_map.find(vi_idx) != is_vertex_fixed_map.end()
       && is_vertex_fixed_map.find(vj_idx) != is_vertex_fixed_map.end())
      {
        constrains_map.set_is_constrained(*eb, true);
      }
    }
  }

  int iteratively_collapse_edges()
  {
    Constrains_map constrains_map;
    init_constraint_map(constrains_map);

    int num_collapses = 0;
    while (true)
    {
//      int cnt = collapse_short_edges();
      int cnt = collapse_edges(constrains_map);
      if (cnt == 0)
      {
        break;
      }
      else
      {
        num_collapses += cnt;
      }
    }
    return num_collapses;
  }

  void compute_incident_angle()
  {
    halfedge_angle.clear();
    int ne = boost::num_edges(polyhedron);
    halfedge_angle.resize(ne, 0);

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
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
        vertex_descriptor vi = boost::source(ed, polyhedron);
        vertex_descriptor vj = boost::target(ed, polyhedron);
        edge_descriptor ed_next = ed->next();
        vertex_descriptor vk = boost::target(ed_next, polyhedron);
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
        }
      }
    }
  }

  Point project_vertex(const vertex_descriptor vs,
                       const vertex_descriptor vt,
                       const vertex_descriptor vk)
  {
    Point ps = vs->point();
    Point pt = vt->point();
    Point pk = vk->point();
    CGAL::internal::Vector vec_st = CGAL::internal::Vector(ps, pt);
    CGAL::internal::Vector vec_sk = CGAL::internal::Vector(ps, pk);

    vec_st.normalize();
    double t = vec_st.dot(vec_sk);
    Point st = Point(vec_st[0] * t, vec_st[1] * t, vec_st[2] * t);
    Point pn = Point(ps[0] + st[0], ps[1] + st[1], ps[2] + st[2]);

    // project the pole
    if (is_medially_centered)
    {
      int sid = boost::get(vertex_id_pmap, vs);
      int tid = boost::get(vertex_id_pmap, vt);
      Point pole_s = cell_dual[poles[sid]];
      Point pole_t = cell_dual[poles[tid]];
      Vector pole_st = pole_t - pole_s;
      Vector p_projector = pole_st / sqrt(pole_st.squared_length());
      Point pole_n = pole_s + p_projector * t;
      poles[vertex_id_count] = cell_dual.size();
      cell_dual.push_back(pole_n);
    }
    return pn;
  }

  int split_flat_triangle()
  {
    int ne = boost::num_edges(polyhedron);
    compute_incident_angle();

    int cnt = 0;
    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
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

      vertex_descriptor vs = boost::source(ei, polyhedron);
      vertex_descriptor vt = boost::target(ei, polyhedron);

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
      vertex_descriptor vk = boost::target(ek, polyhedron);
      Point pn = project_vertex(vs, vt, vk);
      edge_descriptor en = CGAL::internal::mesh_split(polyhedron, ei, pn);
      // set id for new vertex
      boost::put(vertex_id_pmap, en->vertex(), vertex_id_count++);
      cnt++;
    }
    return cnt;
  }

  int iteratively_split_triangles()
  {
    MCFSKEL_DEBUG(std::cerr << "before split\n";)

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
        num_splits += cnt;
      }
    }

    MCFSKEL_DEBUG(std::cerr << "after split\n";)

    return num_splits;
  }

  int update_topology()
  {
    MCFSKEL_DEBUG(std::cerr << "before collapse edges\n";)

    int num_collapses = iteratively_collapse_edges();
    MCFSKEL_INFO(std::cerr << "collapse " << num_collapses << " edges.\n";)

    int num_splits = iteratively_split_triangles();
    MCFSKEL_INFO(std::cerr << "split " << num_splits << " edges.\n";)

    return num_collapses + num_splits;
  }

  bool is_vertex_degenerate(vertex_descriptor root)
  {
    std::set<vertex_descriptor> vertices_in_disk;
    std::set<edge_descriptor> edges_in_disk;
    std::set<Face_handle> faces_in_disk;

    vertices_in_disk.clear();
    search_vertices_in_disk(root, vertices_in_disk);

    typename std::set<vertex_descriptor>::iterator v_iter;
    for (v_iter = vertices_in_disk.begin(); v_iter != vertices_in_disk.end(); ++v_iter)
    {
      vertex_descriptor vd = *v_iter;
      out_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = boost::out_edges(vd, polyhedron); e != e_end; ++e)
      {
        edge_descriptor ed = *e;
        edge_descriptor ed_op = ed->opposite();
        vertex_descriptor target = boost::target(ed, polyhedron);
        if (vertices_in_disk.find(target) != vertices_in_disk.end())
        {
          edges_in_disk.insert(ed);
          edges_in_disk.insert(ed_op);
        }
        Face_handle f = ed->face();
        Halfedge_facet_circulator j = f->facet_begin();
        bool in = true;
        do
        {
          vertex_descriptor v = j->vertex();
          if (vertices_in_disk.find(v) == vertices_in_disk.end())
          {
            in = false;
            break;
          }
        } while (++j != f->facet_begin());

        if (in)
        {
          faces_in_disk.insert(f);
        }
      }
    }

    int V = vertices_in_disk.size();
    int E = edges_in_disk.size() / 2;
    int F = faces_in_disk.size();
    int euler = V + F - E;
    if (euler != 1)
    {
      return true;
    }
    return false;
  }

  void search_vertices_in_disk(vertex_descriptor root,
                               std::set<vertex_descriptor>& vertices_in_disk)
  {
    std::map<vertex_descriptor, bool> vertex_visited;

    std::queue<vertex_descriptor> Q;
    Q.push(root);
    vertices_in_disk.insert(root);
    vertex_visited[root] = true;

    double dist_TH = edgelength_TH;
    while (!Q.empty())
    {
      vertex_descriptor v = Q.front();
      Q.pop();

      out_edge_iterator e, e_end;
      for(boost::tie(e, e_end) = boost::out_edges(v, polyhedron); e != e_end; ++e)
      {
        edge_descriptor ed = *e;

        vertex_descriptor new_v = boost::target(ed, polyhedron);
        if (vertex_visited.find(new_v) == vertex_visited.end())
        {
          double distance = sqrt(squared_distance(new_v->point(), root->point()));
          if (distance < dist_TH)
          {
            vertex_visited[new_v] = true;
            Q.push(new_v);
            vertices_in_disk.insert(new_v);
          }
        }
      }
    }
  }

  int detect_degeneracies_in_disk()
  {
    int num_fixed = 0;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);

      if (is_vertex_fixed_map.find(idx) == is_vertex_fixed_map.end())
      {
        bool willbefixed = is_vertex_degenerate(v);
        if (willbefixed)
        {
          is_vertex_fixed_map[idx] = willbefixed;
          num_fixed++;
        }
      }
    }

    MCFSKEL_INFO(std::cerr << "fixed " << num_fixed << " vertices.\n";)

    return num_fixed;
  }

  int detect_degeneracies()
  {
    int num_fixed = 0;
    double elength_fixed = edgelength_TH;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);
      if (is_vertex_fixed_map.find(idx) == is_vertex_fixed_map.end())
      {
        bool willbefixed = false;
        int bad_counter = 0;

        in_edge_iterator eb, ee;
        for (boost::tie(eb, ee) = boost::in_edges(v, polyhedron); eb != ee; ++eb)
        {
          edge_descriptor edge = *eb;
          vertex_descriptor v0 = boost::source(edge, polyhedron);
          vertex_descriptor v1 = boost::target(edge, polyhedron);
          double length = sqrt(squared_distance(v0->point(), v1->point()));
          if (length < elength_fixed)
          {
            if (!is_collapse_ok(edge))
            {
              bad_counter++;
            }
          }
        }
        willbefixed = (bad_counter >= 2);
        if (willbefixed)
        {
          is_vertex_fixed_map[idx] = willbefixed;
          num_fixed++;
        }
      }
    }

    MCFSKEL_INFO(std::cerr << "fixed " << num_fixed << " vertices.\n";)

    return num_fixed;
  }

  bool is_collapse_ok(edge_descriptor v0v1)
  {
    edge_descriptor v1v0 = v0v1->opposite();
    vertex_descriptor v0 = boost::target(v1v0, polyhedron);
    vertex_descriptor v1 = boost::source(v1v0, polyhedron);

    vertex_descriptor vv, vl, vr;
    edge_descriptor  h1, h2;

    // the edges v1-vl and vl-v0 must not be both boundary edges
    if (!(v0v1->is_border()))
    {
      vl = boost::target(v0v1->next(), polyhedron);
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
      vr = boost::target(v1v0->next(), polyhedron);
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
    for (boost::tie(eb, ee) = boost::in_edges(v0, polyhedron); eb != ee; ++eb)
    {
      vv = boost::source(*eb, polyhedron);
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
    for (boost::tie(eb, ee) = boost::in_edges(vj, polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vv = boost::source(*eb, polyhedron);
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
    for (boost::tie(eb, ee) = boost::in_edges(aV, polyhedron); eb != ee; ++eb)
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

//    detect_degeneracies_in_disk();

    double area = get_surface_area();
    MCFSKEL_INFO(std::cout << "area " << area << "\n";)
  }

  void run_to_converge()
  {
    double last_area = 0;
    int num_iteration = 0;
    while (true)
    {
      MCFSKEL_INFO(std::cout << "iteration " << num_iteration + 1 << "\n";)

      contract_geometry();
      update_topology();
      detect_degeneracies();
//      detect_degeneracies_in_disk();

      double area = get_surface_area();
      double area_ratio = fabs(last_area - area) / original_area;
      MCFSKEL_INFO(std::cout << "area " << area << "\n";)
      MCFSKEL_INFO(std::cout << "|area - last_area| / original_area " << area_ratio << "\n";)
      if (area_ratio < area_TH)
      {
        break;
      }
      last_area = area;

      num_iteration++;
      if (num_iteration >= iteration_TH)
      {
        break;
      }
    }
  }

  void convert_to_skeleton()
  {
    Skeleton skeleton(polyhedron);
    std::vector<std::vector<int> > record;
    skeleton.extract_skeleton(g, points, record);

    skeleton_to_surface.resize(record.size());
    for (size_t i = 0; i < record.size(); ++i)
    {
      for (size_t j = 0; j < record[i].size(); ++j)
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
    for (size_t i = 0; i < skeleton_to_surface.size(); ++i)
    {
      cnt += skeleton_to_surface[i].size();
    }
    MCFSKEL_INFO(std::cout << "tracked " << cnt << " vertices\n";)

//    collapse_vertices_without_correspondence();
  }

  void collapse_vertices_without_correspondence()
  {
    int vertex_removed = 0;
    while (true)
    {
      int cnt = 0;
      for (size_t i = 0; i < skeleton_to_surface.size(); ++i)
      {
        if (skeleton_to_surface[i].size() == 0)
        {
          edge_desc ed;
          int u;
          bool found = false;
          in_edge_iter eb, ee;
          // look for a neighbor having non-zero correspondent vertices
          for (boost::tie(eb, ee) = boost::in_edges(i, g); eb != ee; ++eb)
          {
            ed = *eb;
            u = boost::source(ed, g);
            if (skeleton_to_surface[u].size() != 0)
            {
              found = true;
              break;
            }
          }
          if (found)
          {
            // add new edges
            for (boost::tie(eb, ee) = boost::in_edges(i, g); eb != ee; ++eb)
            {
              edge_desc edge = *eb;
              int v = boost::source(edge, g);
              if (v == u)
              {
                continue;
              }

              bool exist;
              boost::tie(edge, exist) = boost::edge(u, v, g);
              // avoid adding parallel edges
              if (!exist)
              {
                boost::add_edge(u, v, g);
              }
            }
            // remove incident edges
            for (boost::tie(eb, ee) = boost::in_edges(i, g); eb != ee; ++eb)
            {
              edge_desc edge = *eb;
              boost::remove_edge(edge, g);
            }
            cnt++;
          }
        }
      }
      vertex_removed += cnt;
      if (cnt == 0)
      {
        break;
      }
    }

    MCFSKEL_INFO(std::cout << "removed " << vertex_removed << " vertices\n";)

    int new_size = skeleton_to_surface.size() - vertex_removed;
    Graph new_g(new_size);

    std::vector<int> new_id;
    new_id.resize(skeleton_to_surface.size(), -1);

    int id = 0;
    for (size_t i = 0; i < skeleton_to_surface.size(); ++i)
    {
      if (skeleton_to_surface[i].size() != 0)
      {
        new_id[i] = id++;
      }
    }
    edge_iter eb, ee;
    for (boost::tie(eb, ee) = boost::edges(g); eb != ee; ++eb)
    {
      edge_desc ed = *eb;
      int s = boost::source(ed, g);
      int t = boost::target(ed, g);
      int ns = new_id[s];
      int nt = new_id[t];
      if (ns == -1 || nt == -1)
      {
        MCFSKEL_DEBUG(std::cerr << "wrong id\n";)
      }
      boost::add_edge(ns, nt, new_g);
    }

    std::vector<Point> new_points;
    new_points.resize(new_size);
    for (size_t i = 0; i < points.size(); ++i)
    {
      int index = new_id[i];
      if (index != -1)
      {
        new_points[index] = points[i];
      }
    }

    std::vector<std::vector<int> > new_skeleton_to_surface;
    new_skeleton_to_surface.resize(new_size);
    for (size_t i = 0; i < skeleton_to_surface.size(); ++i)
    {
      int index = new_id[i];
      if (index != -1)
      {
        new_skeleton_to_surface[index] = skeleton_to_surface[i];
      }
    }

    g = new_g;
    points = new_points;
    skeleton_to_surface = new_skeleton_to_surface;

    MCFSKEL_INFO(std::cout << "new vertices " << boost::num_vertices(g) << "\n";)
    MCFSKEL_INFO(std::cout << "new edges " << boost::num_edges(g) << "\n";)
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

  void compute_voronoi_pole()
  {
    MCFSKEL_DEBUG(std::cout << "start compute_voronoi_pole\n";)
    compute_vertex_normal();

    std::vector<std::pair<Exact_point, unsigned> > points;
    std::vector<std::vector<int> > point_to_pole;

    points.clear();
    cell_dual.clear();
    point_to_pole.clear();
    point_to_pole.resize(boost::num_vertices(polyhedron));

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      Exact_point tp((v->point()).x(), (v->point()).y(), (v->point()).z());
      points.push_back(std::make_pair(tp, vid));
    }

    Delaunay T(points.begin(), points.end());

    Finite_cells_iterator cit;
    int cell_id = 0;
    for (cit = T.finite_cells_begin(); cit != T.finite_cells_end(); ++cit)
    {
      Cell_handle cell = cit;
      Exact_point point = T.dual(cell);
      Point pt(to_double(point.x()), to_double(point.y()), to_double(point.z()));
      cell_dual.push_back(pt);
      for (int i = 0; i < 4; ++i)
      {
        TriVertex_handle vt = cell->vertex(i);
        int id = vt->info();
        point_to_pole[id].push_back(cell_id);
      }
      cell_id++;
    }

    poles.clear();
    for (size_t i = 0; i < point_to_pole.size(); ++i)
    {
      Point surface_point = Point(to_double(points[i].first.x()),
                                  to_double(points[i].first.y()),
                                  to_double(points[i].first.z()));

      double max_neg_t = 1;
      int max_neg_i = 0;

      for (size_t j = 0; j < point_to_pole[i].size(); ++j)
      {
        int pole_id = point_to_pole[i][j];
        Point cell_point = cell_dual[pole_id];
        Vector vt = cell_point - surface_point;
        Vector n = normals[i];

        double t = vt * n;
        if (t < 0 && t < max_neg_t)
        {
          max_neg_i = pole_id;
          max_neg_t = t;
        }
      }

      poles[i] = max_neg_i;
    }
  }

  void compute_vertex_normal()
  {
    normals.resize(boost::num_vertices(polyhedron));

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      normals[vid] = internal::get_vertex_normal<typename Polyhedron::Vertex,Kernel>(*v);
    }
  }

  void get_poles(std::vector<Point>& max_poles)
  {
    max_poles.resize(boost::num_vertices(polyhedron));
    vertex_iterator vb, ve;
    int cnt = 0;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      max_poles[cnt++] = cell_dual[poles[vid]];
    }
  }
};

} //namespace CGAL

#endif // MEAN_CURVATURE_SKELETON_H
