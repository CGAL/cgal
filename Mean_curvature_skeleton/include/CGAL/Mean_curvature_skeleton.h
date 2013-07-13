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

template <class Polyhedron, class Graph,
          class PolyhedronVertexIndexMap, class PolyhedronEdgeIndexMap>
class CurveSkeleton
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
  typedef typename Polyhedron::Facet_iterator                                  Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator                Halfedge_facet_circulator;

// Data members
private:
  std::vector<std::vector<int> > edge_to_face;
  std::vector<std::vector<int> > edge_to_vertex;
  std::vector<std::vector<int> > vertex_to_edge;
  std::vector<std::vector<int> > face_to_edge;
  std::vector<bool> is_vertex_deleted;
  std::vector<bool> is_edge_deleted;
  std::vector<bool> is_face_deleted;

  PolyhedronVertexIndexMap vertex_id_pmap;
  PolyhedronEdgeIndexMap edge_id_pmap;

  Polyhedron *polyhedron;

  struct Edge
  {
    int id;
    double length;

    bool operator <(const Edge& rhs) const
    {
      return (length == rhs.length) ? (id < rhs.id) : length < rhs.length;
    }
  };

// Public methods
public:
  CurveSkeleton() : polyhedron(NULL)
  {
  }

  CurveSkeleton(Polyhedron* polyhedron) : polyhedron(polyhedron)
  {
  }

// Private methods
//private:
  void init()
  {
    int num_edges = boost::num_edges(*polyhedron) / 2;
    int num_faces = polyhedron->size_of_facets();
    int num_vertices = boost::num_vertices(*polyhedron);
    edge_to_face.clear();
    edge_to_face.resize(num_edges);
    edge_to_vertex.clear();
    edge_to_vertex.resize(num_edges);
    vertex_to_edge.clear();
    vertex_to_edge.resize(num_vertices);
    face_to_edge.clear();
    face_to_edge.resize(num_faces);

    is_vertex_deleted.clear();
    is_vertex_deleted.resize(num_vertices, false);
    is_edge_deleted.clear();
    is_edge_deleted.resize(num_edges, false);
    is_face_deleted.clear();
    is_face_deleted.resize(num_faces, false);

    // assign vertex id
    vertex_iterator vb, ve;
    int idx = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, idx++);
    }
    std::cerr << "vertex num " << idx << "\n";

    // assign edge id
    // the two halfedges representing the same edge get the same id
    edge_iterator eb, ee;
    idx = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, -1);
    }
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ed = *eb;
      int id = boost::get(edge_id_pmap, ed);
      if (id == -1)
      {
        boost::put(edge_id_pmap, ed, idx);
        edge_descriptor ed_opposite = ed->opposite();
        boost::put(edge_id_pmap, ed_opposite, idx);
        idx++;
      }
    }
    std::cerr << "edge num " << idx << "\n";

    // assign face id and compute edge-face connectivity
    int face_id = 0;
    for (Facet_iterator i = polyhedron->facets_begin(); i != polyhedron->facets_end(); ++i)
    {
      Halfedge_facet_circulator j = i->facet_begin();
      // Facets in polyhedral surfaces are at least triangles.
      CGAL_assertion(CGAL::circulator_size(j) >= 3);
      do
      {
        int id = j->id();
        face_to_edge[face_id].push_back(id);
        edge_to_face[id].push_back(face_id);
      } while (++j != i->facet_begin());
      face_id++;
    }
    std::cerr << "face num " << face_id << "\n";

    // compute vertex-edge connectivity
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor vd = *vb;
      int vid = boost::get(vertex_id_pmap, vd);
      in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = boost::in_edges(*vb, *polyhedron); e != e_end; e++)
      {
        edge_descriptor ed = *e;
        int eid = boost::get(edge_id_pmap, ed);
        vertex_to_edge[vid].push_back(eid);
        edge_to_vertex[eid].push_back(vid);
      }
    }
  }

  void collapse()
  {
    std::set<Edge> queue;
    queue.clear();

    // put all the edges into a priority queue
    // shorter edge has higher priority
    edge_iterator eb, ee;
    std::vector<bool> is_edge_inserted;
    is_edge_inserted.clear();
    is_edge_inserted.resize(edge_to_face.size(), false);
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      Edge edge;

      edge_descriptor ed = *eb;
      edge.id = boost::get(edge_id_pmap, ed);

      if (is_edge_inserted[edge.id])
      {
        continue;
      }
      vertex_descriptor v1 = ed->vertex();
      vertex_descriptor v2 = ed->opposite()->vertex();

      Point source = v1->point();
      Point target = v2->point();
      edge.length = sqrtf(squared_distance(source, target));

      queue.insert(edge);
      is_edge_inserted[edge.id] = true;
    }

    // start collapsing edges until all the edges have no incident faces
    while (!queue.empty())
    {
      Edge edge = *(queue.begin());
      queue.erase(queue.begin());

      int eid = edge.id;

      if (is_edge_deleted[eid])
      {
        std::cerr << "edge already deleted\n";
        continue;
      }
      if (edge_to_face[eid].size() == 0)
      {
        std::cerr << "edge has no incident face\n";
        continue;
      }
      else
      {
        for (size_t i = 0; i < edge_to_face[eid].size(); i++)
        {
          int fid = edge_to_face[eid][i];
          if (is_face_deleted[fid])
          {
//            std::cerr << "wtf\n";
            continue;
          }
        }
      }

      // mark the edge and incident faces as deleted
      is_edge_deleted[eid] = true;
//      std::cerr << "delete edge " << eid << "\n";
      std::vector<int> faces(edge_to_face[eid]);
      for (size_t i = 0; i < faces.size(); i++)
      {
        int fid = faces[i];
        is_face_deleted[fid] = true;
        // remove face from the incident edges
        for (size_t j = 0; j < face_to_edge[fid].size(); j++)
        {
          int e = face_to_edge[fid][j];
          for (size_t k = 0; k < edge_to_face[e].size(); k++)
          {
            if (edge_to_face[e][k] == fid)
            {
              edge_to_face[e].erase(edge_to_face[e].begin() + k);
              break;
            }
          }
        }
//        std::cerr << "delete face " << fid << "\n";
      }

      // p1 to be deleted
      int p1 = edge_to_vertex[eid][0];
      int p2 = edge_to_vertex[eid][1];
      is_vertex_deleted[p1] = true;

      // delete the edge from p1 and p2's incident edges
      for (size_t i = 0; i < vertex_to_edge[p1].size(); i++)
      {
        if (vertex_to_edge[p1][i] == eid)
        {
          vertex_to_edge[p1].erase(vertex_to_edge[p1].begin() + i);
          break;
        }
      }
      for (size_t i = 0; i < vertex_to_edge[p2].size(); i++)
      {
        if (vertex_to_edge[p2][i] == eid)
        {
          vertex_to_edge[p2].erase(vertex_to_edge[p2].begin() + i);
          break;
        }
      }

      // add the incident edges of p1 to p2
      for (size_t i = 0; i < vertex_to_edge[p1].size(); i++)
      {
        int edge = vertex_to_edge[p1][i];
        if (is_edge_deleted[edge])
        {
          continue;
        }
        vertex_to_edge[p2].push_back(edge);
        // change the incident vertex to p2
        for (size_t j = 0; j < edge_to_vertex[edge].size(); j++)
        {
          if (edge_to_vertex[edge][j] == p1)
          {
            edge_to_vertex[edge][j] = p2;
          }
        }
      }

      std::vector<int> vertex_to_edge_p2(vertex_to_edge[p2]);
      for (size_t i = 0; i < vertex_to_edge_p2.size(); i++)
      {
        // ei to be removed
        int ei = vertex_to_edge_p2[i];
        for (size_t j = i + 1; j < vertex_to_edge_p2.size(); j++)
        {
          int ej = vertex_to_edge_p2[j];
          if (is_same_edge(ei, ej) || is_edge_deleted[ei])
          {
            // remove ei from p2
            for (size_t k = 0; k < vertex_to_edge[p2].size(); k++)
            {
              int ek = vertex_to_edge[p2][k];
              // migrate faces from ei to ej
              if (ei == ek)
              {
                for (size_t f = 0; f < edge_to_face[ei].size(); f++)
                {
                  int fid = edge_to_face[ei][f];
                  if (!is_face_deleted[fid])
                  {
                    if (std::find(edge_to_face[ej].begin(),
                                  edge_to_face[ej].end(),
                                  fid)
                        == edge_to_face[ej].end())
                    {
                      edge_to_face[ej].push_back(fid);
                      for (size_t e = 0; e < face_to_edge[fid].size(); e++)
                      {
                        if (face_to_edge[fid][e] == ei)
                        {
                          face_to_edge[fid][e] = ej;
                          break;
                        }
                      }
                    }
                  }
                }
                // finally remove ei from p2
                vertex_to_edge[p2].erase(vertex_to_edge[p2].begin() + k);
                // and also remove ei from the other end point
                for (size_t t = 0; t < edge_to_vertex[ei].size(); t++)
                {
                  int vid = edge_to_vertex[ei][t];
                  if (vid != p2)
                  {
                    for (size_t q = 0; q < vertex_to_edge[vid].size(); q++)
                    {
                      if (vertex_to_edge[vid][q] == ei)
                      {
                        vertex_to_edge[vid].erase(vertex_to_edge[vid].begin() + q);
                        break;
                      }
                    }
                  }
                }
                is_edge_deleted[ei] = true;
//                std::cerr << "delete edge " << ei << "\n";
                break;
              }
            }
          }
        }
      }
    }

    // for debugging purpose
    std::cerr << "finish collapse\n";
    print_detail_stat();
//    print_stat();
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ed = *eb;
      int id = boost::get(edge_id_pmap, ed);
      if (!is_edge_deleted[id])
      {
        if (edge_to_face[id].size() > 0)
        {
          std::cerr << "edge should not have faces " << edge_to_face[id].size() << "\n";
        }
      }
    }
  }

  bool is_same_edge(int ei, int ej)
  {
    if (edge_to_vertex[ei][0] == edge_to_vertex[ej][0]
     && edge_to_vertex[ei][1] == edge_to_vertex[ej][1])
    {
      return true;
    }
    if (edge_to_vertex[ei][1] == edge_to_vertex[ej][0]
     && edge_to_vertex[ei][0] == edge_to_vertex[ej][1])
    {
      return true;
    }
    return false;
  }

  void print_stat()
  {
    int cnt = 0;
    for (size_t i = 0; i < is_vertex_deleted.size(); i++)
    {
      if (!is_vertex_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of vertices " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_edge_deleted.size(); i++)
    {
      if (!is_edge_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of edges " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_face_deleted.size(); i++)
    {
      if (!is_face_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of faces " << cnt << "\n";
  }

  void print_detail_stat()
  {
    int cnt = 0;
    for (size_t i = 0; i < is_vertex_deleted.size(); i++)
    {
      if (!is_vertex_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of vertices " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_edge_deleted.size(); i++)
    {
      if (!is_edge_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of edges " << cnt << "\n";

//    edge_iterator eb, ee;
//    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
//    {
//      edge_descriptor ed = *eb;
//      int id = boost::get(edge_id_pmap, ed);
//      if (is_edge_deleted[id])
//      {
//        continue;
//      }
//      std::cerr << "edge to vertex: ";
//      for (int i = 0; i < edge_to_vertex[id].size(); i++)
//      {
//        std::cerr << edge_to_vertex[id][i] << " ";
//      }
//      std::cerr << "\n";
//    }

    cnt = 0;
    for (size_t i = 0; i < is_face_deleted.size(); i++)
    {
      if (!is_face_deleted[i])
      {
        cnt++;
        std::cerr << "face " << i << "\n";
      }
    }
    std::cerr << "num of faces " << cnt << "\n";
  }

  // extract the skeleton to a boost::graph data structure
  void extract_skeleton(Graph& graph, std::vector<Point>& points)
  {
    std::vector<int> new_vertex_id;
    new_vertex_id.clear();
    new_vertex_id.resize(vertex_to_edge.size(), -1);

    int id = 0;
    for (size_t i = 0; i < is_vertex_deleted.size(); i++)
    {
      if (!is_vertex_deleted[i])
      {
        new_vertex_id[i] = id++;
      }
    }
//    std::cerr << "num of vertices " << id << "\n";

    Graph temp_graph(id);

    for (size_t i = 0; i < is_edge_deleted.size(); i++)
    {
      if (!is_edge_deleted[i])
      {
        int p1 = edge_to_vertex[i][0];
        int p2 = edge_to_vertex[i][1];
        int p1_id = new_vertex_id[p1];
        int p2_id = new_vertex_id[p2];
        boost::add_edge(p1_id, p2_id, temp_graph);
      }
    }

    vertex_iterator vb, ve;
    points.resize(id);
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      int new_id = new_vertex_id[id];
      Point pos = vb->point();
      points[new_id] = pos;
    }
    boost::copy_graph(temp_graph, graph);
  }
};

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
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator              edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator           in_edge_iterator;
  typedef typename internal::Cotangent_weight<Polyhedron,
  internal::Cotangent_value_minimum_zero<Polyhedron,
  internal::Cotangent_value_Meyer_secure<Polyhedron> > >                       Weight_calculator;

  // Skeleton types
  typedef typename internal::CurveSkeleton<Polyhedron, Graph,
  PolyhedronVertexIndexMap, PolyhedronEdgeIndexMap>                            Skeleton;

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

  int vertex_id_count;

  std::map<size_t, bool> is_vertex_fixed_map;
  std::vector<double> halfedge_angle;

  Graph g;
  std::vector<Point> points;

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

// Public methods
public:

  // The constructor gets the polyhedron that we will model
  Mean_curvature_skeleton(Polyhedron* P,
                          PolyhedronVertexIndexMap Vertex_index_map,
                          PolyhedronEdgeIndexMap Edge_index_map,
                          double omega_L, double omega_H,
                          double edgelength_TH, double zero_TH,
                          Weight_calculator weight_calculator = Weight_calculator()
                          )
    :polyhedron(P), vertex_id_pmap(Vertex_index_map), edge_id_pmap(Edge_index_map),
      omega_L(omega_L), omega_H(omega_H), edgelength_TH(edgelength_TH), TH_ALPHA(110),
      weight_calculator(weight_calculator), zero_TH(zero_TH)
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
      int i = boost::get(vertex_id_pmap, *vb);
      if (i >= nver)
      {
        std::cerr << "id is too large\n";
      }
      if (i < 0)
      {
        std::cerr << "id is too small\n";
      }
      if (is_vertex_fixed_map.find(i) != is_vertex_fixed_map.end()
          && is_vertex_fixed_map[i])
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
      int i = boost::get(vertex_id_pmap, vi);
      double omega;
      if (is_vertex_fixed_map.find(i) != is_vertex_fixed_map.end()
          && is_vertex_fixed_map[i])
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
//    std::cerr << "fix " << cnt_fix << " vertices in RHS\n";
//    std::cerr << "end RHS\n";
  }

  void update_vertex_id()
  {
//    std::cerr << "start update id\n";
    std::map<size_t, bool> is_vertex_fixed_map_new;
    is_vertex_fixed_map_new.clear();

    vertex_iterator vb, ve;
    int vertex_id = 0;
    int cnt = 0;
    std::vector<int> vertex_id_old;
    vertex_id_old.clear();
    vertex_id_old.resize(boost::num_vertices(*polyhedron));
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_id_old[cnt++] = boost::get(vertex_id_pmap, *vb);
      boost::put(vertex_id_pmap, *vb, vertex_id++);
    }

    edge_iterator eb, ee;
    int edge_id = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, edge_id++);
    }

    cnt = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int old_id = vertex_id_old[cnt++];
      if (is_vertex_fixed_map.find(old_id) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[old_id])
        {
          int new_id = boost::get(vertex_id_pmap, *vb);
          is_vertex_fixed_map_new[new_id] = true;
//          std::cerr << "fix " << new_id << "\n";
        }
      }
    }
    is_vertex_fixed_map = is_vertex_fixed_map_new;

    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[id])
        {
//          std::cerr << "check fix " << id << "\n";
        }
      }
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
      int i = boost::get(vertex_id_pmap, vi);
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
    }

    int edge_id = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, edge_id++);
    }

    // This is a stop predicate (defines when the algorithm terminates).
    // The simplification stops when the length of all edges is greater than the minimum threshold.
    CGAL::internal::Minimum_length_predicate<Polyhedron> stop(edgelength_TH);

    int cnt = 0;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[id])
        {
          cnt++;
        }
      }
    }
//    std::cerr << "before collapse " << cnt << " fixed vertices\n";

    int r = SMS::edge_collapse
                (*polyhedron
                ,stop
                ,CGAL::get_cost     (SMS::Edge_length_cost  <Polyhedron>())
                      .get_placement(SMS::Midpoint_placement<Polyhedron>())
                      .edge_is_border_map(constrains_map)
                );

    cnt = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        if (is_vertex_fixed_map[id])
        {
          cnt++;
        }
      }
    }
//    std::cerr << "after collapse " << cnt << " fixed vertices\n";
    return r;
  }

  int iteratively_collapse_edges()
  {
    int num_collapses = 0;
    int iter = 0;
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

      if (ed->is_border())
      {
        std::cerr << "compute_incident_angle: edge is border\n";
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
      cnt++;
    }
    return cnt;
  }

  int iteratively_split_triangles()
  {
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
    return num_splits;
  }

  void update_topology()
  {
//    std::cerr << "before collapse edges\n";
    int num_collapses = iteratively_collapse_edges();
    std::cerr << "collapse " << num_collapses << " edges.\n";

    int num_splits = iteratively_split_triangles();
    std::cerr << "split " << num_splits << " edges.\n";
  }

  // TODO: check if the local neighborhood is a disk
  int detect_degeneracies()
  {
    int num_fixed = 0;
    double elength_fixed = edgelength_TH * 0.7;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; vb++)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);
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
  }

  void convert_to_skeleton()
  {
    Skeleton skeleton(polyhedron);
    skeleton.init();
    skeleton.collapse();
    skeleton.extract_skeleton(g, points);
  }

  void get_skeleton(Graph& g, std::vector<Point>& points)
  {
    g = this->g;
    points = this->points;
  }
};

} //namespace CGAL

#endif // MEAN_CURVATURE_SKELETON_H
