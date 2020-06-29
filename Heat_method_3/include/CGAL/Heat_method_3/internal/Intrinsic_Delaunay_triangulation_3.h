// Copyright (c) 2018  Carnegie Mellon University (USA), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Christina Vaz, Keenan Crane, Andreas Fabri


#ifndef CGAL_INTRINSIC_DELAUNAY_TRIANGULATION_3_H
#define CGAL_INTRINSIC_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Heat_method_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/property_map.h>
#include <CGAL/double.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Heat_method_3/internal/V2V.h>


#include <boost/iterator/transform_iterator.hpp>
#include <boost/unordered_map.hpp>

#include <set>
#include <stack>
#include <cmath>

#ifndef DOXYGEN_RUNNING

namespace CGAL {
namespace Heat_method_3 {



// forward declaration
template <typename IDT>
struct IDT_vertex_point_property_map;

// forward declaration
template <typename IDT, typename PM>
struct IDT_vertex_distance_property_map;

template <class TriangleMesh>
struct Intrinsic_Delaunay_triangulation_3_vertex_descriptor {
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  halfedge_descriptor hd;

  bool operator<(const Intrinsic_Delaunay_triangulation_3_vertex_descriptor& other) const
  {
    return hd < other.hd;
  }

  Intrinsic_Delaunay_triangulation_3_vertex_descriptor(const halfedge_descriptor& hd)
    : hd(hd)
  {}

  explicit Intrinsic_Delaunay_triangulation_3_vertex_descriptor(const vertex_descriptor vd, const TriangleMesh& tm)
    : hd(halfedge(vd,tm))
  {}
};

template <class TriangleMesh>
struct Intrinsic_Delaunay_triangulation_3_vertex_iterator_functor
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef vertex_descriptor argument_type;
  typedef Intrinsic_Delaunay_triangulation_3_vertex_descriptor<TriangleMesh> result_type;
  const TriangleMesh& tm;

  Intrinsic_Delaunay_triangulation_3_vertex_iterator_functor(const TriangleMesh& tm)
    :tm(tm)
  {}

  result_type
  operator()(vertex_descriptor vd) const
  {
    return result_type(halfedge(vd, tm));
  }
};

/**
 * \ingroup PkgHeatMethod
 *
 * Class `Intrinsic_Delaunay_triangulation_3` is a remeshing algorithm to improve the approximation of the `Surface_mesh_geodesic_distances_3`.
 * It internally makes a copy of the triangle mesh, performs edge flips, and computes 2D vertex coordinates per face
 * which are stored in the halfedge with the vertex as target.
 *
 * The BGL API of this class .....
 *
 *
 * \tparam TriangleMesh a triangulated surface mesh, model of `FaceListGraph` and `HalfedgeListGraph`
 * \tparam Traits a model of `HeatMethodTraits_3`
 *
 * \cgalModels `FaceListGraph`
 */

template <typename TriangleMesh,
          typename Traits = typename Kernel_traits<
                              typename boost::property_traits<
                                typename boost::property_map<TriangleMesh, vertex_point_t>::const_type
                                >::value_type
                            >::Kernel >
class Intrinsic_Delaunay_triangulation_3
{
  typedef Intrinsic_Delaunay_triangulation_3<TriangleMesh, Traits> Self;

  typedef boost::graph_traits<TriangleMesh>                        graph_traits;
  typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
  typedef typename graph_traits::edge_descriptor                edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                face_descriptor;
  typedef typename std::set<vertex_descriptor>::iterator        vertex_iterator;
  typedef typename std::set<edge_descriptor>::iterator            edge_iterator;
  /// Geometric typedefs
  typedef typename Traits::Point_3                                      Point_3;
  typedef typename Traits::FT                                                FT;
  typedef typename Traits::Vector_3                                    Vector_3;

  typedef std::pair<double,double>                                      Point_2;

  typedef int Index;

  typedef CGAL::dynamic_halfedge_property_t<Point_2> Halfedge_coordinate_tag;
  typedef typename boost::property_map<TriangleMesh, Halfedge_coordinate_tag >::type HalfedgeCoordinateMap;

  typedef typename boost::graph_traits<TriangleMesh>::vertices_size_type vertices_size_type;
  typedef typename boost::graph_traits<TriangleMesh>::edges_size_type edges_size_type;
  typedef typename boost::graph_traits<TriangleMesh>::faces_size_type faces_size_type;

  typedef CGAL::dynamic_edge_property_t<Index> Edge_property_tag;
  typedef typename boost::property_map<TriangleMesh, Edge_property_tag >::type Edge_id_map;
  typedef typename std::stack<edge_descriptor, std::list<edge_descriptor> > edge_stack;


private:
  friend struct IDT_vertex_point_property_map<Self>;
  template <class IDT, class VDM> friend struct IDT_vertex_distance_property_map;

public: // for the BGL functions below. They should maybe become friend?
  typedef CGAL::Heat_method_3::IDT_vertex_point_property_map<Self> Vertex_point_map;


  typedef Intrinsic_Delaunay_triangulation_3_vertex_descriptor<TriangleMesh> Vertex_descriptor;
  typedef Intrinsic_Delaunay_triangulation_3_vertex_iterator_functor<TriangleMesh> Vertex_iterator_functor;

public:
  /// Constructor
  /// \param input_tm the triangle mesh
  Intrinsic_Delaunay_triangulation_3(const TriangleMesh& input_tm)
    : m_intrinsic_tm(), m_input_tm(input_tm), m_vpm(*this), hcm(get(Halfedge_coordinate_tag(), m_intrinsic_tm))
  {
    build();
  }

  template <class VertexPointMap>
  Intrinsic_Delaunay_triangulation_3(const TriangleMesh& input_tm, VertexPointMap vpm)
    : m_intrinsic_tm(), m_input_tm(input_tm), m_vpm(*this), hcm(get(Halfedge_coordinate_tag(), m_intrinsic_tm))
  {
    build(vpm);
  }

  typedef TriangleMesh Triangle_mesh;

  const Triangle_mesh&
  triangle_mesh() const
  {
    return m_intrinsic_tm;
  }


  Triangle_mesh&
  triangle_mesh()
  {
    return m_intrinsic_tm;
  }

  const HalfedgeCoordinateMap&
  hcmap() const
  {
    return hcm;
  }

  template <class VertexDistanceMap>
  IDT_vertex_distance_property_map<Self,VertexDistanceMap>
  vertex_distance_map(VertexDistanceMap vdm) const
  {
    return IDT_vertex_distance_property_map<Self,VertexDistanceMap>(*this, vdm);
  }

  Vertex_point_map
  vertex_point_map() const
  {
    return m_vpm;
  }

private:

  double
  get_cotan_weight(edge_descriptor ed)
  {
    double cotan_weight = 0;
    halfedge_descriptor hd = halfedge(ed, m_intrinsic_tm);
    halfedge_descriptor hd2 = next(hd,m_intrinsic_tm);
    halfedge_descriptor hd3 = next(hd2,m_intrinsic_tm);
    Index a_i = get(edge_id_map, ed);
    Index b_i = get(edge_id_map, edge(hd2,m_intrinsic_tm));
    Index c_i = get(edge_id_map, edge(hd3,m_intrinsic_tm));
    double a = edge_lengths[a_i] + 0.0;
    double b = edge_lengths[b_i] + 0.0;
    double c = edge_lengths[c_i] + 0.0;

    double tan2 = CGAL::sqrt(CGAL::abs(((a-b+c)*(a+b-c))/((a+b+c)*(-a+b+c))));
    cotan_weight+=(1-(tan2*tan2))/(2*tan2);

    hd = opposite(hd,m_intrinsic_tm);
    hd2 =next(hd,m_intrinsic_tm);
    hd3 = next(hd2,m_intrinsic_tm);
    b_i = get(edge_id_map, edge(hd2,m_intrinsic_tm));
    c_i = get(edge_id_map, edge(hd3,m_intrinsic_tm));
    b = edge_lengths[b_i] + 0.0;
    c = edge_lengths[c_i] + 0.0;
    tan2 = CGAL::sqrt(CGAL::abs(((a-b+c)*(a+b-c))/((a+b+c)*(-a+b+c))));
    cotan_weight+=(1-(tan2*tan2))/(2*tan2);
    return cotan_weight;
  }


  //returns true if edge is locally Delaunay (opposing angles are less than pi):
  //Two ways of doing this: taking angles directly (not good with virtual edges)
  //OR: taking edge length and using law of cosines,
  //The second way checks cotan weights
  bool
  is_edge_locally_delaunay(edge_descriptor ed)
  {
    return (get_cotan_weight(ed)>=0);
  }


  void
  change_edge_length(Index i, edge_descriptor ed)
  {
    halfedge_descriptor hd = halfedge(ed,m_intrinsic_tm);
    halfedge_descriptor hd2 = next(hd,m_intrinsic_tm);
    halfedge_descriptor hd3 = next(hd2,m_intrinsic_tm);
    Index b_i = get(edge_id_map, edge(hd2,m_intrinsic_tm));
    Index c_i = get(edge_id_map, edge(hd3,m_intrinsic_tm));
    double a = edge_lengths[i];
    double b1 = edge_lengths[b_i];
    double c1 = edge_lengths[c_i];
    double tan2a = CGAL::sqrt(CGAL::abs(((c1-a+b1)*(-b1+a+c1))/((a+b1+c1)*(b1+a-c1))));
    hd = opposite(hd,m_intrinsic_tm);
    hd2 =next(hd,m_intrinsic_tm);
    hd3 = next(hd2,m_intrinsic_tm);
    b_i = get(edge_id_map, edge(hd2,m_intrinsic_tm));
    c_i = get(edge_id_map, edge(hd3,m_intrinsic_tm));
    double b2 = edge_lengths[b_i];
    double c2 = edge_lengths[c_i];
    double tan2d = CGAL::sqrt(CGAL::abs(((-a+b2+c2)*(a+b2-c2))/((a+b2+c2)*(a-b2+c2))));
    double tan2ad = (tan2a + tan2d)/(1-tan2a*tan2d);
    double cosad = (1-tan2ad*tan2ad)/(1+tan2ad*tan2ad);
    double new_length = CGAL::sqrt( CGAL::abs(b1*b1 + c2*c2 - 2*b1*c2*cosad));
    edge_lengths[i] = new_length;
  }


  //Heron's formula
  double
  face_area(double a, double b, double c)
  {
    double S = (a+b+c)/2;
    return CGAL::sqrt(S*(S-a)*(S-b)*(S-c));
  }


  void
  loop_over_edges(edge_stack stack, std::vector<int>& marked_edges)
  {
    int a = 0;
    while(!stack.empty()) {
      edge_descriptor ed = stack.top();
      stack.pop();

      Index edge_i = get(edge_id_map,ed);

      marked_edges[edge_i]=0;
      //if the edge itself is not locally delaunay, go back
      if(!(is_edge_locally_delaunay(ed))) {
        if(!(is_border(ed,m_intrinsic_tm))) {
          a++;
          change_edge_length(edge_i,ed);
          halfedge_descriptor hd = (halfedge(ed, m_intrinsic_tm));
          CGAL::Euler::flip_edge(hd, m_intrinsic_tm);
          edge_descriptor next_edge= edge(next(hd,m_intrinsic_tm),m_intrinsic_tm);
          Index next_edge_i =  get(edge_id_map, next_edge);

          //if edge was already checked, go back and check again
          //for the 4 surrounding edges, since local 'geometry' changed,
          if(!(marked_edges[next_edge_i])) {
            stack.push(next_edge);
            marked_edges[next_edge_i] = 1;
          }
          next_edge = edge(prev(hd,m_intrinsic_tm),m_intrinsic_tm);
          next_edge_i = get(edge_id_map,next_edge);
          if(!(marked_edges[next_edge_i])) {
            stack.push(next_edge);
            marked_edges[next_edge_i] = 1;
          }
          next_edge = edge(next(opposite(hd,m_intrinsic_tm),m_intrinsic_tm),m_intrinsic_tm);
          next_edge_i = get(edge_id_map,next_edge);
          if(!(marked_edges[next_edge_i])) {
            stack.push(next_edge);
            marked_edges[next_edge_i] = 1;
          }
          next_edge = edge(prev(opposite(hd,m_intrinsic_tm),m_intrinsic_tm),m_intrinsic_tm);
          next_edge_i = get(edge_id_map,next_edge);
          if(!(marked_edges[next_edge_i])) {
            stack.push(next_edge);
            marked_edges[next_edge_i] = 1;
          }
        }
        //then go back to top of the stack
      }
    }
  }


  template <class VertexPointMap>
  void
  build(VertexPointMap vpm)
  {
    CGAL_precondition(is_triangle_mesh(m_intrinsic_tm));

    typename Traits::Compute_squared_distance_3 squared_distance = Traits().compute_squared_distance_3_object();

    std::vector<std::pair<vertex_descriptor,
                          vertex_descriptor> > pairs;
    copy_face_graph(m_input_tm, m_intrinsic_tm,
                    parameters::vertex_to_vertex_output_iterator(std::back_inserter(pairs)).
                    vertex_point_map(vpm));

    for(std::size_t i=0; i < pairs.size(); i++) {
      v2v[pairs[i].second] = pairs[i].first;
      vtov[pairs[i].first] = pairs[i].second;
    }

    edge_stack stack;
    std::size_t number_of_edges = num_edges(m_intrinsic_tm);
    edge_lengths.resize(number_of_edges);
    mark_edges.resize(number_of_edges, 1);
    edge_id_map = get(Edge_property_tag(), m_intrinsic_tm);
    Index edge_i = 0;
    VertexPointMap vpm_intrinsic_tm = get(boost::vertex_point,m_intrinsic_tm);

    for(edge_descriptor ed : edges(m_intrinsic_tm)) {
      edge_lengths[edge_i] = CGAL::sqrt(to_double(squared_distance(get(vpm_intrinsic_tm, source(ed,m_intrinsic_tm)),
                                                                   get(vpm_intrinsic_tm, target(ed,m_intrinsic_tm)))));
        //  Polygon_mesh_processing::edge_length(halfedge(ed,m_intrinsic_tm),m_intrinsic_tm);
      put(edge_id_map, ed, edge_i++);
      stack.push(ed);
    }
    loop_over_edges(stack, mark_edges);
    //now that edges are calculated, go through and for each face, calculate the vertex positions around it

    for(face_descriptor f : faces(m_intrinsic_tm)) {
      CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;

      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,m_intrinsic_tm),m_intrinsic_tm);
      halfedge_descriptor hd = halfedge(f,m_intrinsic_tm);
      if(face(hd,m_intrinsic_tm) != f) {
        hd = opposite(hd,m_intrinsic_tm);
      }
      hd = next(hd,m_intrinsic_tm);
      //each 'local' set of coordinates will have 0,0 at the first vertex/halfedge
      Point_2 p11(0,0);
      put(hcm, prev(hd,m_intrinsic_tm),p11);
      edge_descriptor ed1 = edge(hd, m_intrinsic_tm);
      hd = next(hd,m_intrinsic_tm);
      //the second local coordinate will be edge_length(first edge),0
      Point_2 p21(edge_lengths[get(edge_id_map,ed1)], 0);
      put(hcm,prev(hd,m_intrinsic_tm),p21);

      //use basic trigonometry to compute third coordinate
      edge_descriptor ed2 = edge(hd, m_intrinsic_tm);
      hd = next(hd,m_intrinsic_tm);
      edge_descriptor ed3 = edge(hd, m_intrinsic_tm);
      Index e1 = get(edge_id_map, ed1);
      Index e2 = get(edge_id_map, ed2);
      Index e3 = get(edge_id_map, ed3);
      double e1_len = edge_lengths[e1];
      double e2_len = edge_lengths[e2];
      double e3_len = edge_lengths[e3];
      double angle_a = -(e2_len*e2_len) + e3_len*e3_len + e1_len*e1_len;
      angle_a = acos(angle_a/(2*e3_len*e1_len));
      Point_2 p31(e3_len*std::cos(angle_a), e3_len*std::sin(angle_a));
      put(hcm,prev(hd,m_intrinsic_tm),p31);

    }
  }

  void
  build()
  {
    build( get(boost::vertex_point, m_input_tm) );
  }


  //todo:: determine which can be const
  TriangleMesh m_intrinsic_tm; // this is the copy where edges get flipped
  const TriangleMesh& m_input_tm; // this is the reference to the original
  Vertex_point_map m_vpm;
  HalfedgeCoordinateMap hcm;
  Edge_id_map edge_id_map;

  std::vector<double> edge_lengths;
  std::vector<int> mark_edges;
public:

  boost::unordered_map<vertex_descriptor,vertex_descriptor> v2v, vtov;
};

} // namespace Heat_method_3

namespace Heat_method_3 {

template <typename TM,
          typename T>
struct V2V<CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T> >
{
  typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T> Idt;
  const Idt& idt;
  /**
   * Vertex descriptor for iDT
   */
  typedef typename boost::graph_traits<Idt>::vertex_descriptor Idt_vertex_descriptor;
  /**
   * Vertex descriptor for TriangleMesh
   */
  typedef typename boost::graph_traits<TM>::vertex_descriptor TM_vertex_descriptor;


  /**
   * Default constructor
   */
  V2V(const Idt& idt)
    : idt(idt)
  {}

    /**
     * Create iDT vertex descriptor map for vertex bijection between original mesh and iDT mesh for Heat method
     */
  Idt_vertex_descriptor operator()(const TM_vertex_descriptor& vd) const
  {
    return Idt_vertex_descriptor(idt.vtov.at(vd),idt.triangle_mesh());
  }
};

} // namespace Heat_method_3
} // namespace CGAL


namespace boost {

template <typename TM,
          typename T>
struct graph_traits<CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T> > {

  typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3_vertex_descriptor<TM> vertex_descriptor;
  typedef boost::transform_iterator<
    CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3_vertex_iterator_functor<TM>,
    typename boost::graph_traits<TM>::vertex_iterator> vertex_iterator;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::halfedge_iterator halfedge_iterator;
  typedef typename boost::graph_traits<TM>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TM>::edge_iterator edge_iterator;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TM>::face_iterator face_iterator;

  typedef typename boost::graph_traits<TM>::vertices_size_type vertices_size_type;

  static face_descriptor null_face() { return boost::graph_traits<TM>::null_face(); }
};

} // namespace boost


namespace CGAL {
namespace Heat_method_3 {

template <typename TM,
          typename T>
typename boost::graph_traits<TM>::vertices_size_type
num_vertices(const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return num_vertices(idt.triangle_mesh());
}


template <typename TM,
          typename T>
typename boost::graph_traits<TM>::edges_size_type
num_edges(const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return num_edges(idt.triangle_mesh());
}


template <typename TM,
          typename T>
typename boost::graph_traits<TM>::faces_size_type
num_faces(const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return num_faces(idt.triangle_mesh());
}


template <typename TM,
          typename T>
Iterator_range<typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::vertex_iterator>
vertices(const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
 {
   std::pair<typename boost::graph_traits<TM>::vertex_iterator,
             typename boost::graph_traits<TM>::vertex_iterator> p = vertices(idt.triangle_mesh());

  typedef typename Intrinsic_Delaunay_triangulation_3<TM,T>::Vertex_iterator_functor Fct;
  Fct fct(idt.triangle_mesh());
  return make_range(boost::make_transform_iterator(p.first, fct),
                    boost::make_transform_iterator(p.second,fct));
 }


template <typename TM,
          typename T>
Iterator_range<typename boost::graph_traits<TM>::halfedge_iterator>
halfedges(const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return make_range( halfedges(idt.triangle_mesh()) );
}


template <typename TM,
          typename T>
Iterator_range<typename boost::graph_traits<TM>::edge_iterator>
edges(const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return make_range( edges(idt.triangle_mesh()) );
}


template <typename TM,
          typename T>
Iterator_range<typename boost::graph_traits<TM>::face_iterator>
faces(const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
 {
   return make_range( faces(idt.triangle_mesh()) );
 }

template <typename TM,
          typename T>
typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::vertex_descriptor
vertex(typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor hd,
       const Intrinsic_Delaunay_triangulation_3<TM,T>& )
{
  return typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::vertex_descriptor(hd);
}


template <typename TM,
          typename T>
typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::face_descriptor fd,
         const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return halfedge(fd, idt.triangle_mesh());
}


template <typename TM,
          typename T>
typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::edge_descriptor ed,
       const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return halfedge(ed, idt.triangle_mesh());
}


template <typename TM,
          typename T>
typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor
opposite(typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor hd,
         const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return opposite(hd, idt.triangle_mesh());
}


template <typename TM,
          typename T>
typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor
next(typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor hd,
     const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return next(hd, idt.triangle_mesh());
}


template <typename TM,
          typename T>
typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::face_descriptor
face(typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor hd,
     const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return face(hd, idt.triangle_mesh());
}


template <typename TM,
          typename T>
typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::vertex_descriptor
source(typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor hd,
       const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  typedef typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::vertex_descriptor vertex_descriptor;

  return vertex_descriptor(opposite(hd, idt.triangle_mesh()));
}


template <typename TM,
          typename T>
typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::vertex_descriptor
target(typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::halfedge_descriptor hd,
       const Intrinsic_Delaunay_triangulation_3<TM,T>&)
{
  typedef typename boost::graph_traits<Intrinsic_Delaunay_triangulation_3<TM,T> >::vertex_descriptor vertex_descriptor;

  return vertex_descriptor(hd);
}

template <typename IDT>
struct IDT_vertex_point_property_map {
  const IDT& idt;
  typedef typename IDT::Triangle_mesh TM;
  typedef typename boost::graph_traits<IDT>::vertex_descriptor key_type;
  typedef typename IDT::Point_3 value_type;
  typedef typename IDT::Point_2 Point_2;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  /**
   * Default constructor for vertex/point property map
   */
  IDT_vertex_point_property_map(const IDT& idt)
    : idt(idt)
    {}


  /**
   * friend function for Heat method to get vertex descriptor's coordinates in iDT's local coordinate system
   */
  friend value_type get(const IDT_vertex_point_property_map<IDT>& pm,
                        key_type vd)
  {
    const Point_2& p = get(pm.idt.hcmap(), vd.hd);
    return value_type(p.first, p.second, 0);
  }
};

template <typename IDT, typename PM>
struct IDT_vertex_distance_property_map {
  const IDT& idt;
  PM pm;

  typedef typename IDT::Triangle_mesh TM;
  typedef typename IDT::Vertex_descriptor key_type;
  typedef double value_type;
  typedef value_type reference;

  IDT_vertex_distance_property_map(const IDT& idt,
                                   PM pm)
    : idt(idt), pm(pm)
    {}


  // no need for a get()

  friend void put(IDT_vertex_distance_property_map<IDT,PM> idtpm,
                  key_type vd,
                  value_type v)
  {
    typename boost::graph_traits<TM>::vertex_descriptor tm_vd = target(vd.hd, idtpm.idt.triangle_mesh());

    put(idtpm.pm, idtpm.idt.v2v.at(tm_vd), v);
  }
};

} // namespace Heat_method_3
} // namespace CGAL


namespace boost {

template <typename TM,
          typename T>
struct property_map<CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T>,
                    CGAL::vertex_point_t > {
  typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T> IDT;
  typedef CGAL::Heat_method_3::IDT_vertex_point_property_map<IDT> type;
  typedef type const_type;
};

template <typename TM,
          typename T>
struct property_map<CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T>,
                    CGAL::face_index_t > {
  typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T> IDT;
  typedef typename property_map<TM, CGAL::face_index_t>::type type;
  typedef typename property_map<TM, CGAL::face_index_t>::const_type const_type;
};

} // boost


namespace CGAL {
namespace Heat_method_3 {

template <typename TM,
          typename T>
typename boost::property_map<TM,CGAL::face_index_t>::type
get(CGAL::face_index_t fi,
    const CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return get(fi, idt.triangle_mesh());
}

template <typename TM,
          typename T>
CGAL::Heat_method_3::IDT_vertex_point_property_map<CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T> >
get(CGAL::vertex_point_t,
    const CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  return CGAL::Heat_method_3::IDT_vertex_point_property_map<CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T> >(idt);
}

template <typename IDT, typename PM, typename K, typename V>
class IDT_dynamic_vertex_property_map {
  const IDT& idt;
  PM pm;

public:
  typedef IDT_dynamic_vertex_property_map<IDT,PM,K,V> Self;
  typedef typename IDT::Triangle_mesh TM;
  typedef typename boost::graph_traits<TM>::vertex_descriptor TM_vertex_descriptor;


  IDT_dynamic_vertex_property_map(const IDT& idt, PM pm)
    : idt(idt), pm(pm)
  {}


  friend V get(const Self& idpm, const K& k)
  {
    return get(idpm.pm, target(k.hd, idpm.idt.triangle_mesh()));
  }


  friend void put(const Self& idpm, const K& k, const V& v)
  {
    put(idpm.pm, target(k.hd, idpm.idt.triangle_mesh()), v);
  }


  friend V get(const Self& idpm, const TM_vertex_descriptor& k)
  {
    return get(idpm.pm, k);
  }


  friend void put(const Self& idpm, const TM_vertex_descriptor& k, const V& v)
  {
    put(idpm.pm, k, v);
  }

};

} } // CGAL::Heat_method_3

namespace boost {

template <typename TM,
          typename T,
          typename dT>
struct property_map<CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T>,
                    CGAL::dynamic_vertex_property_t<dT> > {
  typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TM,T> IDT;
  typedef CGAL::Heat_method_3::IDT_dynamic_vertex_property_map<IDT,
                                                               typename property_map<TM, CGAL::dynamic_vertex_property_t<dT> >::type,
                                                               typename graph_traits<IDT>::vertex_descriptor,
                                                               dT> type;
  typedef CGAL::Heat_method_3::IDT_dynamic_vertex_property_map<IDT,
                                                               typename property_map<TM, CGAL::dynamic_vertex_property_t<dT> >::const_type,
                                                               typename graph_traits<IDT>::vertex_descriptor,
                                                               dT> const_type;
};

} // namespace boost


namespace CGAL {
namespace Heat_method_3 {

template <typename TM,
          typename T,
          typename dT>
typename boost::property_map<Intrinsic_Delaunay_triangulation_3<TM,T>, CGAL::dynamic_vertex_property_t<dT> >::const_type
get(CGAL::dynamic_vertex_property_t<dT> dvp,
    const Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  typedef Intrinsic_Delaunay_triangulation_3<TM,T> IDT;
  typedef IDT_dynamic_vertex_property_map<IDT,
                                          typename boost::property_map<TM, CGAL::dynamic_vertex_property_t<dT> >::const_type,
                                          typename boost::graph_traits<IDT>::vertex_descriptor,
                                          dT> PM;
  return PM(idt,get(dvp,idt.triangle_mesh()));
}


template <typename TM,
          typename T,
          typename dT>
typename boost::property_map<Intrinsic_Delaunay_triangulation_3<TM,T>, CGAL::dynamic_vertex_property_t<dT> >::type
get(CGAL::dynamic_vertex_property_t<dT> dvp,
    Intrinsic_Delaunay_triangulation_3<TM,T>& idt)
{
  typedef Intrinsic_Delaunay_triangulation_3<TM,T> IDT;
  typedef IDT_dynamic_vertex_property_map<IDT,
                                          typename boost::property_map<TM, CGAL::dynamic_vertex_property_t<dT> >::type,
                                          typename boost::graph_traits<IDT>::vertex_descriptor,
                                          dT> PM;

  return PM(idt, get(dvp,idt.triangle_mesh()));
}

} // namespace Heat_method_3
} // namespace CGAL

#endif // DOXYGEN_RUNNING

#include <CGAL/enable_warnings.h>
#endif // CGAL_INTRINSIC_DELAUNAY_TRIANGULATION_3_H
