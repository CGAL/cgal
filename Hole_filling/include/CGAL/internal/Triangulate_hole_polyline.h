#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/utility.h>
#include <vector>
#include <limits>
#include <CGAL/value_type_traits.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {
namespace internal {

struct Weight {

  std::pair<double,double> w;

  Weight() : w(std::make_pair<double>(0, 0)) { }

  Weight(double angle, double area) : w(angle, area) { }

  template<class Point_3>
  Weight(const std::vector<Point_3>& P, 
         const std::vector<Point_3>& Q, 
         int i, int j, int k)
  {
    CGAL_assertion(i+1 ==j);
    CGAL_assertion(j+1 ==k);
    const Point_3& p = P[i];
    const Point_3& q = P[j];
    const Point_3& r = P[k];
    // The CGAL::dihedral angle is measured between the oriented triangles, that is it goes from [-pi, pi]
    // What we need is the angle between the normals of the triangles between [0, pi]
    double ang1 = 0, ang2 = 0;

    if(!Q.empty()){
      ang1 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(p,q,r,Q[i]));
      ang2 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(q,r,p, Q[j]));
    }
    w = std::make_pair((std::max)(ang1, ang2), std::sqrt(CGAL::squared_area(p,q,r)));
  }

  template<class Point_3>
  Weight(const std::vector<Point_3>& P, 
         const std::vector<Point_3>& Q, 
         int i, int m, int k, 
         const std::vector<int>& lambda)
  {
    CGAL_assertion(i < m);
    CGAL_assertion(m < k);
    int n = P.size() -1; // because the first and last point are equal
    const Point_3& pi = P[i];
    const Point_3& pm = P[m];
    const Point_3& pk = P[k];
    double ang1=0, ang2=0;

    if(!Q.empty()){
      if(lambda[i*n+m] != -1){
        const Point_3& pim = P[lambda[i*n+m]];
        ang1 = 180 - CGAL::abs(CGAL::Mesh_3::dihedral_angle(pi,pm,pk, pim));
      }
      if(lambda[m*n+k] != -1){
        const Point_3& pmk = P[lambda[m*n+k]];
        ang2 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(pm,pk,pi, pmk));
      }
    }
    w = std::make_pair((std::max)(ang1, ang2), std::sqrt(CGAL::squared_area(pi,pm,pk)));
  }

  Weight operator+(const Weight& w2) const
  {
    return Weight((std::max)(w.first, w2.w.first), w.second + w2.w.second);
  }

  bool operator<(const Weight& w2) const
  {
    if(w.first == w2.w.first)
    { return w.second < w2.w.second; }
    return w.first < w2.w.first;
  }

  bool operator==(const Weight& w2) const
  {
    return w.first == w2.w.first && w.second == w2.w.second;
  }

  bool operator!=(const Weight& w2) const 
  {
    return !(*this == w2);
  }
};

template<typename K>
class Triangulate_hole_polyline_DT 
{
public:
  typedef typename K::Point_3 Point_3;
  typedef std::vector<Point_3> Polyline_3;
private:
  typedef Triangulation_vertex_base_with_info_3<int, K>  VB_with_id;
  typedef Triangulation_data_structure_3<VB_with_id>     TDS;
  typedef Delaunay_triangulation_3<K, TDS>               Triangulation;

  typedef typename Triangulation::Finite_edges_iterator  Finite_edges_iterator;
  typedef typename Triangulation::Finite_cells_iterator  Finite_cells_iterator;
  typedef typename Triangulation::Facet_circulator       Facet_circulator;

  typedef typename Triangulation::Cell_handle            Cell_handle;
  typedef typename Triangulation::Vertex_handle          Vertex_handle;
  typedef typename Triangulation::Edge                   Edge;
  typedef typename Triangulation::Facet                  Facet;
  
  template<class Point>
  struct Auto_count {
    typedef std::pair<Point, int> result_type;

    Auto_count() : count(0)  { }

    template<class Point>
    std::pair<Point, int> operator()(const Point& p) const {
      return std::make_pair(p, count++);
    }

    mutable int count;
  };

public:

  template <typename OutputIteratorValueType, typename OutputIterator>
  OutputIterator
  trace(int n,
        const std::vector<int>& lambda, 
        int i, 
        int k, 
        OutputIterator out)
  {
    if(i + 1 == k) { return out; }

    if(i+2 == k){
      *out++ = OutputIteratorValueType(i%n, (i+1)%n, k%n);
    } 
    else {
      int la = lambda[i*n + k];
      out = trace<OutputIteratorValueType>(n, lambda, i, la, out);
      *out++ = OutputIteratorValueType(i%n, la%n, k%n);
      out = trace<OutputIteratorValueType>(n, lambda, la, k, out);
    }
    return out;
  }

  template <typename OutputIteratorValueType, typename OutputIterator>
  OutputIterator 
  triangulate(const Polyline_3& P, 
              const Polyline_3& Q,
              OutputIterator out)
  {
    CGAL_assertion(P.front() == P.back());
    CGAL_assertion(Q.empty() || (Q.front() == Q.back()));
    CGAL_assertion(Q.empty() || (P.size() == Q.size()));

    int n = P.size() - 1; // because the first and last point are equal

    Triangulation T(boost::make_transform_iterator(P.begin(), Auto_count<Point_3>()),
                    boost::make_transform_iterator(--P.end(), Auto_count<Point_3>()));
    T.infinite_vertex()->info() = -1;

    Finite_edges_iterator v0_vn_edge; // v0 vn-1 edge
    std::vector<bool> edge_exist(n, false);
    for(Finite_edges_iterator eb = T.finite_edges_begin(); eb != T.finite_edges_end(); ++eb) {
      Cell_handle  ch = eb->first;
      int v0_id = ch->vertex(eb->second)->info();
      int v1_id = ch->vertex(eb->third)->info();
      if(v0_id > v1_id) { std::swap(v0_id, v1_id); }
      // to start from v0 vn-1 edge
      if(v0_id == 0 && v1_id == n-1) { v0_vn_edge = eb; }
      // check whether the edge is border edge
      if(v0_id + 1 == v1_id) { edge_exist[v0_id] = true; }
      else if(v0_id == 0 && v1_id == P.size() -2) { edge_exist[v1_id] = true; }
    }

    int not_exists = 0;
    for(std::vector<bool>::iterator it = edge_exist.begin(); it != edge_exist.end(); ++it) {
      if(!(*it)) { not_exists++; }
    }

    if(not_exists > 0) {
      CGAL_TRACE_STREAM << "Not all border edges are included in 3D Triangulation!" << std::endl;
      CGAL_TRACE_STREAM << "  Not inside " << not_exists << " of " << (P.size() - 1) << std::endl;
      CGAL_warning(!"Returning no output!");
      return out;
    }

    if(T.dimension() < 2) {
      CGAL_TRACE_STREAM << "Dimension of 3D Triangulation is above 2!" << std::endl;
      CGAL_warning(!"Returning no output!");
      return out;
    }

    if(T.dimension() == 2) {
      // in case of dimension 2 return triangles directly
      for(Finite_cells_iterator cb = T.finite_cells_begin(); cb != T.finite_cells_begin(); ++cb) {
        int v0(cb->vertex(0)->info()), 
          v1(cb->vertex(1)->info()), 
          v2(cb->vertex(2)->info());

        // TODO not sure about orientation
        *out++ = OutputIteratorValueType(v0, v1, v2);
      }
      return out;
    }
    
    std::vector<Weight> W(n*n,Weight(0,0));
    std::vector<int> lambda(n*n,-1);
    bool success = triangulate_DT(P, Q, W, lambda, *v0_vn_edge, T, n);
    CGAL_assertion(success);

    return trace<OutputIteratorValueType>(n, lambda, 0, n-1, out);
  }

  // this also may return infinite vertex
  std::pair<int, int> // <vertex id, index in cell>
  get_facet_remaining_vertex(Cell_handle ch, Facet f, int v0, int v1) 
  {
    int f0(ch->vertex(f.second)->info());

    for(int i = 0; i < 4; ++i) {
      int f3 = ch->vertex(i)->info();
      if(f3 != f0 && f3 != v0 && f3 != v1) {
        return std::make_pair(f3, i); 
      }
    }
    CGAL_assertion(false);
    return std::make_pair(-1, -1);
  }

  int get_vertex_index(Cell_handle ch, int id) {
    for(int i = 0; i < 4; ++i) {
      int v = ch->vertex(i)->info();
      if(v == id) { return i; }
    }
    CGAL_assertion(false);
    return -1;
  }

  bool 
  triangulate_DT(const Polyline_3& P, 
                     const Polyline_3& Q, 
                     std::vector<Weight>& W, 
                     std::vector<int>& lambda, 
                     Edge e,
                     Triangulation& T,
                     int n)
  {
    int v0 = e.first->vertex(e.second)->info();
    int v1 = e.first->vertex(e.third)->info();
    if(v0 > v1) { std::swap(v0, v1); }
    // edge can not be incident to infinite vertex
    CGAL_assertion(v0 != -1);

    // the range is previously processed
    if(W[v0*n + v1].w.first != 0 && W[v0*n + v1].w.second != 0) 
    { return true; }

    // border edge - just return true
    // should not check v0 = 0, v1 = n-1, because it is the initial edge where the algorithm starts
    if(v0 + 1 == v1)
    { return true; }

    // one triangle remains
    if(v0 + 2 == v1){
      // check whether it is included in DT
      Facet_circulator fb(T.incident_facets(e)), end(fb);
      do {
        int f3 = get_facet_remaining_vertex(fb->first, *fb, v0, v1).first;
        if(f3 == v0 + 1) {
          W[v0*n + v1] = Weight(P, Q, v0, v0+1, v1);
          lambda[v0*n + v1] = v0+1;
          return true;
        }
      } while(++fb != end);
      return false;
    }

    int m_min = -1;
    Weight w_min((std::numeric_limits<double>::max)(), 
                 (std::numeric_limits<double>::max)());
    bool found_any = false;
    Facet_circulator fb(T.incident_facets(e)), end(fb);
    do {
      int v2, v2_cell_index;
      boost::tie(v2, v2_cell_index) = get_facet_remaining_vertex(fb->first, *fb, v0, v1);
      
      if(v2 < v0 || v2 > v1) { continue; } // this will also skip infinite vertex

      Edge e0 = Edge(fb->first, get_vertex_index(fb->first, v0) , v2_cell_index);
      Edge e1 = Edge(fb->first, get_vertex_index(fb->first, v1) , v2_cell_index);

      bool found = triangulate_DT(P, Q, W, lambda, e0, T, n) &&
                   triangulate_DT(P, Q, W, lambda, e1, T, n);
      if(!found) { continue; }

      CGAL_assertion(W[v0*n + v2].w.first != -1 && W[v0*n + v2].w.second != -1);
      CGAL_assertion(W[v2*n + v1].w.first != -1 && W[v2*n + v1].w.second != -1);

      Weight w = W[v0*n + v2] + W[v2*n + v1] + Weight(P,Q, v0,v2,v1, lambda);
      if(w < w_min){
        w_min = w;
        m_min = v2;
      }
    } while(++fb != end);

    if(m_min == -1) // which means no triangulation exists between v0 - v1
    { return false; }

    W[v0*n+v1] = w_min;
    lambda[v0*n+v1] = m_min;
    return true;
  }
  
};

template <typename K>
class Triangulate_hole_polyline {
public:
  typedef typename K::Point_3 Point_3;
  typedef std::vector<Point_3> Polyline_3;

private:

  template <typename OutputIteratorValueType, typename OutputIterator>
  OutputIterator
  trace(int n,
        const std::vector<int>& lambda, 
        int i, 
        int k, 
        OutputIterator out)
  {
    if(i + 1 == k) { return out; }

    if(i+2 == k){
      *out++ = OutputIteratorValueType(i%n, (i+1)%n, k%n);
    } 
    else {
      int la = lambda[i*n + k];
      out = trace<OutputIteratorValueType>(n, lambda, i, la, out);
      *out++ = OutputIteratorValueType(i%n, la%n, k%n);
      out = trace<OutputIteratorValueType>(n, lambda, la, k, out);
    }
    return out;
  }
  
public:

  template <typename OutputIteratorValueType, typename OutputIterator>
  OutputIterator 
  triangulate(const Polyline_3& P, 
             const Polyline_3& Q,
             OutputIterator out,
             const std::set<std::pair<int, int> >& existing_edges)
  {
    CGAL_assertion(P.front() == P.back());
    CGAL_assertion(Q.empty() || (Q.front() == Q.back()));
    CGAL_assertion(Q.empty() || (P.size() == Q.size()));
    
    int n = P.size() - 1; // because the first and last point are equal
    std::vector<Weight> W(n*n,Weight(0,0));
    std::vector<int> lambda(n*n,-1);
    
    // calculate weights for boundary triangles
    for(int i=0; i < n-2; ++i){
      if(existing_edges.find(std::make_pair(i, i+2)) != existing_edges.end()) {
        W[i*n + (i+2)] = Weight(-1,-1);
        lambda[i*n + (i+2)] = -1;
      }
      else {
        W[i*n + (i+2)] = Weight(P, Q, i, i+1, i+2);
        lambda[i*n + (i+2)] = i+1;
      }
    }
    
    for(int j = 3; j< n; ++j){ // 3 - 4 - 5 (range)

      for(int i=0; i<n-j; ++i){ // 0-3, 1-4, 2-5 find min
        int k = i+j;
        int m_min = -1;
        Weight w_min((std::numeric_limits<double>::max)(), 
                     (std::numeric_limits<double>::max)());

        for(int m = i+1; m<k; ++m){ 

          if( existing_edges.find(std::make_pair(i,m)) != existing_edges.end() ||
              existing_edges.find(std::make_pair(i,k)) != existing_edges.end() ||
              existing_edges.find(std::make_pair(m,k)) != existing_edges.end() ) {
            // we can not construct i-m-k triangle
            continue;
          }
          // now the regions i-m and m-k might be valid(constructed) patches,
          // if not, we can not construct i-m-k triangle
          if( W[i*n + m] == Weight(-1, -1) || W[m*n + k] == Weight(-1, -1) ) {
            continue;
          }

          Weight w = W[i*n + m] + W[m*n + k] + Weight(P,Q,i,m,k, lambda);
          if(w < w_min){
            w_min = w;
            m_min = m;
          }
        }
        // if any found, update weight
        if(m_min != -1) { W[i*n+k] = w_min; }
        else            { W[i*n+k] = Weight(-1,-1); }
        lambda[i*n+k] = m_min;
      }
    }

    if(lambda[0*n + (n-1)] == -1) {
      CGAL_warning(!"Due to existing edges, no possible triangulation is found!");
      return out;
    }

    return  trace<OutputIteratorValueType>(n, lambda, 0, n-1, out);
  }
};

template <typename OutputIteratorValueType, typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          InputIterator qbegin, InputIterator qend, 
                          OutputIterator out,
                          const std::set<std::pair<int, int> >& existing_edges) 
{
  typedef typename CGAL::Kernel_traits< typename std::iterator_traits<InputIterator>::value_type>::Kernel Kernel;
  typedef CGAL::internal::Triangulate_hole_polyline<Kernel> Fill;
  typename Fill::Polyline_3 P(pbegin, pend);
  typename Fill::Polyline_3 Q(qbegin, qend);
  if(P.front() != P.back()){
    P.push_back(P.front());
    if( !Q.empty() && P.size() > Q.size()) {
      Q.push_back(Q.front());
    }
  }
  Fill fill;
  return fill.template triangulate<OutputIteratorValueType>(P,Q,out,existing_edges);
}

} // namespace internal

/*!
Creates triangles to fill the hole defined by points in the range `(pbegin,pend)`.
The range `(qbegin,qend)` indicate for each pair of consecutive points in the aforementioned range,
the third point of the facet this segment is incident to. Triangles are put into `out`
using the indices of the input points in the range `(pbegin,pend)`.
@tparam OutputIteratorValueType value type of OutputIterator having a constructor `OutputIteratorValueType(int p0, int p1, int p2)` available. 
        It is default to value_type_traits<OutputIterator>::type, and can be omitted when the default is fine
@tparam InputIterator iterator over input points
@tparam OutputIterator iterator over patch triangles
@param pbegin first iterator of the range of points
@param pend past-the-end iterator of the range of points
@param qbegin first iterator of the range of third points
@param qend past-the-end iterator of the range of third points
@param out iterator over output patch triangles
*/
template <typename OutputIteratorValueType, typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          InputIterator qbegin, InputIterator qend, 
                          OutputIterator out)
{
  return internal::triangulate_hole_polyline<OutputIteratorValueType>
    (pbegin, pend, qbegin, qend, out, std::set<std::pair<int, int> >());
}

// overload for OutputIteratorValueType
template <typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          InputIterator qbegin, InputIterator qend, 
                          OutputIterator out)
{
  return triangulate_hole_polyline<typename value_type_traits<OutputIterator>::type>
    (pbegin, pend, qbegin, qend, out);
}


/*!
Creates triangles to fill the hole defined by points in the range `(pbegin,pend)`.
Triangles are put into `out` using the indices of the input points in the range `(pbegin,pend)`.
@tparam OutputIteratorValueType value type of OutputIterator having a constructor `OutputIteratorValueType(int p0, int p1, int p2)` available. 
        It is default to value_type_traits<OutputIterator>::type, and can be omitted when the default is fine
@tparam InputIterator iterator over input points
@tparam OutputIterator iterator over output patch triangles
@param pbegin first iterator of the range of points
@param pend past-the-end iterator of the range of points
@param out iterator over output patch triangles
*/
template <typename OutputIteratorValueType, typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          OutputIterator out)
{

  typename Fill::Polyline_3 Q;
  return triangulate_hole_polyline<OutputIteratorValueType>
    (pbegin, pend, Q.begin(), Q.end(), out);
}

// overload for OutputIteratorValueType
template <typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          OutputIterator out)
{
  return triangulate_hole_polyline<typename value_type_traits<OutputIterator>::type>
    (pbegin, pend, out);
}

} // namespace CGAL

#endif // CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H
