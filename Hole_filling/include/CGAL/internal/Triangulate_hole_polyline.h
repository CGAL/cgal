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
/************************************************************************/
/* Common functionality: Weight + Tracer
/************************************************************************/
struct Weight {

  std::pair<double,double> w;

  Weight() : w(std::make_pair<double>(0, 0)) { }

  Weight(double angle, double area) : w(angle, area) { }

  template<class Point_3>
  Weight(const std::vector<Point_3>& P, 
         const std::vector<Point_3>& Q, 
         int i, int j, int k, 
         const std::vector<int>& lambda)
  {
    CGAL_assertion(i < j);
    CGAL_assertion(j < k);
    int n = P.size() -1; // because the first and last point are equal
    
    // The CGAL::dihedral angle is measured between the oriented triangles, that is it goes from [-pi, pi]
    // What we need is the angle between the normals of the triangles between [0, pi]
    double ang_max = 0;

    if(!Q.empty()){
      // Test each edge
      int vertices[] = {i, j, k};
      for(int e = 0; e < 3; ++e) 
      {
        int v0      = vertices[e];
        int v1      = vertices[(e+1)%3];
        int v_other = vertices[(e+2)%3];
        double angle = 0;
        // check whether the edge is border
        if(v0 + 1 == v1 || v0 == n-1 && v1 == 0) {
          angle = 180 - CGAL::abs( 
            CGAL::Mesh_3::dihedral_angle(P[v0],P[v1],P[v_other],Q[v0]) );
        }
        else {
          if(e == 2) { continue; }
          if(lambda[v0*n+v1] != -1){
            const Point_3& p01 = P[lambda[v0*n+v1]];
            angle = 180 - CGAL::abs( 
              CGAL::Mesh_3::dihedral_angle(P[v0],P[v1],P[v_other],p01) );
          }
        }
        ang_max = (std::max)(ang_max, angle);
      }
    } // if !Q.empty() 
    w = std::make_pair(ang_max, std::sqrt(CGAL::squared_area(P[i],P[j],P[k])));
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

struct Tracer {

  template <typename OutputIteratorValueType, typename OutputIterator>
  OutputIterator
  trace(int n,
        const std::vector<int>& lambda, 
        int i, 
        int k, 
        OutputIterator out)
  {
    CGAL_assertion(i >= 0 && i < n);
    CGAL_assertion(k >= 0 && k < n);

    if(i + 1 == k) { return out; }

    int la = lambda[i*n + k];
    CGAL_assertion(la >= 0 && la < n);

    out = trace<OutputIteratorValueType>(n, lambda, i, la, out);
    *out++ = OutputIteratorValueType(i, la, k);
    out = trace<OutputIteratorValueType>(n, lambda, la, k, out);
    return out;
  }
};

/************************************************************************/
/* Triangulate hole with support of 3D Triangulation
/************************************************************************/

// to support incident_facets(Edge e) function for both dimension 2 and 3
template<unsigned int Dimension, class Triangulator>
struct Incident_facet_circulator;

// Use the fact that an edge can be incident to 2 facets in dimension 2
// and all valid facets (which contains finite + infinite vertices but not the NULL vertex) are
// pointed by index 3 in cells
template<class Triangulator>
struct Incident_facet_circulator<2, Triangulator>
{
  typedef typename Triangulator::Facet         Facet;
  typedef typename Triangulator::Edge          Edge;
  typedef typename Triangulator::Triangulation Triangulation;
  
  Incident_facet_circulator(Edge e, Triangulation* t) 
  {
    f1 = Facet(e.first, 3);
    int remaining_index = 0;
    for( ; remaining_index < 3; ++remaining_index) {
      if(remaining_index != e.second && remaining_index != e.third) {
        break;
      }
    }
    f2 = Facet(e.first->neighbor(remaining_index), 3);
    it = f1;
  }
  Incident_facet_circulator& operator++() {
    it = it == f1 ? f2 : f1;
    return *this;
  }
  operator bool() const { return it != f1; }
  const Facet* operator->() const { return &it; }
  Facet operator*() const { return it; }

  Facet f1, f2, it;
};

// Just a wrapper around Facet_circulator
template<class Triangulator>
struct Incident_facet_circulator<3, Triangulator>
{
  typedef typename Triangulator::Facet            Facet;
  typedef typename Triangulator::Edge             Edge;
  typedef typename Triangulator::Triangulation    Triangulation;
  typedef typename Triangulator::Facet_circulator Facet_circulator;

  Incident_facet_circulator(Edge e, Triangulation* t)
    : it(t->incident_facets(e)), end(it)
  { }
  Incident_facet_circulator& operator++() {
    ++it;
    return *this;
  }
  operator bool() const { return it != end; }
  Facet_circulator operator->() const { return it; }
  Facet operator*() const { return *it; }

  Facet_circulator it;
  Facet_circulator end;
};

template<typename K>
class Triangulate_hole_polyline_DT 
{
public:
  typedef Triangulate_hole_polyline_DT<K> Self;
  typedef typename K::Point_3 Point_3;
  typedef std::vector<Point_3> Polyline_3;

  typedef Triangulation_vertex_base_with_info_3<int, K>  VB_with_id;
  typedef Triangulation_data_structure_3<VB_with_id>     TDS;
  typedef Delaunay_triangulation_3<K, TDS>               Triangulation;

  typedef typename Triangulation::Finite_edges_iterator  Finite_edges_iterator;
  typedef typename Triangulation::Facet_circulator       Facet_circulator;

  typedef typename Triangulation::Cell_handle            Cell_handle;
  typedef typename Triangulation::Vertex_handle          Vertex_handle;
  typedef typename Triangulation::Edge                   Edge;
  typedef typename Triangulation::Facet                  Facet;
  
  template<class Point>
  struct Auto_count {
    typedef std::pair<Point, int> result_type;

    Auto_count() : count(0)  { }
    std::pair<Point, int> operator()(const Point& p) const 
    { return std::make_pair(p, count++); }
    mutable int count;
  };

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

    Triangulation T(boost::make_transform_iterator(P.begin(), Auto_count<Point_3>()),
                    boost::make_transform_iterator(--P.end(), Auto_count<Point_3>()));
    T.infinite_vertex()->info() = -1;

    // Check whether all edges are included in DT, and get v0-vn-1 edge
    int nb_exists = 0;
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
      int border_id = -1;
      if(v0_id + 1 == v1_id)              { border_id = v0_id; }
      else if(v0_id == 0 && v1_id == n-1) { border_id = v1_id; }
      if(border_id != -1 && !edge_exist[border_id]) {
        ++nb_exists;
        edge_exist[border_id] = true;
      }
    }
    CGAL_assertion(n >= nb_exists);

    if(nb_exists != n) {
      CGAL_TRACE_STREAM << "Not all border edges are included in 3D Triangulation!" << std::endl;
      CGAL_TRACE_STREAM << "  Not inside " << (n - nb_exists) << " of " << n << std::endl;
      CGAL_warning(!"Returning no output!");
      return out;
    }

    if(T.dimension() < 2) {
      CGAL_TRACE_STREAM << "Dimension of 3D Triangulation is above 2!" << std::endl;
      CGAL_warning(!"Returning no output!");
      return out;
    }
    
    std::vector<Weight> W(n*n,Weight(0,0));
    std::vector<int> lambda(n*n,-1);
    if(T.dimension() == 3) {
      triangulate_DT<Incident_facet_circulator<3, Self> >
        (P, Q, W, lambda, *v0_vn_edge, T, n, existing_edges);
    }
    else {
      CGAL_assertion(T.dimension() == 2);
      triangulate_DT<Incident_facet_circulator<2, Self> >
        (P, Q, W, lambda, *v0_vn_edge, T, n, existing_edges);
    }
    
    if(lambda[0*n + (n-1)] == -1) {
      CGAL_warning(!"No possible triangulation is found!");
      return out;
    }
    CGAL_TRACE_STREAM << "Triangulation Weight = [max-angle: "
                      << W[n-1].w.first << ", area: " << W[n-1].w.second <<"]" << std::endl;
    return Tracer().trace<OutputIteratorValueType>(n, lambda, 0, n-1, out);
  }

private:
  // Finds other vertex then v0 and v1 in facet f
  // Note that this may return infinite vertex
  std::pair<int, int> // <vertex id(info), index in cell>
  get_facet_remaining_vertex(Facet f, int v0_info, int v1_info) 
  {
    // warning: it should be designed to handle dimension 2 (e.g. f.first->vertex(3)->info() will crash)
    for(int i = 0; i < 4; ++i) {
      if(i == f.second) { continue; }
      int f3 = f.first->vertex(i)->info();
      if(f3 != v0_info && f3 != v1_info) {
        return std::make_pair(f3, i); 
      }
    }
    CGAL_assertion(false);
    return std::make_pair(-1, -1);
  }

  int get_vertex_index(Cell_handle ch, int info) {
    for(int i = 0; i < 4; ++i) {
      int v = ch->vertex(i)->info();
      if(v == info) { return i; }
    }
    CGAL_assertion(false);
    return -1;
  }

  template<class IncidentFacetCirculator>
  void 
  triangulate_DT( const Polyline_3& P, 
                  const Polyline_3& Q, 
                  std::vector<Weight>& W, 
                  std::vector<int>& lambda, 
                  Edge e,
                  Triangulation& T,
                  int n,
                  const std::set<std::pair<int, int> >& existing_edges)
  {
    /**********************************************************************
     *  + Default W value is (0,0), default lambda value is -1.
     *  + W value (0,0) is used to check whether the region (v0-v1) is processed.
     *  + If a range v0-v1 does not contains any possible triangulation, then W[v0,v1] = (-1,-1) and lambda[v0,v1] = -1
     */
    int v0 = e.first->vertex(e.second)->info();
    int v1 = e.first->vertex(e.third)->info();
    if(v0 > v1) { std::swap(v0, v1); }
    
    // edge can not be incident to infinite vertex
    CGAL_assertion(v0 != -1);

    // the range is previously processed
    if( W[v0*n + v1] != Weight(0, 0) ) { return; }

    // border edge - just return
    // should not check v0 = 0, v1 = n-1, because it is the initial edge where the algorithm starts
    if(v0 + 1 == v1) { return; }

    // check whether the edge is valid
    if(existing_edges.find(std::make_pair(v0, v1)) != existing_edges.end()) {
      W[v0*n + v1] = Weight(-1,-1);
      return;
    }

    int m_min = -1;
    Weight w_min((std::numeric_limits<double>::max)(), 
                 (std::numeric_limits<double>::max)());
    bool found_any = false;
    IncidentFacetCirculator fb(e, &T);
    do {
      int v2, v2_cell_index;
      boost::tie(v2, v2_cell_index) = get_facet_remaining_vertex(*fb, v0, v1);

      if(v2 < v0 || v2 > v1) { continue; } // this will also skip infinite vertex

      if( existing_edges.find(std::make_pair(v0,v2)) != existing_edges.end() ||
          existing_edges.find(std::make_pair(v2,v1)) != existing_edges.end() ) {
          // we can not construct i-m-k triangle
          continue;
      }

      Edge e0 = Edge(fb->first, get_vertex_index(fb->first, v0) , v2_cell_index);
      Edge e1 = Edge(fb->first, get_vertex_index(fb->first, v1) , v2_cell_index);

      triangulate_DT<IncidentFacetCirculator>(P, Q, W, lambda, e0, T, n, existing_edges); // v0-v2
      triangulate_DT<IncidentFacetCirculator>(P, Q, W, lambda, e1, T, n, existing_edges); // v2-v1
      if( W[v0*n + v2] == Weight(-1,-1) || W[v2*n + v1] == Weight(-1,-1) )
      { continue; }

      Weight w = W[v0*n + v2] + W[v2*n + v1] + Weight(P,Q, v0,v2,v1, lambda);
      if(w < w_min){
        w_min = w;
        m_min = v2;
      }
    } while(++fb);

    if(m_min != -1) { W[v0*n+v1] = w_min; }
    else            { W[v0*n+v1] = Weight(-1,-1); } // which means no triangulation exists between v0 - v1
    lambda[v0*n+v1] = m_min;
  }
  
};

/************************************************************************/
/* Triangulate hole by using all search space
/************************************************************************/
template <typename K>
class Triangulate_hole_polyline {
public:
  typedef typename K::Point_3 Point_3;
  typedef std::vector<Point_3> Polyline_3;
  
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
    
    for(int j = 2; j< n; ++j){ // determines range (2 - 3 - 4 )
      for(int i=0; i<n-j; ++i){ // iterates over ranges and find min triangulation in those ranges 
        int k = i+j;            // like [0-2, 1-3, 2-4, ...], [0-3, 1-4, 2-5, ...]
        int m_min = -1;
        Weight w_min((std::numeric_limits<double>::max)(), 
                     (std::numeric_limits<double>::max)());

        if(existing_edges.find(std::make_pair(i,k)) == existing_edges.end()) 
        {
          for(int m = i+1; m<k; ++m){ 
            if( existing_edges.find(std::make_pair(i,m)) != existing_edges.end() ||
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
        }
        
        if(m_min != -1) { W[i*n+k] = w_min; }
        else            { W[i*n+k] = Weight(-1,-1); }
        lambda[i*n+k] = m_min; // if not found lambda should be -1
      }
    }

    if(lambda[0*n + (n-1)] == -1) {
      CGAL_warning(!"Due to existing edges, no possible triangulation is found!");
      return out;
    }
    CGAL_TRACE_STREAM << "Triangulation Weight = [max-angle: "
                      << W[n-1].w.first << ", area: " << W[n-1].w.second <<"]" << std::endl;
    return Tracer().trace<OutputIteratorValueType>(n, lambda, 0, n-1, out);
  }
};

template <typename OutputIteratorValueType, typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          InputIterator qbegin, InputIterator qend, 
                          OutputIterator out,
                          const std::set<std::pair<int, int> >& existing_edges,
                          bool use_delaunay_triangulation) 
{
  typedef typename CGAL::Kernel_traits< typename std::iterator_traits<InputIterator>::value_type>::Kernel Kernel;
  typedef CGAL::internal::Triangulate_hole_polyline_DT<Kernel> Fill_DT;
  typedef CGAL::internal::Triangulate_hole_polyline<Kernel>    Fill;

  typename Fill::Polyline_3 P(pbegin, pend);
  typename Fill::Polyline_3 Q(qbegin, qend);
  if(P.front() != P.back()){
    P.push_back(P.front());
    if( !Q.empty() && P.size() > Q.size()) {
      Q.push_back(Q.front());
    }
  }

  return use_delaunay_triangulation ?
         Fill_DT().template triangulate<OutputIteratorValueType>(P,Q,out,existing_edges) :
         Fill().template triangulate<OutputIteratorValueType>(P,Q,out,existing_edges);
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
                          OutputIterator out, bool use_delaunay_triangulation = false)
{
  return internal::triangulate_hole_polyline<OutputIteratorValueType>
    (pbegin, pend, qbegin, qend, out, std::set<std::pair<int, int> >(), use_delaunay_triangulation);
}

// overload for OutputIteratorValueType
template <typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          InputIterator qbegin, InputIterator qend, 
                          OutputIterator out, bool use_delaunay_triangulation = false)
{
  return triangulate_hole_polyline<typename value_type_traits<OutputIterator>::type>
    (pbegin, pend, qbegin, qend, out, use_delaunay_triangulation);
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
                          OutputIterator out, bool use_delaunay_triangulation = false)
{
  typedef typename CGAL::Kernel_traits< typename std::iterator_traits<InputIterator>::value_type>::Kernel Kernel;
  typedef std::vector<typename Kernel::Point_3> Polyline_3;
  Polyline_3 Q;
  return triangulate_hole_polyline<OutputIteratorValueType>
    (pbegin, pend, Q.begin(), Q.end(), out, use_delaunay_triangulation);
}

// overload for OutputIteratorValueType
template <typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          OutputIterator out, bool use_delaunay_triangulation = false)
{
  return triangulate_hole_polyline<typename value_type_traits<OutputIterator>::type>
    (pbegin, pend, out, use_delaunay_triangulation);
}

} // namespace CGAL

#endif // CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H
