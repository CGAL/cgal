#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <vector>
#include <stack>
#include <map>
#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {
namespace internal {

/************************************************************************/
/* Lookup tables
/************************************************************************/
// Wrapper around vector
template<class T>
class Lookup_table {
public:
  Lookup_table(int n, const T& t) : n(n), table(n*n, t) { }
  void put(int i, int j, const T& t) {
    CGAL_assertion(bound_check(i,j));
    table[i*n + j] = t;
  }
  const T& get(int i, int j) const {
    CGAL_assertion(bound_check(i,j));
    return table[i*n + j];
  }

  int n;
private:
  bool bound_check(int i, int j) const {
    CGAL_assertion(i >= 0 && i < n);
    CGAL_assertion(j >= 0 && j < n);
    CGAL_assertion(i < j); 
    // previous implementation was based on directly vector and i supposed to be always smaller than j.
    // this check actually can be removed and i =min(i,j) j = max(i,j) can be used for reflexive access 
    return true;
  }
  std::vector<T> table;
};

// Wrapper around map, where if i,j is not found a default value is returned,
// and if default value inserted i,j erased.
template<class T>
class Lookup_table_map {
public:
  Lookup_table_map(int n, const T& default) : n(n), default(default) { }

  void put(int i, int j, const T& t) {
    CGAL_assertion(bound_check(i,j));
    
    if(t == default) {
      table.erase(std::make_pair(i,j));
      return;
    }
    // table[std::make_pair(i,j)] = t (to not require def constructor)
    table.insert(std::make_pair(std::make_pair(i,j), default)).first->second = t;
  }
  const T& get(int i, int j) const {
    CGAL_assertion(bound_check(i,j));
    typename std::map<std::pair<int,int>, T>::const_iterator ij 
      = table.find(std::make_pair(i,j));
    if(ij != table.end()) {
      return ij->second;
    }
    return default;
  }

  int n;
private:
  bool bound_check(int i, int j) const {
    CGAL_assertion(i >= 0 && i < n);
    CGAL_assertion(j >= 0 && j < n);
    CGAL_assertion(i < j);
    return true;
  }
  std::map<std::pair<int,int>, T> table;
  T default;
};
/************************************************************************/
/* Is_valid classes (to be used in Weight_calculator)
/************************************************************************/
// to be used in existing_edges set
struct Edge_comp {
  bool operator()(const std::pair<int, int>& p0, const std::pair<int, int>& p1) const {
    int p0_min = (std::min)(p0.first, p0.second);
    int p0_max = (std::max)(p0.first, p0.second);
    int p1_min = (std::min)(p1.first, p1.second);
    int p1_max = (std::max)(p1.first, p1.second);
    if(p0_min == p1_min) {
      return p0_max < p1_max;
    }
    return p0_min < p1_min;
  }
};
typedef std::set<std::pair<int, int>, Edge_comp> Edge_set;

struct Is_valid_existing_edges 
{
  Is_valid_existing_edges(const Edge_set& existing_edges) : existing_edges(existing_edges) { }

  template<class Point_3, class LookupTable>
  bool operator()(const std::vector<Point_3>&, const std::vector<Point_3>&, 
                  int v0, int v1, int v2, const LookupTable&) const 
  {
    return existing_edges.find(std::make_pair(v0,v1)) == existing_edges.end() &&
           existing_edges.find(std::make_pair(v1,v2)) == existing_edges.end() &&
           existing_edges.find(std::make_pair(v0,v2)) == existing_edges.end();
  }
  const Edge_set& existing_edges;
};

struct Is_valid_degenerate_triangle 
{
  template<class Point_3, class LookupTable>
  bool operator()(const std::vector<Point_3>& P, const std::vector<Point_3>&, 
                  int v0, int v1, int v2, const LookupTable&) const 
  {
    return !CGAL::collinear(P[v0], P[v1], P[v2]);
  }
};

// Combine above two
struct Is_valid_existing_edges_and_degenerate_triangle 
{
  Is_valid_existing_edges_and_degenerate_triangle(const Edge_set& edges) : is_valid_edges(edges) { }

  template<class Point_3, class LookupTable>
  bool operator()(const std::vector<Point_3>& P, const std::vector<Point_3>& Q, 
                  int v0, int v1, int v2, const LookupTable& t) const 
  {
    return Is_valid_degenerate_triangle()(P,Q,v0,v1,v2,t) && is_valid_edges(P,Q,v0,v1,v2,t);
  }

  Is_valid_existing_edges is_valid_edges;
};

/************************************************************************/
/* Weights
/************************************************************************/

// Weight calculator class is both responsible from calculating weights, and checking validity of triangle
template<class Weight_, class IsValid>
struct Weight_calculator 
{
  typedef Weight_ Weight;
  Weight_calculator(const IsValid& is_valid = IsValid()) : is_valid(is_valid) { }

  template<class Point_3, class LookupTable>
  Weight operator()(const std::vector<Point_3>& P, 
                    const std::vector<Point_3>& Q, 
                    int i, int j, int k, 
                    const LookupTable& lambda) const 
  {
    if( !is_valid(P, Q, i,j,k, lambda) ) 
    { return Weight::NOT_VALID(); }
    return Weight(P, Q, i,j,k, lambda);
  }

  IsValid is_valid;
};

class Weight_min_max_dihedral_and_area 
{
  template<class Weight_, class IsValid>
  friend struct Weight_calculator;

public:
  // these two should not be used (used in test code)
  std::pair<double,double> w;
  Weight_min_max_dihedral_and_area(double angle, double area) : w(angle, area) { }

// below required by Weight concept
private:
  template<class Point_3, class LookupTable>
  Weight_min_max_dihedral_and_area(const std::vector<Point_3>& P, 
                                   const std::vector<Point_3>& Q, 
                                   int i, int j, int k, 
                                   const LookupTable& lambda)
  {
    CGAL_assertion(i < j);
    CGAL_assertion(j < k);
    int n = P.size() -1; // because the first and last point are equal
    
    // The CGAL::dihedral angle is measured between the oriented triangles, that is it goes from [-pi, pi]
    // What we need is the angle between the normals of the triangles between [0, pi]
    double ang_max = 0;

    // Test each edge
    int vertices[] = {i, j, k};
    for(int e = 0; e < 3; ++e) 
    {
      int v0      = vertices[e];
      int v1      = vertices[(e+1)%3];
      int v_other = vertices[(e+2)%3];
      double angle = 0;
      // check whether the edge is border
      if( (v0 + 1 == v1 || v0 == n-1 && v1 == 0) && !Q.empty() ) {
        angle = 180 - CGAL::abs( 
          CGAL::Mesh_3::dihedral_angle(P[v0],P[v1],P[v_other],Q[v0]) );
      }
      else {
        if(e == 2) { continue; }
        if(lambda.get(v0, v1) != -1){
          const Point_3& p01 = P[lambda.get(v0, v1)];
          angle = 180 - CGAL::abs( 
            CGAL::Mesh_3::dihedral_angle(P[v0],P[v1],P[v_other],p01) );
        }
      }
      ang_max = (std::max)(ang_max, angle);
    }
   
    w = std::make_pair(ang_max, std::sqrt(CGAL::squared_area(P[i],P[j],P[k])));
  }

public:
  Weight_min_max_dihedral_and_area operator+(const Weight_min_max_dihedral_and_area& w2) const 
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());

    return Weight_min_max_dihedral_and_area((std::max)(w.first, w2.w.first), w.second + w2.w.second);
  }

  bool operator<(const Weight_min_max_dihedral_and_area& w2) const
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());

    if(w.first == w2.w.first)
    { return w.second < w2.w.second; }
    return w.first < w2.w.first;
  }

  bool operator==(const Weight_min_max_dihedral_and_area& w2) const 
  { return w.first == w2.w.first && w.second == w2.w.second; }

  bool operator!=(const Weight_min_max_dihedral_and_area& w2) const 
  { return !(*this == w2); }

  static const Weight_min_max_dihedral_and_area DEFAULT() // rule: x + DEFAULT() == x
  { return Weight_min_max_dihedral_and_area(0,0); }
  static const Weight_min_max_dihedral_and_area NOT_VALID() 
  { return Weight_min_max_dihedral_and_area(-1,-1); }

  friend std::ostream& operator<<(std::ostream& out, const Weight_min_max_dihedral_and_area& w) {
    out << "Max dihedral: " << w.w.first << ", Total area: " << w.w.second;
    return out;
  }
};

// For proof of concept. Tested weakly.
class Weight_total_edge {
private:
  double total_length;

  Weight_total_edge(double total_length = 0) : total_length(total_length) { }

public:
  template<class Point_3, class LookupTable>
  Weight_total_edge(const std::vector<Point_3>& P, 
                    const std::vector<Point_3>&, 
                    int i, int j, int k, 
                    const LookupTable&)
    : total_length(0)
  {
    CGAL_assertion(i < j);
    CGAL_assertion(j < k);
    int n = P.size() -1; // because the first and last point are equal

    // Test each edge
    int vertices[] = {i, j, k};
    for(int e = 0; e < 3; ++e) 
    {
      int v0      = vertices[e];
      int v1      = vertices[(e+1)%3];
      int v_other = vertices[(e+2)%3];
        
      // check whether the edge is border
      bool border = (v0 + 1 == v1) || (v0 == n-1 && v1 == 0);
      if(!border) {
        total_length += std::sqrt(CGAL::squared_distance(P[v0],P[v1]));
      }
    }
  }

  Weight_total_edge operator+(const Weight_total_edge& w2) const 
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());
    return Weight_total_edge(total_length + w2.total_length);
  }

  bool operator<(const Weight_total_edge& w2) const
  { 
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());
    return total_length < w2.total_length; 
  }

  bool operator==(const Weight_total_edge& w2) const 
  { return total_length == w2.total_length; }
  bool operator!=(const Weight_total_edge& w2) const 
  { return !(*this == w2); }

  static const Weight_total_edge DEFAULT() { return Weight_total_edge(0); } // rule: x + DEFAULT() == x
  static const Weight_total_edge NOT_VALID() { return Weight_total_edge(-1); }
  friend std::ostream& operator<<(std::ostream& out, const Weight_total_edge& w) {
    out << "Total edge length : " << w.total_length;
    return out;
  }
};

/************************************************************************/
/* Tracer
/************************************************************************/
template<class OutputIteratorValueType, class OutputIterator>
struct Tracer_polyline {
  Tracer_polyline(OutputIterator& out) : out(out) 
  { }

  template <class LookupTable>
  void
  operator()(const LookupTable& lambda, 
             int v0, 
             int v1)
  {
    const int n = lambda.n;
    std::stack<std::pair<int, int> > ranges;
    ranges.push(std::make_pair(v0, v1));

    while(!ranges.empty()) {
      std::pair<int, int> r = ranges.top(); 
      ranges.pop();
      CGAL_assertion(r.first >= 0 && r.first < n);
      CGAL_assertion(r.second >= 0 && r.second < n);

      if(r.first + 1 == r.second) { continue; }

      int la = lambda.get(r.first, r.second);
      CGAL_assertion(la >= 0 && la < n);
      CGAL_assertion(r.first < la && r.second > la);
      *out++ = OutputIteratorValueType(r.first, la, r.second);

      ranges.push(std::make_pair(r.first, la));
      ranges.push(std::make_pair(la, r.second));
    }
  }

  OutputIterator& out;
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
  
  Incident_facet_circulator(Edge e, const Triangulation&)
    : f1( Facet(e.first, 3) ),
      f2( Facet(e.first->neighbor(3 - e.second - e.third), 3) ),
      it(f1)
  {
     CGAL_assertion(f1 != f2);
     CGAL_assertion(e.second < 3 && e.third < 3);
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

  Incident_facet_circulator(Edge e, const Triangulation& T)
    : it(T.incident_facets(e)), end(it)
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

// By default Lookup_table_map is used, since Lookup_table requires n*n mem.
// Performance decrease is nearly 2x (for n = 10,000, for larger n Lookup_table just goes out of mem) 
template<
  class Kernel,
  class Tracer,
  class WeightCalculator,
  template <class> class LookupTable = Lookup_table_map
>
class Triangulate_hole_polyline_DT 
{
public:
  typedef Triangulate_hole_polyline_DT                        Self;
  typedef typename WeightCalculator::Weight                   Weight;
  typedef typename Kernel::Point_3                            Point_3;
  typedef std::vector<Point_3>                                Polyline_3;

  typedef Triangulation_vertex_base_with_info_3<int, Kernel>  VB_with_id;
  typedef Triangulation_data_structure_3<VB_with_id>          TDS;
  typedef Delaunay_triangulation_3<Kernel, TDS>               Triangulation;

  typedef typename Triangulation::Finite_edges_iterator       Finite_edges_iterator;
  typedef typename Triangulation::Facet_circulator            Facet_circulator;
  typedef typename Triangulation::Cell_handle                 Cell_handle;
  typedef typename Triangulation::Vertex_handle               Vertex_handle;
  typedef typename Triangulation::Edge                        Edge;
  typedef typename Triangulation::Facet                       Facet;

  Weight operator()(const Polyline_3& P, 
                    const Polyline_3& Q,
                    Tracer& tracer,
                    const WeightCalculator& WC)
  {
    struct Auto_count {
      typedef std::pair<Point_3, int> result_type;

      Auto_count() : count(0)  { }
      std::pair<Point_3, int> operator()(const Point_3& p) const 
      { return std::make_pair(p, count++); }
      mutable int count;
    };

    CGAL_assertion(P.front() == P.back());
    CGAL_assertion(Q.empty() || (Q.front() == Q.back()));
    CGAL_assertion(Q.empty() || (P.size() == Q.size()));

    int n = P.size() - 1; // because the first and last point are equal

    Triangulation T(boost::make_transform_iterator(P.begin(), Auto_count()),
                    boost::make_transform_iterator(--P.end(), Auto_count()));
    T.infinite_vertex()->info() = -1;

    // Check whether all edges are included in DT, and get v0-vn-1 edge
    int nb_exists = 0;
    Finite_edges_iterator v0_vn_edge; // v0 vn-1 edge
    std::vector<bool> edge_exist(n, false);
    for(Finite_edges_iterator eb = T.finite_edges_begin(); eb != T.finite_edges_end(); ++eb) 
    {
      int v0_id = eb->first->vertex(eb->second)->info();
      int v1_id = eb->first->vertex(eb->third )->info();
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
      CGAL_TRACE_STREAM << "Not inside " << (n - nb_exists) << " of " << n << std::endl;
      CGAL_warning(!"Returning no output. Not all border edges are included in 3D Triangulation!");
      return Weight::NOT_VALID();
    }

    if(T.dimension() < 2) {
      CGAL_warning(!"Returning no output. Dimension of 3D Triangulation is below 2!");
      return Weight::NOT_VALID();
    }
    
    LookupTable<Weight> W(n, Weight::DEFAULT()); // do not forget that these default values are not changed for [i, i+1]
    LookupTable<int>    lambda(n,-1);

    if(T.dimension() == 3) {
      triangulate_DT<Incident_facet_circulator<3, Self> >
        (P, Q, W, lambda, *v0_vn_edge, T, WC);
    }
    else {
      CGAL_assertion(T.dimension() == 2);
      triangulate_DT<Incident_facet_circulator<2, Self> >
        (P, Q, W, lambda, *v0_vn_edge, T, WC);
    }
    
    if(lambda.get(0, n-1) == -1) {
      CGAL_warning(!"Returning no output. No possible triangulation is found!");
      return Weight::NOT_VALID();
    }

    tracer(lambda, 0, n-1);
    return W.get(0,n-1);
  }

private:
  // Finds other vertex than v0 and v1 in facet f
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
                  LookupTable<Weight>& W, 
                  LookupTable<int>& lambda, 
                  Edge e,
                  const Triangulation& T,
                  const WeightCalculator& WC)
  {
    /**********************************************************************
     *  + Default W value is Weight::DEFAULT(), default lambda value is -1.
     *  + DEFAULT() is used to check whether the region (v0-v1) is processed.
     *  + If a range v0-v1 does not contains any possible triangulation, then W[v0,v1] = NOT_VALID() and lambda[v0,v1] = -1
     *  + Note that w + DEFAULT() == w must hold
     */
    int v0 = e.first->vertex(e.second)->info();
    int v1 = e.first->vertex(e.third)->info();
    if(v0 > v1) { std::swap(v0, v1); }
    CGAL_assertion(v0 != -1); // edge can not be incident to infinite vertex

    if( W.get(v0, v1) != Weight::DEFAULT() || // the range is previously processed
        v0 + 1 == v1 ) // border edge - should not check v0 = 0, v1 = n-1, because it is the initial edge where the algorithm starts
    { return; }

    int m_min = -1;
    Weight w_min = Weight::NOT_VALID();

    IncidentFacetCirculator fb(e, T);
    do {
      int v2, v2_cell_index;
      boost::tie(v2, v2_cell_index) = get_facet_remaining_vertex(*fb, v0, v1);

      if(v2 < v0 || v2 > v1) { continue; } // this will also skip infinite vertex

      Edge e0 = Edge(fb->first, get_vertex_index(fb->first, v0) , v2_cell_index); // edge v0-v2
      Edge e1 = Edge(fb->first, get_vertex_index(fb->first, v1) , v2_cell_index); // edge v1-v2

      triangulate_DT<IncidentFacetCirculator>(P, Q, W, lambda, e0, T, WC); // region v0-v2
      triangulate_DT<IncidentFacetCirculator>(P, Q, W, lambda, e1, T, WC); // region v2-v1

      Weight w_021 = WC(P,Q, v0,v2,v1, lambda);
      if( W.get(v0, v2) == Weight::NOT_VALID() || W.get(v2, v1) == Weight::NOT_VALID() || w_021 == Weight::NOT_VALID())
      { continue; }

      Weight w = W.get(v0, v2) + W.get(v2, v1) + w_021;
      if(m_min == -1 || w < w_min){
        w_min = w;
        m_min = v2;
      }
    } while(++fb);

    // can be m_min = -1 and w_min = NOT_VALID which means no possible triangulation between v0-v1
    W.put(v0, v1, w_min);
    lambda.put(v0,v1, m_min);
  }
  
};

/************************************************************************/
/* Triangulate hole by using all search space
/************************************************************************/
template<
  class Kernel,
  class Tracer,
  class WeightCalculator,
  template <class> class LookupTable = Lookup_table
>
class Triangulate_hole_polyline {
public:
  typedef typename WeightCalculator::Weight  Weight;
  typedef typename Kernel::Point_3           Point_3;
  typedef std::vector<Point_3>               Polyline_3;

  Weight operator()(const Polyline_3& P, 
                    const Polyline_3& Q,
                    Tracer& tracer,
                    const WeightCalculator& WC)
  {
    CGAL_assertion(P.front() == P.back());
    CGAL_assertion(Q.empty() || (Q.front() == Q.back()));
    CGAL_assertion(Q.empty() || (P.size() == Q.size()));
    
    int n = P.size() - 1;                       // because the first and last point are equal
    LookupTable<Weight> W(n,Weight::DEFAULT()); // do not forget that these default values are not changed for [i, i+1]
    LookupTable<int>    lambda(n,-1);
    
    for(int j = 2; j< n; ++j) {   // determines range (2 - 3 - 4 )
      for(int i=0; i<n-j; ++i) {  // iterates over ranges and find min triangulation in those ranges 
        int k = i+j;              // like [0-2, 1-3, 2-4, ...], [0-3, 1-4, 2-5, ...]

        int m_min = -1;
        Weight w_min = Weight::NOT_VALID();

        for(int m = i+1; m<k; ++m) { 
          // now the regions i-m and m-k might be valid(constructed) patches,
          // if not, we can not construct i-m-k triangle
          if( W.get(i,m) == Weight::NOT_VALID() || W.get(m,k) == Weight::NOT_VALID() ) 
          { continue; }

          Weight w_imk = WC(P,Q,i,m,k, lambda);
          if(w_imk == Weight::NOT_VALID()) 
          { continue; }

          Weight w = W.get(i,m) + W.get(m,k) + w_imk;
          if(m_min == -1 || w < w_min) {
            w_min = w;
            m_min = m;
          }
        }

        // can be m_min = -1 and w_min = NOT_VALID which means no possible triangulation between i-k
        W.put(i,k,w_min);
        lambda.put(i,k, m_min);
      }
    }

    if(lambda.get(0,n-1) == -1) {
      CGAL_warning(!"Returning no output. No possible triangulation is found!");
      return Weight::NOT_VALID();
    }

    tracer(lambda, 0, n-1);
    return W.get(0,n-1);
  }
};

/***********************************************************************************/
/* Internal entry point for both polyline and Polyhedron_3 triangulation functions
/***********************************************************************************/
template <
  typename InputIterator,
  typename Tracer,
  typename WeightCalculator
>
typename WeightCalculator::Weight
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          InputIterator qbegin, InputIterator qend, 
                          Tracer& tracer,
                          const WeightCalculator& WC,
                          bool use_delaunay_triangulation) 
{
  typedef typename CGAL::Kernel_traits< typename std::iterator_traits<InputIterator>::value_type>::Kernel Kernel;
  typedef CGAL::internal::Triangulate_hole_polyline_DT<Kernel, Tracer, WeightCalculator> Fill_DT;
  typedef CGAL::internal::Triangulate_hole_polyline<Kernel, Tracer, WeightCalculator>    Fill;

  typename Fill::Polyline_3 P(pbegin, pend);
  typename Fill::Polyline_3 Q(qbegin, qend);
  if(P.front() != P.back()){
    P.push_back(P.front());
    if( !Q.empty() && P.size() > Q.size()) {
      Q.push_back(Q.front());
    }
  }

  typename WeightCalculator::Weight w = use_delaunay_triangulation ?
    Fill_DT().operator()(P,Q,tracer,WC) :
    Fill().operator()(P,Q,tracer,WC);
  CGAL_TRACE_STREAM << w << std::endl;
  return w;
}

} // namespace internal

/*!
\ingroup PkgHoleFilling
Creates triangles to fill the hole defined by points in the range (@a pbegin, @a pend). Triangles are put into @a out
using the indices of the input points in the range (@a pbegin, @a pend).
Note that no degenerate triangle is allowed during filling. If no possible patch is found, then no triangle is put into @a out.

The optional range (@a qbegin, @a qend) indicate for each pair of consecutive points in the range (@a pbegin, @a pend),
the third point of the facet this segment is incident to. 

Note that the range (@a pbegin, @a pend) and (@a qbegin, @a qend) may or may not contain duplicated first point at the end of sequence.

@tparam OutputIteratorValueType value type of OutputIterator having a constructor `OutputIteratorValueType(int p0, int p1, int p2)` available. 
        It is default to value_type_traits<OutputIterator>::type, and can be omitted when the default is fine
@tparam InputIterator iterator over input points
@tparam OutputIterator iterator over patch triangles
@param pbegin first iterator of the range of points
@param pend past-the-end iterator of the range of points
@param qbegin first iterator of the range of third points, can be omitted
@param qend past-the-end iterator of the range of third points, can be omitted
@param out iterator over output patch triangles
*/
template <typename OutputIteratorValueType, typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
                          InputIterator qbegin, InputIterator qend, 
                          OutputIterator out, bool use_delaunay_triangulation = false)
{
  typedef internal::Weight_calculator<
    internal::Weight_min_max_dihedral_and_area, 
    internal::Is_valid_degenerate_triangle
  > WC;

  internal::Tracer_polyline<OutputIteratorValueType, OutputIterator> tracer(out);
  internal::triangulate_hole_polyline(pbegin, pend, qbegin, qend, tracer, WC(), use_delaunay_triangulation);
  return tracer.out;
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

// overload no (qbegin, qend)
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
