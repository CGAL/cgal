#ifndef _PSURF_OPERATIONS_H_
#define _PSURF_OPERATIONS_H_

using namespace std;


//----------------------------------------------------------------
// some functors 
//----------------------------------------------------------------
struct Edge_length {
  template <class HalfEdge> 
  double operator() (HalfEdge h) {
    return CGAL::sqrt(CGAL::squared_distance(h.prev()->vertex()->point(),
					     h.vertex()->point()));
  }
};


//----------------------------------------------------------------
// operations on hedges, facets etc, handled using proeprty maps
//----------------------------------------------------------------

template <class TPoly , class HEdgePropertyMap> class T_PolyhedralSurf_ops
{
protected:
  //Polyhedron
  typedef typename TPoly::Vertex Vertex;
  typedef typename TPoly::Halfedge Halfedge;
  typedef typename TPoly::Facet Facet;
  
  typedef typename TPoly::Vertex_iterator Vertex_iterator;
  typedef typename TPoly::Halfedge_iterator Halfedge_iterator;
  typedef typename TPoly::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
  
public:  
  
  static void compute_edges_length(TPoly& P, HEdgePropertyMap& hepm);
  static double compute_mean_edges_length_around_vertex(Vertex * v, HEdgePropertyMap& hepm);
};

template <class TPoly , class HEdgePropertyMap>
void T_PolyhedralSurf_ops<TPoly, HEdgePropertyMap>::
compute_edges_length(TPoly& P, HEdgePropertyMap& hepm)
{
 //  std::for_each(this->halfedges_begin(), this->halfedges_end(),Edge_length());
  Halfedge_iterator itb = P.halfedges_begin(), ite = P.halfedges_end();
  for(; itb!=ite; itb++) {
    Halfedge h=*itb;
    //BUG??? put(hepm, h, Edge_length()(h));
  }
}

template <class TPoly , class HEdgePropertyMap>
double T_PolyhedralSurf_ops<TPoly, HEdgePropertyMap>::
compute_mean_edges_length_around_vertex(Vertex * v, HEdgePropertyMap& hepm)
{
  Halfedge_around_vertex_circulator  hedgeb = v->vertex_begin(), hedgee = hedgeb;
  int  count_he = 0;
  double sum = 0.;
  CGAL_For_all(hedgeb, hedgee){
    Halfedge h = *hedgeb;
    double d = get(hepm, h);
    sum += d;
    count_he++;
  }
  return sum/count_he;
}



#endif
