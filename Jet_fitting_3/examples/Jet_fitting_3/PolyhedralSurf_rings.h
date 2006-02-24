#ifndef _PSURF_RINGS_H_
#define _PSURF_RINGS_H_

using namespace std;



template < class TPoly , class VertexPropertyMap> class T_PolyhedralSurf_rings
{
protected:
  //Polyhedron
  typedef typename TPoly::Vertex Vertex;
  typedef typename TPoly::Halfedge Halfedge;
  typedef typename TPoly::Facet Facet;
  typedef typename TPoly::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
  typedef typename TPoly::Vertex_iterator Vertex_iterator;
  
  //vertex property map
//   typedef typename boost::property_traits<VertexPropertyMap>::value_type vpm_value_type;
//   typedef typename boost::property_traits<VertexPropertyMap>::key_type vpm_key_type;
  
public:  

  //collect all neighbours up to ith ring
  static int
  push_neighbours_of(Vertex * start, int ith,
		     std::vector < Vertex * >&nextRing,
		     std::vector < Vertex * >&all,
		     VertexPropertyMap& vpm);

  //from currentRing, include neighbors of the next ring
  static int
  collect_ith_ring_neighbours(int ith,
			      std::vector < Vertex * >&currentRing,
			      std::vector < Vertex * >&nextRing,
			      std::vector < Vertex * >&all,
			      VertexPropertyMap& vpm);


  static void reset_ring_indices(std::vector < Vertex * >&vces,
				 VertexPropertyMap& vpm);

};


template < class TPoly , class VertexPropertyMap>
int T_PolyhedralSurf_rings <TPoly, VertexPropertyMap>::
push_neighbours_of(Vertex * start, int ith,
		   std::vector < Vertex * >&nextRing,
		   std::vector < Vertex * >&all,
		   VertexPropertyMap& vpm)
{
  int nb = 0;
  Vertex *v;
  Halfedge *h;
  
  Halfedge_around_vertex_circulator
    hedgeb = start->vertex_begin(), hedgee = hedgeb;
  
  do{
    h = &(*hedgeb);
    hedgeb++;
    v = &(*h->opposite()->vertex());
    if(get(vpm, v) != -1) continue;//if visited: next
    
    put(vpm, v, ith);
    nb++;
    nextRing.push_back(v);
    all.push_back(v);
  }
  while (hedgeb != hedgee);
  return nb;
}

template <class TPoly, class VertexPropertyMap>
int T_PolyhedralSurf_rings <TPoly, VertexPropertyMap>::
collect_ith_ring_neighbours(int ith, std::vector < Vertex * >&currentRing,
			    std::vector < Vertex * >&nextRing,
			    std::vector < Vertex * >&all,
			    VertexPropertyMap& vpm)
{
  typename std::vector < Vertex * >::iterator 
    itb =    currentRing.begin(), ite = currentRing.end();
  
  CGAL_For_all(itb, ite) push_neighbours_of(*itb, ith, nextRing, all, vpm);
  return nextRing.size();
}

template <class TPoly, class VertexPropertyMap>
void T_PolyhedralSurf_rings <TPoly, VertexPropertyMap>::
reset_ring_indices(std::vector < Vertex * >&vces,
		   VertexPropertyMap& vpm)
{
  typename std::vector < Vertex * >::iterator 
    itb = vces.begin(), ite = vces.end();
  CGAL_For_all(itb, ite)  put(vpm, *itb, -1);
}



#endif
