#ifndef CGAL_PSURF_RINGS_H_
#define CGAL_PSURF_RINGS_H_

#include <cassert>

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

  //vertex indices are initialised to -1
  static void reset_ring_indices(std::vector < Vertex * >&vces,
				 VertexPropertyMap& vpm);

  //i >= 1; from a start vertex on the current i-1 ring, push non-visited neighbors
  //of start in the nextRing and set indices to i. Also add these vertices in all.
  static void push_neighbours_of(Vertex * start, int ith,
			  std::vector < Vertex * >&nextRing,
			  std::vector < Vertex * >&all,
			  VertexPropertyMap& vpm);

  //i >= 1, from a currentRing i-1, collect all neighbors, set indices
  //to i and store them in nextRing and all.
  static void collect_ith_ring(int ith,
			std::vector < Vertex * >&currentRing,
			std::vector < Vertex * >&nextRing,
			std::vector < Vertex * >&all,
			VertexPropertyMap& vpm);

 public:
  //collect i>=1 rings : all neighbours up to the ith ring,
  static void
    collect_i_rings(Vertex* v,
		    int ring_i,
		    std::vector < Vertex * >& all,
		    VertexPropertyMap& vpm);

  //collect enough rings (at least 1), to get at least min_nb of neighbors
  static void
    collect_enough_rings(Vertex* v,
			 unsigned int min_nb,
			 std::vector < Vertex * >& all,
			 VertexPropertyMap& vpm);
};

////IMPLEMENTATION/////////////////////////////////////////////////////////////////////

template < class TPoly , class VertexPropertyMap>
void T_PolyhedralSurf_rings <TPoly, VertexPropertyMap>::
push_neighbours_of(Vertex * start, int ith,
		   std::vector < Vertex * >&nextRing,
		   std::vector < Vertex * >&all,
		   VertexPropertyMap& vpm)
{
  Vertex *v;
  Halfedge_around_vertex_circulator
    hedgeb = start->vertex_begin(), hedgee = hedgeb;

 CGAL_For_all(hedgeb, hedgee)
  {
    v = &*(hedgeb->opposite()->vertex());
    if(get(vpm, v) != -1) continue;//if visited: next

    put(vpm, v, ith);
    nextRing.push_back(v);
    all.push_back(v);
  }
}

template <class TPoly, class VertexPropertyMap>
void T_PolyhedralSurf_rings <TPoly, VertexPropertyMap>::
collect_ith_ring(int ith, std::vector < Vertex * >&currentRing,
		 std::vector < Vertex * >&nextRing,
		 std::vector < Vertex * >&all,
		 VertexPropertyMap& vpm)
{
  typename std::vector < Vertex * >::iterator
    itb =    currentRing.begin(), ite = currentRing.end();

  CGAL_For_all(itb, ite) push_neighbours_of(*itb, ith, nextRing, all, vpm);
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

template <class TPoly, class VertexPropertyMap>
  void T_PolyhedralSurf_rings <TPoly, VertexPropertyMap>::
collect_i_rings(Vertex* v,
		int ring_i,
		std::vector < Vertex * >& all,
		VertexPropertyMap& vpm)
{
  std::vector<Vertex*> current_ring, next_ring;
  std::vector<Vertex*> *p_current_ring, *p_next_ring;
  assert(ring_i >= 1);
  //initialize
  p_current_ring = &current_ring;
  p_next_ring = &next_ring;
  put(vpm, v, 0);
  current_ring.push_back(v);
  all.push_back(v);

  for (int i=1; i<=ring_i; i++)
    {
      collect_ith_ring(i, *p_current_ring, *p_next_ring, all, vpm);
      //next round must be launched from p_nextRing...
      p_current_ring->clear();
      std::swap(p_current_ring, p_next_ring);
    }
  //clean up
  reset_ring_indices(all, vpm);
}

template <class TPoly, class VertexPropertyMap>
  void T_PolyhedralSurf_rings <TPoly, VertexPropertyMap>::
collect_enough_rings(Vertex* v,
		     unsigned int min_nb,
		     std::vector < Vertex * >& all,
		     VertexPropertyMap& vpm)
{
  std::vector<Vertex*> current_ring, next_ring;
  std::vector<Vertex*> *p_current_ring, *p_next_ring;

  //initialize
  p_current_ring = &current_ring;
  p_next_ring = &next_ring;
  put(vpm, v, 0);
  current_ring.push_back(v);
  all.push_back(v);

  int i = 1;

  while ( (all.size() < min_nb) &&  (p_current_ring->size() != 0) )
    {
      collect_ith_ring(i, *p_current_ring, *p_next_ring, all, vpm);
      //next round must be launched from p_nextRing...
      p_current_ring->clear();
      std::swap(p_current_ring, p_next_ring);
      i++;
    }
  //clean up
  reset_ring_indices(all, vpm);
}


#endif
