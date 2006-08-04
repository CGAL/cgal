#ifndef _PSURF_RINGS_H_
#define _PSURF_RINGS_H_

using namespace std;

//---------------------------------------------------------------------------
//T_PolyhedralSurf_rings
//---------------------------------------------------------------------------
template < class TPoly > class T_PolyhedralSurf_rings
{
protected:
  //Polyhedron
  typedef typename TPoly::Vertex Vertex;
  typedef typename TPoly::Halfedge Halfedge;
  typedef typename TPoly::Facet Facet;
  typedef typename TPoly::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
  typedef typename TPoly::Vertex_iterator Vertex_iterator;
  typedef std::map<Vertex*, int> Vertex2int_map_type;
  Vertex2int_map_type ring_index_map;
  
  //vertex indices are initialised to -1
  void reset_ring_indices(std::vector < Vertex * >&vces);

  //i >= 1; from a start vertex on the current i-1 ring, push non-visited neighbors
  //of start in the nextRing and set indices to i. Also add these vertices in all.
  void push_neighbours_of(Vertex * start, int ith,
			  std::vector < Vertex * >&nextRing,
			  std::vector < Vertex * >&all);

  //i >= 1, from a currentRing i-1, collect all neighbors, set indices
  //to i and store them in nextRing and all.
  void collect_ith_ring(int ith,
			std::vector < Vertex * >&currentRing,
			std::vector < Vertex * >&nextRing,
			std::vector < Vertex * >&all);
 
 public:  
  T_PolyhedralSurf_rings(TPoly& P);
  
  //collect i>=1 rings : all neighbours up to the ith ring,
  void collect_i_rings(Vertex* v, 
		       int ring_i, 
		       std::vector < Vertex * >& all);

  //collect enough rings (at least 1), to get at least min_nb of neighbors
  void collect_enough_rings(Vertex* v,
			    unsigned int min_nb,
			    std::vector < Vertex * >& all);
};

////IMPLEMENTATION/////////////////////////////////////////////////////////////////////

template < class TPoly >
T_PolyhedralSurf_rings <TPoly>::
T_PolyhedralSurf_rings(TPoly& P)
{
  //init the ring_index_map
  Vertex_iterator itb = P.vertices_begin(), ite = P.vertices_end();
  for(;itb!=ite;itb++) ring_index_map[&*itb] = -1; 
}

template < class TPoly >
void T_PolyhedralSurf_rings <TPoly>::
push_neighbours_of(Vertex * start, int ith,
		   std::vector < Vertex * >&nextRing,
		   std::vector < Vertex * >&all)
{
  Vertex *v;  
  Halfedge_around_vertex_circulator
    hedgeb = start->vertex_begin(), hedgee = hedgeb;
 
 CGAL_For_all(hedgeb, hedgee)
  {
    v = &*(hedgeb->opposite()->vertex());
    if (ring_index_map[v] != -1)  continue;//if visited: next
    
    ring_index_map[v] = ith;
    nextRing.push_back(v);
    all.push_back(v);
  }
}

template <class TPoly>
void T_PolyhedralSurf_rings <TPoly>::
collect_ith_ring(int ith, std::vector < Vertex * >&currentRing,
		 std::vector < Vertex * >&nextRing,
		 std::vector < Vertex * >&all)
{
  typename std::vector < Vertex * >::iterator 
    itb =    currentRing.begin(), ite = currentRing.end();
  
  CGAL_For_all(itb, ite) push_neighbours_of(*itb, ith, nextRing, all);
}

template <class TPoly>
  void T_PolyhedralSurf_rings <TPoly>::
reset_ring_indices(std::vector < Vertex * >&vces)
{
  typename std::vector < Vertex * >::iterator 
    itb = vces.begin(), ite = vces.end();
  CGAL_For_all(itb, ite)  ring_index_map[*itb] = -1;
}
 
template <class TPoly>
  void T_PolyhedralSurf_rings <TPoly>::
collect_i_rings(Vertex* v, 
		int ring_i, 
		std::vector < Vertex * >& all)
{
  std::vector<Vertex*> current_ring, next_ring; 
  std::vector<Vertex*> *p_current_ring, *p_next_ring;
  assert(ring_i >= 1);
  //initialize
  p_current_ring = &current_ring;
  p_next_ring = &next_ring;
  ring_index_map[v] = 0;
  current_ring.push_back(v);
  all.push_back(v);

  for (int i=1; i<=ring_i; i++)
    {
      collect_ith_ring(i, *p_current_ring, *p_next_ring, all);
      //next round must be launched from p_nextRing...
      p_current_ring->clear();
      std::swap(p_current_ring, p_next_ring);
    }
  //clean up
  reset_ring_indices(all);
}

template <class TPoly>
  void T_PolyhedralSurf_rings <TPoly>::
collect_enough_rings(Vertex* v,
		     unsigned int min_nb,
		     std::vector < Vertex * >& all)
{
  std::vector<Vertex*> current_ring, next_ring; 
  std::vector<Vertex*> *p_current_ring, *p_next_ring;

  //initialize
  p_current_ring = &current_ring;
  p_next_ring = &next_ring;
  ring_index_map[v] = 0;
  current_ring.push_back(v);
  all.push_back(v);

  int i = 1;

  while ( (all.size() < min_nb) &&  (p_current_ring->size() != 0) ) 
    {
      collect_ith_ring(i, *p_current_ring, *p_next_ring, all);
      //next round must be launched from p_nextRing...
      p_current_ring->clear();
      std::swap(p_current_ring, p_next_ring);
      i++;
    }
  //clean up
  reset_ring_indices(all);
}


#endif
