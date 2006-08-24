#ifndef _POLYHEDRALSURF_NEIGHBORS_H_
#define _POLYHEDRALSURF_NEIGHBORS_H_

#include <queue>
#include <algorithm>
#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

//---------------------------------------------------------------------------
//T_Gate : element of the priority queue. A gate is a halfedge and a
//number giving the max distance from v to the vertices of the
//triangle incident to the halfedge.
//---------------------------------------------------------------------------
template < class TPoly > class T_Gate
{  
public:
  typedef typename TPoly::Traits::FT       FT;
  typedef typename TPoly::Traits::Vector_3 Vector_3;
  typedef typename TPoly::Traits::Point_3  Point_3;
  typedef typename TPoly::Vertex_handle    Vertex_handle;
  typedef typename TPoly::Halfedge_handle  Halfedge_handle;
 
  T_Gate( Vertex_handle v, Halfedge_handle he);
  FT& d() { return m_d;}
  const FT d() const { return m_d;}            
  Halfedge_handle he() { return m_he;}

private:
  FT m_d;
  Halfedge_handle m_he;
};

//////////////IMPLEMENTATION//////////////////////////
template < class TPoly > 
T_Gate<TPoly>::T_Gate( Vertex_handle v, Halfedge_handle he)
  : m_he(he)
{
  Point_3 p0 = v->point(),
    p1 = he->vertex()->point(),
    p2 = he->next()->vertex()->point(),
    p3 = he->prev()->vertex()->point();
  Vector_3 p0p1 = p0 - p1,
    p0p2 = p0 - p2,
    p0p3 = p0 - p3;
  FT d1 = p0p1*p0p1,
    d2 = p0p2*p0p2,
    d3 = p0p3*p0p3;
  m_d = CGAL::sqrt( (std::max)( (std::max)(d1,d2), d3) );
}

//---------------------------------------------------------------------------
// functor for priority queue
// order so than the top element is the smallest in the queue
//---------------------------------------------------------------------------
template<class g>
struct compare_gates 
{       
        bool operator()(const g& g1, 
                        const g& g2) const
        {       
                return g1.d() > g2.d();
        }
};

//---------------------------------------------------------------------------
//T_PolyhedralSurf_neighbors : MAIN class for computation, it uses the
//class Gate and the functor compare_gates for the definition of a
//priority queue
//---------------------------------------------------------------------------
template < class TPoly > class T_PolyhedralSurf_neighbors
{
public:
  typedef typename TPoly::Traits::FT        FT;
  typedef typename TPoly::Traits::Vector_3  Vector_3;
  typedef typename TPoly::Traits::Point_3   Point_3;
  typedef typename TPoly::Vertex_handle     Vertex_handle;
  typedef typename TPoly::Halfedge_handle   Halfedge_handle;
  typedef typename TPoly::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;
  typedef typename TPoly::Vertex_iterator   Vertex_iterator;
  typedef T_Gate<TPoly> Gate;

  T_PolyhedralSurf_neighbors(TPoly& P);
  // vertex_neigh stores the vertex v and its 1Ring neighbors contour
  // stores halfedges, oriented CW, following the 1Ring disk border
  // OneRingSize is the max distance from v to its OneRing
  // neighbors. (the tag is_visited is not mofified)
  void compute_one_ring(Vertex_handle v,
			std::vector<Vertex_handle> &vertex_neigh,
			std::list<Halfedge_handle> &contour,
			FT &OneRingSize);
  // call compute_one_ring and expand the contour (circle of halfedges
  // CW), vertex_neigh are vertices on and inside the contour (there
  // tag is_visited is set to true, but reset to false at the end),
  // size is such that gates with distance less than size*OneRingSize
  // are processed
  void compute_neighbors(Vertex_handle v,
			 std::vector<Vertex_handle> &vertex_neigh,
			 std::list<Halfedge_handle> &contour,
			 FT size); 
  //vertex tags is_visited are set to false
  void reset_is_visited_map(std::vector<Vertex_handle> &vces);

 protected:
  //tag to visit vertices
  struct Vertex_cmp{//comparison is wrt vertex addresses
    bool operator()(Vertex_handle a,  Vertex_handle b) const{
      return &*a < &*b;
    }
  };
  typedef std::map<Vertex_handle, bool, Vertex_cmp> Vertex2bool_map_type;
  Vertex2bool_map_type is_visited_map;
};

//////////////IMPLEMENTATION//////////////////////////
template < class TPoly >
T_PolyhedralSurf_neighbors < TPoly >::
T_PolyhedralSurf_neighbors(TPoly& P)
{
  //init the is_visited_map
  Vertex_iterator itb = P.vertices_begin(), ite = P.vertices_end();
  for(;itb!=ite;itb++) is_visited_map[itb] = false; 
}

template < class TPoly >
void T_PolyhedralSurf_neighbors < TPoly >::
compute_one_ring(Vertex_handle v,
		      std::vector<Vertex_handle> &vertex_neigh,
		      std::list<Halfedge_handle> &contour,
		      FT &OneRingSize)
{
  typedef typename std::list<Halfedge_handle>::iterator list_it;
  vertex_neigh.push_back(v);
  Halfedge_around_vertex_circulator he_circ = v->vertex_begin(), 
                                    he_end = he_circ;
  do {
      if ( he_circ->is_border() )//then he and he->next follow the contour CW
	{contour.push_back(he_circ);
	contour.push_back(he_circ->next());}
      else contour.push_back(he_circ->prev()->opposite());//not border, he->prev->opp on contour CW
      vertex_neigh.push_back(he_circ->opposite()->vertex());
      he_circ++;
  } while (he_circ != he_end);

  //compute OneRingSize = distance(v, 1Ring)
  OneRingSize = 0;
  typename std::vector<Vertex_handle>::iterator itb = vertex_neigh.begin(),
    ite = vertex_neigh.end();
  itb++;//the first vertex v is the center to which distances are
	//computed from, for other 1ring neighbors
  Point_3 p0 = v->point(), p;
  Vector_3 p0p;
  FT d = OneRingSize;
  for (; itb != ite; itb++){

    p = (*itb)->point();
    p0p = p0 - p;
    d =  CGAL::sqrt(p0p*p0p);
    if (d > OneRingSize) OneRingSize = d;
  }
}

template < class TPoly >
void T_PolyhedralSurf_neighbors < TPoly >::
compute_neighbors(Vertex_handle v,
		   std::vector<Vertex_handle> &vertex_neigh,
		   std::list<Halfedge_handle> &contour,
		   FT size)  
{
  FT OneRingSize;
  compute_one_ring(v, vertex_neigh, contour, OneRingSize);
  FT d_max = OneRingSize*size;
  std::priority_queue< Gate, std::vector< Gate >, compare_gates< Gate > > GatePQ;
  // tag neighbors 
  typename std::vector<Vertex_handle>::iterator itbv = vertex_neigh.begin(),
    itev = vertex_neigh.end();
  for (; itbv != itev; itbv++) is_visited_map.find(*itbv)->second = true;

  // init GatePQ
  typename std::list<Halfedge_handle>::iterator itb = contour.begin(),
                                       ite = contour.end();
  for (; itb != ite; itb++) {
    if (!( (*itb)->is_border() )) GatePQ.push(Gate(v, *itb));
  }
  // init d_current
  Gate firstGate = GatePQ.top();
  FT d_current = firstGate.d();
  // main loop
  while ( !GatePQ.empty() && d_current <= d_max ) {
    Gate gate = GatePQ.top();
    GatePQ.pop();
    d_current = gate.d();
    Halfedge_handle he = gate.he(), he1, he2;
    Vertex_handle v1;
    // find the gate on the contour
    typename std::list<Halfedge_handle>::iterator pos_he, pos_prev, pos_next, iter;
   
    pos_he = find(contour.begin(), contour.end(), he);
    iter = pos_he;
    /**
       there are different cases to expand the contour : 
       (case 3) he is not on the contour, nothing to do
       (case 2) he is on the contour and either the previous or the next 
       following edge in the triangle is also on the contour, then delete
       these 2 he from the contour and add the third one to the contour 
       and the PQ.
       (case1) the vertex opposite to he is not visited, then the he is removed
       from the contour, the two others are added to the contour and PQ, the 
       vertex is set visited.
    */

    // if the gate is not encountered on the contour (case 3)
    if ( pos_he == contour.end() ) continue;
    // simulate a circulator on the contour: 
    // find the prev and next pos on coutour
    if ( (++iter) != ite ) pos_next = iter;
    else pos_next = contour.begin();
    iter = pos_he;
    if ( iter != contour.begin() ) pos_prev = --iter;
    else pos_prev = --contour.end();

    if ( he->next() == *pos_next )
      {  // case 2a
	//contour
	he1 = he->prev()->opposite();
	contour.insert(pos_he, he1);
	contour.erase(pos_he);
	contour.erase(pos_next);
	//GatePQ
	if ( !(he1->is_border()) ) GatePQ.push(Gate(v, he1));
	continue;
      }
    else if ( he->prev() == (*pos_prev) )
      {  // case 2b
	//contour
	he1 = he->next()->opposite();
	contour.insert(pos_prev, he1);
	contour.erase(pos_prev);
	contour.erase(pos_he);
	//GatePQ
	if ( !(he1->is_border()) ) GatePQ.push(Gate(v, he1));
	continue;
      }
    v1 = he->next()->vertex();
    if ( !is_visited_map.find(v1)->second )
      {  // case 1
	//vertex
	is_visited_map.find(v1)->second = true;
	vertex_neigh.push_back(v1);
	//contour
	he1 = he->prev()->opposite();
	he2 = he->next()->opposite();
	contour.insert(pos_he, he1);
	contour.insert(pos_he, he2);
	contour.erase(pos_he);
	//GatePQ
	if ( !(he1->is_border()) ) GatePQ.push(Gate(v, he1));
	if ( !(he2->is_border()) ) GatePQ.push(Gate(v, he2));
	continue;
      }
    //else case non admissible
    CGAL_postcondition ( "case non admissible" == 0);
  }// end while
  
  reset_is_visited_map(vertex_neigh);
}

template < class TPoly >
void T_PolyhedralSurf_neighbors < TPoly >::
reset_is_visited_map(std::vector<Vertex_handle> &vces)
{
  typename std::vector<Vertex_handle>::iterator 
    itb = vces.begin(), ite = vces.end();
  for (;itb != ite; itb++) is_visited_map[*itb] = false;
}

CGAL_END_NAMESPACE

#endif
