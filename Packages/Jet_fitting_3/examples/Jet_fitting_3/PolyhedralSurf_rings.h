#ifndef _SURPOLYRINGS_H_
#define _SURPOLYRINGS_H_

using namespace std;

//--------------------------------------------------
//we store hedges. exple: 1ring neighbours of 
//v=hedges anchored at v
//--------------------------------------------------
// when storing hedges to a ith ring: do we store hedges to the ring or
// opposite hedges?
enum HTO_HOPPOSITE_ONERING { HTO_ONERING, HOPPOSITE_ONERING };

template < class TPoly > class T_PolyhedralSurf_rings
{
public:
  typedef typename TPoly::Point_3 Point_3;
  typedef typename TPoly::Vertex Vertex;
  typedef typename TPoly::Halfedge Halfedge;
  typedef typename TPoly::Facet Facet;
  typedef typename TPoly::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;

  //debug
  static void check_vertex_indices(Vertex* v); 
  static void check_ring_indices(TPoly& poly);

  //  static void from_into(ptrSet & a, ptrSet & b, char erase_a = 0);


  //---------------VERTEX---------------------------------
  //neighbours using a tag on each vertex ---default is -1, vertex is 0, 
  //  then 1,2,.. for ith rings
  //---------------------------------------------------------
  static int
  push_neighbours_of(Vertex * start, int ith,
		     std::vector < Vertex * >&nextRing,
		     std::vector < Vertex * >&all);

  static int
  collect_ith_ring_neighbours(int ith,
			      std::vector < Vertex * >&currentRing,
			      std::vector < Vertex * >&nextRing,
			      std::vector < Vertex * >&all);


  static void reset_ring_indices(std::vector < Vertex * >&vces);
  static void reset_ring_indices(std::vector< Vertex_handle >&vces);
  //Neighbours collected in ccw order on each ring
  static int collect_ith_ring_ccw(int ith,
				  vector<int>& sizes,
				  vector<Vertex*>& vces);
  static int i_rings_ccw(int i, 
			 Vertex* start, 
			 vector<int>& sizes, 
			 vector<Vertex*>& vces);


  //------------------HEDGES--------------------------------
  //hdges incident to endpoint of h
  //we store opposite hedges, ie those
  //to the next ith ring neighbours
  //---------------------------------------------------------
  static int
  push_hedges_of(Vertex * start, int ith,
		 HTO_HOPPOSITE_ONERING to_opp,
		 std::vector < Halfedge * >&nextRing,
		 std::vector < Halfedge * >&recordedOnce);
	//, char& multiple);

  //grabbing the next ring from currentRing
  static int
  next_ring_hedges(int ith, HTO_HOPPOSITE_ONERING to_opp,
		   std::vector < Halfedge * >&currentRing,
		   std::vector < Halfedge * >&nextRing, 
		   std::vector < Halfedge * >&recordedOnce);
	//, char& multiple);

  //ith ring from scratch ie the center vertex
  static int
  collect_ith_ring_hedges(Vertex * start, int ith,
			  HTO_HOPPOSITE_ONERING to_opp,
			  std::vector < Halfedge * >&recordedOnce,
			  char cleanup);

  //Collect contour edges of the ith ring, these edges oriented CCW
  //  i.e. their incident facet is inside.
  static void
  contour_ith_ring_hedges(Vertex * start, int ith,
			  std::vector < Halfedge * > &contour);

  static int
  contour_ccw_ith_ring_hedges(Vertex * start, int ith,
			      std::vector < Halfedge * > &contourCCW);

  static void reset_ring_indices(std::vector < Halfedge * >&hedges);

  static void reset_ring_indices(std::vector < Facet * >&facets);

  // to report 1st crossing edge with a sphere
  static Point_3 orig_point;
  static void set_orig_point(Point_3 & p);
  static void
  getFarthestNeighbourOfV_ofOrigPoint(Vertex * v,
				      Halfedge * &h_farthest);

};

//--------------------------------------------------
//-- collecting ith ring neighbours
//--------------------------------------------------

template < class TPoly >
void T_PolyhedralSurf_rings <TPoly >::
check_vertex_indices(Vertex* v)
{
  Halfedge_around_vertex_circulator itb = v->vertex_begin(), 
    ite = v->vertex_begin(); 
  assert( ((&(*itb))->opposite()->vertex()->getRingIndex()==-1) );
  itb++;
  for( ; itb != ite ; itb++)
    assert(((&(*itb))->opposite()->vertex()->getRingIndex()==-1) );
} 


template < class TPoly >
void T_PolyhedralSurf_rings < TPoly >::
check_ring_indices(TPoly& poly)
{
  Vertex_iterator itb = poly.vertices_begin (),
    ite = poly.vertices_end();
  for( ; itb != ite ; itb++)
    T_PolyhedralSurf_rings < TPoly >::check_vertex_indices( &(*itb));
}

// //neighbours using a map to make sure we do not collect twice
// //---------------------------------------------------------
// template < class TPoly >
// void T_PolyhedralSurf_rings < TPoly >::
// from_into(ptrSet & a, ptrSet & b, char erase_a)
// {
//   ptrSet_iterator itb = a.begin(), ite = a.end();
//   for (; itb != ite; itb++)
//     b.insert(*itb);

//   if (erase_a)
//     a.erase(a.begin(), ite);
// }

//neighbours using a tag on each vertex ---default is -1, vertex is 0, 
//  then 1,2,.. for ith rings
//-------------------------------------------------------------------

template < class TPoly >
int T_PolyhedralSurf_rings < TPoly >::
push_neighbours_of(Vertex * start, int ith,
		   std::vector < Vertex * >&nextRing,
		   std::vector < Vertex * >&all)
{
  int nb = 0;
  Vertex *v;
  Halfedge *h;

  Halfedge_around_vertex_circulator
    hedgeb = start->vertex_begin(), hedgee = hedgeb;

  do
    {
      h = &(*hedgeb);
      hedgeb++;
      v = &(*h->opposite()->vertex());

      //make sure vertex not visited yet
      if (v->getRingIndex() != -1)
	continue;
      //tag the vertex as belonging to the its ring   
      v->setRingIndex(ith);
      nb++;
      nextRing.push_back(v);
      all.push_back(v);
    }
  while (hedgeb != hedgee);
  return nb;
}

template < class TPoly >
int T_PolyhedralSurf_rings < TPoly >::
collect_ith_ring_neighbours(int ith, std::vector < Vertex * >&currentRing,
			    std::vector < Vertex * >&nextRing,
			    std::vector < Vertex * >&all)
{
  Vertex *v;
  typename std::vector < Vertex * >::iterator itb =
    currentRing.begin(), ite = currentRing.end();

  for (; itb != ite; itb++)
    {
      v = *itb;
      T_PolyhedralSurf_rings < TPoly >::
	push_neighbours_of(v, ith, nextRing, all);
    }
  return nextRing.size();
}

template < class TPoly >
void T_PolyhedralSurf_rings < TPoly >::
reset_ring_indices(std::vector < Vertex * >&vces)
{
  typename std::vector < Vertex * >::iterator itb = vces.begin(), ite =
    vces.end();
  for (; itb != ite; itb++)
    (*itb)->resetRingIndex();
}

//surcharge pour les handles
template < class TPoly >
void T_PolyhedralSurf_rings < TPoly >::
reset_ring_indices(std::vector < Vertex_handle >&vces)
{
  typename std::vector < Vertex_handle >::iterator itb = vces.begin(), ite =
    vces.end();
  for (; itb != ite; itb++)
    (&(*(*itb)))->resetRingIndex();
}

//see also the other method below--------------
//Neighbours collected in ccw order on each ring
//if a border edge is encountered, return 0
//sizes[j] = nb of vertives of the jth ring
// template <class TPoly>
// int T_PolyhedralSurf_rings<TPoly>::
// collect_ith_ring_ccw(int ith,
// 		     vector<int>& sizes,
// 		     vector<Vertex*>& vces)
// {
//   //find a vertex of the ith ring from a vertex of the (i-1)th ring
//   int size=0;
//   Vertex* v=vces[vces.size()-1];
//   Halfedge* h;
//   Halfedge_around_vertex_circulator
// 	hedgeb = v->vertex_begin(), hedgee = hedgeb;
//   h = &(*hedgeb);
//   if (h->is_border_edge()) return 0;
//   while (h->opposite()->vertex()->getRingIndex() != -1)
//     {
//      hedgeb++;
//      h = &(*hedgeb);
//      if (h->is_border_edge()) return 0;
//     }
//   v = &(*(h->opposite()->vertex()));
//   vces.push_back(v);
//   size++;
//   v->setRingIndex(ith);
//   h=&(*h->opposite()->next());
//   if (h->is_border_edge()) return 0;
//   v=&(*h->vertex());
//   //follow ccw the ith ring
//   while (v->getRingIndex() != ith)
//     {
//       while ( (v->getRingIndex() != -1) && (v->getRingIndex() != ith) ) 
// 	{
// 	  h=&(*h->opposite()->next());
// 	  if (h->is_border_edge()) return 0;
// 	  v=&(*h->vertex());
// 	}
//       if (v->getRingIndex() == -1)
// 	{
// 	  vces.push_back(v);
// 	  size++;
// 	  v->setRingIndex(ith);
// 	  h=&(*h->next());
// 	  if (h->is_border_edge()) return 0;
// 	  v=&(*h->vertex());
// 	}
//     }
//   sizes.push_back(size);
//   return 1;
// }
//Other method using contour_ccw_ith_ring_hedges--------------------------
template <class TPoly>
int T_PolyhedralSurf_rings<TPoly>::
collect_ith_ring_ccw(int ith,
		     vector<int>& sizes,
		     vector<Vertex*>& vces)
{
  Vertex* v=vces[0];
  std::vector < Halfedge * > contourCCW;
  //check if no border edge is encountered and set contourCCW of halfedges
  if (T_PolyhedralSurf_rings<TPoly>::
      contour_ccw_ith_ring_hedges(v, ith, contourCCW) == 0) return 0;

  typename std::vector < Halfedge * >::iterator itb =
    contourCCW.begin(), ite = contourCCW.end();

  //retrieve vertices from hedge contour
  for (; itb != ite; itb++)  vces.push_back(&(*((*itb)->vertex())));
  sizes.push_back(contourCCW.size());
  return 1;
}


template <class TPoly>
int T_PolyhedralSurf_rings<TPoly>::
i_rings_ccw(int i, 
	    Vertex* start, 
	    vector<int>& sizes, 
	    vector<Vertex*>& vces)
{
  //set the 0-ring
  start->setRingIndex(0);
  sizes.push_back(1);
  vces.push_back(start);
  int ith;
  //set ith ring, ith from 1 to i
  for(ith=1;ith<=i;ith++) 
    if (T_PolyhedralSurf_rings<TPoly>::
	collect_ith_ring_ccw(ith, sizes, vces) == 0) 
      {
	T_PolyhedralSurf_rings<TPoly>::reset_ring_indices(vces);
	return 0;
      }
  T_PolyhedralSurf_rings<TPoly>::reset_ring_indices(vces);
  return 1;
}

//-----------------------------------------------------
//---------- hedges. 
//-----------------------------------------------------
template < class TPoly > 
int T_PolyhedralSurf_rings < TPoly >:: 
push_hedges_of(Vertex * start, int ith, HTO_HOPPOSITE_ONERING to_opp,
	       std::vector < Halfedge * >&nextRing,
	       std::vector < Halfedge * >&recordedOnce)
{
  // char ok;
  int n = 0;
  Halfedge *h, *hopp, *hpushed;
  Halfedge_around_vertex_circulator
    hedgeb = start->vertex_begin(), hedgee = hedgeb;
  do
    {
      h = &(*hedgeb);
      hopp = &(*h->opposite());
      hedgeb++;
      //store either h or hopp
      hpushed = (to_opp == HTO_ONERING) ? hopp : h;

      //in any case
      nextRing.push_back(hpushed);

      //record or not?
      if (hpushed->getRingIndex() == -1)
	{
	  //tag the vertex as belonging to the its ring   
	  //marc: isn't it: tag the edge and its opposite?
	  h->setRingIndex(ith);
	  hopp->setRingIndex(ith);
	  n++;		//xfc
	  recordedOnce.push_back(hpushed);
	}
      //hedge already recorded
      //else multiple=1;
    }
  while (hedgeb != hedgee);
  return n;
}

//notice that for each helfedge, we need to recover the endpoint, ie
//the vertex on the current ring being explored
template < class TPoly >
 int T_PolyhedralSurf_rings < TPoly >:: 
next_ring_hedges(int ith, HTO_HOPPOSITE_ONERING to_opp,
		 std::vector < Halfedge * >&currentRing, 
		 std::vector < Halfedge * >&nextRing, 
		 std::vector < Halfedge * >&recordedOnce)
{
  int s = 0;
  Vertex *v;
  Halfedge *h;
  typename std::vector < Halfedge * >::iterator itb =
    currentRing.begin(), ite = currentRing.end();
  for (; itb != ite; itb++)
    {
      h = (Halfedge *) (*itb);

      //grab vertex 
      if (to_opp == HTO_ONERING)
	v = &*(h->vertex());
      else
	v = &*(h->opposite()->vertex());

      s += T_PolyhedralSurf_rings < TPoly >::
	push_hedges_of(v, ith, to_opp, nextRing, recordedOnce);	
    }
  return s;
}

//Warning collect hedges from rings 1 to ith, better to call it
//collect_i_rings_hedges ???????
template < class TPoly >
int T_PolyhedralSurf_rings < TPoly >::
collect_ith_ring_hedges(Vertex * start, int ith,
			HTO_HOPPOSITE_ONERING to_opp,
			std::vector < Halfedge * >&recordedOnce,
			char cleanup)
{
  // char multiple=0;
  int i = 1, np;
  std::vector < Halfedge * >baseRing, nextRing;
  // typename std::vector<Halfedge*>::iterator itb, ite;
  std::vector < Halfedge * >*p_baseRing = &baseRing, *p_nextRing =
    &nextRing;

    //collect the one ring
  np = T_PolyhedralSurf_rings < TPoly >::
    push_hedges_of(start, int (1), to_opp, nextRing, recordedOnce);

  while (1)
    {
      //exit?
      i++;
      if (i > ith)
	break;

	//find next ring
      p_baseRing->erase(p_baseRing->begin(), p_baseRing->end());
      std::swap(p_baseRing, p_nextRing);

      np = T_PolyhedralSurf_rings < TPoly >::
	next_ring_hedges(i, HTO_ONERING, *p_baseRing, *p_nextRing, 
			 recordedOnce);
    }
  //clean up hedges tags
  np = recordedOnce.size();	//#triangles
  if (cleanup)
    T_PolyhedralSurf_rings < TPoly >::reset_ring_indices(recordedOnce);
  return np;
}

//Warning contour edges are not given in a cyclic order
template < class TPoly >
void T_PolyhedralSurf_rings < TPoly >::
contour_ith_ring_hedges(Vertex * start, int ith,
			std::vector < Halfedge * > &contour)
{
  std::vector < Halfedge * > ringHedges;
  typename std::vector < Halfedge * >::iterator itb =
    ringHedges.begin(), ite = ringHedges.end();
  T_PolyhedralSurf_rings < TPoly >::
    collect_ith_ring_hedges(start, ith, HTO_ONERING, ringHedges, 0);
  //debug
  // fprintf(stderr, "\n");
  // int iii;
  //  for( iii=0;iii<ringHedges.size();iii++)
  //   fprintf(stderr, "[ ringedge %d , edgeindex %d, vertex %d ,opposite vertex %d] \n",
  // 	  ringHedges[iii],
  // 	  ringHedges[iii]->getRingIndex(),
  // 	  &(*( ringHedges[iii]->vertex() )), 
  // 	   &(*(ringHedges[iii]->opposite()->vertex() )) );
  //debug


  //set ringHedges index to i
  // for (; itb != ite; itb++)   (*itb)->setRingIndex(ith);non! garder
  // les index!!

  //for a hedge h in ringHedges and h index = ith, h->next is on the
  // contour if its opposite is not visited
  //Note that contour edges are then oriented CCW
  itb = ringHedges.begin(), ite = ringHedges.end();
  for (; itb != ite; itb++) 
    {
      if ( ( (*itb)->getRingIndex()==ith )&&
	   ((*itb)->next()->opposite()->getRingIndex() == -1) ) 
	contour.push_back( &(*((*itb)->next())) );
    }
  //debug
  // fprintf(stderr, "\n\n");
  //  for( iii=0;iii<contour.size();iii++)
  //   fprintf(stderr, "[ contour edge %d , edgeindex %d, vertex %d ,opposite vertex %d] \n",
  // 	  contour[iii], contour[iii]->getRingIndex(),
  //   &(*( contour[iii]->vertex() )), 
  // 	   &(*(contour[iii]->opposite()->vertex() )) );
  //  //debug
 
  //reset indices
  T_PolyhedralSurf_rings < TPoly >::reset_ring_indices(ringHedges);
}

//if the contour is a circle, then put edges head to tail, ccw
template < class TPoly >
int T_PolyhedralSurf_rings < TPoly >::
contour_ccw_ith_ring_hedges(Vertex * start, int ith,
			    std::vector < Halfedge * > &contourCCW)
{
  std::vector < Halfedge * > contour;
  T_PolyhedralSurf_rings < TPoly >::
    contour_ith_ring_hedges(start, ith, contour);

  //check if no border edge is encountered
  typename std::vector < Halfedge * >::iterator itb =
    contour.begin(), ite = contour.end();
  for (; itb != ite; itb++)  if ( (*itb)->is_border_edge()) return 0;
 
  int size = contour.size();
  Halfedge* h_cur, *h;
  h_cur = contour[0];
  // fprintf(stderr, "\n\n");
  // int ii;
  //  for( ii=0, itb=contour.begin();ii<contour.size();ii++, itb++){
  //   fprintf(stderr, "[%d  ", contour[ii]->vertex());
  //   fprintf(stderr, " %d] ", (*itb)->vertex());}

  // fprintf(stderr, "\n\n");
  contourCCW.push_back(h_cur);
  do
    {
      for (itb=contour.begin(); itb != ite; itb++)
	{//debug
	  h = &(*( ((*itb)->opposite()) ));
	  //	  cerr << "h_cur->vertex() " << h_cur->vertex() << endl;
	  //   << "h->vertex() "     << h->vertex() << endl;
	  if ( (h_cur->vertex()) == (h->vertex()) ) 
	    {//the edge following the current edge is found
	      contourCCW.push_back(*itb);
	      h_cur = *itb;
	      break;
	    }
	}//debug
    }
  while (contourCCW.size() != size);
 
  //debug
  // int iii;

  // fprintf(stderr, "\n");

  //  for( iii=0;iii<contour.size();iii++)
  //   fprintf(stderr, "[ contour edge %d , edgeindex %d, vertex %d ,opp vertex %d] \n",
  // 	  contour[iii],
  // 	  contour[iii]->getRingIndex(),
  // 	  &(*( contour[iii]->vertex() )), 
  // 	   &(*(contour[iii]->opposite()->vertex() )) );
  //    //debug
  //   //check the contour debug!  //debug
  // fprintf(stderr, "\n");
  //  for( iii=0;iii<contourCCW.size();iii++)
  //   fprintf(stderr, "[ contourCCW edge %d , edgeindex %d, vertex %d ,opp vertex %d] \n",
  // 	  contourCCW[iii],
  // 	  contourCCW[iii]->getRingIndex(),
  // 	  &(*( contourCCW[iii]->vertex() )), 
  // 	   &(*(contourCCW[iii]->opposite()->vertex() )) );
  //    //debug

  for (int i=0; i<contourCCW.size(); i++)
    {
      assert ( &(*(contourCCW[i]->vertex())) == 
	       &(*(contourCCW[(i+1) % size]->opposite()->vertex())) );
    }
  return 1;
}


template < class TPoly >
void T_PolyhedralSurf_rings < TPoly >::
reset_ring_indices(std::vector < Halfedge * >&hedges)
{
  typename std::vector < Halfedge * >::iterator itb =
    hedges.begin(), ite = hedges.end();
  for (; itb != ite; itb++)
    {
      (*itb)->resetRingIndex();
      (*itb)->opposite()->resetRingIndex();
    }
}

template < class TPoly >
void T_PolyhedralSurf_rings < TPoly >::
reset_ring_indices(std::vector < Facet * >&facets)
{
  typename std::vector < Facet * >::iterator itb = facets.begin(), ite =
    facets.end();
  for (; itb != ite; itb++)
    (*itb)->resetRingIndex();
}

// to report 1st crossing edge with a sphere
template < class TPoly >
typename TPoly::Point_3 T_PolyhedralSurf_rings < TPoly >::orig_point;

template < class TPoly >
void T_PolyhedralSurf_rings < TPoly >::set_orig_point(Point_3 & p)
{
  T_PolyhedralSurf_rings < TPoly >::orig_point = p;
}

template < class TPoly >
void T_PolyhedralSurf_rings < TPoly >::
getFarthestNeighbourOfV_ofOrigPoint(Vertex * v, Halfedge * &h_farthest)
{
  double d, dref;
  Halfedge *h;
  Halfedge_around_vertex_circulator
    hedgeb = v->vertex_begin(), hedgee = hedgeb;

  //init
  h = h_farthest = &(*hedgeb);
  hedgeb++;
  dref = squared_distance(h->opposite()->vertex()->point(),
			  T_PolyhedralSurf_rings < TPoly >::orig_point);
  //Visit neighbours  
  do
    {
      h = &(*hedgeb);
      hedgeb++;

      assert(v == (&(*h->vertex())));

      d = squared_distance(h->opposite()->vertex()->point(),
			   T_PolyhedralSurf_rings < TPoly >::orig_point);
      if (d > dref)
	{
	  dref = d;
	  h_farthest = h;
	}
    }
  while (hedgeb != hedgee);
}

#endif
