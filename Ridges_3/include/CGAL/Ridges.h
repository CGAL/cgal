#ifndef _RIDGE_3_H_
#define _RIDGE_3_H_

#include <CGAL/basic.h>
#include <utility>
#include <list>
#include "PolyhedralSurf_neighbors.h"
#include "Umbilic.h"

//note : one has to orient monge normals according to mesh normals to
//define min/max curv

CGAL_BEGIN_NAMESPACE

enum Ridge_type {NONE=0, BLUE_RIDGE, RED_RIDGE, CREST, BE, BH, BC, RE, RH, RC};

//---------------------------------------------------------------------------
//Ridge_line : a connected sequence of edges crossed by a ridge, with
//type and weigths
//--------------------------------------------------------------------------
template < class Poly > class Ridge_line
{
public:
  typedef typename Poly::Traits::FT FT;
  typedef typename Poly::Traits::Vector_3 Vector_3;
  typedef typename Poly::Traits::Point_3 Point_3;
  typedef typename Poly::Vertex_handle Vertex_handle;
  typedef typename Poly::Halfedge_handle Halfedge_handle;
  typedef std::pair< Halfedge_handle, FT> ridge_he;

protected:
  Ridge_type m_line_type;//one of BE, BH, BC, RE, RH or RC
  std::list<ridge_he> m_line;
  FT m_strength;
  FT m_sharpness;

public:
  const Ridge_type line_type() const {return m_line_type;}
  const FT strength() const {return m_strength;}
  const FT sharpness() const {return m_sharpness;}
  std::list<ridge_he>* line() { return &m_line;}
  const std::list<ridge_he>* line() const { return &m_line;}

  //constructor
  //a ridge line begins with a segment in a triangle
  Ridge_line( Halfedge_handle h1, Halfedge_handle h2, Ridge_type r_type);

  //compute the barycentric coordinate of the xing point (blue or red)
  //for he: p->q  coord is st xing_point = coord*p + (1-coord)*q
  FT bary_coord( Halfedge_handle he); 

  //When the line is extended with a he, the bary coord of the
  //crossing point is computed, the pair (he,coord) is added and the
  //weigths are updated 
  void addback( Halfedge_handle he);
  void addfront( Halfedge_handle he);
  
  void dump_4ogl(std::ostream& out_stream);
};

// IMPLEMENTATION OF Ridge_line members
//////////////////////////////////////////////////////////////////////////////

 //constructor
 template < class Poly >
 Ridge_line<Poly>::
 Ridge_line( Halfedge_handle h1, Halfedge_handle h2, Ridge_type r_type)
 : m_line_type(r_type), m_strength(0.), m_sharpness(0.)
    {
      m_line.push_back(ridge_he(h1, bary_coord(h1)));
      addback(h2);
    }

template < class Poly >
void Ridge_line<Poly>::
addback( Halfedge_handle he) 
{
  Halfedge_handle he_cur = ( --(m_line.end()) )->first;
  FT coord_cur = ( --(m_line.end()) )->second;//bary_coord(he_cur);
  FT coord = bary_coord(he);
  Vertex_handle v_p = he->opposite()->vertex(), v_q = he->vertex(),
    v_p_cur = he_cur->opposite()->vertex(), v_q_cur = he->vertex(); // he: p->q
  FT k;//abs value of the ppal curvature at the Xing point on he.
  if ( (m_line_type == BE) || (m_line_type == BH) || (m_line_type == BC) ) {
    k = CGAL::abs(v_p->k1()) * coord + CGAL::abs(v_q->k1()) * (1-coord) ;   
  }
  if ( (m_line_type == RE) || (m_line_type == RH) || (m_line_type == RC) ) {
    k = CGAL::abs(v_p->k2()) * coord + CGAL::abs(v_q->k2()) * (1-coord) ;   
  }
  Vector_3 segment = (v_p->point()-ORIGIN)*coord + (v_q->point()-ORIGIN)*(1-coord) - 
    ((v_p_cur->point()-ORIGIN)*coord_cur + (v_q_cur->point()-ORIGIN)*(1-coord_cur));
  m_strength += k * CGAL::sqrt(segment * segment); 
  //TODO update sharpness
  m_line.push_back( ridge_he(he, coord));
}

template < class Poly >
void Ridge_line<Poly>::
addfront( Halfedge_handle he) 
{
  Halfedge_handle he_cur = ( m_line.begin() )->first;
  FT coord_cur = ( m_line.begin() )->second;
  FT coord = bary_coord(he);
  Vertex_handle v_p = he->opposite()->vertex(), v_q = he->vertex(),
    v_p_cur = he_cur->opposite()->vertex(), v_q_cur = he->vertex(); // he: p->q
  FT k;
  if ( (m_line_type == BE) || (m_line_type == BH) || (m_line_type == BC) ) {
    k = CGAL::abs(v_p->k1()) * coord + CGAL::abs(v_q->k1()) * (1-coord) ;   
  }
  if ( (m_line_type == RE) || (m_line_type == RH) || (m_line_type == RC) ) {
    k = CGAL::abs(v_p->k2()) * coord + CGAL::abs(v_q->k2()) * (1-coord) ;   
  }
  Vector_3 segment = (v_p->point()-ORIGIN)*coord + (v_q->point()-ORIGIN)*(1-coord) - 
    ((v_p_cur->point()-ORIGIN)*coord_cur + (v_q_cur->point()-ORIGIN)*(1-coord_cur));
  m_strength += k * CGAL::sqrt(segment * segment); 
  //TODO update sharpness
  m_line.push_front( ridge_he(he, coord));
}


template < class Poly >
typename Poly::Traits::FT Ridge_line<Poly>::
bary_coord( Halfedge_handle he)
{
  FT b_p, b_q; // extremalities at p and q for he: p->q
  if ( (m_line_type == BE) || (m_line_type == BH) || (m_line_type == BC) ) {
    b_p = he->opposite()->vertex()->b0();
    b_q = he->vertex()->b0();    
  }
  if ( (m_line_type == RE) || (m_line_type == RH) || (m_line_type == RC) ) {
    b_p = he->opposite()->vertex()->b3();
    b_q = he->vertex()->b3();    
  }
  return CGAL::abs(b_q) / ( CGAL::abs(b_q) + CGAL::abs(b_p) );
}

template < class Poly >
void Ridge_line<Poly>::
dump_4ogl(std::ostream& out_stream)
{
  out_stream << line_type() << " "
	     << strength() << " "
	     << sharpness() << " ";

  typename std::list<ridge_he >::iterator
    iter = line()->begin(), 
    ite =  line()->end();
  for (;iter!=ite;iter++){
    //he: p->q, r is the crossing point
    Point_3 p = iter->first->opposite()->vertex()->point(),
      q = iter->first->vertex()->point();
    Vector_3 r = (p-CGAL::ORIGIN)*iter->second +
      (q-CGAL::ORIGIN)*(1-iter->second); 
    out_stream << " " << r ;	
  }
  out_stream  << std::endl;  
}

//---------------------------------------------------------------------------
//Ridge_approximation
//--------------------------------------------------------------------------

template < class Poly, class OutputIt >
class Ridge_approximation
{
public:
  typedef typename Poly::Traits::FT FT;
  typedef typename Poly::Traits::Vector_3 Vector_3;
  typedef typename Poly::Vertex_handle Vertex_handle;
  typedef typename Poly::Halfedge_handle Halfedge_handle;
  typedef typename Poly::Facet_handle Facet_handle;
  typedef typename Poly::Facet_iterator Facet_iterator;
  typedef Ridge_line<Poly> Ridge_line;
  //  typedef T_PolyhedralSurf_neighbors<Poly> Poly_neighbors;//for umbilics
  //are ridges tagged as elliptic or hyperbolic using 3rd or 4th order
  //differential quantitities?
  enum Tag_order {Tag_3 = 3, Tag_4 = 4};

public:
  Ridge_approximation(){};

  OutputIt compute_all_ridges(Poly &P, OutputIt it, Tag_order ord = Tag_3);
  
  //Find BLUE_RIDGE, RED_RIDGE or CREST ridges 
  //iterate on P facets, find a non-visited, regular, 2Xing triangle,
  //follow non-visited, regular, 2Xing triangles in both sens to create
  //a Ridge line.
  //Each time an edge is added the strength of the current line is updated
  // + length(ridge segment in the facet)*|k|
  void compute_ridges(Poly &P, 
		      Ridge_type r_type, 
		      OutputIt ridge_lines_it,
		      Tag_order ord = Tag_3);
  // void compute_umbilics(Poly &P	);//container, class for umbilics?

protected:
  //is a facet crossed by a BLUE, RED or CREST ridge? if so, return
  //the crossed edges and more precise type from BE, BH, BC, RE, RH,
  //RC or NONE
  Ridge_type facet_ridge_type(Facet_handle f, 
			      Halfedge_handle& he1, 
			      Halfedge_handle& he2,
			      Ridge_type r_type,
			      Tag_order ord = Tag_3);
  
  //is an edge crossed by a BLUE/RED ridge? (color is BLUE_RIDGE or
  //RED_RIDGE ).  As we only test edges of regular triangles, the ppal
  //direction at endpoints d_p and d_q cannot be orthogonal. If both
  //extremalities vanish, we consider no crossing occurs. If only one
  //of them vanishes, we consider it as an positive infinitesimal and
  //apply the general rule. The general rule is that for both
  //non-vanishing extremalities, a crossing occurs if their sign
  //differ; Assuming the accute rule to orient the ppal directions,
  //there is a crossing iff d_p.d_q * b_p*b_q < 0
  void xing_on_edge(Halfedge_handle he, 
		    bool& is_crossed, 
		    Ridge_type color);
 
  //for the computation with tag_order = 3
  //for a ridge segment [r1,r2] in a triangle (v1,v2,v3), let r = r2 -
  //r1 and normalize, the projection of a point p on the line (r1,r2)
  //is pp=r1+tr, with t=(p-r1)*r then the vector v starting at p is
  //pointing to the ridge line (r1,r2) if (pp-p)*v >0. Return the sign
  //of b, for a ppal direction pointing to the ridge segment,
  //appearing at least at two vertices of the facet.
  // for color = BLUE_RIDGE, sign = 1 if BE, -1 if BH
  // for color = RED_RIDGE, sign = -1 if RE, 1 if RH
  int b_sign_pointing_to_ridge(Vertex_handle v1, 
			       Vertex_handle v2,
			       Vertex_handle v3,
			       Vector_3 r1, Vector_3 r2, 
			       Ridge_type color);
};


// IMPLEMENTATION OF Ridge_approximation members
/////////////////////////////////////////////////////////////////////////////
template < class Poly, class OutputIt >
OutputIt Ridge_approximation<Poly, OutputIt>::
compute_all_ridges(Poly &P, OutputIt it, Tag_order ord)
{
  compute_ridges(P, BLUE_RIDGE, it, ord);
  compute_ridges(P, RED_RIDGE, it, ord);
  compute_ridges(P, CREST, it, ord);
  return it;
}

template < class Poly, class OutputIt >
void Ridge_approximation<Poly, OutputIt>::
compute_ridges(Poly &P, Ridge_type r_type,
	       OutputIt ridge_lines_it, Tag_order ord)
{
  //set all facets non visited
  Facet_iterator itb = P.facets_begin(), ite = P.facets_end();
  for(;itb!=ite;itb++) itb->reset_is_visited();
 
  itb = P.facets_begin();
  for(;itb!=ite;itb++)
    {
      Facet_handle f= &(*itb);
      if (f->is_visited()) continue;
      f->set_visited(true);
      Halfedge_handle h1, h2, curhe1, curhe2, curhe;
      
      //h1 h2 are the hedges crossed if any, r_type should be
      //BLUE_RIDGE, RED_RIDGE or CREST ; cur_ridge_type should be BE,
      //BH, BC, RE, RH, RC or NONE
      Ridge_type cur_ridge_type = facet_ridge_type(f,h1,h2,r_type, ord);
      if ( cur_ridge_type == NONE ) continue;
      
      //a ridge_line is begining and stored
      //     Ridge_line cur_ridge_line(h1,h2,cur_ridge_type); 
      Ridge_line* cur_ridge_line = new Ridge_line(h1,h2,cur_ridge_type); 
      *ridge_lines_it++ = cur_ridge_line;
      //debug
      //  cur_ridge_line->dump_4ogl(std::cout);
      // std::cout << "??????????????????????????" << endl;

      //next triangle adjacent to h1 (push_front)
      if ( !(h1->is_border_edge()) ) 
	{
	  f = h1->opposite()->facet();
	  curhe = h1;
	  while (cur_ridge_type == facet_ridge_type(f,curhe1,curhe2,
						    r_type, ord))
	    {
	      //follow the ridge from curhe
	      if (f->is_visited()) break;
	      f->set_visited(true);
	      if (curhe->opposite() == curhe1) curhe = curhe2;
	      else curhe = curhe1;//curhe stays at the ridge extremity
	      cur_ridge_line->addfront(curhe);
	      if ( !(curhe->is_border_edge()) ) f =
						  curhe->opposite()->facet();
	      else break;
	    }
	  //exit from the while if
	  //1. border or already visited (this is a ridge loop)
	  //2. not same type, then do not set visisted cause a BE
	  //	  follows a BH
	}

      //next triangle adjacent to h2 (push_back)
      if ( !(h2->is_border_edge()) ) 
	{
	  f = h2->opposite()->facet();
	  curhe = h2;
	  while (cur_ridge_type ==
		 facet_ridge_type(f,curhe1,curhe2,r_type, ord))
	    {
	      //follow the ridge from curhe
	      if (f->is_visited()) break;
	      f->set_visited(true);
	      if (curhe->opposite() == curhe1) curhe = curhe2;
	      else curhe = curhe1;
	      cur_ridge_line->addback(curhe);
	      if ( !(curhe->is_border_edge()) ) f =
						  curhe->opposite()->facet();
	      else break;
	    }
	} 
    }
}

template < class Poly, class OutputIt >
Ridge_type Ridge_approximation<Poly, OutputIt>::
facet_ridge_type(Facet_handle f, Halfedge_handle& he1, Halfedge_handle&
		 he2, Ridge_type r_type, Tag_order ord)
{
  //polyhedral data
  //we have v1--h1-->v2--h2-->v3--h3-->v1
  Halfedge_handle h1 = f->halfedge(), h2, h3;
  Vertex_handle v1, v2, v3;
  v2 = h1->vertex();
  h2 = h1->next();
  v3 = h2->vertex();
  h3 = h2->next();
  v1 = h3->vertex();

  //check for regular facet
  //i.e. if there is a coherent orientation of ppal dir at the facet vertices
  if ( v1->d1()*v2->d1() * v1->d1()*v3->d1() * v2->d1()*v3->d1() < 0 ) 
    return NONE;
   
  //determine potential crest color
  //BC if |sum(k1)|>|sum(k2)| sum over facet vertices vi
  //RC if |sum(k1)|<|sum(k2)|
  Ridge_type crest_color = NONE;
  if (r_type == CREST) 
    {
      if ( CGAL::abs(v1->k1()+v2->k1()+v3->k1()) > CGAL::abs(v1->k2()+v2->k2()+v3->k2()) ) 
	crest_color = BC; 
      if ( CGAL::abs(v1->k1()+v2->k1()+v3->k1()) < CGAL::abs(v1->k2()+v2->k2()+v3->k2()) ) 
	crest_color = RC;
      if ( CGAL::abs(v1->k1()+v2->k1()+v3->k1()) == CGAL::abs(v1->k2()+v2->k2()+v3->k2()) ) 
	return NONE;
    }
  
  //compute Xing on the 3 edges
  bool h1_is_crossed = false, h2_is_crossed = false, h3_is_crossed = false;
  if ( r_type == BLUE_RIDGE || crest_color == BC ) 
    {
      xing_on_edge(h1, h1_is_crossed, BLUE_RIDGE);
      xing_on_edge(h2, h2_is_crossed, BLUE_RIDGE);
      xing_on_edge(h3, h3_is_crossed, BLUE_RIDGE);
    }
  if ( r_type == RED_RIDGE || crest_color == RC ) 
    {
      xing_on_edge(h1, h1_is_crossed, RED_RIDGE);
      xing_on_edge(h2, h2_is_crossed, RED_RIDGE);
      xing_on_edge(h3, h3_is_crossed, RED_RIDGE);
    }

  //there are either 0 or 2 crossed edges
  if ( !h1_is_crossed && !h2_is_crossed && !h3_is_crossed ) 
    return NONE; 
  if (h1_is_crossed && h2_is_crossed && !h3_is_crossed)
    {
      he1 = h1; 
      he2 = h2;
    }
  if (h1_is_crossed && !h2_is_crossed && h3_is_crossed)
    {
      he1 = h1; 
      he2 = h3;
    }
  if (!h1_is_crossed && h2_is_crossed && h3_is_crossed)
    {
      he1 = h2; 
      he2 = h3;
    }
  //check there is no other case (just on edge crossed)
  assert ( !( (h1_is_crossed && !h2_is_crossed && !h3_is_crossed)
	      || (!h1_is_crossed && h2_is_crossed && !h3_is_crossed)
	      || (!h1_is_crossed && h2_is_crossed && !h3_is_crossed)) );

  //There is a ridge segment in the triangle, determine its type
  Vertex_handle v_p1 = he1->opposite()->vertex(), v_q1 = he1->vertex(),
    v_p2 = he2->opposite()->vertex(), v_q2 = he2->vertex(); // he1: p1->q1
 
  if ( r_type == BLUE_RIDGE || crest_color == BC ) {
    FT coord1 = CGAL::abs(v_q1->b0()) / ( CGAL::abs(v_p1->b0()) +
					  CGAL::abs(v_q1->b0()) ), 
      coord2 = CGAL::abs(v_q2->b0()) / ( CGAL::abs(v_p2->b0()) +
					 CGAL::abs(v_q2->b0()) ); 
    if ( ord == Tag_3 ) {
      Vector_3 r1 = (v_p1->point()-ORIGIN)*coord1 +
	(v_q1->point()-ORIGIN)*(1-coord1), 
	r2 = (v_p2->point()-ORIGIN)*coord2 +
	(v_q2->point()-ORIGIN)*(1-coord2); 
      int b_sign = b_sign_pointing_to_ridge(v1, v2, v3, r1, r2, BLUE_RIDGE); 
      if (r_type == CREST) {if (b_sign == 1) return BC; else return NONE;} 
      if (b_sign == 1) return BE; else return BH; 
    }
    else {//ord == Tag_4, check the sign of the meanvalue of the signs
      //      of P1 at the two crossing points
      FT sign_P1 =  v_p1->P1()*coord1 + v_q1->P1()*(1-coord1) 
	+ v_p2->P1()*coord2 + v_q2->P1()*(1-coord2);
      if (r_type == CREST) {if ( sign_P1 < 0 ) return BC; else return NONE;}
      if ( sign_P1 < 0 ) return BE; else return BH;
    }
  }
 
  if ( r_type == RED_RIDGE || crest_color == RC ) {
    FT coord1 = CGAL::abs(v_q1->b3()) / ( CGAL::abs(v_p1->b3()) +
					  CGAL::abs(v_q1->b3()) ), 
      coord2 = CGAL::abs(v_q2->b3()) / ( CGAL::abs(v_p2->b3()) +
					 CGAL::abs(v_q2->b3()) ); 
    if ( ord == Tag_3 ) {
      Vector_3 r1 = (v_p1->point()-ORIGIN)*coord1 +
	(v_q1->point()-ORIGIN)*(1-coord1), 
	r2 = (v_p2->point()-ORIGIN)*coord2 +
	(v_q2->point()-ORIGIN)*(1-coord2); 
      int b_sign = b_sign_pointing_to_ridge(v1, v2, v3, r1, r2, RED_RIDGE);
      if (r_type == CREST) {if (b_sign == -1) return RC; else return NONE;} 
      if (b_sign == -1) return RE; else return RH; 
    } 
    else {//ord == Tag_4, check the sign of the meanvalue of the signs
      //      of P2 at the two crossing points
      FT sign_P2 =  v_p1->P2()*coord1 + v_q1->P2()*(1-coord1) 
	+ v_p2->P2()*coord2 + v_q2->P2()*(1-coord2);
      if (r_type == CREST) {if ( sign_P2 < 0 ) return RC; else return NONE;}
      if ( sign_P2 < 0 ) return RE; else return RH;
    } 
  }
  assert(0);//should return before!
}

template < class Poly, class OutputIt >
void Ridge_approximation<Poly, OutputIt>::
xing_on_edge(Halfedge_handle he, bool& is_crossed, Ridge_type color)
{
  is_crossed = false;
  FT sign;
  FT b_p, b_q; // extremalities at p and q for he: p->q
  Vector_3  d_p = he->opposite()->vertex()->d1(),
    d_q = he->vertex()->d1(); //ppal dir
  if ( color == BLUE_RIDGE ) {
    b_p = he->opposite()->vertex()->b0();
    b_q = he->vertex()->b0();
  }
  else {     
    b_p = he->opposite()->vertex()->b3();
    b_q = he->vertex()->b3();
  }
  if ( b_p == 0 && b_q == 0 ) return;
  if ( b_p == 0 && b_q !=0 ) sign = d_p*d_q * b_q;
  if ( b_p != 0 && b_q ==0 ) sign = d_p*d_q * b_p;
  if ( b_p != 0 && b_q !=0 ) sign = d_p*d_q * b_p * b_q;
  if ( sign < 0 ) is_crossed = true;
}


template < class Poly, class OutputIt >
int Ridge_approximation<Poly, OutputIt>::
b_sign_pointing_to_ridge(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
			 Vector_3 r1, Vector_3 r2, Ridge_type color)
{
  Vector_3 r = r2 - r1, dv1, dv2, dv3;
  FT bv1, bv2, bv3;
  if ( color == BLUE_RIDGE ) {
    bv1 = v1->b0();
    bv2 = v2->b0();
    bv3 = v3->b0();
    dv1 = v1->d1();
    dv2 = v2->d1();
    dv3 = v3->d1();
  }
  else {
    bv1 = v1->b3();
    bv2 = v2->b3();
    bv3 = v3->b3();
    dv1 = v1->d2();
    dv2 = v2->d2();
    dv3 = v3->d2();    
  }
  if ( r != CGAL::NULL_VECTOR ) r = r/CGAL::sqrt(r*r);
  FT sign1, sign2, sign3;
  sign1 = (r1 - (v1->point()-ORIGIN) + (((v1->point()-ORIGIN)-r1)*r)*r )*dv1;
  sign2 = (r1 - (v2->point()-ORIGIN) + (((v2->point()-ORIGIN)-r1)*r)*r )*dv2;
  sign3 = (r1 - (v3->point()-ORIGIN) + (((v3->point()-ORIGIN)-r1)*r)*r )*dv3;
  
  int compt = 0;
  if ( sign1 > 0 ) compt++; else if (sign1 < 0) compt--;
  if ( sign2 > 0 ) compt++; else if (sign2 < 0) compt--;
  if ( sign3 > 0 ) compt++; else if (sign3 < 0) compt--;
  
  if (compt > 0) return 1; else return -1;
}

//---------------------------------------------------------------------------
//Umbilic_approximation
//--------------------------------------------------------------------------

template < class Poly, class OutputIt >
class Umbilic_approximation
{
public:
  typedef typename Poly::Traits::FT FT;
  typedef typename Poly::Traits::Vector_3 Vector_3;
  typedef typename Poly::Vertex_handle Vertex_handle;
  typedef typename Poly::Halfedge_handle Halfedge_handle;
  typedef typename Poly::Facet_handle Facet_handle;
  typedef typename Poly::Facet_iterator Facet_iterator;
  typedef typename Poly::Vertex_iterator Vertex_iterator;
  typedef T_PolyhedralSurf_neighbors<Poly> Poly_neighbors;
  typedef Umbilic<Poly> Umbilic;
  CGAL::Abs<FT> cgal_abs;
  static FT neigh_size;//the size of neighbourhood for umbilic
  //  computation is (neigh_size * OneRingSize)
 
  Umbilic_approximation(){};
  OutputIt compute(Poly &P, OutputIt it, FT size);

};

template < class Poly, class OutputIt >
  OutputIt Umbilic_approximation< Poly, OutputIt >::
compute(Poly &P, OutputIt umbilics_it, FT size)
{
  std::vector<Vertex_handle> vces;
  std::list<Halfedge_handle> contour;
  double umbilicEstimatorVertex, umbilicEstimatorNeigh;
  
  bool is_umbilic = true;

  //MAIN loop on P vertices
  Vertex_iterator itb = P.vertices_begin(), ite = P.vertices_end();
  for (;itb != ite; itb++) {
    Vertex_handle vh = itb;
    umbilicEstimatorVertex = cgal_abs(vh->k1()-vh->k2());
    //reset vector, list and bool
    vces.clear();
    contour.clear();
    is_umbilic = true;
    Poly_neighbors::compute_neighbors(vh, vces, contour, size);
    
    
    // OPTIONAL: avoid umbilics whose contours touch the border
    typename std::list<Halfedge_handle>::iterator itb_cont = contour.begin(),
      ite_cont = contour.end();
    for (; itb_cont != ite_cont; itb_cont++)
      if ( (*itb_cont)->is_border() ) {is_umbilic = false; continue;}
    if (is_umbilic == false) continue;
    
    //is v an umbilic?
    //a priori is_umbilic = true, and it switches to false as soon as a 
    //  neigh vertex has a lower umbilicEstimator value
    typename std::vector<Vertex_handle>::iterator itbv = vces.begin(),
      itev = vces.end();
    assert(*itbv ==  vh);
    itbv++;
    for (; itbv != itev; itbv++)
      {	umbilicEstimatorNeigh = cgal_abs( (*itbv)->k1() - (*itbv)->k2() );
	if ( umbilicEstimatorNeigh < umbilicEstimatorVertex ) 
	  {is_umbilic = false; break;}
      }
    if (is_umbilic == false) continue;
    
    //v is an umbilic, compute the index
    Umbilic*  cur_umbilic = new Umbilic(vh, contour);
    cur_umbilic->compute_type();
    *umbilics_it++ = cur_umbilic;
  }
  return umbilics_it;
}

CGAL_END_NAMESPACE

#endif
