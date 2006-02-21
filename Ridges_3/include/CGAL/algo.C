#ifndef _RIDGE_3_H_
#define _RIDGE_3_H_

#include <pair>


//Orient monge normal according to mesh normals.

CGAL_BEGIN_NAMESPACE


enum Ridge_type {NONE=0, BLUE, RED, CREST, BE, BH, BC, RE, RH, RC};


//---------------------------------------
//Ridge_line : a connected sequence of edges crossed by a ridge, with type and weigths
//---------------------------------------
template < class Poly > class Ridge_line
{
public:
  typedef typename Poly::Traits::FT FT;
  typedef typename Poly::Vertex_handle Vertex_handle;
  typedef typename Poly::Halfedge_handle Halfedge_handle;
  typedef typename Poly::Facet_handle Facet_handle;
  typedef std::pair<Halfedge_handle, FT> ridge_he;

protected:
  Ridge_type line_type
  std::list<ridge_he> line;
  FT strength;
  FT sharpness;

public:
  const FT weight() const {return weight;}
  const FT sharpness() const {return sharpness;}
  std::list<ridge_he>* line() { return &line;}

  //constructor
  Ridge_line( Halfedge_handle h1, Halfedge_handle h2, Ridge_type r_type) :
    line_type(r_type), strength(0.)
  {};

  //compute the barycentric coordinate of the xing point (blue or red)
  //for he: p->q  coord is st xing_point = coord*p + (1-coord)*q
  FT bary_coord(Halfedge_handle he); 
  void compute_weight(char color);
  void compute_sharpness(char color);
  //When the line is extended with a he, the bary coord of the
  //crossing point is computed, the pair (he,coord) is added and the
  //weigths are updated 
  void addback( Halfedge_handle he);
  void addfront( Halfedge_handle he);

};

// IMPLEMENTATION OF Ridge_line members
//////////////////////////////////////////////////////////////////////////////

//constructor
template < class Poly >
Ridge_line( Halfedge_handle h1, Halfedge_handle h2, Ridge_type r_type) :
  line_type(r_type), strength(0.)
{
  line.push_back(h1);
  addback(h2);
}

template < class Poly >
void Ridge_line<Poly>::
addback( Halfedge_handle he) 
{
  Halfedge_handle he_cur = ( --(line.end()) )->first;
  FT coord = bary_coord(he);
  FT coord_cur = bary_coord(he_cur);
  Vertex_handle v_p = he->opposite()->vertex(), v_q = he->vertex(),
    v_p_cur = he_cur->opposite()->vertex(), v_q_cur = he->vertex(); // he: p->q
  FT k;
  if ( (line_type == BE) || (line_type == BH) || (line_type == BC) ) {
    k =( std::fabs(v_p->k1) * coord + std::fabs(v_q->k1) * (1-coord) )/2;   
  }
  if ( (line_type == RE) || (line_type == RH) || (line_type == RC) ) {
    k =( std::fabs(v_p->k2) * coord + std::fabs(v_q->k2) * (1-coord) )/2;   
  }
  Vector_3 segment = (v_p->point()-ORIGIN)*coord + (v_q->point()-ORIGIN)*(1-coord) - 
    ((v_p_cur->point()-ORIGIN)*coord_cur + (v_q_cur->point()-ORIGIN)*(1-coord_cur));
  strength += k * CGAL::sqrt(segment * segment); 
    
  line.push_back( pair(he, coord));
}

template < class Poly >
void Ridge_line<Poly>::
addfront( Halfedge_handle he) 
{
  Halfedge_handle he_cur = ( line.begin() )->first;
  FT coord = bary_coord(he);
  FT coord_cur = bary_coord(he_cur);
  Vertex_handle v_p = he->opposite()->vertex(), v_q = he->vertex(),
    v_p_cur = he_cur->opposite()->vertex(), v_q_cur = he->vertex(); // he: p->q
  FT k;
  if ( (line_type == BE) || (line_type == BH) || (line_type == BC) ) {
    k =( std::fabs(v_p->k1) * coord + std::fabs(v_q->k1) * (1-coord) )/2;   
  }
  if ( (line_type == RE) || (line_type == RH) || (line_type == RC) ) {
    k =( std::fabs(v_p->k2) * coord + std::fabs(v_q->k2) * (1-coord) )/2;   
  }
  Vector_3 segment = (v_p->point()-ORIGIN)*coord + (v_q->point()-ORIGIN)*(1-coord) - 
    ((v_p_cur->point()-ORIGIN)*coord_cur + (v_q_cur->point()-ORIGIN)*(1-coord_cur));
  strength += k * CGAL::sqrt(segment * segment); 
    
  line.push_front( pair(he, coord));
}


template < class Poly >
FT Ridge_line<Poly>::
bary_coord(Halfedge_handle he)
{
  FT b_p, b_q; // extremalities at p and q for he: p->q
  if ( (line_type == BE) || (line_type == BH) || (line_type == BC) ) {
    b_p = he->opposite()->vertex()->b0();
    b_q = he->vertex()->b0();    
  }
  if ( (line_type == RE) || (line_type == RH) || (line_type == RC) ) {
    b_p = he->opposite()->vertex()->b3();
    b_q = he->vertex()->b3();    
  }
  return std::fabs(b_q) / ( std::fabs(b_q) + std::fabs(b_p) );
}


/////////Ridge_approximation//////////////////////////////////////////

//Find BLUE ridges (Elliptic or Hyperbolic) 
//do the same for RED and CREST
//iterate on P facets, find a non-visited, regular, 2BXing triangle,
//follow non-visited, regular, 2BXing triangles in both sens to create
//a Ridge line.
//Each time a edge is added the strength of the current line is updated
// + length(edge)*|k|
void
compute_ridges(Ridge_type r_type)
{
  //set all facets non visited
  

  Facet_iterator itb = P.facets_begin(), ite = P.facets_end();
  for(;itb!=ite;itb++)
    {
      Facet_handle f= &(*itb);
      if (f->is_visited()) continue;
      f->set_visited(true);
      Halfedge_handle h1, h2, curhe1, curhe2, curhe;
      
      //h1 h2 are the hedges crossed if any, r_type should be BLUE,
      //RED or CREST ; cur_ridge_type should be BE, BH, BC, RE, RH, RC or NONE
      Ridge_type cur_ridge_type = facet_ridge_type(f,h1,h2,r_type)
      if ( cur_ridge_type == NONE ) continue;
      
      //When a he is inserted in a Ridge_line, a pair <he,bary_coord>
      //is created and the weight(s) are updated
      Ridge_line *cur_ridge_line = new Ridge_line(h1,h2,cur_ridge_type); 
      *ridge_lines_it++ = cur_ridge_line;

      //next triangle adjacent to h1 (push_front)
      if ( !(h1->is_border_edge()) ) 
	{
	  f = h1->opposite()->facet();
	  curhe = h1;
	  while (cur_ridge_type == facet_ridge_type(f,curhe1,curhe2,r_type))
	    {
	      //follow the ridge from curhe
	      if (f->is_visited()) break;
	      f->set_visited(true);
	      if (curhe->opposite() == curhe1) curhe = curhe2;
	      else curhe = curhe1;
	      cur_ridge_line->addfront(curhe);
	      if ( !(curhe->is_border_edge()) ) f =
						  curhe->opposite()->facet();
	      else break;
	    }
	  //exit from the while if
	  //1. border
	  //2. not same type, then do not set visisted cause a BE
	  //	  follows a BH
	}

      //next triangle adjacent to h2 (push_back)
      if ( !(h2->is_border_edge()) ) 
	{
	  f = h2->opposite()->facet();
	  curhe = h2;
	  while (cur_ridge_type == facet_ridge_type(f,curhe1,curhe2,r_type))
	    {
	      //follow the ridge from curhe
	      if (f->is_visitedB()) break;
	      f->set_visitedB(true);
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


Ridge_type 
facet_ridge_type(Facet_handle f, Halfedge_handle& he1, Halfedge_handle&
		 he2, Ridge_type r_type)
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
  if ( v1->d1()*v2->d1() * v1->d1()*v3->d1() * v2->d1()*v3->d1() < 0 ) return NONE;
   
  //determine potential crest color
  Ridge_type crest_color = NONE;
  if (r_type == CREST) 
    {
      if ( std::fabs(v1->k1()+v2->k1()+v3->k1()) > std::fabs(v1->k2()+v2->k2()+v3->k2()) ) 
	crest_color = BC; 
      if ( std::fabs(v1->k1()+v2->k1()+v3->k1()) < std::fabs(v1->k2()+v2->k2()+v3->k2()) ) 
	crest_color = RC;
      if ( std::fabs(v1->k1()+v2->k1()+v3->k1()) = std::fabs(v1->k2()+v2->k2()+v3->k2()) ) 
	return NONE;
    }
  
  //compute Xing on the 3 edges
  bool h1_is_crossed, h2_is_crossed, h3_is_crossed;
  if ( r_type == BLUE || crest_color == BC ) 
    {
      xing_on_edge(h1, h1_is_crossed, BLUE);
      xing_on_edge(h2, h2_is_crossed, BLUE);
      xing_on_edge(h3, h3_is_crossed, BLUE);
    }
  if ( r_type == RED || crest_color == RC ) 
    {
      xing_on_edge(h1, h1_is_crossed, RED);
      xing_on_edge(h2, h2_is_crossed, RED);
      xing_on_edge(h3, h3_is_crossed, RED);
    }

  //test of 2Xing
  if ( !(h1_is_crossed && h2_is_crossed && !h3_is_crossed) &&
       !(h1_is_crossed && !h2_is_crossed && h3_is_crossed) &&
       (!h1_is_crossed && h2_is_crossed && h3_is_crossed) ) 
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
  
  Vertex_handle v_p1 = he1->opposite()->vertex(), v_q1 = he1->vertex(),
    v_p2 = he2->opposite()->vertex(), v_q2 = he2->vertex(); // he1: p1->q1
 
  if ( r_type == BLUE || crest_color == BC ) 
    {
      FT coord1 = std::fabs(v_q1->b0()) / ( std::fabs(v_p1->b0()) + std::fabs(v_q1->b0()) ),
	coord2 = std::fabs(v_q2->b0()) / ( std::fabs(v_p2->b0()) + std::fabs(v_q2->b0()) );
      Vector_3 r1 = (v_p1->point()-ORIGIN)*coord1 + (v_q1->point()-ORIGIN)*(1-coord1),
	r2 = (v_p2->point()-ORIGIN)*coord2 + (v_q2->point()-ORIGIN)*(1-coord2);
      int b_sign = b_sign_pointing_to_ridge(v1, v2, v3, r1, r2, BLUE);
      if (b_sign == 1) { if (r_type == BLUE) return BE; else return BC;}
      if (b_sign == -1) return BH;
    }
  if ( r_type == RED || crest_color == RC ) 
    {
      FT coord1 = std::fabs(v_q1->b3()) / ( std::fabs(v_p1->b3()) + std::fabs(v_q1->b3()) ),
	coord2 = std::fabs(v_q2->b3()) / ( std::fabs(v_p2->b3()) + std::fabs(v_q2->b3()) );
      Vector_3 r1 = (v_p1->point()-ORIGIN)*coord1 + (v_q1->point()-ORIGIN)*(1-coord1),
	r2 = (v_p2->point()-ORIGIN)*coord2 + (v_q2->point()-ORIGIN)*(1-coord2);
      int b_sign = b_sign_pointing_to_ridge(v1, v2, v3, r1, r2, RED);
      if (b_sign == -1) { if (r_type == RED) return RE; else return RC;}
      if (b_sign == 1) return BH;
    }
  assert(0);
}

void
xing_on_edge(Halfedge_handle he, bool& is_crossed, Ridge_type color)
{
  is_crossed = false;
  FT sign;
  FT b_p, b_q; // extremalities at p and q for he: p->q
  Vector_3  d_p = he->opposite()->vertex()->d1(),
    d_q = he->vertex()->d1(); //ppal dir
  if ( color == BLUE ) {
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



//for a ridge segment in a triangle [r1,r2], let r = r2 - r1 and normalize,
//the projection of a point p on the line (r1,r2) is pp=r1+tr, with t=(p-r1)*r
//then the vector v starting at p is pointing to the ridge line (r1,r2) if 
//(pp-p)*v >0
int 
b_sign_pointing_to_ridge(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
			 Vector_3 r1, Vector_3 r2, Ridge_type color)
{
  Vector_3 r = r2 - r1, dv1, dv2, dv3;
  FT bv1, bv2, bv3;
  if ( color == BLUE ) {
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
  if ( r !=0 ) r = r/CGAL::sqrt(r*r);
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












CGAL_END_NAMESPACE

#endif
