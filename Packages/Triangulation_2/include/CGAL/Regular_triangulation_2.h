// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : Triangulation/include/CGAL/Regular_triangulation_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Frederic Fichel, Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_2_H
#define CGAL_REGULAR_TRIANGULATION_2_H

#include <pair.h>
#include <list.h>
#include <map.h>
#include <assert.h>
#include <stack.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Weighted_point_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>




template < class Gt, class Tds >
class CGAL_Regular_triangulation_2 
  : public CGAL_Triangulation_2<Gt,Tds>
{

public:
  typedef CGAL_Triangulation_2<Gt,Tds> Triangulation;

  typedef typename Gt::Bare_point  Point;
  typedef typename Gt::Weighted_point Weighted_point;
  typedef typename Gt::Weight      Weight;
  
  // a list to memorise temporary the faces around a point
  typedef list<Face_handle > Faces_around_stack; 

  // point list
  typedef list<Weighted_point> Weighted_point_list;

public:
  CGAL_Regular_triangulation_2()
    : Triangulation( )
  {}

  CGAL_Regular_triangulation_2(const Gt& gt) : Triangulation(gt) { }

  CGAL_Regular_triangulation_2(const CGAL_Regular_triangulation_2 &rt )
    : Triangulation( rt )
  {   CGAL_triangulation_postcondition( is_valid() );  }

  CGAL_Regular_triangulation_2(const Vertex_handle&  v, const Gt& gt) 
    : Triangulation(v,gt) 
  {   CGAL_triangulation_postcondition( is_valid() );  }

	
  // Insertion and Deletion

  Vertex_handle push_back(const Weighted_point &p)
  {	
    Locate_type lt;
    return insert(p, lt, NULL);
  }

	#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
	template < class InputIterator >
	int
	insert(InputIterator first, InputIterator last)
	{	int n = number_of_vertices();
		while(first != last)
		{	insert(*first);
			++first;
		}
		return number_of_vertices() - n;
	}
	#else
	#if defined(LIST_H) || defined(__SGI_STL_LIST_H)
	int
	insert(list<Weighted_point>::const_iterator first,
		list<Weighted_point>::const_iterator last)
	{	int n = number_of_vertices();
		while(first != last)
		{	insert(*first);
			++first;
		}
		return number_of_vertices() - n;
	}
	#endif // LIST_H
	#if defined(VECTOR_H) || defined(__SGI_STL_VECTOR_H)
	int
	insert(vector<Weighted_point>::const_iterator first,
		vector<Weighted_point>::const_iterator last)
	{	int n = number_of_vertices();
		while(first != last)
		{	insert(*first);
			++first;
		}
		return number_of_vertices() - n;
	}
	#endif // VECTOR_H
	#ifdef ITERATOR_H
	int
	insert(istream_iterator<Weighted_point, ptrdiff_t> first,
		istream_iterator<Weighted_point, ptrdiff_t> last)
	{	int n = number_of_vertices();
		while(first != last)
		{	insert(*first);
			++first;
		}
		return number_of_vertices() - n;
	}
	#endif // ITERATOR_H
  
	int
	insert(Weighted_point* first, Weighted_point* last)
	{	int n = number_of_vertices();
		while(first != last)
		{	insert(*first);
			++first;
		}
		return number_of_vertices() - n;
	}
	#endif // CGAL_TEMPLATE_MEMBER_FUNCTIONS

public:
  CGAL_Oriented_side
  power_test(const Face_handle& f, const Weighted_point & p) const
  {
    //p is assume to be a finite point
    CGAL_Orientation o;
    if ( ! is_infinite(f) )
      {
	return geom_traits().power_test(f->vertex(0)->point(),
					f->vertex(1)->point(),
					f->vertex(2)->point(),p);
      }

    else if ( f->vertex(0) == infinite_vertex() )
      {	o = geom_traits().orientation(f->vertex(1)->point(),
				 f->vertex(2)->point(),p);
      if (o==CGAL_COLLINEAR)
	{	return geom_traits().power_test(f->vertex(1)->point(),
			                         f->vertex(2)->point(),p);
	}
      }
 
    else if ( f->vertex(1) == infinite_vertex() )
	{	o = geom_traits().orientation(f->vertex(2)->point(),
					 f->vertex(0)->point(),p);
	if (o==CGAL_COLLINEAR)
	  {	return geom_traits().power_test(f->vertex(2)->point(),
					   f->vertex(0)->point(),p);
	  }
	} 

    else if ( f->vertex(2) == infinite_vertex() )
      {	o = geom_traits().orientation(f->vertex(0)->point(),
			                         f->vertex(1)->point(),p);
      if (o==CGAL_COLLINEAR)
	{	
	  return geom_traits().power_test(f->vertex(0)->point(),
					  f->vertex(1)->point(),p);
	}
      }

    return (o == CGAL_NEGATIVE) ? CGAL_ON_NEGATIVE_SIDE :
      (o == CGAL_POSITIVE) ? CGAL_ON_POSITIVE_SIDE : CGAL_ON_ORIENTED_BOUNDARY;
  }

  private :
  // add the point_list of f2 and f3 to the point list of f1
  // for the 3-1 flip 
  void update_hidden_points_3_1(const Face_handle& f1, 
				const Face_handle& f2, 
				const Face_handle& f3 )
  {
     (f1->point_list()).splice(f1->point_list().begin(),f2->point_list());
     (f1->point_list()).splice(f1->point_list().begin(),f3->point_list());
  }


  // the points of the lists of 2 faces are sorted
  // because of a 2-2 flip
  void update_hidden_points_2_2(const Face_handle& f1, const Face_handle& f2)
  {	
    CGAL_triangulation_assertion(f1->has_neighbor(f2));
    
    Weighted_point_list p_list;
    p_list.splice(p_list.begin(),f1->point_list());
    p_list.splice(p_list.begin(),f2->point_list());

    // if both faces are infinite faces
    //if they share the finite edge, any one hide all points
    //if they share an infinite edge : their finite edges
    //should be collinear
    //split hidden points along the midpoint
    if ( is_infinite(f1) && is_infinite(f2) ) {
      int i1 = f1->index(infinite_vertex());
      if (f2 == f1->neighbor(i1)) {
	(f1->point_list()).splice(f1->point_list().begin(),p_list);
	return;
      }
      // else split the hidden points according to finite edges
      // because f2 != f1->neighbor(i1)	
      // the following vertices are finite vertices
      Weighted_point a1 = f1->vertex(f1->index(f2))->point();
      Weighted_point a2 = f2->vertex(f2->index(f1))->point();
      Weighted_point a = f1->vertex(3-(i1+f1->index(f2)))->point();
      while ( ! p_list.empty() ){
	if ( geom_traits().compare_x(a, p_list.front()) == 
	     geom_traits().compare_x(a, a1)  &&
	     geom_traits().compare_y(a, p_list.front()) == 
	     geom_traits().compare_y(a, a1)) {
	  (f1->point_list()).push_back(p_list.front());
	}
	else {
	  (f2->point_list()).push_back(p_list.front());
	}
	 p_list.pop_front();
      }
      return;  
    }

    //otherwise at least one of the face is finite
    // if one is infinite, the finite face hide all the points
    if ( is_infinite(f1)) {
      (f2->point_list()).splice(f2->point_list().begin(),p_list);
      return;
    }
    if ( is_infinite(f2)) {
      (f1->point_list()).splice(f1->point_list().begin(),p_list);
      return;
    }

    // from here f1 and f2 are finite faces
    int idx2;
    idx2=f1->index(f2);
    Vertex_handle v2=f1->vertex(idx2),
      v0=f1->vertex(ccw(idx2)),
      v1=f1->vertex(cw(idx2));
    CGAL_triangulation_assertion((!is_infinite(v0)) && (!is_infinite(v1))); 

    while ( ! p_list.empty() ){
      if (CGAL_orientation(v0->point(), v1->point(), p_list.front()) ==
			     CGAL_COUNTERCLOCKWISE) {
	  (f1->point_list()).push_back(p_list.front());
      }
      else{(f2->point_list()).push_back(p_list.front()) ;}
      p_list.pop_front();
    }
  }
	  
  // The point list of f1 is separated into 3 lists
  // for a 1-3 flip
  void update_hidden_points_1_3(const Face_handle& f1, const Face_handle& f2, 
				const Face_handle& f3) 
  {
    
    CGAL_triangulation_assertion(f1->has_neighbor(f2) &&
				 f2->has_neighbor(f3) &&
				 f3->has_neighbor(f1));

    Weighted_point_list p_list;
    p_list.splice(p_list.begin(),f1->point_list());
    p_list.splice(p_list.begin(),f1->point_list());
    p_list.splice(p_list.begin(),f1->point_list());
    if(p_list.empty()) return;

    // the following does not work if 
    // two of f1,f2 and f3 are twice neighbors
    // but this cannot appear taking the assertion into account;
    int idx2,idx3;
    idx2=f1->index(f2);
    idx3=f1->index(f3);
    Vertex_handle v2=f1->vertex(idx2),
      v3=f1->vertex(idx3),
      v0=f1->vertex(3-(idx2+idx3)),
      v1=f2->vertex(f2->index(f1));

    CGAL_triangulation_assertion(f2->has_vertex(v0) && f1->has_vertex(v0));
    CGAL_triangulation_assertion(f3->has_vertex(v1));
    CGAL_triangulation_assertion( ! is_infinite(v0));

    // if two of f1, f2,and f3 are infinite
    // the list goes entirely to the third finite face
    // no orientation test necessary
    // because the point list of an infinite face
    // is only made of point projecting on its finite edge
    if ( is_infinite(f1 ) && is_infinite(f2)) {
      f3->point_list().splice(f3->point_list().begin(), p_list);
      return;
    }
    if ( is_infinite(f1) && is_infinite(f3)) {
      f2->point_list().splice(f2->point_list().begin(), p_list);
      return;
    }
    if ( is_infinite(f2) && is_infinite(f3)){
      f1->point_list().splice(f1->point_list().begin(), p_list);
      return;
    }
    
    // if here, v1,v2,v3 and v0 are finite vertices
    while(! p_list.empty()) {
      if(CGAL_orientation(v2->point(),v0->point(),p_list.front()) !=
	 CGAL_orientation(v2->point(),v0->point(),v3->point()) )
	{ // not in f1
	if (CGAL_orientation(v1->point(), v0->point(), p_list.front() ) !=
	    CGAL_orientation(v1->point(), v0->point(), v3->point() ) )
	  {// not in f2
	    f3->point_list().push_back(p_list.front());
	  }
	else{
	  f2->point_list().push_back(p_list.front());
	}
	}
      else { f1->point_list().push_back(p_list.front());
      }
      p_list.pop_front();
    }
  }

protected:
  bool is_degree(Vertex_handle v, int i) const
  {
    Face_circulator fic = v->incident_faces(), done(fic);
    fic++;
    int count=1;
    while ( fic != done && count <= i) {
      fic ++; count++;
    }
    return (count==i ? true :false );
  }

bool degree_equal_3(Vertex_handle v) const
{
  Face_handle f = v->face();
    if (f == NULL) 
      {	return false;
      }
    int i = f->index(v);
    Face_handle ptr1 = f->neighbor(ccw(i));
    Face_handle ptr2 = f;
    f = f->neighbor(cw(i));
    
    int count = 2;
    
    while(ptr1 != f)
      {	count++;
      i = ptr1->index(ptr2);
      ptr2 = ptr1;
      ptr1 = ptr1->neighbor(cw(i));
      if (count>3)
	break;
      }
    return (count==3 ? true :false );
  }
  
  bool degree_equal_4(Vertex_handle v) const
  {
    return ( v->degree() == 4);
  }
  
  void hide_vertex(const Face_handle& f,const Weighted_point& p)
    // insert the point among the hidden points list
  {

	 f->point_list().push_back(p);
  }



public:
void 
regularize(Vertex_handle v)
  {
    CGAL_triangulation_precondition( v != infinite_vertex());
    Faces_around_stack faces_around;

    //initialise faces_around
    Face_circulator fit = v->incident_faces(), done(fit);
    do {
      faces_around.push_back(fit++);
    } while(fit != done);

      while( ! faces_around.empty() ){
	stack_flip(v, faces_around);
      }
    return;
  }
   
private :
void 	stack_flip(Vertex_handle v, Faces_around_stack &faces_around)
  {
    Face_handle f=faces_around.front();
    int i = i=f->index(v);
    faces_around.pop_front();

    //test the regularity of edge (f,i)
    Face_handle n = f->neighbor(i);
    if( power_test(n, v->point()) == CGAL_ON_NEGATIVE_SIDE){
      return;
    }
    
    if(is_infinite(f,i)) {
      int j = 3 - ( i + f->index(infinite_vertex()));
      if ( is_degree(f->vertex(j),4)) { 
	stack_flip_4_2(f,i,j,faces_around);
      }
      return;
    }

    // now f and n are both finite faces
    int ni= n->index(f);
    CGAL_Orientation occw = geom_traits().orientation(f->vertex(i)->point(),
						  f->vertex(ccw(i))->point(),
						  n->vertex(ni)->point());
    CGAL_Orientation ocw = geom_traits().orientation(f->vertex(i)->point(),
						  f->vertex(cw(i))->point(),
						  n->vertex(ni)->point());
    if(occw == CGAL_LEFTTURN && ocw  ==  CGAL_RIGHTTURN){
      // quadrilater (f,n) is convex
      stack_flip_2_2(f,i, faces_around);
      return;
    }
    else if (occw == CGAL_RIGHTTURN && is_degree(f->vertex(ccw(i)), 3) ){
      stack_flip_3_1(f,i,ccw(i),faces_around);
      return;
    }
    else if (ocw == CGAL_LEFTTURN && is_degree(f->vertex(cw(i)), 3) ){
      stack_flip_3_1(f,i,cw(i),faces_around);
      return;
    }
    else if (occw == CGAL_COLLINEAR && is_degree(f->vertex(ccw(i)), 4) ){
      stack_flip_4_2(f,i,ccw(i),faces_around);
      return;
    }
    else if (ocw == CGAL_COLLINEAR && is_degree(f->vertex(cw(i)), 4) ){ 
      stack_flip_4_2(f,i,cw(i),faces_around);
      return;
    }
    return;
  }
      
      
void stack_flip_4_2(Face_handle f, int i, int j, 
		    Faces_around_stack & faces_around)
  {
    int k = 3-(i+j);
    Face_handle g=f->neighbor(k);
    if (!faces_around.empty()){
      if (faces_around.front()== g){ faces_around.pop_front();}
      else if ( faces_around.back() == g) { faces_around.pop_back();}
    }
    
    //union f with  g and f->neihgbor(i) with g->f->neihgbor(i)
    Face_handle fn = f->neighbor(i);
    Face_handle gn = g->neighbor(g->index(f->vertex(i)));
    Vertex_handle vq = f->vertex(j);
    
    _tds.flip( &(*f), i); //not using flip because the vertex j is flat.
    update_hidden_points_2_2(f,fn);
    Face_handle h1 = ( f->has_vertex(vq) ? fn : f);
    remove_degree_3(vq,g);
    hide_vertex(h1, vq->point());
    faces_around.push_front(g);
    faces_around.push_front(h1);    
  }

void stack_flip_3_1(Face_handle f, int i, int j,
		    Faces_around_stack & faces_around)
{
  int k = 3-(i+j);
  Face_handle g=f->neighbor(k);
  if (!faces_around.empty()){
    if (faces_around.front()== g){ faces_around.pop_front();}
    else if ( faces_around.back() == g) { faces_around.pop_back();}
  }

  Vertex_handle vq= f->vertex(j);
  remove_degree_3(vq,f);
  hide_vertex(f,vq->point());
  faces_around.push_front(f);
}



void stack_flip_2_2(Face_handle f, int i, Faces_around_stack & faces_around)
  {
    Vertex_handle vq = f->vertex(ccw(i));
    flip(f,i);
    if(f->has_vertex(vq)) {
      faces_around.push_front(f->neighbor(ccw(i)));
      faces_around.push_front(f);
    }
    else { 
      faces_around.push_front(f);
      faces_around.push_front(f->neighbor(cw(i)));
    }
  }
  



public:
   // INSERTION / DELETION
   Vertex_handle
   insert(const Weighted_point  &p,
	  CGAL_Triangulation_2<Gt,Tds>::Locate_type& lt,
	  Face_handle f = Face_handle() ) 
  {
    Vertex_handle v;
   
    if(number_of_vertices() == 0) {
      v = new Vertex(p);
      lt = OUTSIDE;
      insert_first(v);
      return v;
    }
		
    if(number_of_vertices() == 1) {
      if (geom_traits().compare(p,finite_vertex()->point())) {
	remove_first(finite_vertex());
	v= new Vertex(p);
	insert_first(v);
	return v;
      }
      else{
	v = new Vertex(p);
	lt = OUTSIDE;
	insert_second(v);	
	return v;
      }
    }
    
    int li;
    Face_handle loc = locate(p, lt, li, f);
    switch(lt){
    case VERTEX:
      remove(loc->vertex(li));
      v = insert(p,lt,f);
      return v;

    case OUTSIDE:
      v = new Vertex(p);
      insert_outside(v,loc); 
      break;
	
    case FACE:
      if (power_test( loc ,p)==CGAL_ON_NEGATIVE_SIDE) {
	hide_vertex(loc, p );
	return v;
      }
      v = new Vertex(p);
      insert_in_face(v,loc);
      break;
	
    case EDGE:
      if (power_test( loc ,p)==CGAL_ON_NEGATIVE_SIDE){
	hide_vertex(loc, p);
	return v;
      }
      v = new Vertex(p);
      insert_in_edge(v,loc,li);
      break;
    
    case COLLINEAR_OUTSIDE:
      v = new Vertex(p);
      insert_collinear_outside( v, loc,li);
      break;  
   
    default:
      CGAL_triangulation_assertion(false);  // locate step failed
    }
   
    regularize(v);
    return v;
  }



    
public:
  Vertex_handle insert(const Weighted_point &p,
		       Face_handle f = Face_handle() )
{
//cerr<<endl<<"INSERT : "<<p<<endl;
  Locate_type lt;
  return insert(p, lt, f);
}

void insert_in_face(Vertex_handle v, Face_handle f)
{
  CGAL_Triangulation_2<Gt,Tds>::insert_in_face(v,f);
  update_hidden_points_1_3(f, f->neighbor(ccw(f->index(v))), 
			   f->neighbor(cw(f->index(v))) );
  return;
}

void insert_in_edge(Vertex_handle v, Face_handle f, int i)
{
  CGAL_Triangulation_2<Gt,Tds>::insert_in_edge(v,f,i);
  Face_handle g = 
	(v==f->vertex(cw(i)) ? f->neighbor(ccw(i)) : f->neighbor(cw(i)) );
  update_hidden_points_2_2(f,g);
  update_hidden_points_2_2(f->neighbor(i), g->neighbor(i));
}    


void flip(Face_handle f,int i)
{
  Face_handle n = f->neighbor(i);
  CGAL_Triangulation_2<Gt,Tds>::flip(f,i);
  update_hidden_points_2_2(f,n);
  return;
}


void remove_degree_3(Vertex_handle v,Face_handle f = Face_handle()) 
{
  if (f == Face_handle()) f=v->face();
  update_hidden_points_3_1(f, f->neighbor(cw(f->index(v))),
			   f->neighbor(ccw(f->index(v))));
  CGAL_Triangulation_2<Gt,Tds>::remove_degree_3(v,f);
}


protected :
void remove_2D(Vertex_handle v)
      {
	// General case
  
	// remove incident faces
	// set up list of faces neighboring the hole
	// in ccw order around the hole
  
	// should be face_handle instead of void*
	// but it doesn't link (maybe because of too long names
	// generated by template classes imbrication)
	typedef pair<void* , int>  Hole_neighbor;

	typedef list<Hole_neighbor> Hole;
	typedef list<Hole> Hole_list;
  
	Hole hole;
	Hole_list hole_list;
	list<Face_handle> to_delete;

	Face_handle  f, ff, fn;
	int i =0,ii =0, in =0;
	Vertex_handle vv;
		
 	Face_circulator fc = v->incident_faces();
	Face_circulator done(fc);
	do 
	  {
	    f = (*fc).handle(); fc++;
	    i  = f->index(v);
	    fn = f->neighbor(i);
	    vv = f->vertex(f->cw(i));
	    if( vv->face()== f) vv->set_face(fn);
	    vv = f->vertex(f->ccw(i));
	    if( vv->face()== f) vv->set_face(fn);
	    in = fn->index( f );
	    fn->set_neighbor(in, NULL);
	    hole.push_back(Hole_neighbor(&(*fn),in));
	    to_delete.push_back(f);
	  }
	while(fc != done);

	while (! to_delete.empty()){
	 to_delete.front().Delete();
	 to_delete.pop_front();
       }
	
	hole_list.push_front(hole);
  
	while( ! hole_list.empty())
	  {
	    hole = hole_list.front();
	    hole_list.pop_front();
	    Hole::iterator hit = hole.begin();
	    
	    // if the hole has only three edges, create the triangle
	    if (hole.size() == 3)
	      {
		Face_handle  newf = new Face();
		hit = hole.begin();
		for(int j = 0;j<3;j++)
		  {
		    ff = (Face*)(*hit).first;
		    ii = (*hit).second;
		    hit++;
		    ff->set_neighbor(ii,newf);
		    newf->set_neighbor(j,ff);
		    newf->set_vertex(newf->ccw(j),ff->vertex(ff->cw(ii)));
		    
		  }
		continue;
	      }
  
			// else find an edge with two finite vertices
			// on the hole boundary
			// and the new triangle adjacent to that edge
			//  cut the hole and push it back
  
			// first, ensure that a neighboring face
			// whose vertices on the hole boundary are finite
			// is the first of the hole
			bool finite= false;
			while (!finite)
			{
				ff = (Face*)(hole.front()).first;
				ii = (hole.front()).second;
				if ( is_infinite(ff->vertex(cw(ii))) ||
					is_infinite(ff->vertex(ccw(ii)))) 
				{
					hole.push_back(hole.front());
					hole.pop_front();
				}
				else 
					finite=true;
			}
  
			// take the first neighboring face and pop it;
			ff = (Face*)(hole.front()).first;
			ii =(hole.front()).second;
			hole.pop_front();
  
  
			Vertex_handle v0 = ff->vertex(ff->cw(ii)); Weighted_point p0 =v0->point();
			Vertex_handle v1 = ff->vertex(ff->ccw(ii)); Weighted_point p1 =v1->point();
			Vertex_handle v2 = infinite_vertex(); Weighted_point p2;
			Vertex_handle vv; Weighted_point p;
  
			Hole::iterator hdone = hole.end();
			hit =  hole.begin();
			Hole::iterator cut_after(hit);
  
			// if tested vertex is c with respect to the vertex opposite
			// to NULL neighbor,
			// stop at the before last face;
			hdone--;
			while( hit != hdone) 
			{
				fn = (Face*)(*hit).first;
				in = (*hit).second;
				vv = fn->vertex(ccw(in));
				if (is_infinite(vv)) 
				{
					if(is_infinite(v2))
						cut_after = hit;
				}else 
				{	// vv is a finite vertex
					p = vv->point();
					if (  geom_traits().orientation(p0,p1,p) == CGAL_COUNTERCLOCKWISE)
					{	if(is_infinite(v2)) 
						{	v2=vv; p2=p;
							cut_after=hit;
						}else
						{	if( geom_traits().power_test (p0,p1,p2,p) == CGAL_ON_POSITIVE_SIDE)
							{	v2=vv; p2=p; cut_after=hit;
							}
						}
					}
				}
				++hit;
			}
  
  
			// create new triangle and update adjacency relations
			Face_handle  newf = new Face(v0,v1,v2);
			newf->set_neighbor(2,ff);
			ff->set_neighbor(ii, newf);
  
  
			//update the hole and push back in the Hole_List stack
			// if v2 belongs to the neighbor following or preceding *f
			// the hole remain a single hole
			// otherwise it is split in two holes
  
			fn = (Face*)(hole.front()).first;
			in = (hole.front()).second;
			if (fn->has_vertex(v2, i) && i == fn->ccw(in)) 
			{
				newf->set_neighbor(0,fn);
				fn->set_neighbor(in,newf);
				hole.pop_front();
				hole.push_front(Hole_neighbor(&(*newf),1));
				hole_list.push_front(hole);
			}else
			{
				fn = (Face*)(hole.back()).first;
				in = (hole.back()).second;
				if (fn->has_vertex(v2, i) && i== fn->cw(in)) 
				{
					newf->set_neighbor(1,fn);
					fn->set_neighbor(in,newf);
					hole.pop_back();
					hole.push_back(Hole_neighbor(&(*newf),0));
					hole_list.push_front(hole);
				}else
				{
					// split the hole in two holes
					Hole new_hole;
					++cut_after;
					while( hole.begin() != cut_after )
					{
						new_hole.push_back(hole.front());
						hole.pop_front();
					}
  
					hole.push_front(Hole_neighbor(&(*newf),1));
					new_hole.push_front(Hole_neighbor(&(*newf),0));
					hole_list.push_front(hole);
					hole_list.push_front(new_hole);
				}
			}
		}

	}

public:
      void  remove(Vertex_handle v )
    {
	  CGAL_triangulation_precondition(! v.is_null());
	  CGAL_triangulation_precondition( !is_infinite(v));
	  Weighted_point_list p_list;

	  if  (number_of_vertices() == 1){
	    _tds.remove_first(& (*v));
	    v.clear();
	    return;
	  }
 
	  //  take care of finite_vertex data member
	  if (finite_vertex() == v){
	    Face_handle f = v->face();
	    int i=f->index(v);
	    Vertex_handle vv= is_infinite(f->vertex(cw(i))) ? 
		f->vertex(ccw(i)) : f->vertex(cw(i));
	    set_finite_vertex( vv);
	  }
	  
	  // Collect in p_list
	  // the points hidden by the face to be deleted
	  Face_circulator fc = v->incident_faces(),done(fc);
	  do {
	    p_list.splice(p_list.begin(), fc->point_list());
	    fc++;
	  }  
	  while( fc != done);
	    
	  if (number_of_vertices() == 2) {
	    Face_handle f = v->face();
	    Face_handle ff =f->neighbor(0);
	    ff.Delete();
	    f.Delete();
	  }
	  else{
	    if ( dimension() == 1) remove_1D(v);
	    else  remove_2D(v);
	  }
    
	  v.Delete();
	  set_number_of_vertices(number_of_vertices()-1);


	  Weighted_point p;
	  while ( ! p_list.empty() )
	    {	p=p_list.front();
	    p_list.pop_front();
//cerr<<endl<<"re-insertion : "<<p<<endl;
	    insert(p);
	    }

	  return;
//cerr<<"etat APRES remove"<<endl;
//affiche_tout();

	}


	bool is_valid(bool verbose = false, int level = 0) const
	{
	  if(number_of_vertices() <= 1) {
	    return true;
	  }
  
	  bool result = CGAL_Triangulation_2<Gt,Tds>::is_valid();
  		
	  for( Face_iterator it = faces_begin(); it != faces_end() ; it++) {

	    for(int i=0; i<3; i++) {
	    if 	( ! is_infinite( it->vertex(i))) {
	      result = result &&
	       CGAL_ON_POSITIVE_SIDE != 
		power_test( it->neighbor(i), it->vertex(i)->point());
	    }
	    if ( !result) {
	      cerr << "face : " << (void*)&(*it)<< "  " 
		   <<"["<< it->vertex(0)->point()
		   <<"/"<< it->vertex(1)->point()
		   <<"/"<< it->vertex(2)->point()<<"]"	<<endl;
	      cerr << "voisin : " << (void*)&(*(it->neighbor(i)))<< "  "
		   <<"["<<(it->neighbor(i))->vertex(0)->point()
		   <<"/"<<(it->neighbor(i))->vertex(1)->point()
		   <<"/"<<(it->neighbor(i))->vertex(2)->point()<<"]" <<endl;
	    }
	    CGAL_triangulation_assertion( result );
	    }
	    
	    Weighted_point_list::iterator plit= it->point_list().begin(),
	      pldone= it->point_list().end();
	    for( ; plit != pldone ; plit++) {
	      result = result &&
		power_test( it, *plit) == CGAL_ON_NEGATIVE_SIDE ;
	    	
	      if ( !result) {
		cerr << "face : " << (void*)&(*it)<< "  " 
		     <<"["<< it->vertex(0)->point()
		     <<"/"<< it->vertex(1)->point()
		     <<"/"<< it->vertex(2)->point()<<"]"	<<endl;
		cerr << "hidden point : " << *plit <<endl;
	      }
	     CGAL_triangulation_assertion( result ); 
	    }
	
	    return result;
	  }
	}


void affiche_tout()
{
  cerr<< "AFFICHE TOUTE LA TRIANGULATION :"<<endl;
		Face_iterator fi , fi_end=faces_end();

		cerr << endl<<"====> "<<this<<endl;
		fi=faces_begin();
cerr<<"***"<<endl;
		while(fi != fi_end)
		{	
			cerr << "face : "<<(void*)&(*fi)<<" => "<<endl;
			cerr <<"point :"<<(fi->vertex(0)->point())
					<<" / voisin "<<&(*(fi->neighbor(0)))
				<<"["<<(fi->neighbor(0))->vertex(0)->point()
				<<"/"<<(fi->neighbor(0))->vertex(1)->point()
			<<"/"<<(fi->neighbor(0))->vertex(2)->point()<<"]"
					<<endl;
			cerr <<"point :"<<(fi->vertex(1)->point())
				<<" / voisin "<<&(*(fi->neighbor(1)))
				<<"["<<(fi->neighbor(1))->vertex(0)->point()
				<<"/"<<(fi->neighbor(1))->vertex(1)->point()
			<<"/"<<(fi->neighbor(1))->vertex(2)->point()<<"]"
				<<endl;
			cerr <<"point :"<<(fi->vertex(2)->point())
				<<" / voisin "<<&(*(fi->neighbor(2)))
				<<"["<<(fi->neighbor(2))->vertex(0)->point()
				<<"/"<<(fi->neighbor(2))->vertex(1)->point()
			    <<"/"<<(fi->neighbor(2))->vertex(2)->point()<<"]"
				<<endl;

			Weighted_point_list::iterator current;
			cerr << "  +++++>>>    ";
			for (current= fi->point_list().begin() ; 
			     current!= fi->point_list().end() ; current++ )
			  {
					cerr <<"[ "<< (*current) << " ] ,  ";
			  }
			cerr <<endl;
			++fi;
		}
 cerr <<"faces infinies "<<endl;
			
 if ( number_of_vertices() <= 2) {return;}
Face_circulator fc = infinite_vertex()->incident_faces(),fcdone(fc);
 do {	
   cerr<<(void*)&(*fc) <<" = "<< fc->vertex(0)->point()<<" / "
       << fc->vertex(1)->point()<<" / "<< fc->vertex(2)->point()
       <<" / ";
   Weighted_point_list::iterator current;
   cerr << "  +++++>>>    ";
   for (current= fc->point_list().begin() ; 
	current!= fc->point_list().end() ; current++ )
     {
       cerr <<"[ "<< (*current) << " ] ,  ";
     }
   cerr <<endl;
 }while(++fc != fcdone);


 cerr <<endl;

 if (number_of_vertices()>1) {
   cerr << "affichage des sommets de la triangulation reguliere"<<endl;
   Vertex_iterator vi;
   vi=vertices_begin();
   if (number_of_vertices()>1) {
     while ( vi!=vertices_end() ) {
       //cerr << "* "<< &(*vi) <<"  /  ";
       cerr << "* "<< vi->point() <<"  / face associee : "
	    << (void*)(&(*(vi->face())))<<endl;;
       ++vi;
     }
     cerr<<endl;
   }
 }

}


};



#endif
