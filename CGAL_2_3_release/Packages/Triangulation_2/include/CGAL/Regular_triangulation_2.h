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
// file          : include/CGAL/Regular_triangulation_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Frederic Fichel, Mariette Yvinec, Julia Floetotto
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_2_H
#define CGAL_REGULAR_TRIANGULATION_2_H

#include <CGAL/Triangulation_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>

CGAL_BEGIN_NAMESPACE 

template < class Gt, 
           class Tds  = Triangulation_data_structure_using_list_2 <
                        Triangulation_vertex_base_2<Gt>,
		        Regular_triangulation_face_base_2<Gt> > >
class Regular_triangulation_2 : public Triangulation_2<Gt,Tds>
{
public:
  typedef Triangulation_2<Gt,Tds>     Triangulation;
  typedef Gt                          Geom_traits;
  typedef typename Gt::Bare_point     Point;
  typedef typename Gt::Weighted_point Weighted_point;
  typedef typename Gt::Weight         Weight;

  typedef typename Triangulation::Face_handle    Face_handle;
  typedef typename Triangulation::Vertex_handle  Vertex_handle;
  typedef typename Triangulation::Edge           Edge;
  typedef typename Triangulation::Locate_type    Locate_type;
  typedef typename Triangulation::Face_circulator       Face_circulator;
  typedef typename Triangulation::Edge_circulator       Edge_circulator;
  typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Triangulation::Finite_faces_iterator Finite_faces_iterator;
  typedef typename Triangulation::Finite_vertices_iterator 
                                                     Finite_vertices_iterator;
  // a list to memorise temporary the faces around a point
  typedef std::list<Face_handle>      Faces_around_stack; 
  // point list
  typedef std::list<Weighted_point>   Weighted_point_list;

public:
  Regular_triangulation_2(const Gt& gt=Gt()) : Triangulation(gt) {}

  Regular_triangulation_2(const Regular_triangulation_2 &rt)
    : Triangulation(rt)
  {   CGAL_triangulation_postcondition( is_valid() );  }

 
  // CHECK - QUERY
  Oriented_side power_test(const Weighted_point &p,
			   const Weighted_point &q,
			   const Weighted_point &r,
			   const Weighted_point &s) const;
  Oriented_side power_test(const Weighted_point &p,
			   const Weighted_point &q,
			   const Weighted_point &r) const;
  Oriented_side power_test(const Face_handle &f, 
			   const Weighted_point &p) const;
  Oriented_side power_test(const Face_handle& f, int i,
			   const Weighted_point &p) const;
  
  bool is_valid(bool verbose = false, int level = 0) const;
  void affiche_tout();	

  // DUAL
  Point dual (Face_handle f) const;
  Object dual(const Edge &e) const ;
  Object dual(const Edge_circulator& ec) const;
  Object dual(const Finite_edges_iterator& ei) const;
  Point weighted_circumcenter(Face_handle f) const; 
  Point weighted_circumcenter(const Weighted_point& p0, 
			      const Weighted_point& p1, 
			      const Weighted_point& p2) const;

  // Insertion, Deletion and Flip
  Vertex_handle push_back(const Weighted_point &p);
  Vertex_handle insert(const Weighted_point &p, 
		       Face_handle f = Face_handle() );
  Vertex_handle insert(const Weighted_point &p,
	 	       Locate_type  lt,
		       Face_handle loc, int li );
  Vertex_handle insert_in_face(const Weighted_point &p, Face_handle f);
  Vertex_handle insert_in_edge(const Weighted_point &p, Face_handle f, int i);
  void flip(Face_handle f, int i);
  void remove_degree_3(const Vertex_handle v, Face_handle f = Face_handle());
  void remove(Vertex_handle v);


private:
  void regularize(Vertex_handle v);
  void remove_2D(Vertex_handle v);
  void fill_hole_regular(std::list<Edge> & hole);
  void update_hidden_points_3_1(const Face_handle& f1, const Face_handle& f2, 
				const Face_handle& f3);
  void update_hidden_points_2_2(const Face_handle& f1, const Face_handle& f2);
  void update_hidden_points_1_3(const Face_handle& f1, const Face_handle& f2,
				const Face_handle& f3);
  void hide_vertex(const Face_handle& f, const Weighted_point& p);

  void stack_flip(Vertex_handle v, Faces_around_stack &faces_around);
  void stack_flip_4_2(Face_handle f, int i, int j, 
		      Faces_around_stack &faces_around);
  void stack_flip_3_1(Face_handle f, int i, int j,
		      Faces_around_stack &faces_around);
  void stack_flip_2_2(Face_handle f, int i, 
		      Faces_around_stack &faces_around);
  void stack_flip_dim1(Face_handle f, int i);
		       


public:
  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
  {
      int n = number_of_vertices();
      while(first != last)
      {
          insert(*first);
  	  ++first;
      }
      return number_of_vertices() - n;
  }

  template < class Stream>
  Stream& draw_dual(Stream & ps) const
    {
      Finite_edges_iterator eit= finite_edges_begin();
      for (; eit != finite_edges_end(); ++eit) {
	Object o = dual(eit);
	typename Geom_traits::Line_2  l;
	typename Geom_traits::Ray_2   r;
	typename Geom_traits::Segment_2 s;
	if (CGAL::assign(s,o)) ps << s;
	if (CGAL::assign(r,o)) ps << r;
	if (CGAL::assign(l,o)) ps << l;
      }
      return ps;
    }
    
};




template < class Gt, class Tds >
Oriented_side
Regular_triangulation_2<Gt,Tds>::
power_test(const Face_handle &f, const Weighted_point &p) const
{
  // p is supposed to be a finite point
  // if f is a finite face, 
  // return  ON_NEGATIVE_SIDE if p is above f 
  // (p has to be hidden)
  int i;
  if ( ! f->has_vertex(infinite_vertex(), i) )
    return power_test(f->vertex(0)->point(),
		      f->vertex(1)->point(),
		      f->vertex(2)->point(),p);

  Orientation o = orientation(f->vertex(ccw(i))->point(),
			      f->vertex( cw(i))->point(),
			      p);
  if (o==COLLINEAR)
    return power_test(f->vertex(ccw(i))->point(),
		      f->vertex( cw(i))->point(),p);

  return Oriented_side(o);
}

template < class Gt, class Tds >
Oriented_side
Regular_triangulation_2<Gt,Tds>::
power_test(const Face_handle& f, int i,
	   const Weighted_point &p) const
{
  // f is supposed to be a finite edge
  // p is supposed to be on  edge (f,i)
  // return ON_NEGATIVE_SIDE if p is above (f,i)
  // (p has to be hidden)
  CGAL_triangulation_precondition (!is_infinite(f,i) &&
	     orientation(f->vertex(ccw(i))->point(),
			 f->vertex( cw(i))->point(),
			 p) == COLLINEAR);
  return  power_test(f->vertex(ccw(i))->point(),
		     f->vertex( cw(i))->point(),
		     p);
}

template < class Gt, class Tds >
inline
Oriented_side
Regular_triangulation_2<Gt,Tds>::
power_test(const Weighted_point &p,
	   const Weighted_point &q,
	   const Weighted_point &r,
	   const Weighted_point &s) const
{
  return geom_traits().power_test_2_object()(p,q,r,s);
}

template < class Gt, class Tds >
inline
Oriented_side
Regular_triangulation_2<Gt,Tds>::
power_test(const Weighted_point &p,
	   const Weighted_point &q,
	   const Weighted_point &r) const
{
  return geom_traits().power_test_degenerated_2_object()(p,q,r);
}

template < class Gt, class Tds >
bool
Regular_triangulation_2<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  if (number_of_vertices() <= 1)    return true; 
  bool result = Triangulation_2<Gt,Tds>::is_valid(verbose, level);
 		
  for(Finite_faces_iterator it = finite_faces_begin(); 
      it != finite_faces_end(); it++)
  {
    for(int i=0; i<3; i++)
    {
      if (!is_infinite(it->vertex(i)))
	  result = result && ON_POSITIVE_SIDE != 
	           power_test(it->neighbor(i), it->vertex(i)->point());
//       if (!result)
//       {
//         std::cerr << "face : " << (void*)&(*it)<< "  " 
// 		  <<"["<< it->vertex(0)->point()
// 		  <<"/"<< it->vertex(1)->point()
// 		  <<"/"<< it->vertex(2)->point()<<"]"<< std::endl
// 		  << "voisin : " << (void*)&(*(it->neighbor(i)))<< "  "
// 		  <<"["<<(it->neighbor(i))->vertex(0)->point()
// 		  <<"/"<<(it->neighbor(i))->vertex(1)->point()
// 		  <<"/"<<(it->neighbor(i))->vertex(2)->point()<<"]" << std::endl;
//       }
      CGAL_triangulation_assertion(result);
    }
    
    typename Weighted_point_list::iterator plit = it->point_list().begin(),
                                         pldone = it->point_list().end();
    for (; plit != pldone; plit++)
    {
      result = result && power_test(it, *plit) == ON_NEGATIVE_SIDE;
    	
      if (!result)
      {
	std::cerr << "face : " << (void*)&(*it)<< "  " 
		  <<"["<< it->vertex(0)->point()
		  <<"/"<< it->vertex(1)->point()
		  <<"/"<< it->vertex(2)->point()<<"]"	<< std::endl
		  << "hidden point : " << *plit << std::endl;
      }
      CGAL_triangulation_assertion(result); 
    }
  }
  return result;
}


template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
affiche_tout()
{
  std::cerr<< "AFFICHE TOUTE LA TRIANGULATION :"<<std::endl;
  Finite_faces_iterator fi , fi_end=finite_faces_end();

  std::cerr << std::endl<<"====> "<<this<<std::endl;
  fi=finite_faces_begin();
  std::cerr<<"***"<<std::endl;
  while(fi != fi_end)
    {	
      std::cerr << "face : "<<(void*)&(*fi)<<" => "<<std::endl;
      std::cerr <<"point :"<<(fi->vertex(0)->point())
	   <<" / voisin "<<&(*(fi->neighbor(0)))
	   <<"["<<(fi->neighbor(0))->vertex(0)->point()
	   <<"/"<<(fi->neighbor(0))->vertex(1)->point()
	   <<"/"<<(fi->neighbor(0))->vertex(2)->point()<<"]"
	   <<std::endl;
      std::cerr <<"point :"<<(fi->vertex(1)->point())
	   <<" / voisin "<<&(*(fi->neighbor(1)))
	   <<"["<<(fi->neighbor(1))->vertex(0)->point()
	   <<"/"<<(fi->neighbor(1))->vertex(1)->point()
	   <<"/"<<(fi->neighbor(1))->vertex(2)->point()<<"]"
	   <<std::endl;
      std::cerr <<"point :"<<(fi->vertex(2)->point())
	   <<" / voisin "<<&(*(fi->neighbor(2)))
	   <<"["<<(fi->neighbor(2))->vertex(0)->point()
	   <<"/"<<(fi->neighbor(2))->vertex(1)->point()
	   <<"/"<<(fi->neighbor(2))->vertex(2)->point()<<"]"
	   <<std::endl;

      typename Weighted_point_list::iterator current;
      std::cerr << "  +++++>>>    ";
      for (current= fi->point_list().begin() ; 
	   current!= fi->point_list().end() ; current++ )
	{
	  std::cerr <<"[ "<< (*current) << " ] ,  ";
	}
      std::cerr <<std::endl;
      ++fi;
    }
  std::cerr <<"faces infinies "<<std::endl;
			
  if ( number_of_vertices() <= 2) {return;}
  Face_circulator fc = 
    infinite_vertex()->incident_faces(),fcdone(fc);
  do {	
    std::cerr<<(void*)&(*fc) <<" = "<< fc->vertex(0)->point()<<" / "
	<< fc->vertex(1)->point()<<" / "<< fc->vertex(2)->point()
	<<" / ";
    typename Weighted_point_list::iterator current;
    std::cerr << "  +++++>>>    ";
    for (current= fc->point_list().begin() ; 
	 current!= fc->point_list().end() ; current++ )
      {
	std::cerr <<"[ "<< (*current) << " ] ,  ";
      }
    std::cerr <<std::endl;
  }while(++fc != fcdone);

  std::cerr <<std::endl;

  if (number_of_vertices()>1) {
    std::cerr << "affichage des sommets de la triangulation reguliere"
	      <<std::endl;
    Finite_vertices_iterator vi;
    vi=finite_vertices_begin();
    if (number_of_vertices()>1) {
      while ( vi!=vertices_end() ) {
	//std::cerr << "* "<< &(*vi) <<"  /  ";
	std::cerr << "* "<< vi->point() <<"  / face associee : "
	     << (void*)(&(*(vi->face())))<<std::endl;;
	++vi;
      }
      std::cerr<<std::endl;
    }
  }
}


//DUALITY
template < class Gt, class Tds >
inline
Regular_triangulation_2<Gt,Tds>::Point
Regular_triangulation_2<Gt,Tds>::
dual (Face_handle f) const
{
  return weighted_circumcenter(f);
}

template < class Gt, class Tds >
inline
Regular_triangulation_2<Gt,Tds>::Point
Regular_triangulation_2<Gt,Tds>::
weighted_circumcenter(Face_handle f) const
{
  CGAL_triangulation_precondition (dimension()==2 || !is_infinite(f));
  return weighted_circumcenter(f->vertex(0)->point(),
			       f->vertex(1)->point(),
			       f->vertex(2)->point());
}

template<class Gt, class Tds>
inline
Regular_triangulation_2<Gt,Tds>::Point
Regular_triangulation_2<Gt,Tds>::
weighted_circumcenter(const Weighted_point& p0, 
		      const Weighted_point& p1, 
		      const Weighted_point& p2) const
{
  return 
    geom_traits().construct_weighted_circumcenter_2_object()(p0,p1,p2);
}

template < class Gt, class Tds >
inline
Object
Regular_triangulation_2<Gt,Tds>::
dual(const Edge &e) const
{
  typedef typename Geom_traits::Line_2        Line;
  typedef typename Geom_traits::Ray_2         Ray;
  typedef typename Geom_traits::Direction_2   Direction;
  typedef typename Geom_traits::Segment       Segment;
  
  CGAL_triangulation_precondition (! is_infinite(e));
  if( dimension()== 1 ){
    Weighted_point p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point q = (e.first)->vertex(ccw(e.second))->point();
    Line l  = geom_traits().construct_radical_axis_2_object()(p,q);
    return Object(new Wrapper< Line >(l));
  }
  
  // dimension==2
  if( (! is_infinite(e.first)) &&
      (! is_infinite(e.first->neighbor(e.second))) ) {
    Segment s = geom_traits().construct_segment_2_object()
      (dual(e.first),dual(e.first->neighbor(e.second)));
    return CGAL::Object(new CGAL::Wrapper< Segment >(s));
  } 

  // one of the adjacent face is infinite
  Face_handle f; int i;
  if ( is_infinite(e.first)) {
    f=e.first->neighbor(e.second); f->has_neighbor(e.first,i);
  } 
  else {
    f=e.first; i=e.second;
  }
  Weighted_point p = f->vertex( cw(i))->point();
  Weighted_point q = f->vertex( ccw(i))->point();
  Line l  = geom_traits().construct_radical_axis_2_object()(p,q);
  Direction d =
    geom_traits().construct_direction_of_line_2_object()(l);
  Ray r = geom_traits().construct_ray_2_object()(dual(f), d);
  return CGAL::Object(new CGAL::Wrapper< Ray >(r));
}
  

template < class Gt, class Tds >
inline 
Object
Regular_triangulation_2<Gt,Tds>::  
dual(const Edge_circulator& ec) const
{
  return dual(*ec);
}
  
template < class Gt, class Tds >
inline 
Object
Regular_triangulation_2<Gt,Tds>::
dual(const Finite_edges_iterator& ei) const
{
  return dual(*ei);
}

//INSERTION-REMOVAL
template < class Gt, class Tds >
Regular_triangulation_2<Gt,Tds>::Vertex_handle
Regular_triangulation_2<Gt,Tds>::
push_back(const Weighted_point &p)
{	
    return insert(p);
}

template < class Gt, class Tds >
Regular_triangulation_2<Gt,Tds>::Vertex_handle
Regular_triangulation_2<Gt,Tds>::
insert(const Weighted_point &p, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate(p, lt, li, start);
  return insert(p, lt, loc, li);
}

template < class Gt, class Tds >
Regular_triangulation_2<Gt,Tds>::Vertex_handle
Regular_triangulation_2<Gt,Tds>::
insert(const Weighted_point &p, Locate_type lt, Face_handle loc, int li) 
{
  if (number_of_vertices() <= 1) return Triangulation::insert(p);
 
  Vertex_handle v;
  Oriented_side os;
  switch (lt) {
  case VERTEX:
    remove(loc->vertex(li));
    return insert(p);
  case FACE:
    if (power_test(loc,p) == ON_NEGATIVE_SIDE) {
      hide_vertex(loc,p);
      return v;
    }
    v = insert_in_face(p,loc);
    break;
  case EDGE:
    os = dimension() == 1 ?  power_test(loc,li,p) : 
                             power_test(loc,p);
    if (os == ON_NEGATIVE_SIDE) {
      hide_vertex(loc, p);
      return v;
    }
    v = insert_in_edge(p,loc,li);
    break;
  case OUTSIDE_CONVEX_HULL:
    v = insert_outside_convex_hull(p,loc);
    break;
  case OUTSIDE_AFFINE_HULL:
    v =  insert_outside_affine_hull(p);
    break;
  default:
    CGAL_triangulation_assertion_msg(false, "locate step failed");
  }
  regularize(v);
  return v;
}

template < class Gt, class Tds >
Regular_triangulation_2<Gt,Tds>::Vertex_handle
Regular_triangulation_2<Gt,Tds>::
insert_in_face(const Weighted_point &p, Face_handle f)
{
  Vertex_handle v = Triangulation::insert_in_face(p,f);
  update_hidden_points_1_3(f, 
			   f->neighbor(ccw(f->index(v))), 
			   f->neighbor( cw(f->index(v))) );
  return v;
}

template < class Gt, class Tds >
Regular_triangulation_2<Gt,Tds>::Vertex_handle
Regular_triangulation_2<Gt,Tds>::
insert_in_edge(const Weighted_point &p, Face_handle f, const int i)
{
  Vertex_handle v = Triangulation::insert_in_edge(p,f,i);
  if (dimension() == 1) {
    Face_handle g = f->neighbor(1 - f->index(v));
    update_hidden_points_2_2(f,g);
  }
  else { //dimension()==2
    Face_handle g = (v==f->vertex(cw(i))) ? f->neighbor(ccw(i))
      : f->neighbor( cw(i));
    update_hidden_points_2_2(f,g);
    update_hidden_points_2_2(f->neighbor(i), g->neighbor(i));
  }
  return v;
} 

template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
regularize(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != infinite_vertex());
  Faces_around_stack faces_around;

  //initialise faces_around
  if (dimension() == 1) {
    faces_around.push_back(v->face());
    faces_around.push_back(v->face()->neighbor(1- v->face()->index(v)));
  }
  else{ //dimension==2
    Face_circulator fit = v->incident_faces(), done(fit);
    do {
      faces_around.push_back(fit++);
    } while(fit != done);
  }

  while( ! faces_around.empty() )
    stack_flip(v, faces_around);
}


template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
flip(Face_handle f, int i)
{
  Face_handle n = f->neighbor(i);
  Triangulation_2<Gt,Tds>::flip(f,i);
  update_hidden_points_2_2(f,n);
}

template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
remove_degree_3(const Vertex_handle v, Face_handle f) 
{
  if (f == Face_handle())
    f=v->face();
  update_hidden_points_3_1(f, f->neighbor( cw(f->index(v))),
			   f->neighbor(ccw(f->index(v))));
  Triangulation_2<Gt,Tds>::remove_degree_3(v,f);
}

template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
remove(Vertex_handle v )
{
  CGAL_triangulation_precondition(!v.is_null());
  CGAL_triangulation_precondition(!is_infinite(v));

  //if (number_of_vertices() <= 1) Triangulation::remove(v);
   
  // Collect in p_list
  // the points hidden by the face to be deleted
  Weighted_point_list p_list;
  if (dimension() == 1) {
    Face_handle f = v->face();
    Face_handle n = f->neighbor(1 - f->index(v));
    p_list.splice(p_list.begin(), f->point_list());
    p_list.splice(p_list.begin(), n->point_list());
  }
  else if (dimension() == 2 ) {
    Face_circulator fc = v->incident_faces(),done(fc);
    do {
      p_list.splice(p_list.begin(), fc->point_list());
      fc++;
    }  
    while( fc != done);
  }

  if (dimension() <= 1) Triangulation::remove(v);
  else remove_2D(v);
   
  Weighted_point p;
  while (! p_list.empty())
  {
    p=p_list.front();
    p_list.pop_front();
    insert(p);
  }
}

template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
remove_2D(Vertex_handle v)
{
  if (test_dim_down(v)) {  _tds.remove_dim_down(&(*v));  }
  else {
    std::list<Edge> hole;
    make_hole(v, hole);
    fill_hole_regular(hole);
    delete &(*v);
    set_number_of_vertices(number_of_vertices()-1);
  }
  return;   
}


template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
fill_hole_regular(std::list<Edge> & first_hole)
{
  typedef std::list<Edge> Hole;
  typedef std::list<Hole> Hole_list;
  
  Hole hole;
  Hole_list hole_list;
  Face_handle ff, fn;
  int i, ii, in;
	
  hole_list.push_front(first_hole);
  
  while (! hole_list.empty())
    {
      hole = hole_list.front();
      hole_list.pop_front();
      typename Hole::iterator hit = hole.begin();
	    
      // if the hole has only three edges, create the triangle
      if (hole.size() == 3)
	{
	  Face_handle  newf = create_face();
	  hit = hole.begin();
	  for(int j=0; j<3; j++)
	    {
	      ff = (*hit).first;
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
      bool finite = false;
      while (!finite)
	{
	  ff = hole.front().first;
	  ii = hole.front().second;
	  if ( is_infinite(ff->vertex(cw(ii))) ||
	       is_infinite(ff->vertex(ccw(ii))))
	    {
	      hole.push_back(hole.front());
	      hole.pop_front();
	    }
	  else
	    finite = true;
	}
 
      // take the first neighboring face and pop it;
      ff = hole.front().first;
      ii = hole.front().second;
      hole.pop_front();
 
      Vertex_handle  v0 = ff->vertex(ff->cw(ii)); 
      Weighted_point p0 = v0->point();
      Vertex_handle  v1 = ff->vertex(ff->ccw(ii)); 
      Weighted_point p1 = v1->point();
      Vertex_handle  v2 = infinite_vertex(); 
      Weighted_point p2;
      Vertex_handle  vv;
      Weighted_point p;
 
      typename Hole::iterator hdone = hole.end();
      hit = hole.begin();
      typename Hole::iterator cut_after(hit);
 
      // if tested vertex is c with respect to the vertex opposite
      // to NULL neighbor,
      // stop at the before last face;
      hdone--;
      while (hit != hdone) 
	{
	  fn = (*hit).first;
	  in = (*hit).second;
	  vv = fn->vertex(ccw(in));
	  if (is_infinite(vv))
	    {
	      if (is_infinite(v2))
		cut_after = hit;
	    }
	  else 
	    {	// vv is a finite vertex
	      p = vv->point();
	      if (orientation(p0,p1,p) == 
		  COUNTERCLOCKWISE)
		{
		  if (is_infinite(v2))
		    {
		      v2=vv;
		      p2=p;
		      cut_after=hit;
		    }
		  else if (power_test(p0,p1,p2,p) == 
			   ON_POSITIVE_SIDE)
		    {
		      v2=vv;
		      p2=p;
		      cut_after=hit;
		    }
		}
	    }
	  ++hit;
	}
 
      // create new triangle and update adjacency relations
      Face_handle newf = create_face(v0,v1,v2);
      newf->set_neighbor(2,ff);
      ff->set_neighbor(ii, newf);
 
      //update the hole and push back in the Hole_List stack
      // if v2 belongs to the neighbor following or preceding *f
      // the hole remain a single hole
      // otherwise it is split in two holes
 
      fn = hole.front().first;
      in = hole.front().second;
      if (fn->has_vertex(v2, i) && i == (int)fn->ccw(in)) 
	{
	  newf->set_neighbor(0,fn);
	  fn->set_neighbor(in,newf);
	  hole.pop_front();
	  hole.push_front(Edge(&(*newf),1));
	  hole_list.push_front(hole);
	}
      else
	{
	  fn = hole.back().first;
	  in = hole.back().second;
	  if (fn->has_vertex(v2, i) && i == (int)fn->cw(in)) 
	    {
	      newf->set_neighbor(1,fn);
	      fn->set_neighbor(in,newf);
	      hole.pop_back();
	      hole.push_back(Edge(&(*newf),0));
	      hole_list.push_front(hole);
	    }
	  else
	    { // split the hole in two holes
	      Hole new_hole;
	      ++cut_after;
	      while (hole.begin() != cut_after)
		{
		  new_hole.push_back(hole.front());
		  hole.pop_front();
		}
 
	      hole.push_front(Edge(&(*newf),1));
	      new_hole.push_front(Edge(&(*newf),0));
	      hole_list.push_front(hole);
	      hole_list.push_front(new_hole);
	    }
	}
    }
}


// add the point_list of f2 and f3 to the point list of f1
// for the 3-1 flip
template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
update_hidden_points_3_1(const Face_handle& f1, const Face_handle& f2, 
			 const Face_handle& f3)
{
  (f1->point_list()).splice(f1->point_list().begin(),f2->point_list());
  (f1->point_list()).splice(f1->point_list().begin(),f3->point_list());
}


// the points of the lists of 2 faces are sorted
// because of a 2-2 flip
template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
update_hidden_points_2_2(const Face_handle& f1, const Face_handle& f2)
{	
  CGAL_triangulation_assertion(f1->has_neighbor(f2));
    
  Weighted_point_list p_list;
  p_list.splice(p_list.begin(),f1->point_list());
  p_list.splice(p_list.begin(),f2->point_list());

  // if one of the face is infinite, 
  // the other face hide all the points
  if ( is_infinite(f1)) {
    (f2->point_list()).splice(f2->point_list().begin(),p_list);
    return;
  }
  if ( is_infinite(f2)) {
    (f1->point_list()).splice(f1->point_list().begin(),p_list);
    return;
  }

  if (dimension() == 1) {
    Weighted_point a1 = f1->vertex(f1->index(f2))->point();
    Weighted_point a2 = f2->vertex(f2->index(f1))->point();
    Weighted_point a  = f1->vertex(1-f1->index(f2))->point();
    while ( ! p_list.empty() ) {
      if ( compare_x(a, p_list.front()) == 
	   compare_x(a, a1)  &&
	   compare_y(a, p_list.front()) == 
	   compare_y(a, a1))
	(f1->point_list()).push_back(p_list.front());
      else
	(f2->point_list()).push_back(p_list.front());
      p_list.pop_front();
    }
    return;
  }

  // from here f1 and f2 are finite 2-dimensional faces
  int idx2 = f1->index(f2);
  Vertex_handle v0=f1->vertex(ccw(idx2));
  Vertex_handle v1=f1->vertex(cw(idx2));
  CGAL_triangulation_assertion( !is_infinite(v0) && !is_infinite(v1)); 

  while ( ! p_list.empty() )
    {
      if (orientation(v0->point(), v1->point(), p_list.front()) ==
	  COUNTERCLOCKWISE)
	(f1->point_list()).push_back(p_list.front());
      else
	(f2->point_list()).push_back(p_list.front());
      p_list.pop_front();
    }
}
	  
// The point list of f1 is separated into 3 lists
// for a 1-3 flip
template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
update_hidden_points_1_3(const Face_handle& f1, const Face_handle& f2, 
			 const Face_handle& f3)
{
    CGAL_triangulation_assertion(f1->has_neighbor(f2) &&
				 f2->has_neighbor(f3) &&
				 f3->has_neighbor(f1));

    Weighted_point_list p_list;
    p_list.splice(p_list.begin(),f1->point_list());
    p_list.splice(p_list.begin(),f1->point_list());
    p_list.splice(p_list.begin(),f1->point_list());
    if (p_list.empty())
	return;

    // the following does not work if 
    // two of f1,f2 and f3 are twice neighbors
    // but this cannot appear taking the assertion into account;
    int idx2 = f1->index(f2),
        idx3 = f1->index(f3);
    Vertex_handle v2 = f1->vertex(idx2),
                  v3 = f1->vertex(idx3),
                  v0 = f1->vertex(3-(idx2+idx3)),
                  v1 = f2->vertex(f2->index(f1));

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
    while(! p_list.empty())
    {
      if(orientation(v2->point(),v0->point(),p_list.front()) !=
	 orientation(v2->point(),v0->point(),v3->point()) )
      { // not in f1
	if (orientation(v1->point(), v0->point(), p_list.front() ) !=
	    orientation(v1->point(), v0->point(), v3->point() ) )
	  // not in f2
	    f3->point_list().push_back(p_list.front());
	else
	    f2->point_list().push_back(p_list.front());
      }
      else
	  f1->point_list().push_back(p_list.front());
      p_list.pop_front();
    }
}

// insert the point among the hidden points list
template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
hide_vertex(const Face_handle& f, const Weighted_point& p)
{
  f->point_list().push_back(p);
}


template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
stack_flip(Vertex_handle v, Faces_around_stack &faces_around)
{
  Face_handle f=faces_around.front();
  faces_around.pop_front();
  int i = f->index(v);
  Face_handle n = f->neighbor(i);
    
  // TODO : the 1dim-dimensional case
  if (dimension() == 1 ) {
    if ( is_infinite(f)  || is_infinite(n) ) return;
    if ( power_test( v->point(),
		     f->vertex(1-i)->point(),
		     n->vertex(n->index(f))->point()) ==
	 ON_NEGATIVE_SIDE) return;
    stack_flip_dim1(f,i);
    return;
  }  

  // now dimension() == 2
  //test the regularity of edge (f,i)
  if( power_test(n, v->point()) == ON_NEGATIVE_SIDE)
    return;
    
  if(is_infinite(f,i))
    {
      int j = 3 - ( i + f->index(infinite_vertex()));
      if ( f->vertex(j)->degree() == 4)
	stack_flip_4_2(f,i,j,faces_around);
      return;
    }

    // now f and n are both finite faces
    int ni = n->index(f);
    Orientation occw = orientation(f->vertex(i)->point(),
				   f->vertex(ccw(i))->point(),
				   n->vertex(ni)->point());
    Orientation ocw  = orientation(f->vertex(i)->point(),
				   f->vertex(cw(i))->point(),
				   n->vertex(ni)->point());
    if (occw == LEFTTURN && ocw == RIGHTTURN) {
      // quadrilater (f,n) is convex
      stack_flip_2_2(f,i, faces_around);
      return;
    }
    if (occw == RIGHTTURN && f->vertex(ccw(i))->degree() == 3) {
      stack_flip_3_1(f,i,ccw(i),faces_around);
      return;
    }
    if (ocw == LEFTTURN && f->vertex(cw(i))->degree() == 3) {
      stack_flip_3_1(f,i,cw(i),faces_around);
      return;
    }
    if (occw == COLLINEAR && f->vertex(ccw(i))->degree() == 4) {
      stack_flip_4_2(f,i,ccw(i),faces_around);
      return;
    }
    if (ocw == COLLINEAR && f->vertex(cw(i))->degree() == 4)
      stack_flip_4_2(f,i,cw(i),faces_around);
}


template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
stack_flip_4_2(Face_handle f, int i, int j, Faces_around_stack & faces_around)
{
    int k = 3-(i+j);
    Face_handle g=f->neighbor(k);
    if (!faces_around.empty())
    {
      if (faces_around.front() == g)
	  faces_around.pop_front();
      else if (faces_around.back() == g) 
	  faces_around.pop_back();
    }
    
    //union f with  g and f->neihgbor(i) with g->f->neihgbor(i)
    Face_handle fn = f->neighbor(i);
    //Face_handle gn = g->neighbor(g->index(f->vertex(i)));
    Vertex_handle vq = f->vertex(j);
    
    _tds.flip( &(*f), i); //not using flip because the vertex j is flat.
    update_hidden_points_2_2(f,fn);
    Face_handle h1 = ( f->has_vertex(vq) ? fn : f);
    hide_vertex(h1, vq->point());
    remove_degree_3(vq,g);
    faces_around.push_front(g);
    faces_around.push_front(h1);    
}


template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
stack_flip_3_1(Face_handle f, int i, int j, Faces_around_stack & faces_around)
{
  int k = 3-(i+j);
  Face_handle g=f->neighbor(k);
  if (!faces_around.empty())
  {
    if (faces_around.front()== g)
	faces_around.pop_front();
    else if ( faces_around.back() == g)
	faces_around.pop_back();
  }

  Vertex_handle vq= f->vertex(j);
  hide_vertex(f,vq->point());
  remove_degree_3(vq,f);
  faces_around.push_front(f);
}


template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
stack_flip_2_2(Face_handle f, int i, Faces_around_stack & faces_around)
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
  
template < class Gt, class Tds >
void
Regular_triangulation_2<Gt,Tds>::
//stack_flip_dim1(Face_handle f, int i, Faces_around_stack &faces_around)
stack_flip_dim1(Face_handle f, int i)
{
  Face_handle n= f->neighbor(i);
  int in = n->index(f);
  f->set_vertex(1-i, n->vertex(in));
  f->set_neighbor(i, n->neighbor(1-in));
  n->neighbor(1-in)->set_neighbor(n->neighbor(1-in)->index(n), f);
  (f->point_list()).splice(f->point_list().end(),n->point_list());
  delete_face(n);
}
CGAL_END_NAMESPACE 

#endif // CGAL_REGULAR_TRIANGULATION_2_H
