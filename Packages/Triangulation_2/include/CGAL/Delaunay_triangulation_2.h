// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-12 $
// release_date  : $CGAL_Date: 1999/04/28 $
//
// file          : include/CGAL/Delaunay_triangulation_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_DELAUNAY_TRIANGULATION_2_H
#define CGAL_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/Triangulation_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
class Delaunay_triangulation_2 : public Triangulation_2<Gt,Tds>
{
friend std::istream& operator>> CGAL_NULL_TMPL_ARGS
                (std::istream& is, Delaunay_triangulation_2<Gt,Tds> &T);
public:
  typedef Gt Geom_traits;
  typedef typename Geom_traits::Distance Distance;
  typedef typename Geom_traits::Ray Ray;
  typedef typename Geom_traits::Line Line;
  typedef typename Geom_traits::Direction Direction;

  Delaunay_triangulation_2()
    : Triangulation_2<Gt,Tds>() {}
  
  Delaunay_triangulation_2(const Gt& gt)
  : Triangulation_2<Gt,Tds>(gt) {}
  
  
  Delaunay_triangulation_2(const Vertex_handle& v, const Gt& gt=Gt())
    : Triangulation_2<Gt,Tds>(v,gt)
  {   CGAL_triangulation_postcondition( is_valid() );  }
  
  
  // copy constructor duplicates vertices and faces
  // no copy constructed needed :
  // the default will invoke the copy constructor of Triangulation_2 ?
  Delaunay_triangulation_2(
	       const Delaunay_triangulation_2<Gt,Tds> &tr)
      : Triangulation_2<Gt,Tds>(tr)
  {   CGAL_triangulation_postcondition( is_valid() );  }
  
  #ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
  {
      int n = number_of_vertices();
      while(first != last){
          insert(*first);
          ++first;
      }
      return number_of_vertices() - n;
  }
  #else
  #if defined(LIST_H) || defined(__SGI_STL_LIST_H)
  int
  insert(std::list<Point>::const_iterator first,
         std::list<Point>::const_iterator last)
  {
      int n = number_of_vertices();
      while(first != last){
          insert(*first);
          ++first;
      }
      return number_of_vertices() - n;
  }
  #endif // LIST_H
  #if defined(VECTOR_H) || defined(__SGI_STL_VECTOR_H)
  int
  insert(std::vector<Point>::const_iterator first,
         std::vector<Point>::const_iterator last)
  {
      int n = number_of_vertices();
      while(first != last){
          insert(*first);
          ++first;
      }
      return number_of_vertices() - n;
  }
  #endif // VECTOR_H
  #ifdef ITERATOR_H
  int
  insert(std::istream_iterator<Point, std::ptrdiff_t> first,
         std::istream_iterator<Point, std::ptrdiff_t> last)
  {
      int n = number_of_vertices();
      while(first != last){
          insert(*first);
          ++first;
      }
      return number_of_vertices() - n;
  }
  #endif // ITERATOR_H
  
  int
  insert(Point* first,
         Point* last)
  {
      int n = number_of_vertices();
      while(first != last){
          insert(*first);
          ++first;
      }
      return number_of_vertices() - n;
   }
  #endif // TEMPLATE_MEMBER_FUNCTIONS

  Oriented_side
  side_of_oriented_circle(const Face_handle& f, const Point & p) const;
  
  Vertex_handle
  nearest_vertex(const Point& p, const Face_handle& f= Face_handle()) const;
   
  Vertex_handle
  insert(const Point  &p,
         Face_handle f = Face_handle() )
  {
      Locate_type lt;
      return insert(p,lt,f);
  }
  
  Vertex_handle
  push_back(const Point &p)
  {
      Locate_type lt;
      return insert(p, lt, NULL);
  }
  
  Vertex_handle
  insert(const Point  &p,
         Triangulation_2<Gt,Tds>::Locate_type& lt,
         Face_handle f = Face_handle() )
  {
      Vertex_handle v = Triangulation_2<Gt,Tds>::insert(p,lt,f);
      restore_Delaunay(v);
      return(v);
  }

private:
  void restore_Delaunay(Vertex_handle v);
  Vertex_handle nearest_vertex_2D(const Point& p, const Face_handle& f);
  Vertex_handle nearest_vertex_1D(const Point& p);

  public :
  void  remove(Vertex_handle v );
  
  bool is_valid(bool verbose = false, int level = 0) const;
  
  
private:
  void
  look_nearest_neighbor(const Point& p,
                        const Face_handle& f,
                        int i,
                        int& min,
                        Vertex_handle& nn,
                        Distance& closer) const;
 

  void
  propagating_flip(Face_handle& f,int i);
  void 
  remove_2D(Vertex_handle v );
  

  public:
  Point dual (const Face_handle &f) const
  {
    return geom_traits().circumcenter(f->vertex(0)->point(),
                   f->vertex(1)->point(),
                   f->vertex(2)->point());
  }
  
  Object dual(const Edge &e) const ;
  
  
  Object dual(const Edge_circulator& ec) const
  {
      return dual(*ec);
  }
  
  Object dual(const Edge_iterator& ei) const
  {
      return dual(*ei);
  }
  
};

template < class Gt, class Tds >
Oriented_side
Delaunay_triangulation_2<Gt,Tds>::
side_of_oriented_circle(const Face_handle& f, const Point & p) const
{
    //Orientation o;
      if ( ! is_infinite(f) ) {
          return geom_traits().side_of_oriented_circle(f->vertex(0)->point(),
                                                  f->vertex(1)->point(),
                                                  f->vertex(2)->point(),p);
      }

      int i = f->index(infinite_vertex());
      Orientation o =
	geom_traits().orientation(f->vertex(ccw(i))->point(),
				  f->vertex(cw(i))->point(),
				  p);
						     
      return (o == NEGATIVE) ? ON_NEGATIVE_SIDE :
          (o == POSITIVE) ? ON_POSITIVE_SIDE :
          ON_ORIENTED_BOUNDARY;
  }

template < class Gt, class Tds >
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex_2D(const Point& p, const Face_handle& f) const
  {
    CGAL_triangulation_precondition(dimension() == 2);
      Vertex_handle nn;
      if (f== Face_handle()) f = locate(p);
      Distance closer(p,&geom_traits());
      int min;
      int i;
  
      i = ( ! is_infinite(f->vertex(0)) ) ? 0 : 1;
  
      closer.set_point(1,f->vertex(i)->point());
      min = 1;
      nn = f->vertex(i);
      if ( ! is_infinite(f->vertex(ccw(i)))){
         closer.set_point( 3-min, f->vertex(ccw(i))->point() );
          if (  ( (min==1)? LARGER : SMALLER )
                == closer.compare() ) {
              min = 3-min;
              nn=f->vertex(ccw(i));
          }
      }
      if ( ! is_infinite(f->vertex(cw(i)))){
          closer.set_point( 3-min, f->vertex(cw(i))->point() );
          if (  ( (min==1)? LARGER : SMALLER )
                == closer.compare() ) {
              min = 3-min;
              nn=f->vertex(cw(i));
          }
      }
      look_nearest_neighbor(p,f,0,min,nn,closer);
      look_nearest_neighbor(p,f,1,min,nn,closer);
      look_nearest_neighbor(p,f,2,min,nn,closer);
      return nn;
  }

template < class Gt, class Tds >
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex_1D(const Point& p) const
    {
      Vertex_handle nn;
      Distance closer(p,&geom_traits());
      int min;

      Vertex_iterator vit=vertices_begin();
      closer.set_point(1,vit->point());
      min = 1;
      nn = vit->handle();
      do {
	closer.set_point( 3-min, (++vit)->point());
	if (  ( (min==1)? CGAL_LARGER : CGAL_SMALLER )
                == closer.compare() ) {
              min = 3-min;
              nn=vit->handle();
          }
      }while( vit != vertices_end());
      return nn;
    }
  
template < class Gt, class Tds >
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex(const Point  &p, const Face_handle& f) const
  {
    switch (dimension()) {
    case 0:
      if (number_of_vertices() == 0) return NULL;
      if (number_of_vertices() == 1) return finite_vertex();
      break;
    case 1:
      return nearest_vertex_1D(p);
      break;      
    case 2:
      return nearest_vertex_2D(p,f);
      break;
    }
    return NULL;
  }
  


template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
restore_Delaunay(Vertex_handle v)
{
      if(dimension() <= 1) return;

      Face_handle f=v->face();
      Face_handle next;
      int i;
      Face_handle start(f);
      do {
          i = f->index(v);
          next = f->neighbor(ccw(i));  // turn ccw around v
          propagating_flip(f,i);
          f=next;
      } while(next != start);
      return;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove(Vertex_handle v )
{
    //CGAL_triangulation_precondition(v != (CGAL_NULL_TYPE) NULL);
    CGAL_triangulation_precondition(! v.is_null());
     CGAL_triangulation_precondition( !is_infinite(v));
        
     if ( dimension() <= 1) CGAL_Triangulation_2<Gt,Tds>::remove(v);
     else  remove_2D(v);
        
     return;
  }

template < class Gt, class Tds >
bool
Delaunay_triangulation_2<Gt,Tds>::
is_valid(bool verbose = false, int level = 0) const
{
       if(number_of_vertices() <= 1){
            return true;
       }

    bool result = Triangulation_2<Gt,Tds>::is_valid(verbose, level);

    for( Face_iterator it = faces_begin(); it != faces_end() ; it++) {

      for(int i=0; i<3; i++) {
	if ( ! is_infinite( it->vertex(i))) {
	  result = result &&
	    ON_POSITIVE_SIDE != 
	    side_of_oriented_circle( it->neighbor(i), it->vertex(i)->point());
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
    }
    return result;
  }

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
look_nearest_neighbor(const Point& p,
                        const Face_handle& f,
                        int i,
                        int& min,
                        Vertex_handle& nn,
                        Distance& closer) const
 {
      Face_handle  ni=f->neighbor(i);
      if ( ON_POSITIVE_SIDE != side_of_oriented_circle(ni,p) ) {
          return;
      }
      i = ni->index(f);
      if ( ! is_infinite(ni->vertex(i))){
          closer.set_point( 3-min, ni->vertex(i)->point() );
          if (  ( (min==1)? LARGER : SMALLER )
                == closer.compare() ) {
              min = 3-min;
              nn=ni->vertex(i);
          }
      }
      // recursive exploration of triangles whose circumcircle contains p
      look_nearest_neighbor(p, ni, ccw(i), min, nn, closer);
      look_nearest_neighbor(p, ni, cw(i),  min, nn, closer);
  }

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
propagating_flip(Face_handle& f,int i)
{
      Face_handle n = f->neighbor(i);
      
      if ( ON_POSITIVE_SIDE != 
	   side_of_oriented_circle(n,  f->vertex(i)->point()) ) {          
	return;
      }
      flip(f, i);
      propagating_flip(f,i);
      i = n->index(f->vertex(i));
      propagating_flip(n,i);
  }


template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_2D(Vertex_handle v )
{
    // General case
  
    // remove incident faces
    // set up list of faces neighboring the hole
    // in ccw order around the hole
  
    // problem with gcc link
    typedef std::pair<void *, int> Hole_neighbor;
    //typedef pair<Face_handle  , int>  Hole_neighbor;
    typedef std::list<Hole_neighbor> Hole;
    typedef std::list<Hole> Hole_list;
  
    Hole hole;
    Hole_list hole_list;
    std::list<Face_handle> to_delete;
  
    Face_handle  f, ff, fn;
    int i =0,ii =0, in =0;
    Hole_list hole_list;
    Hole hole;
        
    hole_list.push_front(first_hole);
  
    while( ! hole_list.empty())
      {
        hole = hole_list.front();
        hole_list.pop_front();
        Hole::iterator hit = hole.begin();
  
        // if the hole has only three edges, create the triangle
          if (hole.size() == 3) {
	    Face_handle  newf = new Face();
	    hit = hole.begin();
	    for(int j = 0;j<3;j++) {
	      ff = (*hit).first;
	      ii = (*hit).second;
	      hit++;
	      ff->set_neighbor(ii,newf);
	      newf->set_neighbor(j,ff);
	      newf->set_vertex(newf->ccw(j),ff->vertex(ff->cw(ii)));
	    }
	  }
  
        // else find an edge with two finite vertices
        // on the hole boundary
        // and the new triangle adjacent to that edge
       //  cut the hole and push it back
  
          // first, ensure that a neighboring face
          // whose vertices on the hole boundary are finite
          // is the first of the hole
       bool finite= false;
       while (!finite){
          ff = (hole.front()).first;
          ii = (hole.front()).second;
          if ( is_infinite(ff->vertex(cw(ii))) ||
               is_infinite(ff->vertex(ccw(ii)))) {
          hole.push_back(hole.front());
          hole.pop_front();
          }
          else finite=true;
        }
  
        // take the first neighboring face and pop it;
        ff = (hole.front()).first;
        ii =(hole.front()).second;
        hole.pop_front();
  
  
        Vertex_handle v0 = ff->vertex(ff->cw(ii)); Point p0 =v0->point();
        Vertex_handle v1 = ff->vertex(ff->ccw(ii)); Point p1 =v1->point();
        Vertex_handle v2 = infinite_vertex(); Point p2;
        Vertex_handle vv; Point p;
  
        Hole::iterator hdone = hole.end();
        hit =  hole.begin();
        Hole::iterator cut_after(hit);
  
        // if tested vertex is c with respect to the vertex opposite
        // to NULL neighbor,
        // stop at the before last face;
        hdone--;
        while( hit != hdone) {
          fn = (*hit).first;
          in = (*hit).second;
          vv = fn->vertex(ccw(in));
          if (is_infinite(vv)) {
            if(is_infinite(v2)) cut_after = hit;
          }
          else {     // vv is a finite vertex
            p = vv->point();
            if (geom_traits().orientation(p0,p1,p) == COUNTERCLOCKWISE) {
              if(is_infinite(v2)) { v2=vv; p2=p; cut_after=hit;}
              else{
                if( geom_traits().side_of_oriented_circle (p0,p1,p2,p) ==
                    ON_POSITIVE_SIDE){
               v2=vv; p2=p; cut_after=hit;}
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
  
        fn = (hole.front()).first;
        in = (hole.front()).second;
        if (fn->has_vertex(v2, i) && i == (int) fn->ccw(in)) {
          newf->set_neighbor(0,fn);
          fn->set_neighbor(in,newf);
          hole.pop_front();
          hole.push_front(Hole_neighbor( &(*newf),1));
          hole_list.push_front(hole);
        }
        else{
          fn = (hole.back()).first;
          in = (hole.back()).second;
          if (fn->has_vertex(v2, i) && i== (int) fn->cw(in)) {
            newf->set_neighbor(1,fn);
            fn->set_neighbor(in,newf);
            hole.pop_back();
            hole.push_back(Hole_neighbor(&(*newf),0));
            hole_list.push_front(hole);
          }
          else{
            // split the hole in two holes
            Hole new_hole;
            ++cut_after;
            while( hole.begin() != cut_after )
            {
              new_hole.push_back(hole.front());
              hole.pop_front();
            }
  
            hole.push_front(Hole_neighbor( &(*newf),1));
            new_hole.push_front(Hole_neighbor( &(*newf),0));
            hole_list.push_front(hole);
            hole_list.push_front(new_hole);
          }
        }
      }
  }

template < class Gt, class Tds >
Object
Delaunay_triangulation_2<Gt,Tds>::
dual(const Edge &e) const
 {
    if( (!is_infinite(e.first)) 
	&& (!is_infinite(e.first->neighbor(e.second))) ) {
      	Segment s(dual(e.first),dual(e.first->neighbor(e.second)));
      	return Object(new Wrapper< Segment >(s));
    } 
    // at least one of the adjacent face is infinite
    Face_handle f; int i;
    if (is_infinite(e.first)) {
      f=e.first->neighbor(e.second); f->has_neighbor(e.first,i);
    } 
    else {
      f=e.first; i=e.second;
    }
    Line l = geom_traits().bisector(segment(f,i)).opposite();

    if (! is_infinite(f)) {
      Ray r(dual(f),l.direction());
      return Object(new Wrapper< Ray >(r));
    } 
    // both adjacent faces are infinite
    return Object(new Wrapper< Line >(l));
  }
  

template < class Gt, class Tds >
std::ostream&
operator<<(std::ostream& os, const Delaunay_triangulation_2<Gt,Tds> &DT)
{
  return os << (const Triangulation_2<Gt,Tds>&)DT;
}


template < class Gt, class Tds >
std::istream&
operator>>(std::istream& is, Delaunay_triangulation_2<Gt,Tds> &dt)
{
  is >> (Triangulation_2<Gt,Tds>&)dt;
  dt.is_valid();
  return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_2_H
