// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three ki '/u/alcor/0/prisme_util/CGAL/CGAL-1.0/make/makefile_sparc_SunOS-5.5_eg++-egcs-2.90.27_LEDA'
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the CGAL Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the CGAL Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Max-Planck-Institute Saarbrucken
// (Germany), RISC Linz (Austria), and Tel-Aviv University (Israel).
//
// ============================================================================
//
// release       : CGAL-1.0
// date          : 21 Apr 1998
//file          : include/CGAL/Delaunay_triangulation_2.h
// author(s)     : Olivier Devillers
//                 Andreas Fabri
//                 Monique Teillaud
//                 Mariette Yvinec
//
// email         : cgal@cs.uu.nl
//
// ============================================================================


#ifndef CGAL_DELAUNAY_TRIANGULATION_2_H
#define CGAL_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/Triangulation_2.h>

template < class Gt, class Tds>
class CGAL_Delaunay_triangulation_2 : public CGAL_Triangulation_2<Gt,Tds>
{
friend istream& operator>> CGAL_NULL_TMPL_ARGS
                (istream& is, CGAL_Delaunay_triangulation_2<Gt,Tds> &T);
public:
  typedef Gt Geom_traits;
  typedef typename Geom_traits::Distance Distance;

  CGAL_Delaunay_triangulation_2()
    : CGAL_Triangulation_2<Gt,Tds>() {}
  
  CGAL_Delaunay_triangulation_2(const Gt& gt)
  : CGAL_Triangulation_2<Gt,Tds>(gt) {}
  
  
  CGAL_Delaunay_triangulation_2(const Vertex_handle& v, const Gt& gt=Gt())
    : CGAL_Triangulation_2<Gt,Tds>(v,gt)
  {   CGAL_triangulation_postcondition( is_valid() );  }
  
  
  // copy constructor duplicates vertices and faces
  // no copy constructed needed :
  // the default will invoke the copy constructor of CGAL_Triangulation_2 ?
  CGAL_Delaunay_triangulation_2(const CGAL_Delaunay_triangulation_2<Gt,Tds> &tr)
      : CGAL_Triangulation_2<Gt,Tds>(tr)
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
  insert(list<Point>::const_iterator first,
         list<Point>::const_iterator last)
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
  insert(vector<Point>::const_iterator first,
         vector<Point>::const_iterator last)
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
  insert(istream_iterator<Point, ptrdiff_t> first,
         istream_iterator<Point, ptrdiff_t> last)
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
  #endif // CGAL_TEMPLATE_MEMBER_FUNCTIONS

  CGAL_Oriented_side
  side_of_oriented_circle(const Face_handle& f, const Point & p) const
  {
      CGAL_Orientation o;
      if ( ! is_infinite(f) ) {
          return geom_traits().side_of_oriented_circle(f->vertex(0)->point(),
                                                  f->vertex(1)->point(),
                                                  f->vertex(2)->point(),p);
      } else if ( f->vertex(0) == infinite_vertex() ) {
          o = geom_traits().orientation(f->vertex(1)->point(),
                                f->vertex(2)->point(),p);
      } else if ( f->vertex(1) == infinite_vertex() ) {
          o = geom_traits().orientation(f->vertex(2)->point(),
                                f->vertex(0)->point(),p);
      } else if ( f->vertex(2) == infinite_vertex() ) {
          o = geom_traits().orientation(f->vertex(0)->point(),
                                f->vertex(1)->point(),p);
      }
      return (o == CGAL_NEGATIVE) ? CGAL_ON_NEGATIVE_SIDE :
          (o == CGAL_POSITIVE) ? CGAL_ON_POSITIVE_SIDE :
          CGAL_ON_ORIENTED_BOUNDARY;
  }

  Vertex_handle
  nearest_vertex(const Point& p, const Face_handle& f) const
  {
      Vertex_handle nn;
      Distance closer(p,&geom_traits());
      int min;
      int i;
  
      if (number_of_vertices() == 0) return NULL;
      if (number_of_vertices() == 1) return finite_vertex();
  
      i = ( ! is_infinite(f->vertex(0)) ) ? 0 : 1;
  
      closer.set_point(1,f->vertex(i)->point());
      min = 1;
      nn = f->vertex(i);
      if ( ! is_infinite(f->vertex(ccw(i)))){
         closer.set_point( 3-min, f->vertex(ccw(i))->point() );
          if (  ( (min==1)? CGAL_LARGER : CGAL_SMALLER )
                == closer.compare() ) {
              min = 3-min;
              nn=f->vertex(ccw(i));
          }
      }
      if ( ! is_infinite(f->vertex(cw(i)))){
          closer.set_point( 3-min, f->vertex(cw(i))->point() );
          if (  ( (min==1)? CGAL_LARGER : CGAL_SMALLER )
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
  
  inline Vertex_handle
  nearest_vertex(const Point  &p) const
  {
      Face_handle f = locate(p);
      return nearest_vertex(p, f);
  }

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
         CGAL_Triangulation_2<Gt,Tds>::Locate_type& lt,
         Face_handle f = Face_handle() )
  {
      Vertex_handle v = CGAL_Triangulation_2<Gt,Tds>::insert(p,lt,f);

      if(dimension() <= 1) return v;

      f=v->face();
      Face_handle next;
      int i;
      Face_handle start(f);
      do {
          i = f->index(v);
          next = f->neighbor(ccw(i));  // turn ccw around v
          propagating_flip(f,i);
          f=next;
      } while(next != start);
      return v;
  }


  void  remove(Vertex_handle& v )
  {
     CGAL_triangulation_precondition(v != NULL);
     CGAL_triangulation_precondition( !is_infinite(v));
  
       if  (number_of_vertices() == 1) {
        CGAL_Triangulation_2<Gt,Tds>::remove(v);
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
  
     if (number_of_vertices() == 2) {
      Face_handle f = v->face();
      Face_handle ff = f->neighbor(0);
      ff.Delete();
      f.Delete();
    }
    else{
      if ( dimension() == 1) remove_1D(v);
      else  remove_2D(v);
    }
    v.Delete();
    set_number_of_vertices(number_of_vertices()-1);
    return;
   }

  bool is_valid(bool verbose = false, int level = 0) const
  {
       if(number_of_vertices() <= 1){
            return true;
       }

    bool result = CGAL_Triangulation_2<Gt,Tds>::is_valid();

    Face_iterator it = faces_begin(),
                        done = faces_end();
          while(it != done){
               if ( ! is_infinite( it->neighbor(0)) ){
                  result = result &&
                      (CGAL_ON_POSITIVE_SIDE !=
                       side_of_oriented_circle(it->neighbor(0), it->vertex(0)->point()) );
	       }
	       if (!result)
		{
		  	cerr << "face : "<<(void*)&(*it)<<" => "<<endl;
			cerr <<"point :"<<(it->vertex(0)->point())<<" / voisin "<<&(*(it->neighbor(0)))
				<<"["<<(it->neighbor(0))->vertex(0)->point()
				<<"/"<<(it->neighbor(0))->vertex(1)->point()
				<<"/"<<(it->neighbor(0))->vertex(2)->point()<<"]"
					<<endl;
		}
              if ( ! is_infinite(it->neighbor(1)) ) {
                  result = result && (CGAL_ON_POSITIVE_SIDE !=
                    side_of_oriented_circle( it->neighbor(1), it->vertex(1)->point()) );
	      } 
	      if (!result)
		{
		  	cerr << "face : "<<(void*)&(*it)<<" => "<<endl;
			cerr <<"point :"<<(it->vertex(1)->point())<<" / voisin "<<&(*(it->neighbor(1)))
				<<"["<<(it->neighbor(1))->vertex(0)->point()
				<<"/"<<(it->neighbor(1))->vertex(1)->point()
				<<"/"<<(it->neighbor(1))->vertex(2)->point()<<"]"
				<<endl;
		}
              if ( ! is_infinite(it->neighbor(2)) ){
                  result = result && (CGAL_ON_POSITIVE_SIDE !=
                      side_of_oriented_circle( it->neighbor(2), it->vertex(2)->point()));
              }
	      if (!result)
		{
		  cerr << "face : "<<(void*)&(*it)<<" => "<<endl;
		  cerr <<"point :"<<(it->vertex(2)->point())<<" / voisin "<<&(*(it->neighbor(2)))
				<<"["<<(it->neighbor(2))->vertex(0)->point()
				<<"/"<<(it->neighbor(2))->vertex(1)->point()
				<<"/"<<(it->neighbor(2))->vertex(2)->point()<<"]"
				<<endl;
		}
	      CGAL_triangulation_assertion( result );
	      ++it;
          }
      return result;
  }
  
private:
  void
  look_nearest_neighbor(const Point& p,
                        const Face_handle& f,
                        int i,
                        int& min,
                        Vertex_handle& nn,
                        Distance& closer)const
  {
      Face_handle  ni=f->neighbor(i);
      if ( CGAL_ON_POSITIVE_SIDE != side_of_oriented_circle(ni,p) ) {
          return;
      }
      i = ni->index(f);
      if ( ! is_infinite(ni->vertex(i))){
          closer.set_point( 3-min, ni->vertex(i)->point() );
          if (  ( (min==1)? CGAL_LARGER : CGAL_SMALLER )
                == closer.compare() ) {
              min = 3-min;
              nn=ni->vertex(i);
          }
      }
      // recursive exploration of triangles whose circumcircle contains p
      look_nearest_neighbor(p, ni, ccw(i), min, nn, closer);
      look_nearest_neighbor(p, ni, cw(i),  min, nn, closer);
  }

  void
  propagating_flip(Face_handle& f,int i)
  {
      Face_handle n = f->neighbor(i);
      
      if ( CGAL_ON_POSITIVE_SIDE != 
	   side_of_oriented_circle(n,  f->vertex(i)->point()) ) {          
	return;
      }
      flip(f, i);
      propagating_flip(f,i);
      i = n->index(f->vertex(i));
      propagating_flip(n,i);
  }


  void remove_2D(Vertex_handle& v )
  {
    // General case
  
    // remove incident faces
    // set up list of faces neighboring the hole
    // in ccw order around the hole
  
    // problem with gcc link
    typedef pair<void *, int> Hole_neighbor;
    //typedef pair<Face_handle  , int>  Hole_neighbor;
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
     do {
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
          if (hole.size() == 3) {
          Face_handle  newf = new Face();
          hit = hole.begin();
          for(int j = 0;j<3;j++){
            ff = (Face *) ((*hit).first);
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
       while (!finite){
          ff = (Face *) ((hole.front()).first);
          ii = (hole.front()).second;
          if ( is_infinite(ff->vertex(cw(ii))) ||
               is_infinite(ff->vertex(ccw(ii)))) {
          hole.push_back(hole.front());
          hole.pop_front();
          }
          else finite=true;
        }
  
        // take the first neighboring face and pop it;
        ff = (Face *) ((hole.front()).first);
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
          fn = (Face *) ((*hit).first);
          in = (*hit).second;
          vv = fn->vertex(ccw(in));
          if (is_infinite(vv)) {
            if(is_infinite(v2)) cut_after = hit;
          }
          else {     // vv is a finite vertex
            p = vv->point();
            if (  geom_traits().orientation(p0,p1,p) == CGAL_COUNTERCLOCKWISE) {
              if(is_infinite(v2)) { v2=vv; p2=p; cut_after=hit;}
              else{
                if( geom_traits().side_of_oriented_circle (p0,p1,p2,p) ==
                    CGAL_ON_POSITIVE_SIDE){
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
  
        fn = (Face *) ((hole.front()).first);
        in = (hole.front()).second;
        if (fn->has_vertex(v2, i) && i == fn->ccw(in)) {
          newf->set_neighbor(0,fn);
          fn->set_neighbor(in,newf);
          hole.pop_front();
          hole.push_front(Hole_neighbor( &(*newf),1));
          hole_list.push_front(hole);
        }
        else{
          fn = (Face *) ((hole.back()).first);
          in = (hole.back()).second;
          if (fn->has_vertex(v2, i) && i== fn->cw(in)) {
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
  
};


template < class Gt, class Tds >
ostream&
operator<<(ostream& os, const CGAL_Delaunay_triangulation_2<Gt,Tds> &DT)
{
  return os << (const CGAL_Triangulation_2<Gt,Tds>&)DT;
}


template < class Gt, class Tds >
istream&
operator>>(istream& is, CGAL_Delaunay_triangulation_2<Gt,Tds> &dt)
{
  is >> (CGAL_Triangulation_2<Gt,Tds>&)dt;
  dt.is_valid();
  return is;
}

#endif // CGAL_DELAUNAY_TRIANGULATION_2_H
