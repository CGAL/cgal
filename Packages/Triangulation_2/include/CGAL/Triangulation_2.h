#ifndef CGAL_TRIANGULATION_2_H
#define CGAL_TRIANGULATION_2_H

#include <list.h>
#include <vector.h>
#include <map.h> 
#include <algo.h>
#include <pair.h>
#include <CGAL/Pointer.h>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_face_2.h>
#include <CGAL/Triangulation_vertex_2.h>
#include <CGAL/Triangulation_handles_2.h>
#include <CGAL/Triangulation_iterators_2.h>
#include <CGAL/Triangulation_circulators_2.h>




template < class Gt, class Tds>
class CGAL_Triangulation_face_iterator_2;

template < class Gt, class Tds>
class CGAL_Triangulation_vertex_iterator_2;

template < class Gt, class Tds>
class CGAL_Triangulation_edge_iterator_2;


template < class Gt, class Tds >
class CGAL_Triangulation_2
{
  friend istream& operator>> CGAL_NULL_TMPL_ARGS
                (istream& is, CGAL_Triangulation_2<Gt,Tds> &tr);
  friend ostream& operator<< CGAL_NULL_TMPL_ARGS
                (ostream& os, CGAL_Triangulation_2<Gt,Tds> &tr);
  friend CGAL_Triangulation_face_iterator_2<Gt,Tds>;
  friend CGAL_Triangulation_edge_iterator_2<Gt,Tds>;
  friend CGAL_Triangulation_vertex_iterator_2<Gt,Tds>;

public:
  typedef Tds Triangulation_data_structure;
  typedef CGAL_Triangulation_2<Gt,Tds> Triangulation;

  typedef Gt  Geom_traits;
  typedef typename Geom_traits::Point Point;
  typedef typename Geom_traits::Segment Segment;
  typedef typename Geom_traits::Triangle Triangle;

  //  typedef typename Tds::Vertex Ve;
  //  typedef typename Tds::Face Fa;
  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Gt,Tds> Vertex;

  //typedef typename Vertex::Vertex_handle        Vertex_handle;
  //typedef CGAL_Pointer<Vertex>  Vertex_handle;
  //typedef typename Face::Face_handle            Face_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef CGAL_Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef pair<Face_handle, int>                Edge;

  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;

  typedef CGAL_Triangulation_face_iterator_2<Gt,Tds>   Face_iterator;
  typedef CGAL_Triangulation_edge_iterator_2<Gt,Tds>   Edge_iterator;
  typedef CGAL_Triangulation_vertex_iterator_2<Gt,Tds> Vertex_iterator;

  class Line_face_circulator;

  enum Locate_type {VERTEX=0, EDGE, FACE, OUTSIDE, COLLINEAR_OUTSIDE};

private:
  Tds _tds;
  Gt  _gt;

public:

// CONSTRUCTORS
  CGAL_Triangulation_2()
    :_tds()
  {}

  CGAL_Triangulation_2(const Geom_traits& geom_traits) 
    : _tds(), _gt(geom_traits)
  {}

  CGAL_Triangulation_2(const Vertex_handle&  v, 
		       const Geom_traits& geom_traits=Geom_traits())
    : _tds(&(*v)), _gt(geom_traits)
  { }

//   CGAL_Triangulation_2(Vertex_handle&  v, const Geom_traits& geom_traits=Geom_traits())
//     : _tds( &(*v)), _gt(geom_traits)
//   {}

  // copy constructor duplicates vertices and faces
  CGAL_Triangulation_2(const CGAL_Triangulation_2<Gt,Tds> &tr)
    : _tds(tr._tds), _gt(tr._gt)
  {}
  

  CGAL_Triangulation_2(Tds tds, const Geom_traits& geom_traits=Geom_traits())
    : _tds(),_gt(geom_traits)
  {
    _tds.swap(tds);
  }
  
 

  //Assignement
  CGAL_Triangulation_2 &operator=(const CGAL_Triangulation_2 &tr)
  {
    _tds = tr._tds;
    _gt = tr._gt;
     return *this;
  }

  // Helping functions
   // should be perhaps protected or even removed
  void  copy_triangulation(const CGAL_Triangulation_2 &tr)
  {
     _gt = tr._gt;
    _tds.copy_tds(tr._tds);
  }

  void swap(CGAL_Triangulation_2 &tr)
  {
    Geom_traits t = geom_traits();
    _gt = tr.geom_traits();
    tr._gt = t; 
    _tds.swap(tr._tds);
  }

  void clear()
  {
    _tds.clear();
  }


  //STATIC
  static int ccw(int i) {return (i+1) % 3;}
  static int cw(int i) {return (i+2) % 3;}

  //ACCESS FUNCTIONs
  int  dimension() const { return _tds.dimension();}
  int number_of_vertices() const {return _tds.number_of_vertices();}
  int number_of_faces() const {return _tds.number_of_faces();}
  const Geom_traits& geom_traits() const { return _gt;}
  
  const Vertex_handle finite_vertex() const
  {
    //return ((Vertex *) _tds.finite_vertex())->handle();
    return (Vertex *) _tds.finite_vertex();
  }

  const Vertex_handle infinite_vertex() const
  {
    //return ((Vertex *) _tds.infinite_vertex())->handle();
    return (Vertex *) _tds.infinite_vertex();
  }
   
  Face_handle infinite_face() const
  {
    CGAL_triangulation_precondition(infinite_vertex()->face()->
                                       has_vertex(infinite_vertex()));
    return infinite_vertex()->face();
  }

 
  //SETTING
protected:
  void set_finite_vertex(const Vertex_handle&  v)
  {
        _tds.set_finite_vertex (&(*v));
  }

//  void set_finite_vertex(CGAL_NULL_TYPE v)
//   {
//         _tds.set_finite_vert(v); 
//   }
// 
  void set_infinite_vertex(const Vertex_handle&  v)
  {
        _tds.set_infinite_vertex(&(*v));
  }

  void set_number_of_vertices(int n) { _tds.set_number_of_vertices(n);}

public:   
  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const
  {
    return _tds.is_valid();
  }

  // TEST IF INFINITE FEATURES
  bool is_infinite(const Face_handle& f) const {
    return f->has_vertex(infinite_vertex());
  }

  bool is_infinite(const Vertex_handle& v) const {
        return v == infinite_vertex();
  }

  bool is_infinite(const Face_handle& f, int i) const {
    return is_infinite(f->vertex(ccw(i))) ||
      is_infinite(f->vertex(cw(i)));
  }

  bool is_infinite(Edge e) const {
    return is_infinite(e.first,e.second);
  }
   
  bool is_infinite(Edge_circulator& ec) const {
    return is_infinite(*ec);
  }

  bool is_infinite(Edge_iterator& ei) const {
    return is_infinite(*ei);
  }


 // GEOMETRIC FEATURES
  Triangle triangle(const Face_handle& f) const
  {
    CGAL_triangulation_precondition( ! is_infinite(f) );
    return Triangle(f->vertex(0)->point(),
		    f->vertex(1)->point(),
		    f->vertex(2)->point());
  }

  Segment segment(const Face_handle& f, int i) const
  {
    CGAL_triangulation_precondition
      ((! is_infinite(f->vertex(ccw(i)))) &&
       (! is_infinite(f->vertex(cw(i)))) );
    return Segment(f->vertex(ccw(i))->point(),
		   f->vertex(cw(i))->point());
  }

  Segment segment(Edge e) const
  {
    CGAL_triangulation_precondition(! is_infinite(e));
    return Segment(e.first->vertex(ccw(e.second))->point(),
		   e.first->vertex( cw(e.second))->point());
  }

  Segment segment(const Edge_circulator& ec) const
  {
    return segment(*ec);
  }
    
  Segment segment(const Edge_iterator& ei) const
  {
    return segment(*ei);
  }



  //INSERTION - DELETION - Flip
public:
  void flip(Face_handle& f, int i)
  {
    _tds.flip( &(*f), i);
  }

  Vertex_handle insert(const Point& p,
                         Locate_type& lt,
                         Face_handle f = Face_handle() )
  {
    Vertex_handle v;
    if(number_of_vertices() == 0) {
      v = new Vertex(p);
      lt = OUTSIDE;
      _tds.insert_first(&(*v));
      return v;
    }
    if(number_of_vertices() == 1) {
      if (p == finite_vertex()->point()) {
	lt = VERTEX;
	return finite_vertex();
      }
      v = new Vertex(p);
      lt = OUTSIDE;
      _tds.insert_second(&(*v));
      return v;
    }

    int li;
    Face_handle loc = locate(p, lt, li, f);
    
    switch(lt){
    case FACE:
      {
	v = new Vertex(p);
       _tds.insert_in_face( &(*v), &(*loc));
	break;
      }
    case OUTSIDE:
      {
	v = new Vertex(p);
	insert_outside(v,loc);
	break;
      }
    case EDGE:
      {
	v = new Vertex(p);
	_tds.insert_on_edge( &(*v), &(*loc), li);
	break;
      }
    case COLLINEAR_OUTSIDE:
      {
	v = new Vertex(p);
	_tds.insert_colinear_outside( &(*v), &(*loc),li);
      	 break; 
      }
    case VERTEX:
      return loc->vertex(li);
    default:
      CGAL_triangulation_assertion(false);  // locate step failed
    }
    return v;
  }

  Vertex_handle insert(const Point &p,
                         Face_handle f = Face_handle() )
  {
      Locate_type lt;
    return insert(p, lt, f);
  }

    Vertex_handle push_back(const Point &p)
    {
        Locate_type lt;
        return insert(p, lt, NULL);
    }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
    template < class InputIterator >
    int insert(InputIterator first, InputIterator last)
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
    int insert(list<Point>::const_iterator first,
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
    int insert(vector<Point>::const_iterator first,
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
    int insert(istream_iterator<Point, ptrdiff_t> first,
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
    int insert(Point* first,
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

    
protected:
  void insert_outside(const Vertex_handle& vp, Face_handle& f)
  {
    CGAL_triangulation_precondition(is_infinite(f));

    list<Face_handle> ccwlist;
    list<Face_handle> cwlist;

    int li;
    Point p = vp->point();
    Point q,r;

    li = f->index(infinite_vertex());
    q = f->vertex(ccw(li))->point();
    r = f->vertex(cw(li))->point();
    CGAL_triangulation_precondition( geom_traits().orientation(p,q,r) == CGAL_LEFTTURN);
    
    Face_circulator fc = infinite_vertex()->incident_faces(f);
    bool done = false;
    while(! done) {
      fc--;
      li = fc->index(infinite_vertex());
      q = fc->vertex(ccw(li))->point();
      r = fc->vertex(cw(li))->point();
      if(geom_traits().orientation(p,q,r) == CGAL_LEFTTURN ) {
      ccwlist.push_back(&(*fc));
      }
      else {done=true;}
    }

    fc= infinite_vertex()->incident_faces(f);
    done = false;
    while(! done){
      fc++;
      li = fc->index(infinite_vertex());
      q = fc->vertex(ccw(li))->point();
      r = fc->vertex(cw(li))->point();
      if(geom_traits().orientation(p,q,r) == CGAL_LEFTTURN ) {
	cwlist.push_back(&(*fc));
      }
      else {done=true;}
    }

    _tds.insert_in_face( &(*vp), &(*f));

     Face_handle fh;
     while ( ! ccwlist.empty()) {
      fh = ccwlist.front();
      li = ccw(fh->index(infinite_vertex()));
      flip( fh , li);
      ccwlist.pop_front();
    }

    while ( ! cwlist.empty()) {
      fh = cwlist.front();
      li = cw(fh->index(infinite_vertex()));
      flip( fh , li);
      cwlist.pop_front();
    }

    //reset infinite_vertex()->face();
    fc = vp->incident_faces();
    while( ! is_infinite(&(*fc))) {
      fc++;}
    infinite_vertex()->set_face(&(*fc));
  } 
    
    
public:  
 void remove(Vertex_handle  v)
    {
      CGAL_triangulation_precondition(v != NULL);
      CGAL_triangulation_precondition( !is_infinite(v));
    
      if  (number_of_vertices() == 1) {
        _tds.remove_first(& (*v));
	v.clear();
        return;
      }
     
    
      // take care of finite_vertex data member
      if (finite_vertex() == v){
        Face_handle f = v->face();
        int i=f->index(v);
        Vertex_handle vv=is_infinite(f->vertex(cw(i))) ?
                           f->vertex(ccw(i)) : f->vertex(cw(i));
        set_finite_vertex( vv);
      }
    
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
      return;
    }

protected:
    void remove_1D(Vertex_handle& v)
    {
     //deal with one dimensional case
      Face_handle f1, f2, f1n, f2n;
      int i1,i2,i1n,i2n;
    
        switch( v->degree())  {
        case 2 :
          f1= v->face(); i1= f1->index(v);
          f2= f1->neighbor(cw(i1));  i2= f2->index(v);
          f1n = f1->neighbor(i1); i1n = f1n->index(f1);
          f2n = f2->neighbor(i2); i2n = f2n->index(f2);
          f1n->set_neighbor(i1n,f2n); f2n->set_neighbor(i2n,f1n);
	  f1->vertex(cw(i1))->set_face(f1n);
	  f1->vertex(ccw(i1))->set_face(f1n);
          f1.Delete(); f2.Delete();
          return;
    
        case 4 :
          int j1,j2,j1n,j2n;
           f1= v->face(); i1= f1->index(v);
           j1 = (is_infinite(f1->vertex(cw(i1)))) ? ccw(i1) : cw(i1);
           // j1 is the other finite vertex of f1
           // 3-i1-j1 is the index of the infinite_vertex in face f1
           f2 = f1->neighbor(3-i1-j1); i2 = f2->index(v);
           j2 = f2->index(f1->vertex(j1));
           f1n = f1->neighbor(j1); i1n = f1n->index(v); j1n = f1n->index(f1);
           f2n = f2->neighbor(j2); i2n = f2n->index(v); j2n = f2n->index(f2);
    
           // update f1n and f2n
           f1n->set_vertex(i1n,f1->vertex(j1));
           f2n->set_vertex(i2n,f2->vertex(j2));
           if(f1->neighbor(i1) != f2) {
             Face_handle ff1 = f1->neighbor(i1);
             Face_handle ff2 = f2->neighbor(i2);
             f1n->set_neighbor(j1n,ff1); ff1->set_neighbor(ff1->index(f1),f1n);
             f2n->set_neighbor(j2n,ff2); ff2->set_neighbor(ff2->index(f2),f2n);
           }
           else {
             f1n->set_neighbor(j1n,f2n);
             f2n->set_neighbor(j2n,f1n);
           }

	   // update the face() pointer of the remaining vertices of f1 and f2
	   f1->vertex(j1)->set_face(f1n);
	   f2->vertex(j2)->set_face(f2n);
	   infinite_vertex()->set_face(f1n);

           f1.Delete(); f2.Delete();
           return;
        }
     }
    
  void remove_2D(Vertex_handle& v)
    {
      // General case
    
      // remove incident faces
      // set up list of faces neighboring the hole
      // in ccw order around the hole
    
      // problem with gcc link
      typedef pair<void *, int> Hole_neighbor;
      //typedef pair<Face_handle  , int>  Hole_neighbor;
      //typedef pair<Face*, int> Hole_neighbor;
      typedef list<Hole_neighbor> Hole;
      
    
      Hole hole;
      typename Hole::iterator hit;
      list<Face*> to_delete;

      Face_handle  f, ff, fn;
      int i =0,ii =0, in =0;
      Vertex_handle  vv;
      Point p = v->point();
    
      Face_circulator fc = v->incident_faces();
      Face_circulator done(fc);
       do {
	 f = (*fc).handle(); fc++;
	 i = f->index(v);
	 fn = f->neighbor(i);
	 vv = f->vertex(f->cw(i));
	 if( vv->face()==  f) vv->set_face(fn);
	 vv = f->vertex(fc->ccw(i));
	 if( vv->face()== f) vv->set_face(fn);
	 in = fn->index(f);
	 fn->set_neighbor(in, NULL);
	 hole.push_back(Hole_neighbor(&(*fn),in));
	 //hole.push_back(Hole_neighbor(fn,in));
	 // to_delete.push_back(f);
	 to_delete.push_back(&(*f));
       }
      while(fc != done);

       while (! to_delete.empty()){
	 delete to_delete.front();
	 to_delete.pop_front();
       }
    
       Vertex_handle v0, v1, v2;
       Point   p0, p1, p2;
       Triangle t;
       CGAL_Bounded_side side;
       	CGAL_Orientation or;
       int nhole;

       if( hole.size() != 3) {
	 //  find the first edge  v0v1 on the hole boundary
	 //  such that
	 //  v0, v1 and the next vertex v2 are all finite
         //   v0v1v2 is a left turn and
	 //  triangle v0v1v2 does not contain the removed point
    
	 // if found create face v0v1v2
	 // stop when no more faces can be created that way
                  
         nhole= hole.size(); 
	 // nhole decount the number of hole edges passed
	 // from the last created edges
         while (hole.size()>3 && nhole>0) {           
            ff  = (Face *) ( (hole.front()).first);
	    //ff =  hole.front().first
	    ii  = (hole.front()).second;
	    hole.pop_front();

            v0 = ff->vertex(cw(ii));
            v1 = ff->vertex(ccw(ii));
            if( !is_infinite(v0) && !is_infinite(v1)) {
	      fn  = (Face *) ( (hole.front()).first);
	      //fn =  hole.front().first
	      in  = (hole.front()).second;
            
	      v2 = fn->vertex(ccw(in));
	      if( !is_infinite(v2)) {
		p0 = v0->point();
		p1 = v1->point();
		p2 = v2->point();
		or = geom_traits().orientation(p0,p1,p2);
		if ( or  == CGAL_LEFTTURN) {
		  side = bounded_side(p0,p1,p2, p);
		  if( side == CGAL_ON_UNBOUNDED_SIDE || 
		      (side == CGAL_ON_BOUNDARY && collinear_between(p0, p, p2)) ) {
		    //create face
		    Face_handle  newf = new Face(v0,v1,v2);
		    newf->set_neighbor(2,ff);
		    newf->set_neighbor(0,fn);
		    ff->set_neighbor(ii, newf);
		    fn->set_neighbor(in,newf);
		    hole.pop_front();
		    hole.push_front(Hole_neighbor(&(*newf),1));
		    nhole = hole.size();
		    continue;
		  }
		}
	      }
	    }

	    // not possible to create face v0,v1,v2;
	    hole.push_back(Hole_neighbor(&(*ff),ii));
	    nhole--;
	 }
       }
    
      // either the hole has only three edges
      // or all its finite vertices are reflex or flat
      // except may be one vertex whose corresponding ear 
      // includes the vertex being removed

      // deal with the last leftturn if any
      if(hole.size() != 3) {
	nhole = hole.size();
	while ( nhole>0) {
	  ff = (Face *) ((hole.front()).first);
	  // ff = ((hole.front()).first)
	  ii = (hole.front()).second;
	  hole.push_back(hole.front());
	  hole.pop_front();
	  nhole--;
	  
	  v0 = ff->vertex(cw(ii));
	  v1 = ff->vertex(ccw(ii));
	  if(is_infinite(v0) || is_infinite(v1))  continue;

	  fn = (Face *) ((hole.front()).first);
	  // fn = ((hole.front()).first)
	  in = (hole.front()).second;
	  v2 = fn->vertex(ccw(in));
          if( is_infinite(v2) ) continue;
	  p0 = v0->point();
	  p1 = v1->point();
	  p2 = v2->point();
	  CGAL_Orientation or = geom_traits().orientation(p0,p1,p2);
          if ( or  == CGAL_LEFTTURN) {
	    Face_handle  newf = new Face(v0,v1,v2);
	    newf->set_neighbor(2,ff);
	    newf->set_neighbor(0,fn);
	    ff->set_neighbor(ii, newf);
	    fn->set_neighbor(in,newf);
	    hole.pop_back();
	    hole.pop_front();
	    hole.push_front(Hole_neighbor(&(*newf),1));
	    break;
	  }
	}
      }


      if(hole.size() != 3) {
	// look for infinite vertex
	ff = (Face *) ((hole.front()).first);
	//ff = (hole.front()).first;
	//ff = ((hole.front()).first)->handle();
	ii = (hole.front()).second;
	while ( ! is_infinite(ff->vertex(cw(ii)))){
	  hole.push_back(hole.front());
	  hole.pop_front();
	  ff = (Face *)((hole.front()).first);
	  //ff = (hole.front()).first;
	  //ff = ((hole.front()).first)->handle();
	  ii = (hole.front()).second;
          }
	//create faces
          while(hole.size() != 3){
            ff = (Face *)((hole.front()).first);
	    //ff = (hole.front()).first;
	    //ff = ((hole.front()).first)->handle();
            ii = (hole.front()).second;
            hole.pop_front();
            fn = (Face *)((hole.front()).first);
	    //fn = (hole.front()).first;
	    //fn = ((hole.front()).first)->handle();
            in = (hole.front()).second;
            hole.pop_front();
            Face_handle  newf = new Face(infinite_vertex(),fn->vertex(cw(in)),
                                                     fn->vertex(ccw(in)));
            ff->set_neighbor(ii,newf);
            fn->set_neighbor(in,newf);
            newf->set_neighbor(0,fn);
            newf->set_neighbor(2,ff);
            hole.push_front(Hole_neighbor(&(*newf),1));
	    //hole.push_front(Hole_neighbor(newf,1));
            }
      }
    
 // now hole has three edges
          Face_handle  newf = new Face();
          hit = hole.begin();
          for(int j = 0;j<3;j++) {
            ff = (Face *)((*hit).first);
	    //ff = (*hit).first;
	    //ff = ((*hit).first)->handle();
            ii = (*hit).second;
            hit++;
            ff->set_neighbor(ii,newf);
            newf->set_neighbor(j,ff);
            newf->set_vertex(newf->ccw(j),ff->vertex(ff->cw(ii)));
          }
    }

    

  // POINT LOCATION
public:
    class Line_face_circulator
      :   public CGAL_Bidirectional_circulator_base<Face,ptrdiff_t,size_t>,
          public Face::Face_handle
    {
        public:
            //typedef Vertex::Face_circulator
            //                       Face_circulator ;
            
            typedef Face           value_type;
            typedef Face   &       reference;
            typedef const Face   & const_reference;
            typedef unsigned int   size_type;
            typedef int            distance_type;
            
            enum State {undefined = -1,
                        vertex_vertex,
                        vertex_edge,
                        edge_vertex,
                        edge_edge};
            
            Line_face_circulator()
                : Face::Face_handle(), _tr(NULL), s(undefined), i(-1)
            {}
            
            Line_face_circulator(const Line_face_circulator& lfc)
                : Face::Face_handle(& (*lfc)), _tr(lfc._tr), s(lfc.s), i(lfc.i),  p(lfc.p), q(lfc.q)
            {}
            
            ~Line_face_circulator()
            {}
            
            
            Line_face_circulator(const Face_handle& face,
                                 int index,
                                 State state,
                                 CGAL_Triangulation_2<Gt,Tds> * t,
                                 const Point& pp,
                                 const Point& qq)
                : Face::Face_handle(face), _tr(t), s(state), i(index),  p(pp), q(qq)
            {
                CGAL_triangulation_precondition(p != q);
            }
            
            
            Line_face_circulator&
            operator=(const Line_face_circulator& lfc)
            {
                ptr() = lfc.ptr();
                i = lfc.i;
                s = lfc.s;
                _tr = lfc._tr;
                p = lfc.p;
                q = lfc.q;
                return *this;
            }

            Line_face_circulator(Vertex_handle v,
                                 CGAL_Triangulation_2<Gt,Tds>* tr,
                                 const Point& dir)
                : _tr(tr)
            {
                CGAL_triangulation_precondition((! _tr->is_infinite(v)) &&
                                                (_tr->dimension() == 2) &&
                                                (v->point() != dir));

		p=v->point();
		q=dir;

		cerr << " p " << p << " q " << q << endl;
		
		Face_circulator fc = v->incident_faces();
		Face_circulator done = fc;

		cerr  << "(" << fc->vertex(0)->point() << ", "
		      <<fc->vertex(1)->point() << ", "
		      << fc->vertex(2)->point() << ")" << endl ;

		int ic = fc->index(v);
		Vertex_handle  vt= fc->vertex(ccw(ic));
		CGAL_Orientation ptq;
		if (! _tr->is_infinite(vt)) 
		  ptq = _tr->geom_traits().orientation(p, vt->point(), q);

		
		while( _tr->is_infinite(vt) || ptq == CGAL_RIGHTTURN) {
		  ++fc;
		  if (fc == done) {
		    // no edge on the left of pq , pq is a supporting line
		    // set ptr() to the right infinite face
		    while ( ! _tr->is_infinite(fc)) 
		      {  ++fc;}
		    ic = fc->index(_tr->infinite_vertex());
		    if( _tr->geom_traits().orientation(
			       fc->vertex( cw(i))->point(),
			       fc->vertex( ccw(i))->point(),
							q) != CGAL_LEFTTURN) {  ++fc;}
		    ptr() = &(*fc);
		    i = fc->index(_tr->infinite_vertex());
		    s = vertex_vertex;
		    return;
		  }

		  ic = fc->index(v);
		  vt= fc->vertex(ccw(ic));
		  if (! _tr->is_infinite(vt)) 
		    ptq = _tr->geom_traits().orientation(p,  vt->point(), q);
		}
		
		cerr << "t " << vt->point() << endl;


		// now vt is a finite vertex and ptq is COLLINEAR or LEFTTURN
		Vertex_handle vr = fc-> vertex(cw(ic));
		CGAL_Orientation prq;
		if (! _tr->is_infinite(vr))
		  prq = _tr->geom_traits().orientation(p, vr->point(), q);

		while ( (!_tr->is_infinite(vr)) && (!(prq == CGAL_RIGHTTURN ))){
		  ++fc;
		  ic = fc->index(v);
		  vr = fc-> vertex(cw(ic));
		  if (! _tr->is_infinite(vr))
		    prq = _tr->geom_traits().orientation(p, vr->point(), q);
		}

	
		ptr() = &(*fc);
		// reset vt, vt is finite and ptq is still COLLINEAR or LEFTTURN
		ic = fc->index(v);
		vt= fc->vertex(ccw(ic));
		ptq = _tr->geom_traits().orientation(p,  vt->point(), q);

		cerr << "t " << vt->point() << endl;
		cerr << "r " << vr->point() << endl;
		cerr << "face" << "(" << fc->vertex(0)->point() << ", "
		      <<fc->vertex(1)->point() << ", "
		      << fc->vertex(2)->point() << ")" << endl ;

		if (_tr->is_infinite(vr)) {		  
		  s = vertex_vertex;
		  if (ptq == CGAL_LEFTTURN) {
		    i = fc->index(vr); 
		  }
		  else {// ptq == CGAL_COLLINEAR
		    i= fc->index(vt);
		  }
		}
		else{ // vr is a finite vertex}
		  if (ptq == CGAL_LEFTTURN) {
		    s = vertex_edge;
		    i = ic;
		  }
		  else { // ptq == CGAL_COLLINEAR
		    s = vertex_vertex;
		    i = fc->index(vt);
		  }
		}

	
	    }

		    

            
//                 ptr() = &(*v->face());
//             
//                 while(_tr->is_infinite(ptr()->handle())){
//                     i = ptr()->index(v);
//                     ptr() = &(*ptr()->neighbor(cw(i)));
//                 }
//                 i = ptr()->index(v);
//                 p = v->point();
//                 q = dir;
//             
//                 Point t = ptr()->vertex(ccw(i))->point();
//                 Point r = ptr()->vertex(cw(i))->point();
//             
//                 CGAL_Orientation ptq = _tr->geom_traits().orientation(p, t, q);
//                 CGAL_Orientation prq = _tr->geom_traits().orientation(p, r, q);
//             
//                 if( ptq == CGAL_LEFTTURN  && prq == CGAL_RIGHTTURN ) {
//                     s = vertex_edge;
//                 } else if( (ptq == CGAL_COLLINEAR) && (prq == CGAL_RIGHTTURN) ) {
//                     i= ccw(i);
//                     s = vertex_vertex;
//                 } else if( (prq == CGAL_COLLINEAR) && (ptq == CGAL_LEFTTURN) ) {
//                     Face_handle n = ptr()->neighbor(ccw(i));
//                     i = n->index(ptr()->vertex(cw(i)));
//                     ptr() = &(*n);
//                     s = vertex_vertex;
//                 } else if( ptq == CGAL_LEFTTURN ) {
//                     do{
//                         Face_handle n = ptr()->neighbor(ccw(i));
//                         i = n->index(ptr()->handle());
//                         ptr() = &(*n);
//                         if(_tr->is_infinite(n)){
//                             s = vertex_vertex;
//                             return;
//                         }
//                         r = n->vertex(i)->point();
//                         i = ccw(i);
//                     }while((prq = _tr->geom_traits().orientation(p, r, q)) == CGAL_LEFTTURN);
//             
//                     if(prq == CGAL_COLLINEAR){
//                         Face_handle n = ptr()->neighbor(ccw(i));
//                         i = cw(n->index(ptr()->handle()));
//                         ptr() = &(*n);
//                         s = vertex_vertex;
//                     } else {
//                         s = vertex_edge;
//                     }
//                 } else 
// {
//                     do{
//                         Face_handle n = ptr()->neighbor(cw(i));
//                         i = n->index(ptr()->handle());
//                         ptr() = &(*n);
//                         if(_tr->is_infinite(n)){
//                             s = vertex_vertex;
//                             return;
//                         }
//                         t = n->vertex(i)->point();
//                         i = cw(i);
//                     }while((ptq = _tr->geom_traits().orientation(p, t, q)) == CGAL_RIGHTTURN);
//             
//                     if(ptq == CGAL_COLLINEAR){
//                         s = vertex_vertex;
//                         i = ccw(i);
//                     } else {
//                         s = vertex_edge;
//                     }
//                 }
//            }

            Line_face_circulator(const Point& pp,
                                 const Point& qq,
                                 CGAL_Triangulation_2<Gt,Tds> * t)
                : _tr(t), s(undefined), p(pp), q(qq)
            {
                Vertex_handle inf = _tr->infinite_vertex();
                 Face_circulator fc = inf->incident_faces(),
                    done(fc);
                i = fc->index(inf);
            
                Point l = fc->vertex(cw(i))->point(),
                    r = fc->vertex(ccw(i))->point();
            
                CGAL_Orientation pql = _tr->geom_traits().orientation(p, q, l),
                    pqr = _tr->geom_traits().orientation(p, q, r);
            
            
                do{
                    if( (pql == CGAL_LEFTTURN) && (pqr == CGAL_RIGHTTURN) ){
                        *this = ++Line_face_circulator( fc->handle() ,
                                                       i,
                                                       vertex_edge,
                                                       t,
                                                       p,
                                                       q);
                        return;
                    } else if ( (pql == CGAL_LEFTTURN) && (pqr == CGAL_COLLINEAR) ){
                        *this = ++Line_face_circulator( fc->handle() ,
                                                       ccw(i),
                                                       vertex_vertex,
                                                       t,
                                                       p,
                                                       q);
                        return;
                    } else if( (pql == CGAL_COLLINEAR) && (pqr == CGAL_COLLINEAR) ){
                        Face_handle n = fc->neighbor(i);
                        int ni  = n->index( fc->handle() );
                        Vertex_handle vn = n->vertex(ni);
                        if(_tr->geom_traits().orientation(p, q, vn->point()) == CGAL_LEFTTURN){
                            // the entire triangulation is to the left of line (p,q).
                            // There might be further collinear edges, so we might have
                            // to walk back on the hull.
                            while(1){
                                ++fc;
                                i = fc->index(inf);
                                l = fc->vertex(cw(i))->point();
                                if(_tr->geom_traits().orientation(p, q, l) == CGAL_COLLINEAR){
                                    continue;
                                } else {
                                    // we went one step to far back
                                    --fc;
                                    i = fc->index(inf);
                                    ptr() = &(*fc->neighbor(i));
                                    i = cw(ptr()->index( fc->handle() ));
                                    s = vertex_vertex;
                                    return;
                                }
                            }
                        } else {
                            // the entire triangulation is to the right of line (p,q).
                            // here are no faces to traverse, so we give the circulator
                            // a singular value
                            return;
                        }
                    } else {
                        --fc;
                        l = r;
                        pql = pqr;
                        i = fc->index(inf);
                        r = fc->vertex(ccw(i))->point();
                        pqr = _tr->geom_traits().orientation(p, q, r);
                    }
                }while(fc != done);
                // if line (p,q) does not intersect the convex hull in an edge
                // the circulator has a singular value
            }

            Line_face_circulator(const Point& pp,
                                 const Point& qq,
                                 const Face_handle& ff,
                                 CGAL_Triangulation_2<Gt,Tds>* t)
                : Face::Face_handle(ff), _tr(t), s(undefined), p(pp), q(qq)
            {
                //CGAL_triangulation_precondition(_tr->is_infinite(f) ||
                //                         _tr->oriented_side(f,p) != CGAL_ON_NEGATIVE_SIDE);
            
                int j;
                if(_tr->is_infinite(ptr()->handle())){
                    *this  = Line_face_circulator(p, q, t);
                    return;
                }
            
                // Test whether p lies on a vertex
                for(j = 0; j < 3; j++){
                    if(ptr()->vertex(j)->point() == p){
                        *this = Line_face_circulator( ptr()->vertex(j), t, q);
                        return;
                    }
                }
            
                // Test whether p lies on an edge
                for(j = 0; j < 3; j++){
                    if(_tr->geom_traits().orientation(ptr()->vertex(j)->point(),
                                             ptr()->vertex(ccw(j))->point(),
                                             p) == CGAL_COLLINEAR){
                        CGAL_Orientation jpq =
                            _tr->geom_traits().orientation(ptr()->vertex(j)->point(),
                                                  p,
                                                  q);
                        CGAL_Orientation p_cwj_q =
                            _tr->geom_traits().orientation(p,
                                                  ptr()->vertex(cw(j))->point(),
                                                  q);
                        switch(jpq){
                        case CGAL_COLLINEAR:
                            if(p_cwj_q == CGAL_RIGHTTURN){
                                s = vertex_vertex;
                                i = ccw(j);
                                return;
                            } else if(! _tr->is_infinite(ptr()->neighbor(cw(j)))){
                                Face_handle n = ptr()->neighbor(cw(j));
                                i = cw(n->index(ptr()->handle()));
                                ptr() = &(*n);
                                s = vertex_vertex;
                                return;
                            } else {
                                // singular value
                                return;
                            }
                        case CGAL_RIGHTTURN:
                            i = cw(j);
                            s = (p_cwj_q == CGAL_COLLINEAR) ? vertex_edge :  edge_edge;
                            break;
                        default: //  CGAL_LEFTTURN
                            switch(p_cwj_q){
                            case CGAL_COLLINEAR:
                                s = edge_vertex;
                                i = cw(j);
                                return;
                            case CGAL_RIGHTTURN:
                                s = edge_edge;
                                i = j;
                                return;
                            default:
                                s = edge_edge;
                                i = ccw(j);
                                return;
                            }
                        }
                    }
                }
            
                // p lies in the interior of the face
                CGAL_Orientation or[3];
                for(j=0; j<3; j++){
                    or[j] = _tr->geom_traits().orientation(p,q,ptr()->vertex(j)->point());
                }
                for(j=0; j<3; j++){
                    if(or[j] == CGAL_COLLINEAR){
                        i = j;
                        s = (or[ccw(j)] == CGAL_LEFTTURN) ? edge_vertex : vertex_edge;
                        return;
                    }
                }
                s = edge_edge;
                for(j=0; j<3; j++){
                    if(or[j] == CGAL_RIGHTTURN){
                        i = (or[ccw(j)] == CGAL_RIGHTTURN) ? j : cw(j);
                        return;
                    }
                }
            }

            void increment()
            {
                CGAL_triangulation_precondition(s != undefined);
                if(s == vertex_vertex || s == edge_vertex){
                 CGAL_Orientation o;
                    Point r;
                    do{
                        Face_handle n = ptr()->neighbor(cw(i));
                        i = n->index(ptr()->handle());
                        ptr() = &(*n);
                        if (n->vertex(i) == _tr->infinite_vertex()){
                            o = CGAL_COLLINEAR;
                            i = cw(i);
                            break;
                        }
                        r = n->vertex(i)->point();
                        i = cw(i);
                    }while((o = _tr->geom_traits().orientation(p, q, r)) == CGAL_LEFTTURN);
            
                    if(o == CGAL_COLLINEAR){
                        s = vertex_vertex;
                        i = ccw(i);
                    } else {
                        s = vertex_edge;
                    }
                } else {
                    Face_handle n = ptr()->neighbor(i);
                    int ni = n->index(ptr()->handle());
                    ptr() = &(*n);
                    CGAL_Orientation o = _tr->is_infinite(ptr()->vertex(ni)) ?
                        CGAL_COLLINEAR :
                        _tr->geom_traits().orientation(p, q, ptr()->vertex(ni)->point());
            
                    switch(o){
                    case CGAL_LEFTTURN:
                        s = edge_edge;
                        i = ccw(ni);
                        break;
                    case CGAL_RIGHTTURN:
                        s = edge_edge;
                        i = cw(ni);
                        break;
                    default:
                        s = edge_vertex;
                        i = ni;
                    }
                }
            }
            
            
            void decrement()
            {
                CGAL_triangulation_precondition(s != undefined);
                if(s == vertex_vertex || s == vertex_edge){
                    if(s == vertex_vertex){
                        i = cw(i);
                    }
                    CGAL_Orientation o;
                    Point r;
                    do{
                        Face_handle n = ptr()->neighbor(ccw(i));
                        i = n->index(ptr()->handle());
                        ptr() = &(*n);
                        if (n->vertex(i) == _tr->infinite_vertex()){
                            o = CGAL_COLLINEAR;
                            i = ccw(i);
                            break;
                        }
                        r = n->vertex(i)->point();
                        i = ccw(i);
                    }while((o = _tr->geom_traits().orientation(p, q, r)) == CGAL_LEFTTURN);
            
                    s = (o == CGAL_COLLINEAR) ? vertex_vertex : edge_vertex;
            
                } else { // s == edge_edge  ||  s == edge_vertex
                    // the following is not nice. A better solution is to say
                    // that index i is at the vertex that is alone on one side of l(p,q)
                    if(s == edge_edge){
                        i = (_tr->geom_traits().orientation
                                            (p, q,
                                            ptr()->vertex(i)->point()) == CGAL_LEFTTURN)
                            ? cw(i) : ccw(i);
                    }
                    Face_handle n = ptr()->neighbor(i);
                    i = n->index(ptr()->handle());
                    ptr() = &(*n);
                    CGAL_Orientation o = _tr->is_infinite(ptr()->vertex(i)) ?
                        CGAL_COLLINEAR :
                        _tr->geom_traits().orientation(p, q, ptr()->vertex(i)->point());
            
                    s = (o == CGAL_COLLINEAR) ? vertex_edge : edge_edge;
                }
            }

            bool
            locate(const Point& t,
                   Locate_type &lt,
                   int &li)
            {
                switch(s){
            
                case edge_edge:
                case vertex_edge:
                    {
                        CGAL_Orientation o =
                            _tr->geom_traits().orientation(ptr()->vertex(ccw(i))->point(),
                                                  ptr()->vertex(cw(i))->point(),
                                                  t);
                        if(o == CGAL_RIGHTTURN){
                            return false;
                        }
                        if(o == CGAL_COLLINEAR){
                            lt = EDGE;
                            li = i;
                            return true;
                        }
                        lt = FACE;
                        return true;
                    }
                case vertex_vertex:
                    {
                        if(_tr->is_infinite(ptr()->vertex(i))){
                            CGAL_triangulation_assertion(_tr->geom_traits().orientation(
                                                         ptr()->vertex(cw(i))->point(),
                                                         ptr()->vertex(ccw(i))->point(),
                                                         t) != CGAL_LEFTTURN);
                            lt = OUTSIDE;
                            li = i;
                            return true;
                        }
            
                        Point u = ptr()->vertex(cw(i))->point();
                        Point v = ptr()->vertex(i)->point();
                        // u == t  was detected earlier
                        if(_tr->geom_traits().compare_x(v,t)==CGAL_EQUAL && 
			   _tr->geom_traits().compare_y(v,t)==CGAL_EQUAL){
                            lt = VERTEX;
                            li = i;
                            return true;
                        }
                        if(_tr->collinear_between(u, t, v)){
                            lt = EDGE;
                            li = ccw(i);
                            return true;
                        }
                        return false;
                    }
                default: // edge_vertex
                    {
                       if(_tr->is_infinite(ptr()->vertex(i))){
                            lt = OUTSIDE;
                            li = i;
                            return true;
		       }
		       if(_tr->geom_traits().compare_x(t,ptr()->vertex(i)->point())==CGAL_EQUAL &&
			  _tr->geom_traits().compare_y(t,ptr()->vertex(i)->point())==CGAL_EQUAL ){
                            li = i;
                            lt = VERTEX;
                            return true;
		       }
		       if(_tr->collinear_between(p, t, ptr()->vertex(i)->point())){
                            lt = FACE;
                            return true;
		       }
		       return false;
                    }
                }
            }

            Line_face_circulator&
            operator++()
            {
                if (ptr()==NULL) {
                   return *this; // circulator has singular value
                }
                //xfc while (_tr->is_infinite(ptr()->handle()))
                //    increment();
                //we are looking for a finite face but stop if back to the origin
                Face_handle origin = ptr()->handle();
                do
                 increment();
                while ((_tr->is_infinite(ptr()->handle())) && (ptr()->handle() != origin));
            
                return *this;
            }
            
            
            Line_face_circulator&
            operator--()
            {
                if (ptr()==NULL) {
                   return *this; // circulator has singular value
                }
                while (_tr->is_infinite(ptr()->handle()))
                       decrement();
                return *this;
            }
            
            
            Line_face_circulator
            operator++(int)
            {
                Line_face_circulator tmp(*this);
                ++(*this);
                return tmp;
            }
            
            
            Line_face_circulator
            operator--(int)
            {
                Line_face_circulator tmp(*this);
                --(*this);
                return tmp;
            }

            bool
            operator==(const Line_face_circulator& lfc) const
            {
                CGAL_triangulation_precondition
                          ( ptr() != NULL  &&  lfc.ptr() != NULL );
                return ptr() == lfc.ptr();
            }
            
            bool
            operator!=(const Line_face_circulator& lfc) const
            {
                CGAL_triangulation_precondition
                            ( ptr() != NULL  &&  lfc.ptr() != NULL );
                return ptr() != lfc.ptr();
            }
            
            
            inline bool
            operator==(CGAL_NULL_TYPE n) const
            {
                CGAL_triangulation_assertion( n == NULL);
                return s == undefined;
            }
            
            inline bool
            operator!=(CGAL_NULL_TYPE n) const
            {
                return !(*this == n);
            }
            
            bool
            collinear_outside() const
            {
            
                return (_tr->is_infinite(ptr()->handle()))
                        && (s == vertex_vertex)
                        && (! _tr->is_infinite(ptr()->vertex(i)));
            }
            
      //    private:
      // o debug
    public:
            CGAL_Triangulation_2<Gt, Tds>* _tr;
            State s;
            int i;
            Point p, q;
            
        };


 Face_handle
    march_locate_1D(const Point& t,
                    Locate_type& lt,
                    int& li) const
    
    {
        Face_handle i_f = infinite_face();
        int i = i_f->index(infinite_vertex());
        Point p = i_f->vertex(ccw(i))->point();
        Point q = i_f->vertex(cw(i))->point();
        CGAL_Orientation pqt = geom_traits().orientation(p, q, t);
        if(pqt == CGAL_RIGHTTURN){
            lt = OUTSIDE;
            return i_f->neighbor(i);
        }
        if(pqt == CGAL_LEFTTURN){
            lt = OUTSIDE;
            return i_f;
        }
        Vertex_handle  vi = infinite_vertex();
        Face_circulator fc = vi->incident_faces();
        i = fc->index(vi);
        Vertex_handle  u = fc->vertex(ccw(i));
        Vertex_handle  v = fc->vertex(cw(i));
        Vertex_handle  w;
    
        while(1){
            if(geom_traits().compare_x(t,v->point()) == CGAL_EQUAL && 
	       geom_traits().compare_y(t,v->point()) == CGAL_EQUAL){
                lt = VERTEX;
                li = cw(i);
                return fc->handle() ;
            }
            if(collinear_between(u->point(), t, v->point())){
                lt = EDGE;
                li = i;
                return fc->handle() ;
            }
            ++fc;
            i = fc->index(vi);
            w = fc->vertex(cw(i));
            if( (u == w) && collinear_between(u->point(), v->point(), t)){
                lt = COLLINEAR_OUTSIDE;
                li = ccw(i);
                return fc->handle() ;
            }
            u = v;
            v = w;
        }
    }
    
  Face_handle
    march_locate_2D(const Face_handle& start,
                    const Point& t,
                    Locate_type& lt,
                    int& li) const
    {
        CGAL_triangulation_precondition( ! is_infinite(start) );
        CGAL_Triangulation_2 *ncthis = (CGAL_Triangulation_2 *)this;
    
        Point p(start->vertex(0)->point());
        if(geom_traits().compare_x(t,p) == CGAL_EQUAL &&  
	   geom_traits().compare_y(t,p) == CGAL_EQUAL){
            lt = VERTEX;
            li = 0;
            return start;
        }

        Line_face_circulator lfc(start->vertex(0),
                                 ncthis,
                                 t);
	
        if(lfc.collinear_outside()){
            // point t lies outside or on the convex hull
            // we walk clockwise on the hull to decide
            int i = lfc->index(infinite_vertex());
            p = lfc->vertex(ccw(i))->point();
            if(geom_traits().compare_x(t,p) == CGAL_EQUAL &&  
	       geom_traits().compare_y(t,p) == CGAL_EQUAL){
                lt = VERTEX;
                li = ccw(i);
                return lfc;
            }
         Point q(lfc->vertex(cw(i))->point());
	 CGAL_Orientation pqt;
         Face_handle f(lfc);
         while(1){
           if(geom_traits().compare_x(t,q) == CGAL_EQUAL &&  
	      geom_traits().compare_y(t,q) == CGAL_EQUAL){
             lt = VERTEX;
             li = cw(i);
             return f;
           }
	   pqt = geom_traits().orientation(p,q,t);
	   if (pqt == CGAL_COLLINEAR && collinear_between(p, t, q)){
	     lt = EDGE;
	     li = i;
	     return f;
	   }
	   if (pqt == CGAL_LEFTTURN){
	     lt = OUTSIDE;
	     return f ;
	   }
	   	       
           // go to the next face
           f = f->neighbor(ccw(i));
           i = f->index(infinite_vertex());
           p = q;
           q = f->vertex(cw(i))->point();
     	 }
        }

        while(! lfc.locate(t, lt, li) ){
	  lfc.increment();
        }
        return lfc;
    }
    

    inline Face_handle
    locate(const Point& p,
           Locate_type& lt,
           int& li,
           Face_handle start = Face_handle()) const
    {
        if(number_of_vertices() < 2) {
            if(number_of_vertices() == 0) {
                lt = OUTSIDE;
            } else { // number_of_vertices() == 1
                lt = (p == finite_vertex()->point()) ? VERTEX : OUTSIDE;
            }
            return NULL;
        }
        if(dimension() == 1){
            return march_locate_1D(p, lt, li);
        }
    
        if(start == NULL){
            start = infinite_face()->
                        neighbor(infinite_face()->index(infinite_vertex()));
        }else if(is_infinite(start)){
            start = start->neighbor(start->index(infinite_vertex()));
        }
    
        return march_locate_2D(start, p, lt, li);
    }
    
    inline Face_handle
    locate(const Point &p,
           Face_handle start = Face_handle()) const
    {
        Locate_type lt;
        int li;
        return locate(p, lt, li, start);
    }
    



  //TRAVERSING : ITERATORS AND CIRCULATORS
 Face_iterator faces_begin() const
    {
        CGAL_Triangulation_2<Gt, Tds>* ncthis = (CGAL_Triangulation_2<Gt, Tds> *)this;
        return Face_iterator(ncthis);
    }
    Face_iterator faces_end() const
    {
        CGAL_Triangulation_2<Gt, Tds>* ncthis = (CGAL_Triangulation_2<Gt, Tds> *)this;
        return Face_iterator(ncthis, 1);
    }
    Vertex_iterator vertices_begin() const
    {
        CGAL_Triangulation_2<Gt, Tds>* ncthis = (CGAL_Triangulation_2<Gt, Tds>*)this;
        return Vertex_iterator(ncthis);
    }
    Vertex_iterator vertices_end() const
    {
        CGAL_Triangulation_2<Gt, Tds>* ncthis = (CGAL_Triangulation_2<Gt, Tds>*)this;
        return Vertex_iterator(ncthis,1);
    }
    Edge_iterator edges_begin() const
    {
        CGAL_Triangulation_2<Gt, Tds>* ncthis = (CGAL_Triangulation_2<Gt, Tds>*)this;
        return Edge_iterator(ncthis);
    }
    Edge_iterator edges_end() const
    {
        CGAL_Triangulation_2<Gt, Tds>* ncthis = (CGAL_Triangulation_2<Gt, Tds>*)this;
        return Edge_iterator(ncthis,1);
    }


    Line_face_circulator
    line_walk(const Point& p,
              const Point& q,
              Face_handle f = Face_handle())
    {
        CGAL_triangulation_precondition( (dimension() == 2) && (p != q) );
    
        Line_face_circulator lfc = (f == NULL)
                                   ? Line_face_circulator(p, q, this)
                                   : Line_face_circulator(p, q, f, this);
    
        if((! (lfc == NULL)) && is_infinite( lfc )){
            return Line_face_circulator();
        }
        return lfc;
    }
    
// not documented
 CGAL_Oriented_side
    oriented_side(const Point &p0, const Point &p1,
                  const Point &p2, const Point &p) const
    {
        // depends on the orientation of the vertices
        CGAL_Orientation o1 = geom_traits().orientation(p0, p1, p),
                         o2 = geom_traits().orientation(p1, p2, p),
                         o3 = geom_traits().orientation(p2, p0, p),
                         ot = geom_traits().orientation(p0, p1, p2);

        if (o1 == CGAL_COLLINEAR ||
            o2 == CGAL_COLLINEAR ||
            o3 == CGAL_COLLINEAR)  {
            if ((o1 == CGAL_COLLINEAR &&
                 collinear_between(p0, p, p1)) ||
                (o2 == CGAL_COLLINEAR &&
                 collinear_between(p1, p, p2)) ||
                (o3 == CGAL_COLLINEAR &&
                 collinear_between(p2, p, p0)))
                {
                return  CGAL_ON_ORIENTED_BOUNDARY;
            }
                // for ot == CGAL_ON_ORIENTED_BOUNDARY the following also
                // does the right thing:
                return (ot == CGAL_LEFTTURN) ? CGAL_ON_POSITIVE_SIDE
                                             : CGAL_ON_NEGATIVE_SIDE;
            }
        if (ot == CGAL_RIGHTTURN)
            {
                if (o1 == CGAL_RIGHTTURN &&
                    o2 == CGAL_RIGHTTURN &&
                    o3 == CGAL_RIGHTTURN)
                    {
                        return CGAL_ON_NEGATIVE_SIDE;
                    }
                return CGAL_ON_POSITIVE_SIDE;
            }
        if (o1 == CGAL_LEFTTURN &&
            o2 == CGAL_LEFTTURN &&
            o3 == CGAL_LEFTTURN)
            {
                return CGAL_ON_POSITIVE_SIDE;
            }
        return CGAL_ON_NEGATIVE_SIDE;
    }


    CGAL_Bounded_side
    bounded_side(const Point &p0, const Point &p1,
                 const Point &p2, const Point &p) const
    {
      CGAL_Orientation o1 = geom_traits().orientation(p0, p1, p),
        o2 = geom_traits().orientation(p1, p2, p),
        o3 = geom_traits().orientation(p2, p0, p),
        ot = geom_traits().orientation(p0, p1, p2);

      if(o1 == CGAL_COLLINEAR ||
         o2 == CGAL_COLLINEAR ||
         o3 == CGAL_COLLINEAR)
        {
          if ((o1 == CGAL_COLLINEAR &&
               collinear_between(p0, p, p1)) ||
              (o2 == CGAL_COLLINEAR &&
               collinear_between(p1, p, p2)) ||
              (o3 == CGAL_COLLINEAR &&
               collinear_between(p2, p, p0)))
            {
              return  CGAL_ON_BOUNDARY;
            }
          return CGAL_ON_UNBOUNDED_SIDE;
        }
      if (ot == CGAL_RIGHTTURN)
        {
          if(o1==CGAL_RIGHTTURN &&
             o2==CGAL_RIGHTTURN &&
             o3==CGAL_RIGHTTURN)
            {
              return CGAL_ON_BOUNDED_SIDE;
            }
          return CGAL_ON_UNBOUNDED_SIDE;
        }
      if (o1 == CGAL_LEFTTURN &&
          o2 == CGAL_LEFTTURN &&
          o3 == CGAL_LEFTTURN)
        {
          return CGAL_ON_BOUNDED_SIDE;
        }
      return CGAL_ON_UNBOUNDED_SIDE;
    }

    CGAL_Oriented_side
    oriented_side(const Face_handle& f, const Point &p) const
    {
        return oriented_side(f->vertex(0)->point(),
                             f->vertex(1)->point(),
                             f->vertex(2)->point(),
                             p);
    }

 
  
    
  // NOT DOCUMENTED
    bool collinear_between(const Point& p, const Point& q, const Point& r) const
    {
        CGAL_Comparison_result c_pr = geom_traits().compare_x(p, r);
        CGAL_Comparison_result c_pq;
        CGAL_Comparison_result c_qr;
        if(c_pr == CGAL_EQUAL) {
            c_pr = geom_traits().compare_y(p, r);
            c_pq = geom_traits().compare_y(p, q);
            c_qr = geom_traits().compare_y(q, r);
        } else {
            c_pq = geom_traits().compare_x(p, q);
            c_qr = geom_traits().compare_x(q, r);
        }
        return ( (c_pq == CGAL_SMALLER) && (c_qr == CGAL_SMALLER) ) ||
            ( (c_qr == CGAL_LARGER) && (c_pq == CGAL_LARGER) );
    
    }
    
};



template <class Gt, class Tds >
ostream&
operator<<(ostream& os, const CGAL_Triangulation_2<Gt, Tds> &tr)
{
  // to debug
  //operator<<(os, tr._tds);
  //os << tr._tds;

  map< void*, int, less<void*> > V;
  map< void*, int, less<void*> > F;
  typename CGAL_Triangulation_2<Gt, Tds>::Vertex_handle  v;

    int n = tr.number_of_vertices() + 1;
    int m = tr.number_of_faces();
    if(CGAL_is_ascii(os)){
        os << n << ' ' << m << ' ' << tr.dimension() << endl;
    } else {
        os << n << m << tr.dimension();
    }

    // write the vertex at infinity
    int i = 0;
    v = tr.infinite_vertex();
    V[&(*v)] = i;
    os << v->point();
    if(CGAL_is_ascii(os)){
        os << ' ';
    }
    if(n == 1){
        return os;
    }

    // write the finite vertices
    {
        typename CGAL_Triangulation_2<Gt, Tds>::Vertex_iterator
          it = tr.vertices_begin();

        while(it != tr.vertices_end()){
            V[&(*it)] = ++i;
            os << it->point();
            if(CGAL_is_ascii(os)){
                os << ' ';
            }
            ++it;
        }
    }
    CGAL_triangulation_assertion( (i+1) == n );
    if(CGAL_is_ascii(os)){ os << "\n";}

    if(n == 2){
        return os;
    }

    i = 0;
    // vertices of the finite faces
    {
        typename CGAL_Triangulation_2<Gt, Tds>::Face_iterator
          it = tr.faces_begin();

        while(it != tr.faces_end()){
	  F[&(*it)] = i++;
	  for(int j = 0; j < 3; j++){
	    os << V[&(*it->vertex(j))];
	    if(CGAL_is_ascii(os)){
	      if(j==2) {
		os << "\n";
	      } else {
		os <<  ' ';
	      }
	    }
	  }
	  ++it;
        }
    }

    // vertices of the infinite faces
    {
        typename CGAL_Triangulation_2<Gt, Tds>::Face_circulator
            fc = tr.infinite_vertex()->incident_faces(),
            done(fc);

        do{
            F[&(*fc)] = i++;
            for(int j = 0; j < 3; j++){
                os << V[&(*fc->vertex(j))];
                if(CGAL_is_ascii(os)){
                    if(j==2) {
                        os << "\n";
                    } else {
                        os <<  ' ';
                    }
                }
            }
        }while(++fc != done);
    }
    CGAL_triangulation_assertion( i == m );

    // neighbor pointers of the finite faces
    {
        typename CGAL_Triangulation_2<Gt, Tds>::Face_iterator
            it = tr.faces_begin();
        while(it != tr.faces_end()){
            for(int j = 0; j < 3; j++){
                os << F[&(*it->neighbor(j))];
                if(CGAL_is_ascii(os)){
                    if(j==2) {
                        os << "\n";
                    } else {
                        os <<  ' ';
                    }
                }
            }
            ++it;
        }
    }

    // neighbor pointers of the infinite faces
    {
        typename CGAL_Triangulation_2<Gt, Tds>::Face_circulator
            fc = tr.infinite_vertex()->incident_faces(),
            done(fc);

        do{
            for(int j = 0; j < 3; j++){
                os << F[&(*fc->neighbor(j))];
                if(CGAL_is_ascii(os)){
                    if(j==2) {
                        os << "\n";
                    } else {
                        os <<  ' ';
                    }
                }
            }
        }while(++fc != done);
    }

    
    return os ;
}



template < class Gt, class Tds >
istream&
operator>>(istream& is, CGAL_Triangulation_2<Gt, Tds> &tr)
{

  
  return operator>>(is, tr._tds);
}
    

#endif CGAL_TRIANGULATION_2_H

