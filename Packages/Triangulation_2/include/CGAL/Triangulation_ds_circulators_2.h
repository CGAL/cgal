 #ifndef CGAL_TRIANGULATION_DS_CIRCULATORS_2_H
#define CGAL_TRIANGULATION_DS_CIRCULATORS_2_H



#include <pair.h>
#include <iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>

template < class Vertex, class Face >
class CGAL_Triangulation_ds_face_circulator_2
    : public CGAL_Bidirectional_circulator_base<Face,ptrdiff_t,size_t>
{
public:
  
  typedef CGAL_Triangulation_ds_face_circulator_2<Vertex,Face> Face_circulator;

   static int ccw(int i)
    {
        return (i+1) % 3;
    }


    static int cw(int i)
    {
        return (i+2) % 3;
    }


  CGAL_Triangulation_ds_face_circulator_2()
    : pos(NULL), _v(NULL)
  {}
        
  CGAL_Triangulation_ds_face_circulator_2(Vertex* v)
    : pos(v->face()), _v(v)
  {}

  CGAL_Triangulation_ds_face_circulator_2(Vertex* v,
				     Face* f)
          : pos(f),_v(v)
  {}
        
        
  CGAL_Triangulation_ds_face_circulator_2(const Face_circulator &fc)
    : pos(fc.pos),_v(fc._v)
  {}
        
   
  Face_circulator& operator++()
  {
    CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
    int i = pos->index(_v);
    pos = pos->neighbor(ccw(i));
    if(pos == NULL){
      pos = _v->face();
    }
    return *this;
  }
        
  Face_circulator operator++(int)
  {
    CGAL_triangulation_precondition( (pos  != NULL) && (_v != NULL) );
    Face_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Face_circulator& operator--()
  {
    CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
    int i = pos->index(_v);
    Face* f = pos->neighbor(cw(i));
    if(f == NULL) {
      Face* n;
      do{
	i = pos->index(_v);
	//	n = pos->neighbor(cw(i));
	n = pos->neighbor(ccw(i));
	if(n != NULL){
	  pos = n;
	}
      }while(n != NULL);
    }
    else{
      pos = f;
    }
    return *this;
  }
        
        
  Face_circulator operator--(int)
  {
    CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
    Face_circulator tmp(*this);
    --(*this);
    return tmp;
  }

  inline Face& operator*() const
  {
    return *pos;
  }

  inline Face* operator->() const
  {
    return pos;
  }
  
  bool operator==(const Face_circulator &fc) const
  {
    return (_v == fc._v) &&  (pos == fc.pos);
  }
        
        
  bool operator!=(const Face_circulator &fc) const
  {
    return ! (*this == fc);
  }
        
  bool
  operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_v == NULL) && (pos == NULL);
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return ! (*this == NULL);
  }
        
private:
  Vertex* _v;
  Face* pos;

};


template < class Vertex, class Face >
class CGAL_Triangulation_ds_vertex_circulator_2 :
    public CGAL_Bidirectional_circulator_base<Vertex,ptrdiff_t,size_t>
{
public:
  
  typedef CGAL_Triangulation_ds_vertex_circulator_2<Vertex, Face> Vertex_circulator;
    
  static int ccw(int i)
  {
        return (i+1) % 3;
  }

  static int cw(int i)
  {
    return (i+2) % 3;
  }

  CGAL_Triangulation_ds_vertex_circulator_2()
    :  _v(), _f(NULL)
  {}
                
  CGAL_Triangulation_ds_vertex_circulator_2(Vertex* v,
				       Face*   f)
    : _v( v ), _f(f)
  {
    if( _f != NULL ) {
      int i = _f->index( _v );
      _ri = ccw(i);
    }
  }
        
  CGAL_Triangulation_ds_vertex_circulator_2(Vertex* v)
				      
    : _v( v ), _f(v->face())
  {
    if( _f != NULL ) {
      int i = _f->index( _v );
      _ri = ccw(i);
    }
  }

  CGAL_Triangulation_ds_vertex_circulator_2(const Vertex_circulator &vc)
    : _ri(vc._ri), _v(vc._v), _f(vc._f)
  {}
        
              
  Vertex_circulator &operator=(const Vertex_circulator &vc)
  {
    _v = vc._v;
    _ri = vc._ri;
    _f = vc._f;
    return *this;
  }   

  inline Vertex& operator*() const
  {
    return *(_f->vertex(_ri));
  }

  inline Vertex* operator->() const
  {
    return _f->vertex(_ri);
  }

  Vertex_circulator& operator++()
  {
    Face* n = _f;
    int i = _f->index(_v);
    n = _f->neighbor(ccw(i));
    if(n == NULL){
      if(_ri == ccw(i)){
	_ri = cw(i);
      } else {
	_f = _v->face();
	i = _f->index(_v);
	_ri = ccw(i);
      }
    } 
    else {
      _f = n;
      i = _f->index(_v);
      _ri = ccw(i);
    }
    return *this;
  }
        
  Vertex_circulator operator++(int)
  {
    Vertex_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
 
  Vertex_circulator& operator--()
  {
    Face* n = _f;
    int i = _f->index(_v);
    n = _f->neighbor(cw(i));
    if(n == NULL) {
      do{
	n = _f->neighbor(ccw(i));
	if(n != NULL){
	  _f = n;
	  i = _f->index(_v);
	}
      }while(n != NULL);
      i = _f->index(_v);
      _ri = cw(i);
    } 
    else {
      if(_ri == cw(i)){
	_ri = ccw(i);
      } 
      else {
	_f = n;
	i = _f->index(_v);
	_ri = ccw(i);
      }
    }
    return *this;
  }
        
        
  Vertex_circulator operator--(int)
  {
    Vertex_circulator tmp(*this);
    --(*this);
    return tmp;
  }
        

  bool operator==(const Vertex_circulator &vc) const
  {
    return (_v == vc._v) &&  (_ri == vc._ri) && (_f == vc._f);
  }
               
  bool operator!=(const Vertex_circulator &vc) const
  {
    return ! (*this == vc);
  }
               
  bool operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_v == NULL) && (_f == NULL);
  }
        
        
  bool operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return !(*this == NULL);
  }
        
private:
  int _ri;
  Vertex* _v;
  Face* _f;
        
};



template < class Vertex, class Face >
class CGAL_Triangulation_ds_edge_circulator_2 :
    public CGAL_Bidirectional_circulator_base<
                    pair<Face*,int>,ptrdiff_t,size_t>
{
public:
  typedef CGAL_Triangulation_ds_edge_circulator_2<Vertex,Face> Edge_circulator;
  typedef pair<Face*, int>         Edge;

  static int ccw(int i)
  {
    return (i+1) % 3;
  }

  static int cw(int i)
  {
    return (i+2) % 3;
  }

  CGAL_Triangulation_ds_edge_circulator_2()
    : _v(NULL), _f(NULL)
  {}
            
   CGAL_Triangulation_ds_edge_circulator_2( Vertex* v)
				       
    : _v(v), _f(v->face())
  {
    if( _f != NULL ){
      int i = _f->index(_v);
      _ri = ccw(i);
    }
  }

  CGAL_Triangulation_ds_edge_circulator_2( Vertex* v,
				        Face*   f)
    : _v(v), _f(f)
  {
    if( _f != NULL ){
      int i = _f->index(_v);
      _ri = ccw(i);
    }
  }
        
  CGAL_Triangulation_ds_edge_circulator_2(const Edge_circulator &vc)
    : _ri(vc._ri), _v(vc._v), _f(vc._f)
  {}
        
  Edge_circulator &operator=(const Edge_circulator &vc)
  {
    _v = vc._v;
    _ri = vc._ri;
    _f = vc._f;
    return *this;
  }

  Edge&
  operator*()
  {
    if( _f == NULL) {
      return make_pair(_f, 0);
    }
    return make_pair(_f, _ri);
  }
        
  Edge_circulator& operator++()
  {
    Face* n = _f;
    int i = _f->index(_v);
    n = _f->neighbor(ccw(i));
    if(n == NULL){
      if(_ri == ccw(i)){
	_ri = cw(i);
      } 
      else {
	_f = _v->face();
	i = _f->index(_v);
	_ri = ccw(i);
      }
    } 
    else {
      _f = n;
      i = _f->index(_v);
      _ri = ccw(i);
    }
    return *this;
  }
        
  Edge_circulator operator++(int)
  {
    Edge_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Edge_circulator& operator--()
  {
    Face* n = _f;
    int i = _f->index(_v);
    n = _f->neighbor(cw(i));
    if(n == NULL){
      do{
	n = _f->neighbor(ccw(i));
	if(n != NULL){
	  _f = n;
	  i = _f->index(_v);
	}
      }while(n != NULL);
      i = _f->index(_v);
      _ri = cw(i);
    } else {
      if(_ri == cw(i)){
	_ri = ccw(i);
      } else {
	_f = n;
	i = _f->index(_v);
	_ri = ccw(i);
      }
    }
    return *this;
  }
        
  Edge_circulator operator--(int)
  {
    Edge_circulator tmp(*this);
    --(*this);
    return tmp;
  }
        
   bool operator==(const Edge_circulator &vc) const
  {
    return (_v == vc._v) &&  (_ri == vc._ri) && (_f == vc._f);
  }
                
  bool operator!=(const Edge_circulator &vc) const
  {
    return ! (*this == vc);
  }

  bool operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_v == NULL) && (_f == NULL);
  }
               
  bool operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return !(*this == NULL);
  }
     
private:
  int _ri;
  Vertex* _v;
  Face* _f;
        
};
        

#endif CGAL_TRIANGULATION_DS_CIRCULATORS_2_H
