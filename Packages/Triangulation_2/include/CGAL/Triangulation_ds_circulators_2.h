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
// file          : Triangulation/include/CGAL/Triangulation_ds_circulators_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_CIRCULATORS_2_H
#define CGAL_TRIANGULATION_DS_CIRCULATORS_2_H



#include <utility>
#include <iterator>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE

template < class Vertex, class Face >
class Triangulation_ds_face_circulator_2
  : public Bidirectional_circulator_base<Face,
                                 CGAL_STD::ptrdiff_t,
				 CGAL_STD::size_t>,
    public Triangulation_cw_ccw_2
      
{
public:
    typedef Triangulation_ds_face_circulator_2<Vertex,Face> Face_circulator;

  // typedef Face                           value_type;
//   typedef std:: ptrdiff_t                difference_type;
//   typedef std::size_t                    size_type;
//   typedef Face*                          pointer;
//   typedef const Face*                    const_pointer;
//   typedef Face&                          reference;
//   typedef const Face&                    const_reference;
//   typedef CGAL::Bidirectional_circulator_tag iterator_category;

private: 
  Vertex* _v;
  Face* pos;

public:
  Triangulation_ds_face_circulator_2()
    : _v(NULL), pos(NULL)
  {}
        
  Triangulation_ds_face_circulator_2(Vertex* v)
    :  _v(v), pos(v->face())
  {}


  Triangulation_ds_face_circulator_2(Vertex* v,    Face* f)
          :_v(v), pos(f)
  {
    CGAL_triangulation_precondition( f->has_vertex(v));
  }

        
  Triangulation_ds_face_circulator_2(const Face_circulator &fc)
    : _v(fc._v), pos(fc.pos)
  {}
        
   
  Face_circulator& operator++()
  {
    CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
    //then dimension() cannot be 0
    int i = pos->index(_v);

    if(pos->dimension() == 1) {pos = pos->neighbor(1-i);} //for one dim case
    else{pos = pos->neighbor(ccw(i));}
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
    if(pos->dimension() == 1) {pos = pos->neighbor(1-i);} //for one dim case
    else {pos = pos->neighbor(cw(i));}
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
   
  bool is_empty() const
  {
    return ((_v == NULL) || (pos == NULL));
  }

  bool
  operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_v == NULL || pos == NULL);
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return ! (*this == NULL);
  }
        

};


template < class Vertex, class Face >
class Triangulation_ds_vertex_circulator_2 :
  public Bidirectional_circulator_base<Vertex, 
                       CGAL_STD::ptrdiff_t,CGAL_STD::size_t>,
  public  Triangulation_cw_ccw_2
{
public:
  // typedef Vertex                           value_type;
//   typedef std:: ptrdiff_t                difference_type;
//   typedef std::size_t                    size_type;
//   typedef Vertex*                          pointer;
//   typedef const Vertex*                    const_pointer;
//   typedef Vertex&                          reference;
//   typedef const Vertex&                    const_reference;
//   typedef CGAL::Bidirectional_circulator_tag iterator_category;
  typedef Triangulation_ds_vertex_circulator_2<Vertex, Face> 
                                          Vertex_circulator;
    

private:
  Vertex* _v;
  Face* _f;
  int _ri;
  int _dim;


public:
  Triangulation_ds_vertex_circulator_2()
    :  _v(NULL), _f(NULL), _dim(0)
  {}
                
  Triangulation_ds_vertex_circulator_2( Vertex* v,  Face* f)
    : _v( v ), _f(f)
  {
    if( _f != NULL ) {
      int i = _f->index( _v );
      if (_f->dimension() == 2) {_ri = ccw(i);}
      else {_ri = 1-i;}
    }
  }

  Triangulation_ds_vertex_circulator_2(Vertex* v)
				      
    : _v( v ), _f(v->face())
  {
    if( _f != NULL ) {
      int i = _f->index( _v );
      if (_f->dimension() == 2) {_ri = ccw(i);}
      else {_ri = 1-i;}
    }
  }
        
  Triangulation_ds_vertex_circulator_2(const Vertex_circulator &vc)
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
     CGAL_triangulation_precondition( (_f != NULL) && (_v != NULL) );
     int i = _f->index(_v);
    
    if (_f->dimension() == 1) { 
      _f = _f->neighbor(1-i);
      _ri = 1 - _f->index(_v);
    }
    else{
      _f = _f->neighbor(ccw(i));
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
    CGAL_triangulation_precondition( (_f != NULL) && (_v != NULL) );
     int i = _f->index(_v);
    
    if (_f->dimension() == 1) { 
      _f = _f->neighbor(1-i);
      _ri = 1 - _f->index(_v);
    }
    else{
      _f = _f->neighbor(cw(i));
      i = _f->index(_v);
      _ri = ccw(i);
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
               
  bool is_empty() const
  {
    return ((_v == NULL) || (_f == NULL));
  }


 bool operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_v == NULL) || (_f == NULL);
  }
        
        
  bool operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return !(*this == NULL);
  }
        
};



template < class Vertex, class Face >
class Triangulation_ds_edge_circulator_2 :
  public Bidirectional_circulator_base
             <CGAL_STD::pair<Face*,int>, CGAL_STD::ptrdiff_t,std::size_t>,
  public Triangulation_cw_ccw_2
{
public:
  typedef Triangulation_ds_edge_circulator_2<Vertex,Face> Edge_circulator;
  typedef std::pair<Face*, int>         Edge;

  // typedef Edge                           value_type;
//   typedef std:: ptrdiff_t                difference_type;
//   typedef std::size_t                    size_type;
//   typedef Edge*                          pointer;
//   typedef const Edge*                    const_pointer;
//   typedef Edge&                          reference;
//   typedef const Edge&                    const_reference;
//   typedef CGAL::Bidirectional_circulator_tag iterator_category;

private:
  int _ri;
  Vertex* _v;
  Face* _f;

public:
  Triangulation_ds_edge_circulator_2()
    : _v(NULL), _f(NULL)
  {}
            
   Triangulation_ds_edge_circulator_2( Vertex* v)
				       
    : _v(v), _f(v->face())
  {
    if( _f != NULL ){
      int i = _f->index(_v);
      if (_f->dimension() == 2) {_ri = ccw(i);}
      else {_ri = 1-i;}
    }
  }

  Triangulation_ds_edge_circulator_2( Vertex* v, Face* f)
    : _v(v), _f(f)
  {
    if( _f != NULL ){
      int i = _f->index(_v);
      if (_f->dimension() == 2) {_ri = ccw(i);}
      else {_ri = 1-i;}
    }
  }
        
  Triangulation_ds_edge_circulator_2(const Edge_circulator &vc)
    : _ri(vc._ri), _v(vc._v), _f(vc._f)
  {}
        
  Edge_circulator &operator=(const Edge_circulator &vc)
  {
    _v = vc._v;
    _ri = vc._ri;
    _f = vc._f;
    return *this;
  }


  Edge
  operator*() const
  {
    if( _f == NULL) {
      return std::make_pair(_f, 0);
    }
     return std::make_pair(_f, _ri);
  }
        
  Edge_circulator& operator++()
  {
    CGAL_triangulation_precondition( (_f != NULL) && (_v != NULL) );
    int i = _f->index(_v);
    if (_f->dimension() == 1) { 
      _f = _f->neighbor(1-i);
      _ri = 1 - _f->index(_v);
      return *this;
    }
    else{
     _f = _f->neighbor(ccw(i));
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
    CGAL_triangulation_precondition( (_f != NULL) && (_v != NULL) );
    int i = _f->index(_v);

    if (_f->dimension() == 1) { 
      _f = _f->neighbor(1-i);
      _ri = 1 - _f->index(_v);
      return *this;
    }
    else{
     _f = _f->neighbor(cw(i));
      i = _f->index(_v);
      _ri = ccw(i);
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

   bool is_empty() const
  {
    return ((_v == NULL) || (_f == NULL));
  }

bool operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_v == NULL) || (_f == NULL);
  }
               
  bool operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return !(*this == NULL);
  }
     

        
};
        
CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DS_CIRCULATORS_2_H
