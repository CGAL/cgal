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

private: 
  const Vertex* _v;
  const Face* pos;

public:
  Triangulation_ds_face_circulator_2()
    : _v(NULL), pos(NULL)
  {}
  
  Triangulation_ds_face_circulator_2(const Vertex* v, const Face* f=NULL);
        
  Face_circulator& operator++();
  Face_circulator operator++(int);
  Face_circulator& operator--();
  Face_circulator operator--(int);
  Face& operator*() const  ;
  Face* operator->() const ;
   
  bool operator==(const Face_circulator &fc) const ;
  bool operator!=(const Face_circulator &fc) const;
  bool is_empty() const;
  bool operator==(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const;
  bool operator!=(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const;
  
  
};


template < class Vertex, class Face >
class Triangulation_ds_vertex_circulator_2 :
  public Bidirectional_circulator_base<Vertex, 
                                       CGAL_STD::ptrdiff_t,
                                       CGAL_STD::size_t>,
  public  Triangulation_cw_ccw_2
{
public:
  typedef Triangulation_ds_vertex_circulator_2<Vertex, Face> 
                                          Vertex_circulator;
private:
  const Vertex* _v;
  const Face* pos;
  int _ri;
  
public:
  Triangulation_ds_vertex_circulator_2()
    :  _v(NULL), pos(NULL)
  {}
                
  Triangulation_ds_vertex_circulator_2(const Vertex* v,const Face* f = NULL);
       
  Vertex_circulator& operator++();
  Vertex_circulator  operator++(int);
  Vertex_circulator& operator--();
  Vertex_circulator  operator--(int);
  Vertex& operator*() const ;
  Vertex* operator->() const; 
 
  bool operator==(const Vertex_circulator &vc) const;
  bool operator!=(const Vertex_circulator &vc) const;
  bool is_empty() const;
  bool operator==(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const;
  bool operator!=(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const;
};



template < class Vertex, class Face >
class Triangulation_ds_edge_circulator_2 :
  public Bidirectional_circulator_base  <CGAL_STD::pair<Face*,int>, 
                                         CGAL_STD::ptrdiff_t,
                                         CGAL_STD::size_t>,
  public Triangulation_cw_ccw_2
{
public:
  typedef Triangulation_ds_edge_circulator_2<Vertex,Face> Edge_circulator;
  typedef std::pair<Face*, int>         Edge;

private:
  int _ri;
  const Vertex* _v;
  const Face* pos;

public:
  Triangulation_ds_edge_circulator_2()
    : _v(NULL), pos(NULL)
  {}
            
   Triangulation_ds_edge_circulator_2( const Vertex* v, const Face* f=NULL);

  Edge operator*() const ;
  Edge_circulator& operator++();
  Edge_circulator operator++(int);
  Edge_circulator& operator--();
  Edge_circulator operator--(int);
 
  bool operator==(const Edge_circulator &vc) const;
  bool operator!=(const Edge_circulator &vc) const;
  bool is_empty() const;
  bool operator==(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const;
  bool operator!=(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const;
 };


template < class Vertex, class Face >
Triangulation_ds_face_circulator_2<Vertex,Face> ::
Triangulation_ds_face_circulator_2(const Vertex* v, const Face* f)
  : _v(v), pos(f)
{
  if (_v == NULL) pos = NULL;
  else if ( pos == NULL) pos = v->face();

  if (pos == NULL || pos->dimension() < 2) { _v = NULL; pos = NULL; return;}
  else CGAL_triangulation_precondition( pos->has_vertex(v));
}
  
    
template < class Vertex, class Face >
Triangulation_ds_face_circulator_2<Vertex,Face>&
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator++()
{
  CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
  int i = pos->index(_v);
  pos = pos->neighbor(ccw(i));
  return *this;
}

template < class Vertex, class Face >
Triangulation_ds_face_circulator_2<Vertex,Face>
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator++(int)
{
  CGAL_triangulation_precondition( (pos  != NULL) && (_v != NULL) );
  Face_circulator tmp(*this);
  ++(*this);
  return tmp;
}

template < class Vertex, class Face >
Triangulation_ds_face_circulator_2<Vertex,Face>&
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator--()
{
   CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
   int i = pos->index(_v);
   pos = pos->neighbor(cw(i));
   return *this;
}

template < class Vertex, class Face >
Triangulation_ds_face_circulator_2<Vertex,Face>
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator--(int)
{
  CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
  Face_circulator tmp(*this);
  --(*this);
  return tmp;
}

template < class Vertex, class Face >
inline Face&
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator*() const
{
  CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
  return const_cast<Face&>(*pos);
}

template < class Vertex, class Face >
inline Face*
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator->() const
{
  CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
  return const_cast<Face *>(pos);
}


template < class Vertex, class Face >
inline bool
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator==(const Face_circulator &fc) const
{    
  return (_v == fc._v) &&  (pos == fc.pos);  
}

template < class Vertex, class Face >
inline bool
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator!=(const Face_circulator &fc) const
{    
return ! (*this == fc);  
}
   
template < class Vertex, class Face >
inline bool
Triangulation_ds_face_circulator_2<Vertex,Face> ::
is_empty() const
{    
return ((_v == NULL) || (pos == NULL));  
}

template < class Vertex, class Face >
inline bool
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator==(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const
{
  CGAL_triangulation_assertion( n == NULL);
  return (_v == NULL || pos == NULL);
}
        
template < class Vertex, class Face >
inline bool
Triangulation_ds_face_circulator_2<Vertex,Face> ::
operator!=(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const
{
  CGAL_triangulation_assertion( n == NULL);
  return ! (*this == NULL);
}

template < class Vertex, class Face >
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
Triangulation_ds_vertex_circulator_2 (const Vertex* v,  const Face* f)
  : _v( v ), pos(f)
{
  if (_v == NULL) { pos = NULL;}
  else if (pos==NULL) {pos = v->face();}

  if (pos == NULL || pos->dimension() < 1){_v = NULL; pos = NULL;return;}
  int i = pos->index( _v );
  if (pos->dimension() == 2) {_ri = ccw(i);}
  else {_ri = 1-i;}
  return;
}


template < class Vertex, class Face >
Triangulation_ds_vertex_circulator_2<Vertex,Face>&
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
operator++()
{
  CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
  int i = pos->index(_v);
    
  if (pos->dimension() == 1) { 
    pos = pos->neighbor(1-i);
    _ri = 1 - pos->index(_v);
  }
  else{
    pos = pos->neighbor(ccw(i));
    i = pos->index(_v);
    _ri = ccw(i);
  }
  return *this;
}
        
template < class Vertex, class Face >
Triangulation_ds_vertex_circulator_2<Vertex,Face>
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
operator++(int)
{
  Vertex_circulator tmp(*this);
  ++(*this);
  return tmp;
}

template < class Vertex, class Face >
Triangulation_ds_vertex_circulator_2<Vertex,Face>&
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
operator--()
{
  CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
  int i = pos->index(_v);
    
  if (pos->dimension() == 1) { 
    pos = pos->neighbor(1-i);
    _ri = 1 - pos->index(_v);
  }
  else{
    pos = pos->neighbor(cw(i));
    i = pos->index(_v);
    _ri = ccw(i);
  }
  return *this;
}
        
template < class Vertex, class Face >
Triangulation_ds_vertex_circulator_2<Vertex,Face>
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::     
operator--(int)
{
  Vertex_circulator tmp(*this);
  --(*this);
  return tmp;
}


template < class Vertex, class Face >
inline Vertex&
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
operator*() const
{
   CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
   return *(pos->vertex(_ri));
}

template < class Vertex, class Face >
inline Vertex*
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
operator->() const
{
   CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
   return pos->vertex(_ri);
}

template < class Vertex, class Face >
inline bool
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
operator==(const Vertex_circulator &vc) const
{
  return (_v == vc._v) &&  (_ri == vc._ri) && (pos == vc.pos);
}
               
template < class Vertex, class Face >
inline bool
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
operator!=(const Vertex_circulator &vc) const
{
  return ! (*this == vc);
}
               
template < class Vertex, class Face >
inline bool
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
is_empty() const
{
  return ((_v == NULL) || (pos == NULL));
}

template < class Vertex, class Face >
inline bool
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::
operator==(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const
{
  CGAL_triangulation_assertion( n == NULL);
  return (_v == NULL) || (pos == NULL);
}
        
template < class Vertex, class Face >
inline bool
Triangulation_ds_vertex_circulator_2<Vertex,Face> ::        
operator!=(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const
{
  CGAL_triangulation_assertion( n == NULL);
  return !(*this == NULL);
}
        

template < class Vertex, class Face >
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
Triangulation_ds_edge_circulator_2(const Vertex* v, const Face* f)
    : _v(v), pos(f)
{
  if (_v == NULL) pos = NULL;
  else if ( pos == NULL) pos = v->face();
 
  if (pos == NULL || pos->dimension() < 1){_v = NULL; pos = NULL;return;}
  int i = pos->index( _v );
  if (pos->dimension() == 2) {_ri = ccw(i);}
  else {_ri = 2;}
  return;
}


template < class Vertex, class Face >
inline std::pair<Face*, int> 
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator*() const
{
   CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
   return std::make_pair(const_cast<Face*>(pos), _ri);
}


template < class Vertex, class Face >
Triangulation_ds_edge_circulator_2<Vertex,Face>&
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator++()
{
  CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
  int i = pos->index(_v);
  if (pos->dimension() == 1) { 
    pos = pos->neighbor(1-i);
    return *this;
  }
  else{
    pos = pos->neighbor(ccw(i));
    i = pos->index(_v);
    _ri = ccw(i);
  }    
  return *this;
}

template < class Vertex, class Face >
Triangulation_ds_edge_circulator_2<Vertex,Face>
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator++(int)
{
  Edge_circulator tmp(*this);
  ++(*this);
  return tmp;
}

template < class Vertex, class Face >
Triangulation_ds_edge_circulator_2<Vertex,Face>&
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator--()
{
  CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
  int i = pos->index(_v);

  if (pos->dimension() == 1) { 
    pos = pos->neighbor(1-i);
    return *this;
  }
  else{
    pos = pos->neighbor(cw(i));
    i = pos->index(_v);
    _ri = ccw(i);
  }    
  return *this;
}
   
template < class Vertex, class Face >
Triangulation_ds_edge_circulator_2<Vertex,Face>
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator--(int)
{
  Edge_circulator tmp(*this);
  --(*this);
  return tmp;
}

template < class Vertex, class Face >
inline bool
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator==(const Edge_circulator &vc) const
{
  return (_v == vc._v) &&  (_ri == vc._ri) && (pos == vc.pos);
}
                
template < class Vertex, class Face >
inline bool
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator!=(const Edge_circulator &vc) const
{
  return ! (*this == vc);
}

template < class Vertex, class Face >
inline bool
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
is_empty() const
{
  return ((_v == NULL) || (pos == NULL));
}

template < class Vertex, class Face >
inline bool
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator==(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const
{
  CGAL_triangulation_assertion( n == NULL);
  return (_v == NULL) || (pos == NULL);
}
               
template < class Vertex, class Face >
inline bool
Triangulation_ds_edge_circulator_2<Vertex,Face> ::
operator!=(CGAL_NULL_TYPE CGAL_triangulation_assertion_code(n)) const
{
  CGAL_triangulation_assertion( n == NULL);
  return !(*this == NULL);
}
       

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DS_CIRCULATORS_2_H
