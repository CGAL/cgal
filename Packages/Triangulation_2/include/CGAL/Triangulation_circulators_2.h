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
// file          : Triangulation/include/CGAL/Triangulation_circulators_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_CIRCULATORS_2_H
#define CGAL_TRIANGULATION_CIRCULATORS_2_H

#include <pair.h>
#include <iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>

template <class Gt, class Tds >
class CGAL_Triangulation_2;

template <class Gt,  class Tds >
class CGAL_Triangulation_face_2;

template <class Gt,  class Tds >
class CGAL_Triangulation_vertex_2;

template < class Gt, class Tds >
class CGAL_Triangulation_face_circulator_2;

template < class Gt, class Tds >
class CGAL_Triangulation_vertex_circulator_2;

template < class Gt, class Tds >
class CGAL_Triangulation_edge_circulator_2;

template < class Gt, class Tds>
class CGAL_Triangulation_face_circulator_2
    : public 
    CGAL_Bidirectional_circulator_base<CGAL_Triangulation_face_2<Gt,Tds>,
					ptrdiff_t,size_t>
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Face_circulator  Base_face_circulator;

  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Gt, Tds> Vertex;
  typedef CGAL_Triangulation_vertex_handle_2<Gt, Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt, Tds>  Face_handle;
  typedef pair<Face_handle, int>     Edge;

  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;


  static int ccw(int i)
  {
        return (i+1) % 3;
  }

  static int cw(int i)
  {
    return (i+2) % 3;
  }

  CGAL_Triangulation_face_circulator_2()
    : _bfc()
  {}

  CGAL_Triangulation_face_circulator_2(const Vertex_handle& v)
    : _bfc(&(*v))
  {
    if (v->face()==NULL  ||
        v->face()->dimension() == 0 || 
	v->face()->dimension() == 1) { _bfc = Base_face_circulator();}
  }

  CGAL_Triangulation_face_circulator_2(const Vertex_handle& v,
				       const Face_handle& f)
    : _bfc( &(*v), &(*f))
  {
    if (v->face()==NULL  ||
        v->face()->dimension() == 0 || 
	v->face()->dimension() == 1) { _bfc = Base_face_circulator();}
  }

 
   
   CGAL_Triangulation_face_circulator_2(const Face_circulator &fc)
    : _bfc(fc._bfc)
  {}
   
  CGAL_Triangulation_face_circulator_2 &operator=(const Face_circulator &fc)
  {
    _bfc = fc._bfc;
    return *this;
  }
  
  Face_circulator& operator++()
  {
    _bfc++;
      return *this;
  }
        
  Face_circulator operator++(int)
  {
    Face_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Face_circulator& operator--()
  {
    _bfc--;
    return *this;
  }
        
  Face_circulator operator--(int)
  {
    Face_circulator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Face_circulator &fc) const
  {
            return (_bfc == fc._bfc) ;
  }
        
        
  bool operator!=(const Face_circulator &fc) const
  {
            return ! (*this == fc);
  }
   
  bool is_empty() const
  {
    return( _bfc.is_empty());
  }


  bool  operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_bfc == NULL);
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return ! (*this == NULL);
  }

  
        
  inline Face& operator*() const
  {
    return (Face&  )(*_bfc);
  }

  inline Face* operator->() const
  {
    return   (Face*)( & (*_bfc));
  }
  
  
private:
  Base_face_circulator  _bfc;
};


template < class Gt, class Tds>
class CGAL_Triangulation_vertex_circulator_2
    : public CGAL_Bidirectional_circulator_base<CGAL_Triangulation_vertex_2<Gt,Tds>,ptrdiff_t,size_t>
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Vertex_circulator  Base_vertex_circulator;

  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef CGAL_Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef pair<Face_handle, int>     Edge;

  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
   typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;


  static int ccw(int i)
  {
        return (i+1) % 3;
  }

  static int cw(int i)
  {
    return (i+2) % 3;
  }

  CGAL_Triangulation_vertex_circulator_2()
    : _bvc()
  {}

  CGAL_Triangulation_vertex_circulator_2(const Vertex_handle& v)
    : _bvc(&(*v))
  {
    if (v->face()==NULL  ||
        v->face()->dimension() == 0 )	 { _bvc = Base_vertex_circulator();}
  }

  CGAL_Triangulation_vertex_circulator_2(const Vertex_handle& v,
				       const Face_handle& f)
    : _bvc( &(*v), &(*f))
  {
    if (v->face()==NULL  ||
        v->face()->dimension() == 0 )	 { _bvc = Base_vertex_circulator();}
  }
   
   CGAL_Triangulation_vertex_circulator_2(const Vertex_circulator &vc)
    : _bvc(vc._bvc)
  {}
   
  CGAL_Triangulation_vertex_circulator_2 &
	operator=(const Vertex_circulator &vc)
  {
    _bvc = vc._bvc;
    return *this;
  }
  
  Vertex_circulator& operator++()
  {
    _bvc++;
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
    _bvc--;
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
            return (_bvc == vc._bvc) ;
  }
        
        
  bool operator!=(const Vertex_circulator &vc) const
  {
            return ! (*this == vc);
  }
  
  bool is_empty() const
  {
    return( _bvc.is_empty());
  }


      
  bool  operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_bvc == (CGAL_NULL_TYPE) NULL);
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return ! (*this == (CGAL_NULL_TYPE) NULL);
  }
        
  inline Vertex& operator*() const
  {
    return (Vertex &)(*_bvc);
  }

  inline Vertex* operator->() const
  {
    return   (Vertex*)( & (*_bvc));
  }
  
  
private:
  Base_vertex_circulator  _bvc;
};


template < class Gt, class Tds>
class CGAL_Triangulation_edge_circulator_2
  : public 
   CGAL_Bidirectional_circulator_base< 
          typename CGAL_Triangulation_2<Gt, Tds>::Edge,ptrdiff_t,size_t>
{
public:
  typedef Tds Triangulation_data_structure;

  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Edge_circulator  Base_edge_circulator;

  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef CGAL_Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef pair<Face_handle, int>     Edge;

  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;

  static int ccw(int i)
  {
        return (i+1) % 3;
  }

  static int cw(int i)
  {
    return (i+2) % 3;
  }

  CGAL_Triangulation_edge_circulator_2()
    : _bec()
  {}

  CGAL_Triangulation_edge_circulator_2(const Vertex_handle& v)
    : _bec(&(*v))
  {
    if (v->face()==NULL  ||
        v->face()->dimension() == 0 )	 { _bec = Base_face_circulator();}
  }

  CGAL_Triangulation_edge_circulator_2(const Vertex_handle& v,
				       const Face_handle& f)
    : _bec( &(*v), &(*f))
  {
    if (v->face()==NULL  ||
        v->face()->dimension() == 0 )	 { _bec = Base_face_circulator();}
  }
   
   CGAL_Triangulation_edge_circulator_2(const Edge_circulator &ec)
    : _bec(ec._bec)
  {}
   
  CGAL_Triangulation_edge_circulator_2 &operator=(const Edge_circulator &ec)
  {
    _bec = ec._bec;
    return *this;
  }
  
  Edge_circulator& operator++()
  {
    _bec++;
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
    _bec--;
    return *this;
  }
        
  Edge_circulator operator--(int)
  {
    Edge_circulator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Edge_circulator &ec) const
  {
            return (_bec == ec._bec) ;
  }
        
        
  bool operator!=(const Edge_circulator &ec) const
  {
            return ! (*this == ec);
  }
        
   bool is_empty() const
  {
    return( _bec.is_empty());
  }

  bool  operator==(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return (_bec == NULL);
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return ! (*this == NULL);
  }
        
  inline Edge operator*() const
  {
    Face_handle fh = (Face *) ((*_bec).first);
    return make_pair( fh  , (*_bec).second );
  }

    
private:
  Base_edge_circulator  _bec;
};


#endif CGAL_TRIANGULATION_CIRCULATORS_2_H
