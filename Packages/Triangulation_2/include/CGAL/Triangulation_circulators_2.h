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

#include <utility>
#include <iterator>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>

CGAL_BEGIN_NAMESPACE 

template <class Gt, class Tds >
class Triangulation_2;

template <class Gt,  class Tds >
class Triangulation_face_2;

template <class Gt,  class Tds >
class Triangulation_vertex_2;

template < class Gt, class Tds >
class Triangulation_face_circulator_2;

template < class Gt, class Tds >
class Triangulation_vertex_circulator_2;

template < class Gt, class Tds >
class Triangulation_edge_circulator_2;

template < class Gt, class Tds>
class Triangulation_face_circulator_2
    : public 
    Bidirectional_circulator_base<Triangulation_face_2<Gt,Tds>,
				   CGAL_STD::ptrdiff_t,CGAL_STD::size_t>,
      public Triangulation_cw_ccw_2
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Face_circulator  Base_face_circulator;

  typedef Triangulation_face_2<Gt,Tds> Face;
  typedef Triangulation_vertex_2<Gt, Tds> Vertex;
  typedef Triangulation_vertex_handle_2<Gt, Tds> Vertex_handle;
  typedef Triangulation_face_handle_2<Gt, Tds>  Face_handle;
  typedef std::pair<Face_handle, int>     Edge;

  typedef Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;


  Triangulation_face_circulator_2()
    : _bfc()
  {}

  Triangulation_face_circulator_2(const Vertex_handle& v)
    : _bfc(&(*v))
  {
    if (v->face()==NULL  ||
        v->face()->dimension() == 0 || 
	v->face()->dimension() == 1) { _bfc = Base_face_circulator();}
  }

  Triangulation_face_circulator_2(const Vertex_handle& v,
				       const Face_handle& f)
    : _bfc( &(*v), &(*f))
  {
    if ( v->face().is_null()  ||
        v->face()->dimension() == 0 || 
	v->face()->dimension() == 1) { _bfc = Base_face_circulator();}
  }

 
   
   Triangulation_face_circulator_2(const Face_circulator &fc)
    : _bfc(fc._bfc)
  {}
   
  Triangulation_face_circulator_2 &operator=(const Face_circulator &fc)
  {
    _bfc = fc._bfc;
    return *this;
  }
  
  Face_circulator& operator++()
  {
    ++_bfc;
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
    --_bfc;
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
    return (_bfc.is_empty());
  }
        
  bool
  operator!=(CGAL_NULL_TYPE n) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return  (! _bfc.is_empty());
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
class Triangulation_vertex_circulator_2
  : public Bidirectional_circulator_base<Triangulation_vertex_2<Gt,Tds>,
                                    CGAL_STD::ptrdiff_t,CGAL_STD::size_t>,
    public Triangulation_cw_ccw_2
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Vertex_circulator  Base_vertex_circulator;

  typedef Triangulation_face_2<Gt,Tds> Face;
  typedef Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef std::pair<Face_handle, int>     Edge;

  typedef Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;

 
  Triangulation_vertex_circulator_2()
    : _bvc()
  {}

  Triangulation_vertex_circulator_2(const Vertex_handle& v)
    : _bvc(&(*v))
  {
    if (v->face()==NULL  ||
        v->face()->dimension() == 0 )	 { _bvc = Base_vertex_circulator();}
  }

  Triangulation_vertex_circulator_2(const Vertex_handle& v,
				       const Face_handle& f)
    : _bvc( &(*v), &(*f))
  {
    if ( v->face().is_null()  ||
        v->face()->dimension() == 0 )	 { _bvc = Base_vertex_circulator();}
  }
   
   Triangulation_vertex_circulator_2(const Vertex_circulator &vc)
    : _bvc(vc._bvc)
  {}
   
  Triangulation_vertex_circulator_2 &
	operator=(const Vertex_circulator &vc)
  {
    _bvc = vc._bvc;
    return *this;
  }
  
  Vertex_circulator& operator++()
  {
    ++_bvc;
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
    --_bvc;
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
class Triangulation_edge_circulator_2
  : public 
   Bidirectional_circulator_base< 
          typename Triangulation_2<Gt, Tds>::Edge, 
           CGAL_STD::ptrdiff_t, CGAL_STD::size_t>,
    public Triangulation_cw_ccw_2
{
public:
  typedef Tds Triangulation_data_structure;

  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Edge_circulator  Base_edge_circulator;

  typedef Triangulation_face_2<Gt,Tds> Face;
  typedef Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef std::pair<Face_handle, int>     Edge;

  typedef Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;
  
  Triangulation_edge_circulator_2()
    : _bec()
  {}

  Triangulation_edge_circulator_2(const Vertex_handle& v)
    : _bec(&(*v))
  {
    if (v->face().is_null()  ||
        v->face()->dimension() == 0 )	 { _bec = Base_face_circulator();}
  }

  Triangulation_edge_circulator_2(const Vertex_handle& v,
				  const Face_handle& f)
    : _bec( &(*v), &(*f))
  {
    if (v->face().is_null()  ||
        v->face()->dimension() == 0 ) { _bec = Base_edge_circulator();}
  }
   
   Triangulation_edge_circulator_2(const Edge_circulator &ec)
    : _bec(ec._bec)
  {}
   
  Triangulation_edge_circulator_2 &operator=(const Edge_circulator &ec)
  {
    _bec = ec._bec;
    return *this;
  }
  
  Edge_circulator& operator++()
  {
    ++_bec;
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
    --_bec;
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
    return std::make_pair( fh  , (*_bec).second );
  }

    
private:
  Base_edge_circulator  _bec;
};

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_CIRCULATORS_2_H
