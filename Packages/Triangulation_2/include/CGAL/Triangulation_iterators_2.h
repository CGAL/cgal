#ifndef CGAL_TRIANGULATION_ITERATORS_2_H
#define CGAL_TRIANGULATION_ITERATORS_2_H



#include <pair.h>
#include <iterator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_ds_iterators_2.h>

template < class Gt, class Tds >
class CGAL_Triangulation_face_2;

template < class Gt, class Tds >
class CGAL_Triangulation_vertex_2;

template < class Gt, class Tds>
class CGAL_Triangulation_face_iterator_2;

template < class Gt, class Tds>
class CGAL_Triangulation_vertex_iterator_2;

template < class Gt, class Tds>
class CGAL_Triangulation_edge_iterator_2;


template < class Gt, class Tds>
class CGAL_Triangulation_face_iterator_2
    : public bidirectional_iterator<CGAL_Triangulation_face_2<Gt,Tds>,ptrdiff_t>
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Face_iterator  Iterator_base;

  typedef CGAL_Triangulation_2<Gt,Tds> Triangulation;
  typedef typename Triangulation::Face Face;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Face_handle Face_handle;
  typedef typename Triangulation::Edge     Edge;

  typedef CGAL_Triangulation_face_iterator_2<Gt,Tds>      Face_iterator;
  typedef CGAL_Triangulation_edge_iterator_2<Gt,Tds>      Edge_iterator;
  typedef CGAL_Triangulation_vertex_iterator_2<Gt,Tds>    Vertex_iterator;


private:
  Iterator_base   _ib;
  Triangulation* _tr;
 
public:
  CGAL_Triangulation_face_iterator_2()
    : _ib(), _tr(NULL)
  {}
        
  CGAL_Triangulation_face_iterator_2(CGAL_Triangulation_2<Gt,Tds> *tr)
            : _ib( &(tr->_tds)), _tr(tr)
  { 
    if (tr->dimension() == 0 ||tr->dimension() ==1){
      _ib = Iterator_base(&(tr->_tds),1);
      return;
    }
    while ( _ib != Iterator_base(&(tr->_tds),1) && _tr->is_infinite((Face *) &(*_ib))){
      ++_ib;
    }
    return;
  }
        
  CGAL_Triangulation_face_iterator_2(CGAL_Triangulation_2<Gt,Tds> *tr, int i)
    : _ib( &(tr->_tds), i), _tr(tr)
  { }
       
  CGAL_Triangulation_face_iterator_2(const Face_iterator& fi)
          : _ib(fi._ib), _tr(fi._tr)
  {}
        
        
   Face_iterator&
        operator=(const Face_iterator& fi)
  { 
    _ib = fi._ib;
    _tr = fi._tr;
    return *this;
  }
  
  bool
  operator==(const Face_iterator& fi) const
  {
    return (  _tr == fi._tr && _ib == fi._ib);
  }

  bool
  operator!=(const Face_iterator& fi)
  {
    return !(*this == fi);
  }

  Face_iterator&
  operator++()
  {
    ++_ib;
    while ( _ib != Iterator_base(&(_tr->_tds),1) && _tr->is_infinite((Face *)&( *_ib))){
      ++_ib;
    }
    return *this;   
  }

  Face_iterator&
  operator--()
  {
     --ib;
     while ( _ib != Iterator_base(&(_tr->_tds),1) && _tr->is_infinite((Face *) &(*_ib))){
      --_ib;
    }
    return *this;   
  }

  Face_iterator
  operator++(int)
  {
    Face_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
        
  Face_iterator
  operator--(int)
  {
    Face_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  inline Face& operator*() const
  {
    return (Face &)(*_ib);
  }

  inline Face* operator->() const
  {
    return   (Face*)( & (*_ib));
  }
     
};


template < class Gt, class Tds>
class CGAL_Triangulation_vertex_iterator_2
    : public bidirectional_iterator<CGAL_Triangulation_vertex_2<Gt, Tds>,ptrdiff_t>
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Vertex_iterator  Iterator_base;

  typedef CGAL_Triangulation_2<Gt,Tds> Triangulation;
  typedef typename Triangulation::Face Face;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Face_handle Face_handle;
  typedef typename Triangulation::Edge     Edge;

  typedef CGAL_Triangulation_face_iterator_2<Gt,Tds>      Face_iterator;
  typedef CGAL_Triangulation_edge_iterator_2<Gt,Tds>      Edge_iterator;
  typedef CGAL_Triangulation_vertex_iterator_2<Gt,Tds>    Vertex_iterator;

private:
  Iterator_base   _ib;
  Triangulation* _tr; 

public:
  CGAL_Triangulation_vertex_iterator_2()
    : _ib(),_tr(NULL) 
  {}
        
  CGAL_Triangulation_vertex_iterator_2(CGAL_Triangulation_2<Gt,Tds> *tr)
            : _ib( &(tr->_tds)), _tr(tr)
  { 
    if (_tr->number_of_vertices() == 0) { _ib = Iterator_base(&(tr->_tds),1);}
    //else if ( _tr->is_infinite( (Vertex &) *_ib) ) { ++_ib;}
    else if ( _tr->is_infinite( (Vertex *) &(*_ib)) ){ ++_ib;}
    return;
  }
        
  CGAL_Triangulation_vertex_iterator_2(CGAL_Triangulation_2<Gt,Tds> *tr, int i)
    : _ib( &(tr->_tds), i), _tr(tr)
  {   }
       
  CGAL_Triangulation_vertex_iterator_2(const Vertex_iterator& vi)
          : _ib(vi._ib), _tr(vi._tr)
  {}
        
        
   Vertex_iterator&
        operator=(const Vertex_iterator& vi)
  { 
    _tr = vi._tr;
    _ib = vi._ib;
    return *this;
  }
  
  bool
  operator==(const Vertex_iterator& vi) const
  {
    return ( _tr == vi._tr && _ib == vi._ib);
  }

  bool
  operator!=(const Vertex_iterator& vi)
  {
    return !(*this == vi);
  }

  Vertex_iterator&
  operator++()
  {
     ++_ib;
     while ( _ib != Iterator_base(&(_tr->_tds),1) && _tr->is_infinite((Vertex *) &(*_ib))){
      ++_ib;
    }
    return *this;   
  }

  Vertex_iterator&
  operator--()
  {
    --_ib;
    while ( _ib != Iterator_base(&(_tr->_tds),1) && _tr->is_infinite((Vertex *) &(*_ib))){
      --_ib;
    }
    return *this;   
  }

  Vertex_iterator
  operator++(int)
  {
    Vertex_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
        
  Vertex_iterator
  operator--(int)
  {
    Vertex_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  inline Vertex& operator*() const
  {
    return (Vertex &)(*_ib);
  }

  inline Vertex* operator->() const
  {
    return   (Vertex*)( & (*_ib));
  }
     
};


template < class Gt, class Tds>
class CGAL_Triangulation_edge_iterator_2
    : public bidirectional_iterator<typename CGAL_Triangulation_2<Gt,Tds>::Edge ,ptrdiff_t>
{
public:
  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;
  typedef typename Tds::Edge_iterator  Iterator_base;

  typedef CGAL_Triangulation_2<Gt,Tds> Triangulation;
  typedef typename Triangulation::Face Face;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Face_handle Face_handle;
  typedef typename Triangulation::Edge     Edge;


  typedef CGAL_Triangulation_face_iterator_2<Gt,Tds>      Face_iterator;
  typedef CGAL_Triangulation_edge_iterator_2<Gt,Tds>      Edge_iterator;
  typedef CGAL_Triangulation_vertex_iterator_2<Gt,Tds>    Vertex_iterator;

private:
   Iterator_base   _ib;
  Triangulation* _tr;

public:  
  CGAL_Triangulation_edge_iterator_2()
    : _ib(), _tr(NULL)
  {}
        
  CGAL_Triangulation_edge_iterator_2(CGAL_Triangulation_2<Gt,Tds> *tr)
            : _ib( &(tr->_tds)), _tr(tr)
  { 
    if (_tr->dimension() == 0 ) {
      _ib = Iterator_base(&(tr->_tds),1);
      return;
    }
    while ( _ib != Iterator_base(&(tr->_tds),1) && _tr->is_infinite((Face *)((*_ib).first), (*_ib).second)){
      ++_ib;
    }
    return;
  }
        
  CGAL_Triangulation_edge_iterator_2(CGAL_Triangulation_2<Gt,Tds> *tr, int i)
    : _ib( &(tr->_tds), i),  _tr(tr)
  { }
       
  CGAL_Triangulation_edge_iterator_2(const Edge_iterator& ei)
          : _ib(ei._ib), _tr(ei._tr)
  {}
        
        
   Edge_iterator&
        operator=(const Edge_iterator& ei)
  { 
    _tr = ei._tr;
    _ib = ei._ib;
    return *this;
  }
  
  bool
  operator==(const Edge_iterator& ei) const
  {
    return ( _tr == ei._tr && _ib == ei._ib);
  }

  bool
  operator!=(const Edge_iterator& ei)
  {
    return !(*this == ei);
  }

  Edge_iterator&
  operator++()
  {
    ++_ib;
    while ( _ib != Iterator_base(&(_tr->_tds),1) && _tr->is_infinite((Face *)((*_ib).first), (*_ib).second)){
      ++_ib;
    }
     return *this;   
  }

  Edge_iterator&
  operator--()
  {
    --_ib;
    while ( _ib != Iterator_base(&(_tr->_tds),1) && _tr->is_infinite((Face *)((*_ib).first), *_ib.second)){
      --_ib;
     }
     return *this;   
  }

  Edge_iterator
  operator++(int)
  {
    Edge_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
        
  Edge_iterator
  operator--(int)
  {
    Edge_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  inline Edge  operator*() const
  {
    Face_handle fh = (Face *)((*_ib).first);
    return make_pair( fh  , (*_ib).second );
  }

};




#endif CGAL_TRIANGULATION_ITERATORS_2_H
