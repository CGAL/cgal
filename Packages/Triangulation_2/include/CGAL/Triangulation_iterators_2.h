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
// file          : Triangulation/include/CGAL/Triangulation_iterators_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_ITERATORS_2_H
#define CGAL_TRIANGULATION_ITERATORS_2_H



#include <utility>
#include <iterator>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_ds_iterators_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
class Triangulation_all_faces_iterator_2 
  : public Tds::Face_iterator
{
public:
  typedef Triangulation_face_2<Gt,Tds>       value_type;
  typedef Triangulation_face_2<Gt,Tds> *     pointer;
  typedef Triangulation_face_2<Gt,Tds> &     reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

  typedef typename Tds::Face_iterator  Base;
  typedef Triangulation_2<Gt,Tds> Triangulation;
  typedef typename Triangulation::Face Face;
  typedef Triangulation_all_faces_iterator_2<Gt,Tds> All_faces_iterator;

public:
  Triangulation_all_faces_iterator_2()
    : Base()
    {}
        
  Triangulation_all_faces_iterator_2(const Triangulation_2<Gt,Tds> *tr)
    : Base(&(tr->_tds))
    {}
 
  Triangulation_all_faces_iterator_2(const Triangulation_2<Gt,Tds> *tr, int i)
    : Base(&(tr->_tds),1) 
    {}
       

  All_faces_iterator&  operator++()
  {
    Base::operator++();
    return *this;   
  }

  All_faces_iterator&  operator--()
  {
     Base::operator--();
     return *this;   
  }

  All_faces_iterator  operator++(int)
  {
    All_faces_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  All_faces_iterator  operator--(int)
  {
    All_faces_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  Face& operator*() const
  {
    return static_cast<Face &>(Base::operator*());
  }

  Face* operator->() const
  {
    return static_cast<Face *>(Base::operator->());
  }
     
};


template < class Gt, class Tds>
class Triangulation_finite_faces_iterator_2
 : public Triangulation_all_faces_iterator_2<Gt,Tds>
{
public:
  typedef Triangulation_2<Gt,Tds>  Triangulation;
  typedef Triangulation_all_faces_iterator_2<Gt,Tds>     All;
  typedef Triangulation_finite_faces_iterator_2<Gt,Tds>  Finite_faces_iterator;
  
private:
   const Triangulation*  _tr;
 
public:
  Triangulation_finite_faces_iterator_2()
    : All(), _tr()
  {}
        
  Triangulation_finite_faces_iterator_2(const Triangulation* tr)
    : All(tr), _tr(tr)
  { 
    while ( *this != All(tr,1) && ( _tr->is_infinite(& **this))) 
      All::operator++();
    return;
  }

  Triangulation_finite_faces_iterator_2(const Triangulation* tr, int i)
    : All(tr,i), _tr(tr)
  { }
       

  Finite_faces_iterator&  operator++()
  {
    do All::operator++();
    while ( *this != All(_tr,1)  && _tr->is_infinite(*this));
    return *this;   
  }

  Finite_faces_iterator&  operator--()
  {
    do All::operator--();
    while ( *this !=  All(_tr,1) &&  _tr->is_infinite(*this));
    return *this;   
  }

  Finite_faces_iterator  operator++(int)
  {
    Finite_faces_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
        
  Finite_faces_iterator  operator--(int)
  {
    Finite_faces_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
};


template < class Gt, class Tds>
class Triangulation_all_vertices_iterator_2
  : public Tds::Vertex_iterator
{
public:
  typedef Triangulation_vertex_2<Gt,Tds>       value_type;
  typedef Triangulation_vertex_2<Gt,Tds> *     pointer;
  typedef Triangulation_vertex_2<Gt,Tds> &     reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

  typedef typename Tds::Vertex_iterator  Base;
  typedef Triangulation_2<Gt,Tds> Triangulation;
  typedef typename Triangulation::Vertex Vertex;
  typedef Triangulation_all_vertices_iterator_2<Gt,Tds> All_vertices_iterator;

public:
  Triangulation_all_vertices_iterator_2()
    : Base() 
  {}
        
  Triangulation_all_vertices_iterator_2(const Triangulation * tr)
    : Base( &(tr->_tds))
  { }

  Triangulation_all_vertices_iterator_2(const Triangulation *tr, int i)
    : Base( &(tr->_tds),i)
  { }

  
  All_vertices_iterator&   operator++()
  {
    Base::operator++();
    return *this;   
  }

  All_vertices_iterator&   operator--()
  {
    Base::operator--();
    return *this;   
  }

  All_vertices_iterator   operator++(int)
  {
    All_vertices_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
         
  All_vertices_iterator   operator--(int)
  {
    All_vertices_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  Vertex& operator*() const
  {
    return static_cast<Vertex&>(Base::operator*());
  }

  Vertex* operator->() const
  {
    return static_cast<Vertex*>(Base::operator->());   
  }
     
};


template < class Gt, class Tds>
class Triangulation_finite_vertices_iterator_2
 : public Triangulation_all_vertices_iterator_2<Gt,Tds>
{
public:
  typedef Triangulation_all_vertices_iterator_2<Gt,Tds> All;
  typedef Triangulation_finite_vertices_iterator_2<Gt,Tds> 
                                            Finite_vertices_iterator;
private:
  const Triangulation* _tr; 

public:
  Triangulation_finite_vertices_iterator_2()
    : All(),_tr(NULL) 
  {}
        
  Triangulation_finite_vertices_iterator_2(const Triangulation *tr)
            : All(tr), _tr(tr)
  { 
    while ( *this != All(_tr,1) && _tr->is_infinite(*this)) 
      All::operator++();
    return;
  }

  Triangulation_finite_vertices_iterator_2(const Triangulation *tr, int i)
    : All(tr,i), _tr(tr)
  { }


  Finite_vertices_iterator&   operator++()
  {
    do All::operator++();
    while ( *this != All(_tr,1)  && _tr->is_infinite(*this));
    return *this;  
  }

  Finite_vertices_iterator&  operator--()
  {
    do All::operator--();
    while ( *this != All(_tr,1)  && _tr->is_infinite(*this));
    return *this;   
  }

  Finite_vertices_iterator
  operator++(int)
  {
    Finite_vertices_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
        
  Finite_vertices_iterator
  operator--(int)
  {
    Finite_vertices_iterator tmp(*this);
    --(*this);
    return tmp;
  }
};

template < class Gt, class Tds>
class Triangulation_all_edges_iterator_2
  : public Tds::Edge_iterator
{
public:
  typedef Triangulation_2<Gt,Tds>     Triangulation;
  typedef typename Triangulation::Edge         Edge;
  typedef typename Triangulation::Face        Face;
  typedef typename Triangulation::Face_handle  Face_handle;
  typedef typename Tds::Edge_iterator Base;
  typedef Triangulation_all_edges_iterator_2<Gt,Tds> All_edges_iterator;

  typedef Edge       value_type;
  typedef Edge *     pointer;
  typedef Edge &     reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;
  


  Triangulation_all_edges_iterator_2()
    : Base()
    {}
        
  Triangulation_all_edges_iterator_2(const Triangulation *tr)
    : Base(&(tr->_tds))
    {} 
  
  Triangulation_all_edges_iterator_2(const Triangulation *tr, int i)
    : Base(&(tr->_tds),i)
  { }
       
  All_edges_iterator&  operator++()
  {
    Base::operator++();
    return *this;
  }

  All_edges_iterator&
  operator--()
  {
    Base::operator--();
    return *this;   
  }

  All_edges_iterator
  operator++(int)
  {
    All_edges_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
        
  All_edges_iterator
  operator--(int)
  {
    All_edges_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  Edge  operator*() const
  {
    Face_handle fh = static_cast<Face *>(Base::operator*().first);
    return std::make_pair( fh  , Base::operator*().second );
  }
};

template < class Gt, class Tds>
class Triangulation_finite_edges_iterator_2
 : public Triangulation_all_edges_iterator_2<Gt,Tds>
{
public:
  
  typedef Triangulation_all_edges_iterator_2<Gt,Tds> All;
  typedef Triangulation_finite_edges_iterator_2<Gt,Tds> Finite_edges_iterator;
  typedef Triangulation_finite_edges_iterator_2<Gt,Tds> Finite;
private:
  const Triangulation* _tr;

public:  
  Triangulation_finite_edges_iterator_2()
    : All(), _tr(NULL)
  {}
        
  Triangulation_finite_edges_iterator_2(const Triangulation *tr)
    : All(tr),_tr(tr) 
  { 
    while ( *this != All(tr,1) &&   _tr->is_infinite(*this) ) {
      All::operator++();
    }
    return;
  }

  Triangulation_finite_edges_iterator_2(const Triangulation *tr, int i)
    : All(tr,i), _tr(tr)
  { }
       
  Finite_edges_iterator&   operator++()
  {
    do {
      All::operator++();
    }  while ( *this != All(_tr,1) &&  _tr->is_infinite(*this));
    return *this;   
  }

  Finite_edges_iterator&
  operator--()
  {
    do All::operator--();
    while ( *this != All(_tr,1) &&  _tr->is_infinite(*this));
    return *this;  
  }

  Finite_edges_iterator
  operator++(int)
  {
    Finite_edges_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Finite_edges_iterator
  operator--(int)
  {
    Finite_edges_iterator tmp(*this);
    --(*this);
    return tmp;
  }

};

CGAL_END_NAMESPACE


#endif //CGAL_TRIANGULATION_ITERATORS_2_H
