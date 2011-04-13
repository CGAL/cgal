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
#include <CGAL/Triangulation_face_2.h>
#include <CGAL/Triangulation_vertex_2.h>
#include <CGAL/Triangulation_handles_2.h>
#include <CGAL/Triangulation_ds_iterators_2.h>

CGAL_BEGIN_NAMESPACE

template <class Gt, class Tds> class Triangulation_2;

template < class Gt, class Tds>
class Triangulation_all_faces_iterator_2 
  : public Tds::Face_iterator
{
public:
  typedef typename Tds::Face_iterator  Base;
  typedef Triangulation_2<Gt,Tds> Triangulation;
  typedef typename Triangulation::Face Face;
  typedef Triangulation_all_faces_iterator_2<Gt,Tds>  All_faces;

  typedef Face      value_type;
  typedef Face*     pointer;
  typedef Face&     reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

public:
  Triangulation_all_faces_iterator_2()     : Base()    {}
  Triangulation_all_faces_iterator_2(const Triangulation_2<Gt,Tds> *tr)
    : Base(&(tr->_tds))
    {}
   Triangulation_all_faces_iterator_2(const Triangulation_2<Gt,Tds> *tr, int)
    : Base(&(tr->_tds),1) 
    {}
  
  All_faces&  operator++();
  All_faces&  operator--();
  All_faces  operator++(int);
  All_faces  operator--(int);
  Face& operator*() const
    {return static_cast<Face &>(All_faces::Base::operator*());}
  Face* operator->() const
    {return static_cast<Face *>(All_faces::Base::operator->());}
};


template < class Gt, class Tds>
class Triangulation_finite_faces_iterator_2
 : public Triangulation_all_faces_iterator_2<Gt,Tds>
{
public:
  typedef Triangulation_all_faces_iterator_2<Gt,Tds>     All_faces;
  typedef Triangulation_finite_faces_iterator_2<Gt,Tds>  Finite_faces;
  typedef typename All_faces::Triangulation              Triangulation;

private:
   const Triangulation*  _tr;
 
public:
  Triangulation_finite_faces_iterator_2()    : All_faces(), _tr()  {}
        
  Triangulation_finite_faces_iterator_2(const Triangulation* tr);

  Triangulation_finite_faces_iterator_2(const Triangulation* tr, int)
    : All_faces(tr,1), _tr(tr)
  { }
  
  Finite_faces&  operator++();
  Finite_faces&  operator--();
  Finite_faces   operator++(int);
  Finite_faces   operator--(int);
};


template < class Gt, class Tds>
class Triangulation_all_vertices_iterator_2
  : public Tds::Vertex_iterator
{
public:
  typedef typename Tds::Vertex_iterator  Base;
  typedef Triangulation_2<Gt,Tds> Triangulation;
  typedef typename Triangulation::Vertex Vertex;
  typedef Triangulation_all_vertices_iterator_2<Gt,Tds> All_vertices;

  typedef Vertex       value_type;
  typedef const Vertex *     pointer;
  typedef const Vertex &     reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;


public:
  Triangulation_all_vertices_iterator_2()
    : Base() 
  {}
        
  Triangulation_all_vertices_iterator_2(const Triangulation * tr)
    : Base( &(tr->_tds))
  { }

  Triangulation_all_vertices_iterator_2(const Triangulation *tr, int)
    : Base( &(tr->_tds),1)
  { }
  
  All_vertices&   operator++();
  All_vertices&   operator--();
  All_vertices    operator++(int);
  All_vertices    operator--(int);
  Vertex& operator*() const;
  Vertex* operator->() const;
};


template < class Gt, class Tds>
class Triangulation_finite_vertices_iterator_2
 : public Triangulation_all_vertices_iterator_2<Gt,Tds>
{
public:
  typedef Triangulation_all_vertices_iterator_2<Gt,Tds>    All_vertices;
  typedef Triangulation_finite_vertices_iterator_2<Gt,Tds> Finite_vertices;
  // typedef Triangulation_2<Gt,Tds>                          Triangulation;
  typedef typename All_vertices::Triangulation             Triangulation;

private:
  const Triangulation* _tr; 

public:
  Triangulation_finite_vertices_iterator_2()    : All_vertices(),_tr(NULL)   {}
          
  Triangulation_finite_vertices_iterator_2(const Triangulation *tr);
    
  Triangulation_finite_vertices_iterator_2(const Triangulation *tr, int )
    : All_vertices(tr,1), _tr(tr)
  { }

  Finite_vertices&  operator++();
  Finite_vertices&  operator--();
  Finite_vertices   operator++(int);
  Finite_vertices   operator--(int);
};

template < class Gt, class Tds>
class Triangulation_all_edges_iterator_2
  : public Tds::Edge_iterator
{
public:
  typedef Triangulation_2<Gt,Tds>     Triangulation;
  typedef typename Triangulation::Edge  Edge;
  typedef typename Triangulation::Face        Face;
  typedef typename Triangulation::Face_handle  Face_handle;
  typedef typename Tds::Edge_iterator Base;
  typedef Triangulation_all_edges_iterator_2<Gt,Tds> All_edges;

  typedef Edge       value_type;
  typedef Edge *     pointer;
  typedef Edge &     reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;
  

  Triangulation_all_edges_iterator_2()    : Base()    {}
        
  Triangulation_all_edges_iterator_2(const Triangulation *tr)
    : Base(&(tr->_tds))
    {} 
  
  Triangulation_all_edges_iterator_2(const Triangulation *tr, int )
    : Base(&(tr->_tds),1)
  { }
       
  All_edges&  operator++();
  All_edges&  operator--();
  All_edges   operator++(int);
  All_edges   operator--(int);
  Edge  operator*() const;
};

template < class Gt, class Tds>
class Triangulation_finite_edges_iterator_2
 : public Triangulation_all_edges_iterator_2<Gt,Tds>
{
public:
  typedef Triangulation_all_edges_iterator_2<Gt,Tds>    All_edges;
  typedef Triangulation_finite_edges_iterator_2<Gt,Tds> Finite_edges;
  typedef typename All_edges::Triangulation             Triangulation;
  
private:
  const Triangulation* _tr;

public:  
  Triangulation_finite_edges_iterator_2()  : All_edges(), _tr(NULL)  {}
        
  Triangulation_finite_edges_iterator_2(const Triangulation *tr);
  
  Triangulation_finite_edges_iterator_2(const Triangulation *tr, int )
    : All_edges(tr,1), _tr(tr)
  { }
       
  Finite_edges& operator++();
  Finite_edges& operator--();
  Finite_edges  operator++(int);
  Finite_edges  operator--(int);
};

  //Helping classes
 //  to be used as adaptators from iterators value_type
 //   Edge of the triangulation
 //  to an iterator with Tds::Edge as value type
  template<class It, class TdsEdge>
  class To_tds_edge_iterator : public It {
  public:
    typedef TdsEdge  Tds_Edge;
    typedef typename std::iterator_traits<It>::value_type Edge;
    To_tds_edge_iterator() {}
    To_tds_edge_iterator(It i) : It(i) {} 
    Tds_Edge  operator*() {
      Edge e = It::operator*();
      return Tds_Edge( &*(e.first), e.second);
    }
  };

 //  to be used as adaptators from iterators with Face_handle value_type
 //   to an iterator with Tds::Face* as value type
  template<class It, class TdsFace>
  class To_tds_face_iterator : public It {
  public:
    typedef TdsFace  Tds_Face;
    To_tds_face_iterator() {}
    To_tds_face_iterator(It i) : It(i) {} 
    Tds_Face* operator*() { return  &*(It::operator*() ); }
  };









template < class Gt, class Tds>
inline
Triangulation_all_faces_iterator_2<Gt,Tds>&
Triangulation_all_faces_iterator_2<Gt,Tds>::
operator++()
{
  Base::operator++();
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_all_faces_iterator_2<Gt,Tds>&
Triangulation_all_faces_iterator_2<Gt,Tds>::
operator--()
{
  Base::operator--();
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_all_faces_iterator_2<Gt,Tds>
Triangulation_all_faces_iterator_2<Gt,Tds>::
operator++(int)
{
  All_faces tmp(*this);
  ++(*this);
  return tmp;
}
        
template < class Gt, class Tds>
inline
Triangulation_all_faces_iterator_2<Gt,Tds>
Triangulation_all_faces_iterator_2<Gt,Tds>::
operator--(int)
{
  All_faces tmp(*this);
  --(*this);
  return tmp;
}
        
// template < class Gt, class Tds>
// inline
// typename Triangulation_2<Gt,Tds>::Face &
// Triangulation_all_faces_iterator_2<Gt,Tds>::
// operator*() const
// {
//   return static_cast<Face &>(Base::operator*());
// }

// template < class Gt, class Tds>
// inline
// typename Triangulation_2<Gt,Tds>::Face *
// Triangulation_all_faces_iterator_2<Gt,Tds>::
// operator->() const
// {
//   return static_cast<Face *>(Base::operator->());
// }
     
template < class Gt, class Tds>
inline
Triangulation_finite_faces_iterator_2<Gt,Tds>::
Triangulation_finite_faces_iterator_2(const Triangulation* tr)
    : All_faces(tr), _tr(tr)
{ 
  while ( *this != All_faces(tr,1) && ( _tr->is_infinite(& **this))) 
    All_faces::operator++();
  return;
}

template < class Gt, class Tds>
inline
Triangulation_finite_faces_iterator_2<Gt,Tds>&
Triangulation_finite_faces_iterator_2<Gt,Tds>::
operator++()
{
  do { All_faces::operator++();}
  while ( *this != All_faces(_tr,1)  && _tr->is_infinite(*this));
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_finite_faces_iterator_2<Gt,Tds>&
Triangulation_finite_faces_iterator_2<Gt,Tds>::
operator--()
{
  do {All_faces::operator--();}
  while ( *this !=  All_faces(_tr,1) &&  _tr->is_infinite(*this));
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_finite_faces_iterator_2<Gt,Tds>
Triangulation_finite_faces_iterator_2<Gt,Tds>::
operator++(int)
{
  Finite_faces tmp(*this);
  ++(*this);
  return tmp;
}
        
template < class Gt, class Tds>
inline
Triangulation_finite_faces_iterator_2<Gt,Tds>
Triangulation_finite_faces_iterator_2<Gt,Tds>::        
operator--(int)
{
  Finite_faces tmp(*this);
  --(*this);
  return tmp;
}

template < class Gt, class Tds>
inline
Triangulation_all_vertices_iterator_2<Gt,Tds>&
Triangulation_all_vertices_iterator_2<Gt,Tds>::
operator++()
{
  Base::operator++();
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_all_vertices_iterator_2<Gt,Tds>&
Triangulation_all_vertices_iterator_2<Gt,Tds>::
operator--()
{
  Base::operator--();
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_all_vertices_iterator_2<Gt,Tds>
Triangulation_all_vertices_iterator_2<Gt,Tds>::
operator++(int)
{
  All_vertices tmp(*this);
  ++(*this);
  return tmp;
}
         
template < class Gt, class Tds>
inline
Triangulation_all_vertices_iterator_2<Gt,Tds>
Triangulation_all_vertices_iterator_2<Gt,Tds>::
operator--(int)
{
  All_vertices tmp(*this);
  --(*this);
  return tmp;
}
        
template < class Gt, class Tds>
inline
Triangulation_all_vertices_iterator_2<Gt,Tds>::Vertex &
Triangulation_all_vertices_iterator_2<Gt,Tds>::
operator*() const
{
  return static_cast<Vertex&>(Base::operator*());
}

template < class Gt, class Tds>
inline
Triangulation_all_vertices_iterator_2<Gt,Tds>::Vertex*
Triangulation_all_vertices_iterator_2<Gt,Tds>::
operator->() const
{
  return static_cast<Vertex*>(Base::operator->());   
}

template < class Gt, class Tds>
inline
Triangulation_finite_vertices_iterator_2<Gt,Tds> ::
Triangulation_finite_vertices_iterator_2(const Triangulation *tr)
  : All_vertices(tr), _tr(tr)
{ 
  while ( *this != All_vertices(_tr,1) && _tr->is_infinite(*this)) 
    All_vertices::operator++();
  return;
} 

template < class Gt, class Tds>
inline
Triangulation_finite_vertices_iterator_2<Gt,Tds>&
Triangulation_finite_vertices_iterator_2<Gt,Tds> ::
operator++()
{
  do All_vertices::operator++();
  while ( *this != All_vertices(_tr,1)  && _tr->is_infinite(*this));
  return *this;  
}

template < class Gt, class Tds>
inline
Triangulation_finite_vertices_iterator_2<Gt,Tds>&
Triangulation_finite_vertices_iterator_2<Gt,Tds> ::
operator--()
{
  do All_vertices::operator--();
  while ( *this != All_vertices(_tr,1)  && _tr->is_infinite(*this));
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_finite_vertices_iterator_2<Gt,Tds>
Triangulation_finite_vertices_iterator_2<Gt,Tds> ::
operator++(int)
{
  Finite_vertices tmp(*this);
  ++(*this);
  return tmp;
}
        
template < class Gt, class Tds>
inline
Triangulation_finite_vertices_iterator_2<Gt,Tds>
Triangulation_finite_vertices_iterator_2<Gt,Tds> ::        
operator--(int)
{
  Finite_vertices tmp(*this);
  --(*this);
  return tmp;
}

template < class Gt, class Tds>
inline
Triangulation_all_edges_iterator_2<Gt,Tds>&
Triangulation_all_edges_iterator_2<Gt,Tds>::
operator++()
{
  Base::operator++();
  return *this;
}

template < class Gt, class Tds>
inline
Triangulation_all_edges_iterator_2<Gt,Tds>&
Triangulation_all_edges_iterator_2<Gt,Tds>::
operator--()
{
  Base::operator--();
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_all_edges_iterator_2<Gt,Tds>
Triangulation_all_edges_iterator_2<Gt,Tds>::
operator++(int)
{
  All_edges tmp(*this);
  ++(*this);
  return tmp;
}
        
template < class Gt, class Tds>
inline
Triangulation_all_edges_iterator_2<Gt,Tds>
Triangulation_all_edges_iterator_2<Gt,Tds>::
operator--(int)
{
  All_edges tmp(*this);
  --(*this);
  return tmp;
}
        
template < class Gt, class Tds>
inline
Triangulation_all_edges_iterator_2<Gt,Tds>::Edge
Triangulation_all_edges_iterator_2<Gt,Tds>::
operator*() const
{
  Face_handle fh = static_cast<Face *>(Base::operator*().first);
  return std::make_pair( fh  , Base::operator*().second );
}

template < class Gt, class Tds>
inline
Triangulation_finite_edges_iterator_2<Gt,Tds>::
Triangulation_finite_edges_iterator_2(const Triangulation *tr)
  : All_edges(tr),_tr(tr) 
{ 
  while ( *this != All_edges(tr,1) &&   _tr->is_infinite(*this) ) {
    All_edges::operator++();
  }
  return;
}

template < class Gt, class Tds>
inline
Triangulation_finite_edges_iterator_2<Gt,Tds>&
Triangulation_finite_edges_iterator_2<Gt,Tds>::
operator++()
{
  do {
    All_edges::operator++();
  }  while ( *this != All_edges(_tr,1) &&  _tr->is_infinite(*this));
  return *this;   
}

template < class Gt, class Tds>
inline
Triangulation_finite_edges_iterator_2<Gt,Tds>&
Triangulation_finite_edges_iterator_2<Gt,Tds>::
operator--()
{
  do {
    All_edges::operator--();
  } while ( *this != All_edges(_tr,1) &&  _tr->is_infinite(*this));
  return *this;  
}

template < class Gt, class Tds>
inline
Triangulation_finite_edges_iterator_2<Gt,Tds>
Triangulation_finite_edges_iterator_2<Gt,Tds>:: 
operator++(int)
{
  Finite_edges tmp(*this);
  ++(*this);
  return tmp;
}
        
template < class Gt, class Tds>
inline
Triangulation_finite_edges_iterator_2<Gt,Tds>
Triangulation_finite_edges_iterator_2<Gt,Tds>:: 
operator--(int)
{
  Finite_edges tmp(*this);
  --(*this);
  return tmp;
}

CGAL_END_NAMESPACE


#endif //CGAL_TRIANGULATION_ITERATORS_2_H
