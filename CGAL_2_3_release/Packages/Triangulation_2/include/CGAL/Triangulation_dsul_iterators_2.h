// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Triangulation_dsul_iterators_2.h
// package       : Triangulation 
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_TRIANGULATION_DSUL_ITERATORS_2_H
#define CGAL_TRIANGULATION_DSUL_ITERATORS_2_H



#include <iterator>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>


CGAL_BEGIN_NAMESPACE

// The  iterator_base visits all the full dimensional faces 
// of the Tds
// whatever may be the dimension of the Tds

template <class Tds>
class Triangulation_dsul_iterator_base_2
{
public:
  typedef typename Tds::Face       value_type;
  typedef typename Tds::Face *     pointer;
  typedef typename Tds::Face &     reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

  typedef Triangulation_dsul_iterator_base_2<Tds> Iterator_base;
  typedef typename Tds::Face       Face;

  Triangulation_dsul_iterator_base_2()
     : _tds(NULL), pos(NULL)
    {}

  Triangulation_dsul_iterator_base_2(const Tds* tds)
     : _tds(tds), pos(tds->dummy().next())
    {}

  Triangulation_dsul_iterator_base_2(const Tds* tds, int)
     :  _tds(tds), pos(tds->past_end())
    {}

  Triangulation_dsul_iterator_base_2(const Tds* tds, const Face* f)
    : _tds(tds), pos(f)
    {}
   
protected:
  const Tds*  _tds;
  const Face* pos;

  void increment() {pos=pos->next();}
  void decrement() {pos=pos->previous();}
  
public:
  Iterator_base&  operator=(const Iterator_base& fi);
  bool   operator==(const Iterator_base& fi) const;
  bool   operator!=(const Iterator_base& fi) const;
  
  Iterator_base&  operator++();
  Iterator_base&  operator--();
  Iterator_base   operator++(int);
  Iterator_base   operator--(int);

  Face& operator*() const
    {
      CGAL_triangulation_precondition(pos != NULL && pos != _tds->past_end());
      return const_cast<Face&>(*pos);
    }

  Face*  operator->() const
    {
      CGAL_triangulation_precondition(pos != NULL);
      return const_cast<Face*>(pos);
    }

};


// The following iterator visit the 2-faces

template<class Tds>
class Triangulation_dsul_face_iterator_2
  : public Triangulation_dsul_iterator_base_2<Tds>
{
public:
  typedef Triangulation_dsul_iterator_base_2<Tds> Iterator_base;
  typedef Triangulation_dsul_face_iterator_2<Tds> Face_iterator;

  Triangulation_dsul_face_iterator_2() : Iterator_base() {}
   
  Triangulation_dsul_face_iterator_2(const Tds * tds) : Iterator_base(tds)
    {
      if (tds->dimension() <2) pos= tds->past_end();
    }

  Triangulation_dsul_face_iterator_2(const Tds* tds, int i)
    : Iterator_base(tds,i)
    {}

  Triangulation_dsul_face_iterator_2(const Face_iterator& fi)
    : Iterator_base(fi._tds, fi.pos)
    {}
        
  Face_iterator&   operator=(const Face_iterator& fi);

  Face_iterator&  operator++();
  Face_iterator&  operator--();
  Face_iterator   operator++(int);
  Face_iterator   operator--(int);

};


template < class Tds>
class Triangulation_dsul_vertex_iterator_2
 : public Triangulation_dsul_iterator_base_2<Tds>
{
public:
  typedef typename Tds::Vertex       Vertex;
  
  typedef Vertex       value_type;
  typedef Vertex *     pointer;
  typedef Vertex &     reference;

  typedef Triangulation_dsul_iterator_base_2<Tds> Iterator_base;
  typedef Triangulation_dsul_vertex_iterator_2<Tds> Vertex_iterator;

private :
  int index;

public:
  Triangulation_dsul_vertex_iterator_2()
    : Iterator_base(), index(0)
    {}
    
  Triangulation_dsul_vertex_iterator_2(const Tds* tds);
   
  Triangulation_dsul_vertex_iterator_2(const Tds* tds, int i)
    : Iterator_base(tds,i), index(0)
    {}

  Vertex_iterator&  operator++();
  Vertex_iterator&  operator--();
  Vertex_iterator operator++(int);
  Vertex_iterator operator--(int);
  bool operator==(const Vertex_iterator& fi) const;
  bool operator!=(const Vertex_iterator& fi) const;

  Vertex& operator*() const  
    {
      CGAL_triangulation_assertion( pos != NULL && pos != _tds->past_end());
      return *(pos->vertex(index));
    }

  Vertex*  operator->() const
    {
      CGAL_triangulation_assertion( pos != NULL && pos != _tds->past_end());
      return pos->vertex(index);
    }

private:
  void increment();
  void decrement() ;
  bool associated_vertex();

 };


template <class Tds>
class Triangulation_dsul_edge_iterator_2
 : public Triangulation_dsul_iterator_base_2<Tds>
{
public:
  typedef typename Tds::Edge       Edge;
  
  typedef Edge            value_type;
  typedef Edge *          pointer;
  typedef Edge &          reference;

  typedef Triangulation_dsul_iterator_base_2<Tds> Iterator_base;
  typedef Triangulation_dsul_edge_iterator_2<Tds> Edge_iterator;
  typedef typename Triangulation_dsul_iterator_base_2<Tds>::Face Face;

private:
int index;

public:
  Triangulation_dsul_edge_iterator_2()
            : Iterator_base(),index(0)
        {}
 
  Triangulation_dsul_edge_iterator_2(const Tds * tds);
   
     
  Triangulation_dsul_edge_iterator_2(const Tds* tds, int i)
            : Iterator_base(tds,i),index(0)
    {
       if (_tds->dimension() == 1) {index = 2;}
    }

  bool  operator==(const Edge_iterator& fi) const ;
  bool  operator!=(const Edge_iterator& fi) const {return !(*this== fi);}
  Edge_iterator&    operator++();
  Edge_iterator&    operator--();
  Edge_iterator    operator++(int);
  Edge_iterator    operator--(int);

  Edge      operator*() const
    {return std::make_pair(const_cast<Face*>(pos), index);}

private: 
  void increment();
  void decrement();
  bool associated_edge();
};


// Iterator base implementation





template<class Tds>
inline
Triangulation_dsul_iterator_base_2<Tds>&
Triangulation_dsul_iterator_base_2<Tds> ::
operator=(const Iterator_base& fi)
{
  pos = fi.pos;
  _tds = fi._tds;
  return *this;
}

template<class Tds>
inline bool
Triangulation_dsul_iterator_base_2<Tds> ::
operator==(const Iterator_base& fi) const
{
  return ((pos == fi.pos )&&(_tds==fi._tds));
}
        
template<class Tds>
inline bool
Triangulation_dsul_iterator_base_2<Tds> ::        
operator!=(const Iterator_base& fi) const
{
  return !(*this == fi);
}

template<class Tds>
inline
Triangulation_dsul_iterator_base_2<Tds>&
Triangulation_dsul_iterator_base_2<Tds> ::
operator++()
{
  CGAL_triangulation_precondition(pos != NULL && pos != _tds->past_end());
  increment();
  return *this;           
}

template<class Tds>
inline
Triangulation_dsul_iterator_base_2<Tds>&
Triangulation_dsul_iterator_base_2<Tds> ::
operator--()
{
  CGAL_triangulation_precondition( pos != NULL && 
				   pos != _tds->dummy().next()); 
  decrement();
  return *this;
}
        
template<class Tds>
inline
Triangulation_dsul_iterator_base_2<Tds>
Triangulation_dsul_iterator_base_2<Tds> ::
operator++(int)
{
  Iterator_base tmp(*this);
  ++(*this);
  return tmp;
}
        
template<class Tds>
inline
Triangulation_dsul_iterator_base_2<Tds>
Triangulation_dsul_iterator_base_2<Tds> ::        
operator--(int)
{
  Iterator_base tmp(*this);
  --(*this);
  return tmp;
}
        

// Face iterator implementation
template<class Tds>
inline 
Triangulation_dsul_face_iterator_2<Tds>&
Triangulation_dsul_face_iterator_2<Tds>::
operator=(const Face_iterator& fi)
{
  Iterator_base::operator=(fi);
  return *this;
}

template<class Tds>
inline 
Triangulation_dsul_face_iterator_2<Tds>&
Triangulation_dsul_face_iterator_2<Tds>::
operator++()
{
  Iterator_base::operator++();    
  return *this;           
}

template<class Tds>
inline 
Triangulation_dsul_face_iterator_2<Tds>&
Triangulation_dsul_face_iterator_2<Tds>::
operator--()
{
  Iterator_base::operator--();    
  return *this;           
}

template<class Tds>
inline 
Triangulation_dsul_face_iterator_2<Tds>
Triangulation_dsul_face_iterator_2<Tds>::
operator++(int)
{
  Face_iterator tmp(*this);
  ++(*this);
  return tmp;
}
        
template<class Tds>
inline 
Triangulation_dsul_face_iterator_2<Tds>
Triangulation_dsul_face_iterator_2<Tds>::
operator--(int)
{
  Face_iterator tmp(*this);
  --(*this);
  return tmp;
}


// Vertex iterator implementation

template<class Tds>
Triangulation_dsul_vertex_iterator_2<Tds> ::
Triangulation_dsul_vertex_iterator_2(const Tds * tds)
  : Iterator_base(tds), index(0)
{
  while ( pos!= tds->past_end() && !associated_vertex()){
      increment();
  }
  return;
}

template<class Tds>
inline
void  
Triangulation_dsul_vertex_iterator_2<Tds> ::
increment()
{
  if (_tds->dimension()== -1) {// there is only one vertex
    CGAL_triangulation_assertion(pos == _tds->past_end()->next() && 
				 index==0);
    pos = _tds->past_end();
    return;
  }

  if ( index==_tds->dimension()) {Iterator_base::increment(); index = 0;}
  else { index +=1; }
  return;
}

template<class Tds>
inline
void   
Triangulation_dsul_vertex_iterator_2<Tds> ::
decrement()
{
  if (_tds->dimension()== -1) {// there is only one vertex
     CGAL_triangulation_assertion(pos == _tds->past_end() && index==0);
     pos = pos->next();
     return;
  }
       
  if (index == 0) { 
    Iterator_base::decrement();
    index = _tds->dimension();
  }
  else {index -= 1;}
  return;
}

template<class Tds>
inline bool 
Triangulation_dsul_vertex_iterator_2<Tds> ::
associated_vertex()
{
  return ( pos->vertex(index)->face() == pos ); // marche en toute dimension
}
   
template<class Tds>
inline
Triangulation_dsul_vertex_iterator_2<Tds>&
Triangulation_dsul_vertex_iterator_2<Tds> ::
operator++()
{
   CGAL_triangulation_precondition(pos != NULL && pos != _tds->past_end());
   do   increment();
   while (pos != _tds->past_end() && !associated_vertex());
   return *this;
}
    
template<class Tds>
inline
Triangulation_dsul_vertex_iterator_2<Tds>&
Triangulation_dsul_vertex_iterator_2<Tds> ::    
operator--()
{
   CGAL_triangulation_assertion(pos != NULL &&  
				*this != Vertex_iterator(_tds));
   do  decrement();
   while ( pos != _tds->past_end() && !associated_vertex()); 
   return *this;
}
    
template<class Tds>
inline
Triangulation_dsul_vertex_iterator_2<Tds>
Triangulation_dsul_vertex_iterator_2<Tds> ::      
operator++(int)
{
  Vertex_iterator tmp(*this);
  ++(*this);
  return tmp;
}

template<class Tds>
inline
Triangulation_dsul_vertex_iterator_2<Tds>
Triangulation_dsul_vertex_iterator_2<Tds> ::      
operator--(int)
{
  Vertex_iterator tmp(*this);
  --(*this);
  return tmp;
}
    
template<class Tds>
inline
bool
Triangulation_dsul_vertex_iterator_2<Tds> :: 
operator==(const Vertex_iterator& fi) const
{
  if (pos==NULL || pos==_tds->past_end())
    return (pos==fi.pos) && (_tds==fi._tds) ;
  return (pos==fi.pos) && (_tds==fi._tds) &&  (index == fi.index);
}
 
template<class Tds>
inline bool
Triangulation_dsul_vertex_iterator_2<Tds> ::     
operator!=(const Vertex_iterator& fi) const
{
  return !(*this == fi);
}
    

// Edge iterator implementation

template<class Tds>
Triangulation_dsul_edge_iterator_2<Tds> ::
Triangulation_dsul_edge_iterator_2(const Tds * tds)
 :  Iterator_base(tds), index(0) 
{
  if (_tds->dimension()<= 0) {
    pos = _tds->past_end();                  // there is no edge
    return;
  }
  if (_tds->dimension() == 1) index = 2;
  while ( pos != _tds->past_end() && !associated_edge() ) increment();
}


template<class Tds>
inline
bool
Triangulation_dsul_edge_iterator_2<Tds> ::
operator==(const Edge_iterator& fi) const
{
  if (pos==NULL || pos==_tds->past_end())
    return (pos==fi.pos) && (_tds==fi._tds) ;
  return (pos==fi.pos) && (_tds==fi._tds) &&  (index == fi.index);
}

template<class Tds>
inline
void
Triangulation_dsul_edge_iterator_2<Tds> ::
increment()
{
  CGAL_triangulation_precondition(_tds->dimension() >= 1);
  if (_tds->dimension() == 1) {Iterator_base::increment();}
  else {
    if (index == 2) {
      Iterator_base::increment();
      index =  0;
    }
    else{
      index +=1;
    }
  }
  return;
}

template<class Tds>
inline
void
Triangulation_dsul_edge_iterator_2<Tds> ::
decrement()
{
  CGAL_triangulation_precondition(_tds->dimension() >= 1);
  if (_tds->dimension() == 1) {Iterator_base::decrement();}
  else {
    if (index == 0) {
      Iterator_base::decrement();
      index = 2;
    }
    else {index -= 1;}
  }
  return;
}

template<class Tds>
inline
bool
Triangulation_dsul_edge_iterator_2<Tds> ::
associated_edge()
{
  if (_tds->dimension() == 1) {return true;}
  return std::less<const Face*>()(pos, pos->neighbor(index));
}

template<class Tds>
inline
Triangulation_dsul_edge_iterator_2<Tds>&
Triangulation_dsul_edge_iterator_2<Tds> ::
operator++()
{
  CGAL_triangulation_precondition(pos != NULL && pos != _tds->past_end());
  do     increment();
  while( pos != _tds->past_end() && !associated_edge());
  return *this;
}
    

template<class Tds>
inline
Triangulation_dsul_edge_iterator_2<Tds>&
Triangulation_dsul_edge_iterator_2<Tds> ::
operator--()
{
  CGAL_triangulation_precondition(pos != NULL && *this != Edge_iterator(_tds));
  do      decrement();
  while ( !associated_edge() && *this != Edge_iterator(_tds) ); 
  return *this;
}

    
template<class Tds>
inline
Triangulation_dsul_edge_iterator_2<Tds>
Triangulation_dsul_edge_iterator_2<Tds> ::    
operator++(int)
{
  Edge_iterator tmp(*this);
  ++(*this);
  return tmp;
}
    
template<class Tds>
inline
Triangulation_dsul_edge_iterator_2<Tds>
Triangulation_dsul_edge_iterator_2<Tds> ::     
operator--(int)
{
  Edge_iterator tmp(*this);
  --(*this);
  return tmp;
}
    

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DSUL_ITERATORS_2_H

