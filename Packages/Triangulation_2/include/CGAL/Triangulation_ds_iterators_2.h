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
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Triangulation_ds_iterators_2.h
// package       : Triangulation (3.7)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_TRIANGULATION_DS_ITERATORS_2_H
#define CGAL_TRIANGULATION_DS_ITERATORS_2_H



#include <utility>
#include <iterator>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE

// The  iterator_base visits all the full dimensional faces 
// of the Tds
// whatever may be the dimension of the Tds

template <class Tds>
class Triangulation_ds_iterator_base_2
  : public Triangulation_cw_ccw_2
{
public:
  typedef typename Tds::Face       value_type;
  typedef typename Tds::Face *     pointer;
  typedef typename Tds::Face &     reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

  typedef typename Tds::Geom_traits Geom_traits;
  typedef typename Geom_traits::Compare_y_2   Compare_y;
  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Face  Face;
  typedef typename Tds::Edge Edge;

  typedef Triangulation_ds_iterator_base_2<Tds> Iterator_base;

  Triangulation_ds_iterator_base_2()
     : _tds(NULL), pos(NULL)
        {}

  Triangulation_ds_iterator_base_2(const Tds* tds)
     : _tds(tds), pos(NULL)
  {
    if(_tds->dimension() < 0) { //no faces
      return;
    }
    pos = _tds->infinite_vertex()->face();
  }

 Triangulation_ds_iterator_base_2(const Tds* tds, int)
     :  _tds(tds), pos(NULL)
 {}

  Triangulation_ds_iterator_base_2(const Tds* tds, const Face* f)
    : _tds(tds), pos(f)
  {}
   
protected:
  const Tds*  _tds;
  const Face* pos;

  void increment();
  void decrement();
  int  maximum(const Face* f) const;

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
      CGAL_triangulation_precondition(pos != NULL);
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
class Triangulation_ds_face_iterator_2
  : public Triangulation_ds_iterator_base_2<Tds>
{
public:
  typedef Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef Triangulation_ds_face_iterator_2<Tds> Face_iterator;

  Triangulation_ds_face_iterator_2() : Iterator_base() {}
   
  Triangulation_ds_face_iterator_2(const Tds * tds) : Iterator_base(tds)
    {
      if (tds->dimension() <2) pos=NULL;
    }

  Triangulation_ds_face_iterator_2(const Tds* tds, int i)
    : Iterator_base(tds,i)
    {}

  Triangulation_ds_face_iterator_2(const Face_iterator& fi)
    : Iterator_base(fi._tds, fi.pos)
    {}
        
  Face_iterator&   operator=(const Face_iterator& fi);

  Face_iterator&  operator++();
  Face_iterator&  operator--();
  Face_iterator   operator++(int);
  Face_iterator   operator--(int);

};


template < class Tds>
class Triangulation_ds_vertex_iterator_2
: public Triangulation_ds_iterator_base_2<Tds>
{
public:
  typedef typename Tds::Vertex       Vertex;
  
  typedef Vertex       value_type;
  typedef Vertex *     pointer;
  typedef Vertex &     reference;

  typedef Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef Triangulation_ds_vertex_iterator_2<Tds> Vertex_iterator;
  typedef typename Triangulation_ds_iterator_base_2<Tds>::Face Face;

private :
  int index;

public:
  Triangulation_ds_vertex_iterator_2()
    : Iterator_base(), index(0)
    {}
    
  Triangulation_ds_vertex_iterator_2(const Tds* tds);
   
  Triangulation_ds_vertex_iterator_2(const Tds* tds, int i)
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
      CGAL_triangulation_assertion( pos != NULL);
      if (pos == (Face*)1) {// only one vertex;
	return *(_tds->infinite_vertex());
      }
      return *(pos->vertex(index));
    }

  Vertex*  operator->() const
    {
      CGAL_triangulation_assertion( pos != NULL);
      if (pos == (Face*)1) {// only one vertex;
	return _tds->infinite_vertex();
      }
      return pos->vertex(index);
    }

private:
  void increment();
  void decrement() ;
  bool associated_vertex();

 };


template <class Tds>
class Triangulation_ds_edge_iterator_2
 : public Triangulation_ds_iterator_base_2<Tds>
{
public:
  typedef typename Tds::Edge       Edge;
  typedef typename Tds::Geom_traits Geom_traits;

  typedef Edge            value_type;
  typedef Edge *          pointer;
  typedef Edge &          reference;

  typedef Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef Triangulation_ds_edge_iterator_2<Tds> Edge_iterator;
  typedef typename Triangulation_ds_iterator_base_2<Tds>::Face Face;

private:
int index;

public:
  Triangulation_ds_edge_iterator_2()
            : Iterator_base(),index(0)
        {}
 
  Triangulation_ds_edge_iterator_2(const Tds * tds);
   
     
  Triangulation_ds_edge_iterator_2(const Tds* tds, int i)
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
void 
Triangulation_ds_iterator_base_2<Tds>::
increment()
{
  CGAL_triangulation_precondition(_tds->dimension() >= 0);
  if (_tds->dimension() == 1 || _tds->dimension() == 0){
    pos = pos->neighbor(0);
    return;
  }

  int max = maximum(pos);
  Face* next=pos->neighbor(max);         // tentative first child
  Face* parent;
  int max2 = maximum(next);
  if ( next->neighbor(cw(max2))== pos){
    // next is the first child of pos
    pos = next;
    return;
  }
  // look for the second child of an ancestor of pos
  next=pos->neighbor(ccw(max));          // tentatative second child
  while (1){
    max2 = maximum(next);
    if ( next->neighbor(cw(max2))== pos) 
      // next is the second child of pos
      { pos = next; return;}
    while (1){
      parent = pos->neighbor(cw(max)); // go to parent
      max = maximum(parent);
      next=parent->neighbor(ccw(max));// tentatative second child
      if (next==pos)   // already coming back from this child
	{ pos = parent; continue; }
      else
	{ pos = parent; break; }
    }
  }
}
        

template<class Tds>
void 
Triangulation_ds_iterator_base_2<Tds>::
decrement()
{
  CGAL_triangulation_precondition(_tds->dimension() >= 0);
  if(_tds->dimension() == 0){
    pos = pos->neighbor(0);
    return;
  }

  if (_tds->dimension() == 1){
    pos = pos->neighbor(1);
    return;
  }

  // dimension() ==2
  int max = maximum(pos);
  Face* next=pos->neighbor(cw(max));     // parent of pos

  int max2 = maximum(next);
  if ( next->neighbor(max2) == pos)  
    // pos is the first child of next
    { pos = next; return;}
  pos = next->neighbor(max2);        // tentative first child
  max = maximum(pos);
  if ( pos->neighbor(cw(max))!=next) 
    // pos is not the first child of next
    { pos = next; return;}
  // look for "last" node in first subtree of next
  while (1){
    next = pos->neighbor(ccw(max));// tentatative second child
    max2 = maximum(next);
    if (next->neighbor(cw(max2))!= pos){
      //next is not the second child of pos
      next=pos->neighbor(max);         // tentative first child
      max2 = maximum(next);
      if ( next->neighbor(cw(max2))!= pos){
	//next is not first child of pos
	return;
      }
    }
    pos=next;
    max=max2;
  }
}

template<class Tds>
int 
Triangulation_ds_iterator_base_2<Tds>::
maximum(const Face* f) const
{
  Compare_y compare_y = _tds->geom_traits().compare_y_2_object();
  
  if ( _tds->is_infinite(f) ){
    return f->index( _tds->infinite_vertex() );
  }
  if(compare_y(f->vertex(0)->point(),
	       f->vertex(1)->point()) == CGAL::SMALLER)
    //  v0 < v1
    if(compare_y(f->vertex(2)->point(),
		 f->vertex(1)->point())==CGAL::SMALLER)
      //  v0,v2 < v1
      { return 1; }
    else
      //  v0 < v1 <= v2
      { return 2; }
  else
    //  v1 <= v0
        
    if(compare_y(f->vertex(1)->point(),
		 f->vertex(2)->point())!=CGAL::SMALLER)
      //  v2 <= v1 <= v0
      if(compare_y(f->vertex(0)->point(),
		   f->vertex(1)->point()) == CGAL::EQUAL)
	//  v2 <= v1 == v0
	{ return 1; }
      else
	//  v2 <= v1 < v0
	{ return 0; }
    else
      //  v1<=v0, v1<v2
      if(compare_y(f->vertex(0)->point(),
		   f->vertex(2)->point())
	 ==CGAL::SMALLER)
	//  v1 <= v0 < v2
	{ return 2; }
      else
	//  v1 < v2 <= v0
	{ return 0; }
        
}

template<class Tds>
inline
Triangulation_ds_iterator_base_2<Tds>&
Triangulation_ds_iterator_base_2<Tds> ::
operator=(const Iterator_base& fi)
{
  pos = fi.pos;
  _tds = fi._tds;
  return *this;
}

template<class Tds>
inline bool
Triangulation_ds_iterator_base_2<Tds> ::
operator==(const Iterator_base& fi) const
{
  return ((pos == fi.pos )&&(_tds==fi._tds));
}
        
template<class Tds>
inline bool
Triangulation_ds_iterator_base_2<Tds> ::        
operator!=(const Iterator_base& fi) const
{
  return !(*this == fi);
}

template<class Tds>
Triangulation_ds_iterator_base_2<Tds>&
Triangulation_ds_iterator_base_2<Tds> ::
operator++()
{
  CGAL_triangulation_precondition(pos != NULL);
  increment();
  if ( pos == _tds->infinite_face() ){ // complete tour
    pos = NULL;  
  }  
  return *this;           
}

template<class Tds>
Triangulation_ds_iterator_base_2<Tds>&
Triangulation_ds_iterator_base_2<Tds> ::
operator--()
{
  CGAL_triangulation_precondition( _tds != NULL && 
			       pos !=  _tds->infinite_face());
  // can't decrement first
  // to decrement past_the_end go to first and decrement
  if(pos == NULL)   pos = _tds->infinite_face(); 
  decrement(); 
  return *this;
}
        
template<class Tds>
Triangulation_ds_iterator_base_2<Tds>
Triangulation_ds_iterator_base_2<Tds> ::
operator++(int)
{
  Iterator_base tmp(*this);
  ++(*this);
  return tmp;
}
        
template<class Tds>
Triangulation_ds_iterator_base_2<Tds>
Triangulation_ds_iterator_base_2<Tds> ::        
operator--(int)
{
  Iterator_base tmp(*this);
  --(*this);
  return tmp;
}
        
// template<class Tds>
// inline 
// typename Tds::Face& 
// Triangulation_ds_iterator_base_2<Tds> ::
// operator*() const
// {
//   CGAL_triangulation_precondition(pos != NULL);
//   return const_cast<Face&>(*pos);
// }
    
// template<class Tds>
// inline 
// typename Tds::Face* 
// Triangulation_ds_iterator_base_2<Tds> ::
// operator->() const
// {
//   CGAL_triangulation_precondition(pos != NULL);
//   return const_cast<Face*>(pos);
// }

// Face iterator implementation
template<class Tds>
inline 
Triangulation_ds_face_iterator_2<Tds>&
Triangulation_ds_face_iterator_2<Tds>::
operator=(const Face_iterator& fi)
{
  Iterator_base::operator=(fi);
  return *this;
}

template<class Tds>
inline 
Triangulation_ds_face_iterator_2<Tds>&
Triangulation_ds_face_iterator_2<Tds>::
operator++()
{
  CGAL_triangulation_precondition(_tds != NULL &&_tds->dimension()==2);
  Iterator_base::operator++();    
  return *this;           
}

template<class Tds>
inline 
Triangulation_ds_face_iterator_2<Tds>&
Triangulation_ds_face_iterator_2<Tds>::
operator--()
{
  CGAL_triangulation_precondition(_tds != NULL &&_tds->dimension()==2);
  Iterator_base::operator--();    
  return *this;           
}

template<class Tds>
inline 
Triangulation_ds_face_iterator_2<Tds>
Triangulation_ds_face_iterator_2<Tds>::
operator++(int)
{
  Face_iterator tmp(*this);
  ++(*this);
  return tmp;
}
        
template<class Tds>
inline 
Triangulation_ds_face_iterator_2<Tds>
Triangulation_ds_face_iterator_2<Tds>::
operator--(int)
{
  Face_iterator tmp(*this);
  --(*this);
  return tmp;
}


// Vertex iterator implementation

template<class Tds>
Triangulation_ds_vertex_iterator_2<Tds> ::
Triangulation_ds_vertex_iterator_2(const Tds * tds)
  : Iterator_base(tds), index(0)
{
  switch( tds->number_of_vertices() ){
  case 0: // past-the-end
    pos = NULL;
    return;
  case 1:
    pos = (Face*)1 ; // different from any pointer;
    return;         // the iterator points to the only vertex 
	
  default:
    pos = tds->infinite_face();
    while ( *this!=Vertex_iterator(tds,1) && !associated_vertex()){
      increment();
    }
    return;
  }
}

template<class Tds>
void   
Triangulation_ds_vertex_iterator_2<Tds> ::
increment()
{
  CGAL_triangulation_precondition(_tds->dimension() >=0);
  if ( index==_tds->dimension()) {Iterator_base::increment(); index = 0;}
  else { index +=1; }
  return;
}

template<class Tds>
void   
Triangulation_ds_vertex_iterator_2<Tds> ::
decrement()
{
  CGAL_triangulation_precondition(_tds->dimension() >=0);
  if (index == 0) { 
    Iterator_base::decrement();
    index = _tds->dimension();
  }
  else {index -= 1;}
  return;
}

template<class Tds>
inline bool 
Triangulation_ds_vertex_iterator_2<Tds> ::
associated_vertex()
{
  return ( pos->vertex(index)->face() == pos ); // marche en toute dimension
}
   
template<class Tds>
Triangulation_ds_vertex_iterator_2<Tds>&
Triangulation_ds_vertex_iterator_2<Tds> ::
operator++()
{
   CGAL_triangulation_precondition(_tds != NULL );
   CGAL_triangulation_precondition(*this != Vertex_iterator(_tds,1));
   if (_tds->dimension()== -1) { //single vertex
     *this = Vertex_iterator(_tds,1);
     return *this;
   }
   do{
     increment();
     if (*this == Vertex_iterator(_tds)) *this= Vertex_iterator(_tds,1);
   }while ( *this!=Vertex_iterator(_tds,1) &&  !associated_vertex() );
   return *this;
}
    
template<class Tds>
Triangulation_ds_vertex_iterator_2<Tds>&
Triangulation_ds_vertex_iterator_2<Tds> ::    
operator--()
{
   CGAL_triangulation_precondition(_tds != NULL );
   CGAL_triangulation_assertion( *this != Vertex_iterator(_tds));
   if (_tds->dimension()== -1) { //single vertex
     *this = Vertex_iterator(_tds);
     return *this;
   }
   if( *this == Vertex_iterator(_tds,1)){
     *this=Vertex_iterator(_tds) ;
   }
   do   
   decrement();
   while ( *this != Vertex_iterator(_tds) && !associated_vertex()); 
   return *this;
}
    
template<class Tds>
Triangulation_ds_vertex_iterator_2<Tds>
Triangulation_ds_vertex_iterator_2<Tds> ::      
operator++(int)
{
  Vertex_iterator tmp(*this);
  ++(*this);
  return tmp;
}

template<class Tds>
Triangulation_ds_vertex_iterator_2<Tds>
Triangulation_ds_vertex_iterator_2<Tds> ::      
operator--(int)
{
  Vertex_iterator tmp(*this);
  --(*this);
  return tmp;
}
    
template<class Tds>
bool
Triangulation_ds_vertex_iterator_2<Tds> :: 
operator==(const Vertex_iterator& fi) const
{
  if (pos==NULL || pos==(Face*)1)
    return (pos==fi.pos) && (_tds==fi._tds) ;
  return (pos==fi.pos) && (_tds==fi._tds) &&  (index == fi.index);
}
 
template<class Tds>
inline bool
Triangulation_ds_vertex_iterator_2<Tds> ::     
operator!=(const Vertex_iterator& fi) const
{
  return !(*this == fi);
}
    
// template<class Tds>
// inline 
// typename Tds::Vertex&
// Triangulation_ds_vertex_iterator_2<Tds> :: 
// operator*() const
// {
//   CGAL_triangulation_assertion( pos != NULL);
//   if (pos == (Face*)1) {// only one vertex;
//       return *(_tds->infinite_vertex());
//   }
//   return *(pos->vertex(index));
// }
    
// template<class Tds>
// inline 
// typename Tds::Vertex*
// Triangulation_ds_vertex_iterator_2<Tds> :: 
// operator->() const
// {
//   CGAL_triangulation_assertion( pos != NULL);
//   if (pos == (Face*)1) {// only one vertex;
//       return _tds->infinite_vertex();
//   }
//   return pos->vertex(index);
// }


// Edge iterator implementation

template<class Tds>
Triangulation_ds_edge_iterator_2<Tds> ::
Triangulation_ds_edge_iterator_2(const Tds * tds)
 :  Iterator_base(tds), index(0) 
{
  if (_tds->dimension()<= 0){
    pos = NULL;                  // there is no edge
    return;
  }
  pos = _tds->infinite_face();
  index = (_tds->dimension() == 1) ? 2 : 0 ;
  while ( *this!= Edge_iterator(_tds,1) && !associated_edge() ) increment();
}


template<class Tds>
bool
Triangulation_ds_edge_iterator_2<Tds> ::
operator==(const Edge_iterator& fi) const
{
  return ( _tds==fi._tds &&
	   ( (pos==NULL && fi.pos==NULL) || 
	     (pos==fi.pos && index==fi.index)) );
}

template<class Tds>
void
Triangulation_ds_edge_iterator_2<Tds> ::
increment()
{
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
void
Triangulation_ds_edge_iterator_2<Tds> ::
decrement()
{
  if (_tds->dimension() == 1) {Iterator_base::decrement();}
  else {
    if (index == 0) {
      Iterator_base::decrement();
      index = 2;
    }
    else {index -= 1;}
    return;
  }
}

template<class Tds>
bool
Triangulation_ds_edge_iterator_2<Tds> ::
associated_edge()
{
  if (_tds->dimension() == 1) {return true;}
  return std::less<const Face*>()(pos, pos->neighbor(index));
}

template<class Tds>
Triangulation_ds_edge_iterator_2<Tds>&
Triangulation_ds_edge_iterator_2<Tds> ::
operator++()
{
  CGAL_triangulation_precondition(_tds != NULL  && _tds->dimension()>=1);
  CGAL_triangulation_precondition(*this != Edge_iterator(_tds,1));
  do {
    increment();
    if (*this == Edge_iterator(_tds)) *this= Edge_iterator(_tds,1);
  } while( *this != Edge_iterator(_tds,1) && !associated_edge());
  return *this;
}
    

template<class Tds>
Triangulation_ds_edge_iterator_2<Tds>&
Triangulation_ds_edge_iterator_2<Tds> ::
operator--()
{
  CGAL_triangulation_precondition(_tds != NULL  && _tds->dimension()>=1);
  CGAL_triangulation_assertion(*this != Edge_iterator(_tds));
  if( *this == Edge_iterator(_tds,1)){
     *this=Edge_iterator(_tds) ;
   }
  do   
   decrement();
  while ( *this != Edge_iterator(_tds) && !associated_edge()); 
  return *this;
}

    
template<class Tds>
Triangulation_ds_edge_iterator_2<Tds>
Triangulation_ds_edge_iterator_2<Tds> ::    
operator++(int)
{
  Edge_iterator tmp(*this);
  ++(*this);
  return tmp;
}
    
template<class Tds>
Triangulation_ds_edge_iterator_2<Tds>
Triangulation_ds_edge_iterator_2<Tds> ::     
operator--(int)
{
  Edge_iterator tmp(*this);
  --(*this);
  return tmp;
}
    
// template<class Tds>
// typename Tds::Edge
// Triangulation_ds_edge_iterator_2<Tds> ::     
// operator*() const
// {
//   return std::make_pair(const_cast<Face*>(pos), index);
// }


CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DS_ITERATORS_2_H

