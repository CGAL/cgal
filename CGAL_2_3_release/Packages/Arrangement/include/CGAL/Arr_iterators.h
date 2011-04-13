// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : ?/include/CGAL/Arr_iterators.h
// package       : arr (1.03)
// author(s)     : Iddo Hanniel
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
// ============================================================================
//
// ============================================================================
#ifndef CGAL_ARR_ITERATORS_H
#define CGAL_ARR_ITERATORS_H

//this part is taken from Lutz's polyhedron, with minor modifications for
//the planar map name convention

#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif

CGAL_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////////////////
//   ARR
////////////////////////////////////////////////////////////////////////



template < class Node, class It, class Ctg>
class _Arr_face_circ : public It {
  // Ptr      nt;    // The internal node ptr inherited from It.
public:
  typedef  It  Base;
  typedef  _Arr_face_circ<Node,It,Ctg> Self;

  typedef  Ctg                iterator_category;
  typedef  Node               value_type;
  typedef  ptrdiff_t          difference_type;
  typedef  size_t             size_type;
  typedef  value_type&        reference;
  typedef  const value_type&  const_reference;
  typedef  value_type*        pointer;
  typedef  const value_type*  const_pointer;


// CREATION
// --------

  //    _Arr_face_circ() : It(0) {}
  _Arr_face_circ() : It() {}
  //_Polyhedron_facet_circ( pointer p) : It(p) {}
  _Arr_face_circ( It i) : It(i) {}

  // OPERATIONS Forward Category
  // ---------------------------

  pointer  ptr() const { return & It::operator*();}

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == NULL);
    return It::operator==( It(NULL));
  }
  bool operator!=( CGAL_NULL_TYPE p) const {
    return !(*this == p);
  }
  bool operator==( const Self& i) const {
    return  It::operator==(i);
  }
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  Self& operator++() {
    nt = (*nt).next_halfedge();
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

  /*  When implemented change it
    Self& operator--() {
        nt = (*nt).previous_halfedge();
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    */

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  friend inline  iterator_category
  iterator_category( const Self&) {
    return iterator_category();
  }
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Circulator_tag
  query_circulator_or_iterator( const Self&) {
    return Circulator_tag();
  }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};


template < class Node, class It, class Ctg>
class _Arr_face_const_circ : public It {
  // Ptr      nt;    // The internal node ptr inherited from It.
public:
  typedef  It  Base;
  typedef  _Arr_face_const_circ<Node,It,Ctg> Self;

  typedef  Ctg                iterator_category;
  typedef  Node               value_type;
  typedef  ptrdiff_t          difference_type;
  typedef  size_t             size_type;
  typedef  value_type&        reference;
  typedef  const value_type&  const_reference;
  typedef  value_type*        pointer;
  typedef  const value_type*  const_pointer;

  // CREATION
  // --------

  _Arr_face_const_circ() : It(0) {}
  _Arr_face_const_circ( const_pointer p) : It(p) {}
  _Arr_face_const_circ( It i) : It(i) {}

  // OPERATIONS Forward Category
  // ---------------------------

  const_pointer  ptr() const { return & It::operator*();}

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == NULL);
    return It::operator==( It(NULL));
  }
  bool operator!=( CGAL_NULL_TYPE p) const {
    return !(*this == p);
  }
  bool operator==( const Self& i) const {
    return  It::operator==(i);
  }
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  Self& operator++() {
    nt = (*nt).next_halfedge();
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------
  /*
    Self& operator--() {
    nt = (*nt).previous_halfedge();
    return *this;
    }
    Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
    }
  */

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  friend inline  iterator_category
  iterator_category( const Self&) {
    return iterator_category();
  }
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Circulator_tag
  query_circulator_or_iterator( const Self&) {
    return Circulator_tag();
  }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};
template < class Node, class It, class Ctg>
class _Arr_vertex_circ : public It {
  // Ptr      nt;    // The internal node ptr inherited from It.
public:
  typedef  It  Base;
  typedef  _Arr_vertex_circ<Node,It,Ctg> Self;

  typedef  Ctg                iterator_category;
  typedef  Node               value_type;
  typedef  ptrdiff_t          difference_type;
  typedef  size_t             size_type;
  typedef  value_type&        reference;
  typedef  const value_type&  const_reference;
  typedef  value_type*        pointer;
  typedef  const value_type*  const_pointer;

  // CREATION
  // --------

  //    _Arr_vertex_circ() : It(0) {}
  _Arr_vertex_circ() : It() {}
  //_Polyhedron_vertex_circ( pointer p) : It(p) {}
  _Arr_vertex_circ( It i) : It(i) {}

  // OPERATIONS Forward Category
  // ---------------------------

  pointer  ptr() const { return & It::operator*();}

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == NULL);
    return It::operator==( It(NULL));
  }
  bool operator!=( CGAL_NULL_TYPE p) const {
    return !(*this == p);
  }
  bool operator==( const Self& i) const {
    return  It::operator==(i);
  }
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  Self& operator++() {
    nt = (*nt).next_halfedge()->twin();
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

  /*
    Self& operator--() {
    nt = (*nt).twin()->previous_halfedge();
    return *this;
    }
    Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
    }
  */

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  friend inline  iterator_category
  iterator_category( const Self&) {
    return iterator_category();
  }
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Circulator_tag
  query_circulator_or_iterator( const Self&) {
    return Circulator_tag();
  }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};


template < class Node, class It, class Ctg>
class _Arr_vertex_const_circ : public It {
  // Ptr      nt;    // The internal node ptr inherited from It.
public:
  typedef  It  Base;
  typedef  _Arr_vertex_const_circ<Node,It,Ctg> Self;

  typedef  Ctg                iterator_category;
  typedef  Node               value_type;
  typedef  ptrdiff_t          difference_type;
  typedef  size_t             size_type;
  typedef  value_type&        reference;
  typedef  const value_type&  const_reference;
  typedef  value_type*        pointer;
  typedef  const value_type*  const_pointer;

  // CREATION
  // --------

  _Arr_vertex_const_circ() : It(0) {}
  _Arr_vertex_const_circ( const_pointer p) : It(p) {}
  _Arr_vertex_const_circ( It i) : It(i) {}

  // OPERATIONS Forward Category
  // ---------------------------

  const_pointer  ptr() const { return & It::operator*();}

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == NULL);
    return It::operator==( It(NULL));
  }
  bool operator!=( CGAL_NULL_TYPE p) const {
    return !(*this == p);
  }
  bool operator==( const Self& i) const {
    return  It::operator==(i);
  }
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  Self& operator++() {
    nt = (*nt).next_halfedge()->twin();
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------
  /*
    Self& operator--() {
    nt = (*nt).twin()->previous_halfedge();
    return *this;
    }
    Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
    }
  */

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  friend inline  iterator_category
  iterator_category( const Self&) {
    return iterator_category();
  }
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Circulator_tag
  query_circulator_or_iterator( const Self&) {
    return Circulator_tag();
  }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};



///////////////////////////////////////////////////////////////////////
// ARR OVERLAP CIRCULATOR
// bidirectional circulator over overlapping edges - like above 

template < class Node, class It, class Ctg>
class Arr_overlap_circulator : public It {
  // Ptr      nt;    // The internal node ptr inherited from It.
public:
  typedef  It  Base;
  typedef  Arr_overlap_circulator<Node,It,Ctg> Self;

  typedef  Ctg                iterator_category;
  typedef  Node               value_type;
  typedef  ptrdiff_t          difference_type;
  typedef  size_t             size_type;
  typedef  value_type&        reference;
  typedef  const value_type&  const_reference;
  typedef  value_type*        pointer;
  typedef  const value_type*  const_pointer;

// CREATION
// --------

  Arr_overlap_circulator() : It() {}

  //try this
  Arr_overlap_circulator( pointer p) : It(p) {}
  Arr_overlap_circulator( It i) : It(i) {}

  // OPERATIONS Forward Category
  // ---------------------------

  pointer  ptr() const { return & It::operator*();}

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == NULL);
    return It::operator==( It(NULL));
  }
  bool operator!=( CGAL_NULL_TYPE p) const {
    return !(*this == p);
  }
  bool operator==( const Self& i) const {
    return  It::operator==(i);
  }
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  Self& operator++() {
    //for this Subcurve_node needs to make the circ a friend 
    nt = (*nt).past_end_child; 
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

  Self& operator--() {
    nt = (*nt).begin_child;
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  friend inline  iterator_category
  iterator_category( const Self&) {
    return iterator_category();
  }
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Circulator_tag
  query_circulator_or_iterator( const Self&) {
    return Circulator_tag();
  }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};


template < class Node, class It, class Ctg>
class Arr_overlap_const_circulator : public It {
  // Ptr      nt;    // The internal node ptr inherited from It.
public:
  typedef  It  Base;
  typedef  Arr_overlap_const_circulator<Node,It,Ctg> Self;

  typedef  Ctg                iterator_category;
  typedef  Node               value_type;
  typedef  ptrdiff_t          difference_type;
  typedef  size_t             size_type;
  typedef  value_type&        reference;
  typedef  const value_type&  const_reference;
  typedef  value_type*        pointer;
  typedef  const value_type*  const_pointer;

  // CREATION
  // --------

  Arr_overlap_const_circulator() : It(0) {}
  Arr_overlap_const_circulator( const_pointer p) : It(p) {}
  Arr_overlap_const_circulator( It i) : It(i) {}

  // OPERATIONS Forward Category
  // ---------------------------

  const_pointer  ptr() const { return & It::operator*();}

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == NULL);
    return It::operator==( It(NULL));
  }
  bool operator!=( CGAL_NULL_TYPE p) const {
    return !(*this == p);
  }
  bool operator==( const Self& i) const {
    return  It::operator==(i);
  }
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  Self& operator++() {
    nt = (*nt).past_end_child; //needs to be friend
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------
  Self& operator--() {
    nt = It((*nt).begin_child);
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  friend inline  iterator_category
  iterator_category( const Self&) {
    return iterator_category();
  }
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Circulator_tag
  query_circulator_or_iterator( const Self&) {
    return Circulator_tag();
  }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};






#ifdef BOGI 
//old implementation

#ifdef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class C, class Fct /*, class Ref, class Ptr*/>
#else
template < class C,
  class Fct>
//           class Ref = typename C::reference,
//           class Ptr = typename C::pointer>
#endif
class Arr_overlap_circulator {
protected:
  C        nt;    // The internal circulator.
public:
  typedef C  Circulator;
  //typedef Circulator_project<C,Fct,Ref,Ptr> Self;
  typedef Arr_overlap_circulator<C,Fct/*,Ref,Ptr*/> Self;

  //typedef typename  C::iterator_category   iterator_category;
  typedef std::bidirectional_iterator_tag   iterator_category;
  //typedef typename  Fct::result_type       value_type;
  typedef Fct       value_type;
  //typedef typename  C::difference_type     difference_type;
  typedef ptrdiff_t     difference_type;
  //typedef typename  C::size_type           size_type;
  typedef size_t           size_type;
  //typedef           Ref                    reference;
  //just for this special case (since Fct is an iterator) :
  typedef           Fct                    reference; 
  //typedef           Ptr                    pointer;
  typedef           C                    pointer;  //just for this special case

  // CREATION
  // --------

  Arr_overlap_circulator() {}
  Arr_overlap_circulator( Circulator j) : nt(j) {}

  // OPERATIONS Forward Category
  // ---------------------------

  Circulator  current_circulator() const { return nt;}
  //Ptr ptr() const {
  pointer ptr() const {
    //Fct fct;
    return nt;
    //return &(fct(*nt));
  }

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == NULL);
    return ( nt == NULL);
  }
  bool  operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
  bool  operator==( const Self& i) const { return ( nt == i.nt); }
  bool  operator!=( const Self& i) const { return !(*this == i); }
  //Ref   operator*()  const { return *ptr(); }
  //Ptr   operator->() const { return ptr(); }
  reference   operator*()  const { return Fct(nt); }
  pointer   operator->() const { return nt; } //this is a bit meaningless
  Self& operator++() {
    //++nt;
    nt=nt->past_end_child; //for this Edge_node needs to make the circ a friend
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

  Self& operator--() {
    //--nt;
    nt=nt->begin_child; //for this Edge_node needs to make the circ a friend
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }

  /*
    // OPERATIONS Random Access Category
    // ---------------------------------

    Self  min_circulator() const {
    return Self( nt.min_circulator());
    }
    Self& operator+=( difference_type n) {
    nt += n;
    return *this;
    }
    Self  operator+( difference_type n) const {
    Self tmp = *this;
    return tmp += n;
    }
    Self& operator-=( difference_type n) {
    return operator+=( -n);
    }
    Self  operator-( difference_type n) const {
    Self tmp = *this;
    return tmp += -n;
    }
    difference_type  operator-( const Self& i) const {
    return nt - i.nt;
    }
    Ref  operator[]( difference_type n) const {
    Self tmp = *this;
    tmp += n;
    return tmp.operator*();
    }
  */
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  friend inline  iterator_category
  iterator_category( const Self&) { return iterator_category(); }
  friend inline  value_type*
  value_type( const Self&) { return (value_type*)(0); }
  friend inline  difference_type*
  distance_type( const Self&) { return (difference_type*)(0); }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class C, class Fct/*, class Ref, class Ptr*/>
inline  Circulator_tag
query_circulator_or_iterator(const Arr_overlap_circulator<C,Fct/*,Ref,Ptr*/>&)
{
  return Circulator_tag();
}
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //

#endif //BOGI





CGAL_END_NAMESPACE


#endif 
// EOF //











