// ======================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
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
// file          : include/CGAL/Kernel_d/Tuple_d.h
// package       : Kernel_d
// maintainer    : Michael Seel <Michael.Seel@mpi-sb.mpg.de>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel
// coordinator   : MPI Saarbruecken (Susan.Hert@mpi-sb.mpg.de)
//
// ======================================================================
#ifndef CGAL_TUPLE_D_H
#define CGAL_TUPLE_D_H

#ifndef NOCGALINCL
#include <strstream>
#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#endif

CGAL_BEGIN_NAMESPACE
#define PointCd PointCd2

template <typename NT, typename LA> class PointHd;
template <typename NT, typename LA> class VectorHd;
template <typename NT, typename LA> class DirectionHd;
template <typename NT, typename LA> class HyperplaneHd;
template <typename NT, typename LA> class Aff_transformationHd;
template <typename FT, typename LA> class PointCd;
template <typename NT, typename LA> class VectorCd;
template <typename FT, typename LA> class DirectionCd;
template <typename FT, typename LA> class HyperplaneCd;
template <typename NT, typename LA> class Aff_transformationCd;

template <typename NT, typename LA>
class Tuple_d : public Ref_counted {
  typedef Tuple_d<NT,LA> Self;
  typedef typename LA::Vector Vector;
  Vector v;
public:
  typedef typename Vector::const_iterator const_iterator;

  struct Cartesian_const_iterator {
    typedef Cartesian_const_iterator self;
    typedef std::random_access_iterator_tag iterator_category;
    typedef CGAL::Quotient<NT>              value_type;
    typedef ptrdiff_t                       difference_type;
    typedef const value_type*               pointer;
    typedef const value_type&               reference;

  Cartesian_const_iterator() : _it(0), _w(0) {}
  Cartesian_const_iterator(const_iterator it, const_iterator w = 0) 
    : _it(it), _w(w) {}
    
  value_type operator*() const { return value_type(*_it,*_w); }

  self& operator++() { ++_it; return *this; }
  self  operator++(int) { self tmp = *this; ++_it; return tmp; }
  self& operator--() { --_it; return *this; }
  self  operator--(int) { self tmp = *this; --_it; return tmp; }

  self& operator+=(difference_type i) { _it+=i; return *this; }
  self& operator-=(difference_type i) { _it-=i; return *this; }
  self operator+(difference_type i) const 
  { self tmp=*this; return tmp += i; }
  self operator-(difference_type i) const 
  { self tmp=*this; return tmp -= i; }

  difference_type operator-(self x) const { return _it-x._it; }
  reference operator[](difference_type i) { return *(*this + i); }

  bool operator==(const self& x) const { return _it==x._it; }
  bool operator!=(const self& x) const { return ! (*this==x); }
  bool operator<(self x) const { (x - *this) > 0; }

  private:
    const_iterator _it, _w;  
  }; // Cartesian_const_iterator

  struct Homogeneous_const_iterator {
    typedef Homogeneous_const_iterator self;
    typedef std::random_access_iterator_tag iterator_category;
    typedef NT                              value_type;
    typedef ptrdiff_t                       difference_type;
    typedef const value_type*               pointer;
    typedef const value_type&               reference;

  Homogeneous_const_iterator() : _it(0), _w(0) {}
  Homogeneous_const_iterator(const_iterator it, const_iterator w = 0) 
    : _it(it), _w(w) {}
    
  value_type operator*() const 
  { if (_it == _w) return value_type(1); else return *_it; }

  self& operator++() { ++_it; return *this; }
  self  operator++(int) { self tmp = *this; ++_it; return tmp; }
  self& operator--() { --_it; return *this; }
  self  operator--(int) { self tmp = *this; --_it; return tmp; }

  self& operator+=(difference_type i) { _it+=i; return *this; }
  self& operator-=(difference_type i) { _it-=i; return *this; }
  self operator+(difference_type i) const 
  { self tmp=*this; return tmp += i; }
  self operator-(difference_type i) const 
  { self tmp=*this; return tmp -= i; }

  difference_type operator-(self x) const { return _it-x._it; }
  reference operator[](difference_type i) { return *(*this + i); }

  bool operator==(const self& x) const { return _it==x._it; }
  bool operator!=(const self& x) const { return ! (*this==x); }
  bool operator<(self x) const { (x - *this) > 0; }

  private:
    const_iterator _it, _w;  
  }; // Homogeneous_const_iterator


  Tuple_d(int d) : v(d) {}
  Tuple_d(const NT& a, const NT& b) : v(2)
  { v[0]=a; v[1]=b; }
  Tuple_d(const NT& a, const NT& b, const NT& c) : v(3)
  { v[0]=a; v[1]=b; v[2]=c; }
  Tuple_d(const NT& a, const NT& b, const NT& c, const NT& d) : v(4)
  { v[0]=a; v[1]=b; v[2]=c; v[3]=d; }

#ifndef CGAL_SIMPLE_INTERFACE

  template <typename I>
  Tuple_d(int d, I& start, I end) : v(d) 
  { int i(0); 
    while ( i < d && start != end ) v[i++] = *start++; 
  } 
  /* this constructor returns the final position of start 
     to offer access to a possible common denominator as
     part of the tuple range */

  template <typename I>
  Tuple_d(int d, I start, I end, NT D) : v(d) 
  { int i(0); 
    while ( i < d && start != end ) v[i++] = *start++; 
    v[d-1] = D; 
  }

#else // provide instantiated constructors:
#define FIXTUPLE(I) \
Tuple_d(int d, I& start, I end) : v(d) \
{ int i(0); while ( i < d && start != end ) v[i++] = *start++; } \
Tuple_d(int d, I start, I end, NT D) : v(d) \
{ int i(0); while ( i < d && start != end ) v[i++] = *start++; v[d-1] = D; }

FIXTUPLE(int*)
FIXTUPLE(const int*)
FIXTUPLE(NT*)
FIXTUPLE(const NT*)

#undef FIXTUPLE
#endif

  int size() const { return v.dimension(); }
  const_iterator begin() const { return v.begin(); }
  const_iterator last() const { return v.end()-1; }
  const_iterator end() const { return v.end(); }
  const_iterator beyondend() const { return v.end()+1; }


  void invert()
  { for (int i=0; i<size(); ++i) v[i]=-v[i]; }
  void invert(int d)
  { for (int i=0; i<d; ++i) v[i]=-v[i]; }

  void print(std::ostream& out, char*) const;
  void read(std::istream& in);
  void homogeneous_add(const Self* a, const Self* b)
  { register int d = a->size()-1;
    if ( d < 0 ) return;
    CGAL_assertion_msg((d == b->size()-1),"dimensions disagree."); 
    CGAL_assertion_msg((d == size()-1),"dimensions disagree."); 
    NT aw = a->v[d], bw = b->v[d]; 
    for (int i = 0; i < d; ++i) { 
      v[i] = a->v[i]*bw + b->v[i]*aw; 
    }
    v[d] = aw*bw; 
  }

  void homogeneous_sub(const Self* a, const Self* b)
  { register int d = a->size()-1; 
    if ( d < 0 ) return;
    CGAL_assertion_msg((d == b->size()-1),"dimensions disagree."); 
    CGAL_assertion_msg((d == size()-1),"dimensions disagree."); 
    NT aw = a->v[d], bw = b->v[d]; 
    for (int i = 0; i < d; ++i) { 
      v[i] = a->v[i]*bw - b->v[i]*aw; 
    }
    v[d] = aw*bw; 
  }

  void cartesian_add(const Self* a, const Self* b)
  { v = a->v + b->v; }
  void cartesian_sub(const Self* a, const Self* b)
  { v = a->v - b->v; }

  friend class PointHd<NT,LA>;
  friend class VectorHd<NT,LA>;
  friend class DirectionHd<NT,LA>;
  friend class HyperplaneHd<NT,LA>;
  friend class PointCd<NT,LA>;
  friend class VectorCd<NT,LA>;
  friend class DirectionCd<NT,LA>;
  friend class HyperplaneCd<NT,LA>;
  
}; // Tuple_d


template <class NT, class LA> 
class Compare_homogeneously
{ 
public:
Comparison_result operator()(
  const typename LA::Vector& v1, const typename LA::Vector& v2)
{ 
  CGAL_assertion_msg((v1.dimension() == v2.dimension()),
    "Compare_homogeneously: dimensions disagree.");
  NT aw = v1[v1.dimension()-1]; 
  NT bw = v2[v2.dimension()-1]; 
  CGAL_assertion(aw>0 && bw>0);
  for (int i = 0; i < v1.dimension()-1; i++ ) { 
    NT aibw = v1[i]*bw; 
    NT biaw = v2[i]*aw; 
    Comparison_result S = (aibw<biaw ? SMALLER : 
                          (biaw<aibw ? LARGER : EQUAL));
    if (S != EQUAL) return S;
  }
  return EQUAL; 
}
}; // Compare_homogeneously

template <class NT, class LA> 
class Compare_componentwise
{ public:
Comparison_result operator()(
  const typename LA::Vector& v1, const typename LA::Vector& v2)
{ 
  CGAL_assertion_msg((v1.dimension() == v2.dimension()),
  "Compare_coefficientwise: dimensions disagree.");
  for (int i = 0; i < v1.dimension(); i++ ) { 
    Comparison_result S = (v1[i]<v2[i] ? SMALLER : 
                          (v2[i]<v1[i] ? LARGER : EQUAL));
    if (S != EQUAL) return S;
  }
  return EQUAL;
}
}; // Compare_coefficientwise


template <typename NT, typename LA>
void Tuple_d<NT,LA>::print(std::ostream& os, char* l) const
{ int i;
  switch( os.iword(CGAL::IO::mode) ) {
    case CGAL::IO::ASCII :
      os << size() << " ";
      for (i = 0; i < size(); ++i) 
        os << v[i] << " "; break;
    case CGAL::IO::BINARY :
      CGAL::write(os, size()); 
      for (i = 0; i < size(); ++i) 
        CGAL::write(os, v[i]); break;
    default :
      os << l << "(" << size() << ", ";
      for (i = 0; i < size(); ++i) {
        os << v[i]; if (i!=size()-1) os<<", "; else os<<")";
      }
  }
}


template <typename NT, typename LA>
void Tuple_d<NT,LA>::read(std::istream& is)
{ int i = 0, d;
  switch( is.iword(CGAL::IO::mode) ) {
    case CGAL::IO::ASCII :
      is >> d; v = Vector(d);
      while (i < d && is >> v[i] ) ++i; 
      break;
    case CGAL::IO::BINARY :
      CGAL::read(is, d); v = Vector(d);
      while (i < d) { CGAL::read(is, v[i]); ++i; } break;
    default:
    CGAL_assertion_msg(0,"\nStream must be in ascii or binary mode\n"); 
  }
}

template <class ForwardIterator>
void tuple_dim_check(ForwardIterator first, ForwardIterator last, 
                     char* file, int line, char* op)
{ if (first==last) return;
  int d = first->dimension(); ++first;
  for (; first!=last; ++first) 
    if (first->dimension() != d) {
      std::ostrstream os;
      os << "Tuple Dimension Error " << 
            "File " << file << "Line " << line << "Operation " << op << '\0';
      CGAL_assertion_msg(0,os.str()); 
      os.freeze(0);
    }
}

#define TUPLE_DIM_CHECK(i1,i2,op) tuple_dim_check(i1,i2,__FILE__,__LINE__,#op)

template <class InputIterator, class OutputIterator>
int copy_and_count(InputIterator first, InputIterator last,
                   OutputIterator result) 
{ int n=0;
  while (first != last) { ++n; *result++ = *first++; }
  return n;
}

#undef PointCd
CGAL_END_NAMESPACE
#endif //CGAL_TUPLE_D_H

