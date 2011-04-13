// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/vector_MSC.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: A simple and constrained STL vector model for MSC
// ============================================================================

#ifndef CGAL_VECTOR_MSC_H
#define CGAL_VECTOR_MSC_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template <class I>
class vector_MSC {
  typedef vector_MSC<I> Self;
public:
  typedef I        value_type;
  typedef size_t   size_type;
  typedef I*       iterator;
  typedef const I* const_iterator;
protected:
  iterator s_,e_,a_;
  void init_by_value(iterator s, iterator e, const I& i=I())
  { for (; s!=e; ++s) *s=i; }
  void allocate(size_type n)
  { s_ = new I[n]; e_=a_=(s_+n); CGAL_assertion(s_!=0); }
public:
  template <class Iterator>
  iterator init_by_range(Iterator s, Iterator e)
  // copies [s,e) into range starting at s_
  { CGAL_assertion(s_!=0); iterator it;
    for (it=s_; s!=e; ++it,++s) *it=*s; return it; }

  vector_MSC() : s_(0),e_(0),a_(0) {}

  vector_MSC(size_type n) { allocate(n); }

  vector_MSC(size_type n, const I& i) 
  { allocate(n); init_by_value(s_,e_,i); }

  vector_MSC(const_iterator s, const_iterator e)
  { allocate(e-s); init_by_range(s,e); }

  vector_MSC(const Self& V) 
  { if (this == &V) return;
    allocate(V.size()); init_by_range(V.begin(),V.end()); }
  Self& operator=(const Self& V)
  { if (this == &V) return *this;
    allocate(V.size()); init_by_range(V.begin(),V.end());
    return *this;  }
  
  ~vector_MSC() { if (s_!=0) delete[] s_; }

  size_type size() const { return e_-s_; }  
  size_type capacity() const { return a_-e_; }  

  const I& back() const { CGAL_assertion(s_!=e_); return *(e_-1); }
  I& back() { CGAL_assertion(s_!=e_); return *(e_-1); }

  iterator begin() { return s_; }
  iterator end() { return e_; }
  const_iterator begin() const { return s_; }
  const_iterator end() const { return e_; }

  void check_index(int i) const { CGAL_assertion(0<=i&&i<int(size())); }
  void check_index(int i) { CGAL_assertion(0<=i&&i<int(size())); }

  const I& operator[](int i) const { check_index(i); return *(s_+i); }
  I& operator[](int i) { check_index(i); return *(s_+i); }

  void resize()
  { int n = (size()==0 ? 1 : 2*size());
    iterator so = s_, eo = e_;
    allocate(n); e_=init_by_range(so,eo);
    init_by_value(e_,a_,I());
    delete [] so;
  }

  void push_back(const I& i)
  { if ( capacity()==0 ) resize(); CGAL_assertion( capacity()>0 );
    ++e_; back()=i; }
  void pop_back()  
  { CGAL_assertion( size()>0 );
    if (a_>e_) *e_=I(); --e_; } 

  void print(std::ostream& os = std::cout) const
  { for(unsigned i=0; i<size(); ++i) os << *(s_+i) << ' '; }
  void read(std::istream& is) 
  { I i; while (is >> i) push_back(i); }

  void dump() const
  { std::cerr << s_ << " " << e_ << " " << a_ << " " << 
      size() << " " << capacity() << " = ";
    print(std::cerr); std::cerr << std::endl;}
  

};

CGAL_END_NAMESPACE

#endif //CGAL_VECTOR_MSC_H
