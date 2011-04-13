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
// file          : include/CGAL/Union_find.h
// package       : Union_find
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Union-find implemented with path compression
// ============================================================================


#ifndef CGAL_UNION_FIND_H
#define CGAL_UNION_FIND_H

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <cassert>

#ifdef _MSC_VER
#define CGAL_UNION_FIND_NO_ALLOCATOR
#endif

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class UFP_, class V_, class R_, class P_>
class UF_forward_iterator {
  UFP_ p_;
public:
  typedef UF_forward_iterator<UFP_, V_, R_, P_> Self;
  typedef V_  value_type;
  typedef R_ reference;
  typedef P_ pointer;
  typedef std::forward_iterator_tag iterator_category;

  UF_forward_iterator() : p_(0) {}
  UF_forward_iterator(UFP_ p) : p_(p) {}

  bool  operator==( const Self& x) const
  { return p_ == x.p_; }
  bool  operator!=( const Self& x) const
  { return !(*this == x); }

  reference operator*()  const { return p_->value_; }
  pointer operator->() const { return &(p_->value_); }

  Self& operator++() { CGAL_assertion(p_); p_ = p_->next_; return *this; }
  Self  operator++(int) { Self tmp = *this; ++*this; return tmp; }

  UFP_ get_pointer() { return p_; }

}; // UF_forward_iterator

} // CGALi

template <typename T, typename A = CGAL_ALLOCATOR(T) > 
class Union_find;

/*{\Manpage {Union_find}{T,A}{Union-Find with path compression}{P}}*/

template <typename T, typename A>
class Union_find {

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
partition of values of type |T| into disjoint sets.  The type |A| has
to be a model of the allocator concept as defined in the C++
standard.}*/

class Union_find_struct {
  typedef Union_find_struct* pointer;
  friend class Union_find<T,A>;
  friend class CGALi::UF_forward_iterator<pointer,T,T&,T*>;
  friend class CGALi::UF_forward_iterator<pointer,T,const T&,const T*>;
protected:
  T              value_;
  size_t         size_;
  pointer next_, up_;
  Union_find_struct(pointer n, const T& x) : 
    value_(x), size_(1), next_(n), up_(0) {}
};


public:
/*{\Mtypes 6}*/ 

typedef Union_find<T,A> Self;
typedef Union_find_struct* pointer;

typedef T value_type;
/*{\Mtypemember values stored in items (equals |T|).}*/

typedef CGALi::UF_forward_iterator<pointer,T,T&,T*> handle;
/*{\Mtypemember handle to values.}*/

typedef CGALi::UF_forward_iterator<pointer,T,T&,T*> iterator;
/*{\Mtypemember iterator over values.}*/

typedef CGALi::UF_forward_iterator<pointer,T,const T&,const T*> 
  const_handle;

typedef CGALi::UF_forward_iterator<pointer,T,const T&,const T*> 
  const_iterator;

/*{\Mtext There are also constant versions |const_handle| and
 |const_iterator|.}*/

#ifndef CGAL_UNION_FIND_NO_ALLOCATOR

typedef typename A::template rebind<Union_find_struct>::other allocator;
/*{\Mtypemember allocator.}*/

#endif

/*{\Mcreation 5}*/

Union_find() : first_(0), sets_(0), values_(0) {}
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| and
initializes it to the empty partition.}*/

~Union_find() { clear(); }

private:
Union_find(const Self&) {}
Self& operator=(const Self&) { return *this; }

pointer find(pointer p) 
{ assert(p); pointer r = p;
  while (r->up_) { r = r->up_; } // now root is r;
  while (p->up_) { pointer u = p->up_; p->up_=r; p=u; }
  return r;
}

bool is_valid(handle v) const
{ return v != handle(0); }
bool is_valid(const_handle v) const
{ return v != const_handle(0); }

public:

/*{\Moperations 2 4}*/

#ifndef CGAL_UNION_FIND_NO_ALLOCATOR

allocator get_allocator() const { return a_; }
/*{\Mop the allocator of |\Mvar|.}*/

#endif

size_t number_of_sets() const { return sets_; }
/*{\Mop returns the number of disjoint sets of |\Mvar|.}*/

size_t size() const { return values_; }
/*{\Mop returns the number of values of |\Mvar|.}*/

size_t bytes() const { return values_*sizeof(Union_find_struct); }
/*{\Mop returns the memory consumed by |\Mvar|.}*/

size_t size(handle p)
/*{\Mop returns the size of the set containing $p$.}*/
{ return find(p).get_pointer()->size_; }
size_t size(const_handle p) const
{ return find(p).get_pointer()->size_; }

void clear();
/*{\Mopl reinitializes |\Mvar| to an empty partition.}*/

handle make_set(const T& x);
/*{\Mop creates a new singleton set containing |x| and returns a
handle to it.}*/

handle push_back(const T& x) { return make_set(x); }
/*{\Mop same as |make_set(x)|.}*/

template <class Forward_iterator>
void insert(Forward_iterator first, Forward_iterator beyond)
/*{\Mop insert the range of values referenced by |[first,beyond)|.
\precond value type of |Forward_iterator| is |T|.}*/
{ while (first != beyond) { push_back(*first++); } }

handle find(handle p) { return find(p.get_pointer()); }
/*{\Mop returns a canonical handle of the set that contains |p|,
i.e., |P.same_set(p,q)| iff  |P.find(p)| and |P.find(q)| return
the same handle. \precond |p| is a handle in |\Mvar|.}*/

const_handle find(const_handle p) const 
{ return const_cast<Self&>(*this)->find(p.get_pointer()); }

void unify_sets(handle p, handle q);
/*{\Mop unites the sets of partition |\Mvar| containing $p$ and
$q$. \precond $p$ and $q$ are in |\Mvar|.}*/

bool same_set(handle p, handle q) { return find(p)==find(q); }
/*{\Mop returns true iff $p$ and $q$ belong to the same set of
|\Mvar|. \precond $p$ and $q$ are in |\Mvar|.}*/

bool same_set(const_handle p, const_handle q) const
{ return find(p)==find(q); }

iterator begin() { return iterator(first_); }
/*{\Mop returns an iterator pointing to the first value of |\Mvar|.}*/
iterator end()   { return iterator(0); }
/*{\Mop returns an iterator pointing beyond the last value of |\Mvar|.}*/

const_iterator begin() const { return const_iterator(first_); }
const_iterator end()   const { return const_iterator(0); }

protected:
  pointer first_;
  size_t sets_;
  size_t values_;

#ifndef CGAL_UNION_FIND_NO_ALLOCATOR
  static allocator a_;
#endif

/*{\Mimplementation |\Mname| is implemented with union by rank and
path compression.  The running time for $m$ set operations on $n$
elements is $O(n \alpha(m,n))$ where $\alpha(m,n)$ is the extremly
slow growing inverse of Ackermann's function.}*/

};

#ifndef CGAL_UNION_FIND_NO_ALLOCATOR
template <typename T, typename A>
typename Union_find<T,A>::allocator Union_find<T,A>::a_;
#endif

template <typename T, typename A>
typename Union_find<T,A>::handle 
Union_find<T,A>::make_set(const T& x) 
{ pointer tmp = first_;
#ifndef CGAL_UNION_FIND_NO_ALLOCATOR
  first_ = a_.allocate(1);
  new ((void*)first_) Union_find_struct(tmp,x);
#else
  first_ = new Union_find_struct(tmp,x);
#endif
  ++sets_; ++values_;
  return handle(first_);
}

template <typename T, typename A>
void Union_find<T,A>::clear()
{ while (first_) { 
     pointer n = first_->next_;
#ifndef CGAL_UNION_FIND_NO_ALLOCATOR
     a_.destroy(first_); a_.deallocate(first_,1);
#else
     delete first_;
#endif
     first_ = n;
  }
  sets_ = values_ = 0;
}


template <typename T, typename A>
void
Union_find<T,A>::unify_sets(handle p, handle q)
{ CGAL_assertion(is_valid(p)&&is_valid(q));
  pointer pit = find(p.get_pointer());
  pointer qit = find(q.get_pointer());
  if (pit == qit) return;
  size_t sp = pit->size_, sq = qit->size_;
  if (sp > sq) std::swap(pit,qit); // now sp <= sq
  pit->up_ = qit;  // linking roots
  qit->size_ += pit->size_; // updating size
  sets_++;
}



CGAL_END_NAMESPACE

#endif // CGAL_UNION_FIND_H

