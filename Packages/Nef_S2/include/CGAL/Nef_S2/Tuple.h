#ifndef CGAL_TUPLE_H
#define CGAL_TUPLE_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template <class T, unsigned n> 
class Tuple {
  typedef Tuple<T,n> Self;
  T object_[n];
public:
  Tuple() { for (unsigned i=0; i<n; ++i) object_[i]=T(); }
  Tuple(const T& t1, const T& t2)
  { CGAL_nef_assertion(n>1); object_[0]=t1; object_[1]=t2; }
  Tuple(const T& t1, const T& t2, const T& t3)
  { CGAL_nef_assertion(n>2); object_[0]=t1; object_[1]=t2; object_[2]=t3; }
  
  Tuple(const Self& t) 
  { for (unsigned i=0; i<n; ++i) object_[i] = t.object_[i]; }
  Self& operator=(const Self& t) 
  { for (unsigned i=0; i<n; ++i) object_[i] = t.object_[i]; 
    return *this; }
  
  const T& operator[](unsigned i) const 
  { CGAL_nef_assertion(i<n); return object_[i]; }
  T& operator[](unsigned i) 
  { CGAL_nef_assertion(i<n); return object_[i]; }  

};

CGAL_END_NAMESPACE
#endif //CGAL_TUPLE_H
