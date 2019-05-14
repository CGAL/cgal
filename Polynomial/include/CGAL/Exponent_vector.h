// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Hemmer
//
// ============================================================================

#ifndef CGAL_EXPONENT_VECTOR_H
#define CGAL_EXPONENT_VECTOR_H

#include <deque>
#include <iterator>
#include <algorithm>
#include <vector>
#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <CGAL/int.h>

#include <boost/operators.hpp>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {


class Exponent_vector :  
  public boost::less_than_comparable1< Exponent_vector >,
  public boost::equality_comparable1< Exponent_vector >
{
  std::vector<int> v; 
public:
  typedef  Exponent_vector Self; 

  Exponent_vector(){};

  Exponent_vector(int e0): v(1) {
    v[0]=e0;
  };
  Exponent_vector(int e0, int e1): v(2) {
    v[0]=e0; v[1]=e1; 
  };
  Exponent_vector(int e0, int e1, int e2): v(3) {
    v[0]=e0; v[1]=e1; v[2]=e2;
  };
  Exponent_vector(int e0, int e1, int e2, int e3): v(4) {
    v[0]=e0; v[1]=e1; v[2]=e2; v[3]=e3;
  };
    
  Exponent_vector(const std::vector<int>& v_): v(v_){};
  Exponent_vector(const Exponent_vector& ev): v(ev.v){};

  template <class InputIterator>
  Exponent_vector(InputIterator begin , InputIterator end)
    :v(begin,end){
    typedef typename std::iterator_traits<InputIterator>::value_type value_type;
    CGAL_USE_TYPE(value_type);
    CGAL_static_assertion(( ::boost::is_same<value_type, int>::value));
  }


  // mirror vector functions 
  typedef  std::vector<int>::value_type value_type; 
  typedef  std::vector<int>::pointer pointer;
  typedef  std::vector<int>::const_pointer const_pointer; 
  typedef  std::vector<int>::reference reference; 
  typedef  std::vector<int>::const_reference const_reference;
  typedef  std::vector<int>::size_type size_type;
  typedef  std::vector<int>::difference_type difference_type;
  typedef  std::vector<int>::iterator iterator;
  typedef  std::vector<int>::const_iterator const_iterator;
  typedef  std::vector<int>::reverse_iterator reverse_iterator;
  typedef  std::vector<int>::const_reverse_iterator const_reverse_iterator;

  iterator begin(){return v.begin();}
  iterator end(){return v.end();}
  const_iterator begin() const {return v.begin();}
  const_iterator end() const {return v.end();}
  reverse_iterator rbegin() {return v.rbegin();}
  reverse_iterator rend(){return v.rend();} 
  const_reverse_iterator rbegin() const {return v.rbegin();}
  const_reverse_iterator rend() const {return v.rend();} 
  size_type size() const {return v.size();}
  size_type max_size() const {return v.max_size();}
  size_type capacity() const {return v.capacity();}
  bool empty() const {return v.empty();}
  reference operator[](size_type n) { return v[n]; }
  const_reference operator[](size_type n) const {return v[n];}
  // vector& operator=(const vector&)
  
  void reserve(size_t s){v.reserve(s);}
  reference front(){return v.front();}
  const_reference front() const {return v.front();}
  reference back() {return v.back();}
  const_reference back() const {return v.back();}
  void push_back(const int& x) { v.push_back(x);}
  void pop_back() {v.pop_back();}
  void swap(Self& ev) {v.swap(ev.v);}
  iterator insert(iterator pos, const int& x){return v.insert(pos,x);}
  
  template <class InputIterator>
  void insert(iterator pos,InputIterator f, InputIterator l){
    v.insert(pos,f,l);
  }
  void insert(iterator pos, size_type n, const int& x){
    v.insert(pos,n,x);
  }
  iterator erase(iterator pos){return v.erase(pos);}
  iterator erase(iterator first, iterator last){return v.erase(first,last);}
  void clear(){v.clear();}
  void resize(size_type n, int t = 0){v.resize(n,t);}
  bool operator==(const Self& ev) const { return v == ev.v; }

  // this is the actual change 
  bool operator<( const Exponent_vector& ev ) const {
    return std::lexicographical_compare (
        this->rbegin(), this->rend(),  ev.rbegin(), ev.rend()); 
}
  
  void output_benchmark( std::ostream& os ) const {
    os << "( ";
    for( unsigned i = 0; i < size(); ++i ) {
      if( i != 0 )
        os << ", ";
      os << v.at(i); 
    }
    os << " )";
  }
};


inline bool is_valid(const Exponent_vector& ev) {
  Exponent_vector::const_iterator it;
  for(it = ev.begin(); it != ev.end();it++){
    if (CGAL::is_negative(*it)) return false;
  }
  return true; 
}

inline std::ostream& operator << (std::ostream& os, const Exponent_vector& ev) {
  Exponent_vector::const_iterator it;
  os << "(" ;
  for(it = ev.begin(); it != ev.end();it++){
    if (it == ev.begin()) {
      os << *it ;
    }
    else{
      os <<"," << *it ;
    }
  }
  os << ")" ;
  return os;
}

} //namespace CGAL

namespace std{
template <> inline 
void swap(CGAL::Exponent_vector& ev1, CGAL::Exponent_vector& ev2){
  ev1.swap(ev2);
}

}


#endif // CGAL_EXPONENT_VECTOR_H
