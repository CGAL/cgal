// Copyright (c) 2007 (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Kinetic_data_structures/include/CGAL/Tools/Log.h $
// $Id: Log.h 39016 2007-06-11 00:08:13Z drussel $
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>


#ifndef CGAL_UTILITY_MACROS_H
#define CGAL_UTILITY_MACROS_H

#define CGAL_SUBSCRIPT(type, expr) type& operator[](unsigned int i){ expr;}\
  const type& operator[](unsigned int i) const { expr;}

#define CGAL_COPY_CONSTRUCTOR(TC) TC(const TC &o){copy_from(o);}\
  TC& operator=(const TC &o) {copy_from(o); return *this;}

#define CGAL_GET(type, name, expr) const type &name() const{expr;}
#define CGAL_GETNR(type, name, expr) const type name() const{expr;}

#define CGAL_IS(name, expr) bool is_##name() const {expr;}

#define CGAL_SET(type, name, expr) void set_##name(const type &k) {expr;}

#define CGAL_FIELDRW(type, name, var) \
  const type &name() const {return var;}\
  void set_##name(const type &k) {var=k;}


#define CGAL_OUTPUT(type)\
  inline std::ostream& operator<<(std::ostream&out, const type &t){	\
    return t.write(out);						\
  }

#define CGAL_OUTPUT1(type)			\
  template <class A>						\
  inline std::ostream& operator<<(std::ostream&out, const type<A> &t){ \
    return t.write(out);						\
  }

#define CGAL_OUTPUT2(T)			\
  template <class A, class B>						\
  inline std::ostream& operator<<(std::ostream&out, const T<A,B> &t){ \
    return t.write(out);						\
  }

#define CGAL_ITERATOR(uc_name, lc_name, it_type, bexpr, eexpr)	\
  typedef it_type uc_name##_iterator;					\
  uc_name##_iterator lc_name##s_begin() {bexpr;}			\
  uc_name##_iterator lc_name##s_end() {eexpr;}

#define CGAL_CONST_ITERATOR(uc_name, lc_name, it_type, bexpr, eexpr) \
  typedef it_type uc_name##_const_iterator;				\
  uc_name##_const_iterator lc_name##s_begin() const {bexpr;}		\
  uc_name##_const_iterator lc_name##s_end() const {eexpr;}

#define CGAL_FIND(ucname, expr)				\
  ucname##_const_iterator find(ucname##_key k) const {expr;}	\
  ucname##_iterator find(ucname##_key k) {expr;}

#define CGAL_INSERT(ucname, expr)			\
  ucname##_iterator insert(ucname##_key k, const ucname &m) {expr;}

#define CGAL_INSERTNK(ucname, expr)			\
  ucname##_iterator insert(const ucname &m) {expr;}

#define CGAL_SIZE(lcname, expr)		\
  size_t number_of_##lcname() const {expr;}

#define CGAL_SWAP(type)			\
  inline void swap(type &a, type &b) {		\
    a.swap_with(b);				\
  }

#define CGAL_ISWAP(name)			\
  std::swap(name, o.name)


#define CGAL_IFNONEQUAL(a,b,cmp) if (a cmp b) return true;	\
  else if (b cmp a) return false;

#define CGAL_COMPARISONS1(ucname, field) bool operator==(const ucname &o) const { \
    return (field== o.field);						\
  }									\
  bool operator!=(const ucname &o) const {				\
    return (field!= o.field);						\
  }									\
  bool operator<(const ucname &o) const {				\
    return (field< o.field);						\
  }									\
  bool operator>(const ucname &o) const {				\
    return (field> o.field);						\
  }									\
  bool operator>=(const ucname &o) const {				\
    return (field>= o.field);						\
  }									\
  bool operator<=(const ucname &o) const {				\
    return (field<= o.field);						\
  }

#define CGAL_COMPARISONS2(ucname, a, b) bool operator==(const ucname &o) const { \
    return (a== o.a && b== o.b);					\
  }									\
  bool operator!=(const ucname &o) const {				\
    return (a!= o.a || b != o.b);					\
  }									\
  bool operator<(const ucname &o) const {				\
    if (a< o.a ) return true;						\
    else if (a > o.a) return false;					\
    else return b < o.b;						\
  }									\
  bool operator>(const ucname &o) const {				\
    if (a> o.a ) return true;						\
    else if (a < o.a) return false;					\
    else return b > o.b;						\
  }									\
  bool operator>=(const ucname &o) const {				\
    return !operator<(o);						\
  }									\
  bool operator<=(const ucname &o) const {				\
    return !operator>(o);						\
  }


#endif
