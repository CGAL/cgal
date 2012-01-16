// Copyright (c) 1997-2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_VECTOR___H
#define CGAL_VECTOR___H

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <CGAL/Kernel_d/debug.h>

#include <cmath>
#include <memory>
#include <new>
#include <iostream>
#include <vector>
#include <iterator>

namespace CGAL {
namespace Linear_Algebra {

template <class NT_, class AL_> class Vector_;
template <class NT_, class AL_> class Matrix_;

/*{\Msubst
<NT_,AL_>#
<NT,AL>#
Vector_#Vector
Matrix_#Matrix
}*/
/*{\Moptions print_title=yes}*/
/*{\Moptions outfile=Vector.man}*/

/*{\Xtext \headerline{Common Notation}
The following data types use the concept of iterator ranges as an
abstraction of tuples and sets. For an iterator range |[first,last)|
we define |S = set [first,last)| as the ordered tuple $(|S[0]|,|S[1]|,
\ldots |S[d-1]|)$ where $|S[i]| = |*| |++|^{(i)}|first|$ (the element
obtained by forwarding the iterator by operator |++| $i$ times and
then dereferencing it to get the value to which it points). We write
|d = size [first,last)|.  This extends the syntax of random access
iterators to input iterators.  If we index the tuple as above then we
require that $|++|^{(d)}|first == last|$ (note that |last| points
beyond the last element to be accepted).}*/

/*{\Manpage {Vector}{}{Vectors with NT Entries}{v}}*/

template <class NT_, class AL_> 
class Vector_
{
/*{\Mdefinition An instance of data type |Vector_| is a vector of
variables of number type |NT|.  Together with the type |Matrix_| it
realizes the basic operations of linear algebra.}*/

public:

/*{\Mtypes 5.5}*/
typedef NT_*       pointer;
typedef const NT_* const_pointer;

typedef NT_    NT;
/*{\Mtypemember the ring type of the components.}*/ 
typedef pointer iterator;
/*{\Mtypemember the iterator type for accessing components.}*/ 
typedef const_pointer const_iterator;
/*{\Mtypemember the const iterator type for accessing components.}*/ 

typedef AL_ allocator_type;
/*{\Xtypemember the allocator type.}*/ 

protected:
  friend class Matrix_<NT_,AL_>;
  NT* v_; int d_;
  static allocator_type MM;

  inline void allocate_vec_space(NT*& vi, int di)
  {
  /* We use this procedure to allocate memory. We first get an appropriate 
     piece of memory from the allocator and then initialize each cell 
     by an inplace new. */

    vi = MM.allocate(di);
    NT* p = vi + di - 1;
    while (p >= vi) { new (p) NT(0);  p--; }   
  }

  inline void deallocate_vec_space(NT*& vi, int di)
  {
  /* We use this procedure to deallocate memory. We have to free it by
     the allocator scheme. We first call the destructor for type NT for each
     cell of the array and then return the piece of memory to the memory
     manager. */

    NT* p = vi + di - 1;
    while (p >= vi)  { MM.destroy(p); p--; }  //af:  as proposed by sylvain
    MM.deallocate(vi, di);
    vi = (NT*)0;
  }

inline void 
check_dimensions(const Vector_<NT_,AL_>& vec) const
{ 
  CGAL_assertion_msg((d_ == vec.d_), 
    "Vector_::check_dimensions: object dimensions disagree.");
}

public:

/*{\Mcreation v 3}*/

Vector_() : v_(0),d_(0) {}
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|.}*/ 

Vector_(int d) 
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|. 
|\Mvar| is initialized to a vector of dimension $d$.}*/ 
{ CGAL_assertion_msg( d >= 0 , 
    "Vector_::constructor: negative dimension.");
  d_ = d; 
  v_ = (NT*)0;
  if (d_ > 0){ 
    allocate_vec_space(v_,d_);
    while (d--) v_[d] = NT(0);
  }
}

Vector_(int d, const NT& x) 
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|. 
|\Mvar| is initialized to a vector of dimension $d$ with entries |x|.}*/ 
{ 
  CGAL_assertion_msg( d >= 0 , 
    "Vector_::constructor: negative dimension.");
  d_ = d; v_ = (NT*)0;
  if (d_ > 0){ 
    allocate_vec_space(v_,d_);
    while (d--) v_[d] = x;
  }
}

template <class Forward_iterator>
Vector_(Forward_iterator first, Forward_iterator last)
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|; 
|\Mvar| is initialized to the vector with entries
|set [first,last)|. \require |Forward_iterator| has value type |NT|.}*/
{ 
#if defined _MSC_VER && _MSC_VER == 1300
  d_ = 0;
  Forward_iterator fit = first;
  while(fit++!=last) d_++;
#else
  d_ = static_cast<int>(std::distance(first, last));
#endif
  allocate_vec_space(v_,d_);
  iterator it = begin();
  while (first != last) { *it = *first; ++it; ++first; }
}

Vector_(const Vector_<NT_,AL_>& p)
{ d_ = p.d_;
  if (d_ > 0) allocate_vec_space(v_,d_);
  else v_ = (NT*)0;
  for(int i=0; i<d_; i++) { v_[i] = p.v_[i]; }
}


Vector_<NT_,AL_>& operator=(const Vector_<NT_,AL_>& vec)
{ 
  if (&vec == this)
    return *this;

  int n = vec.d_;
  if (n != d_) {
    if (d_ > 0) deallocate_vec_space(v_,d_);
    d_=n;
    if (n > 0) allocate_vec_space(v_,n);
    else v_ = (NT*)0;
  }

  while (n--) v_[n] = vec.v_[n];
  return *this;
}

~Vector_() 
{ if (d_ > 0) deallocate_vec_space(v_,d_); }

/*{\Moperations 3 4}*/

int  dimension() const { return d_; }
/*{\Mop returns the dimension of |\Mvar|.}*/ 

bool is_zero() const 
/*{\Mop returns true iff |\Mvar| is the zero vector.}*/ 
{ for(int i=0; i<d_; ++i) if (v_[i]!=NT(0)) return false; 
  return true; }
  
NT& operator[](int i)
/*{\Marrop returns $i$-th component of |\Mvar|.\\
           \precond $0\le i \le |v.dimension()-1|$. }*/
{ CGAL_assertion_msg((0<=i && i<d_), 
    "Vector_::operator[]: index out of range.");
  return v_[i];
}
  
const NT& operator[](int i) const
{ CGAL_assertion_msg((0<=i && i<d_), 
    "Vector_::operator[]: index out of range.");
  return v_[i];
}

iterator begin() { return v_; }
/*{\Mop iterator to the first component.}*/
iterator end() { return v_+d_; }
/*{\Mop iterator beyond the last component.}*/

/*{\Mtext The same operations |begin()|, |end()| exist for 
|const_iterator|.}*/

const_iterator begin() const { return v_; }
const_iterator end() const { return v_+d_; }

Vector_<NT_,AL_>  operator+(const Vector_<NT_,AL_>& v1) const;
/*{\Mbinop Addition. \precond\\ |v.dimension() == v1.dimension()|.}*/

Vector_<NT_,AL_>  operator-(const Vector_<NT_,AL_>& v1) const;
/*{\Mbinop Subtraction. \precond\\ |v.dimension() = v1.dimension()|.}*/

NT operator*(const Vector_<NT_,AL_>& v1) const;
/*{\Mbinop Inner Product. \precond\\ |v.dimension() = v1.dimension()|.}*/

Vector_<NT_,AL_> compmul(const NT& r) const;

Vector_<NT_,AL_>  operator-() const;
/*{\Munop Negation.}*/

Vector_<NT_,AL_>& operator+=(const Vector_<NT_,AL_>& v1);
/*{\Mbinop Addition plus assignment. \precond\\
|v.dimension() == v1.dimension()|.}*/

Vector_<NT_,AL_>& operator-=(const Vector_<NT_,AL_>& v1);
/*{\Mbinop Subtraction plus assignment. \precond\\ |v.dimension() == v1.dimension()|.}*/

Vector_<NT_,AL_>& operator*=(const NT& s);
/*{\Mbinop Scalar multiplication plus assignment.}*/

Vector_<NT_,AL_>& operator/=(const NT& s);
/*{\Mbinop Scalar division plus assignment.}*/
 

bool     operator==(const Vector_<NT_,AL_>& w) const;
bool     operator!=(const Vector_<NT_,AL_>& w) const 
{ return !(*this == w); }

static int  compare(const Vector_<NT_,AL_>&, 
                    const Vector_<NT_,AL_>&);

};


template <class NT, class AL> 

inline Vector_<NT,AL> operator*(const NT& r, const Vector_<NT,AL>& v)
/*{\Mbinopfunc Componentwise multiplication with number $r$.}*/
{ return v.compmul(r); }

template <class NT, class AL> 

inline Vector_<NT,AL> operator*(const Vector_<NT,AL>& v, const NT& r)
/*{\Mbinopfunc Componentwise multiplication with number $r$.}*/
{ return v.compmul(r); }

template <class NT_, class AL_> 
inline Vector_<NT_,AL_>& Vector_<NT_,AL_>::
operator+=(const Vector_<NT_,AL_>& vec)
{ 
  check_dimensions(vec);
  int n = d_;
  while (n--) v_[n] += vec.v_[n];
  return *this;
}

template <class NT_, class AL_> 
inline Vector_<NT_,AL_>& Vector_<NT_,AL_>::
operator-=(const Vector_<NT_,AL_>& vec)
{ 
  check_dimensions(vec);
  int n = d_;
  while (n--) v_[n] -= vec.v_[n];
  return *this;
}

template <class NT_, class AL_> 
inline Vector_<NT_,AL_>& Vector_<NT_,AL_>::
operator*=(const NT& s)
{ int n = d_;
  while (n--) v_[n] *= s;
  return *this;
}

template <class NT_, class AL_>
inline Vector_<NT_,AL_>& Vector_<NT_,AL_>::
operator/=(const NT& s)
{ int n = d_;
  while (n--) v_[n] /= s;
  return *this;
}

template <class NT_, class AL_> 
inline Vector_<NT_,AL_> Vector_<NT_,AL_>::
operator+(const Vector_<NT_,AL_>& vec) const
{ 
  check_dimensions(vec);
  int n = d_;
  Vector_<NT_,AL_> result(n);
  while (n--) result.v_[n] = v_[n]+vec.v_[n];
  return result;
}

template <class NT_, class AL_> 
inline Vector_<NT_,AL_> Vector_<NT_,AL_>::
operator-(const Vector_<NT_,AL_>& vec) const
{ 
  check_dimensions(vec);
  int n = d_;
  Vector_<NT_,AL_> result(n);
  while (n--) result.v_[n] = v_[n]-vec.v_[n];
  return result;
}

template <class NT_, class AL_> 
inline Vector_<NT_,AL_> Vector_<NT_,AL_>::
operator-() const  // unary minus
{ 
  int n = d_;
  Vector_<NT_,AL_> result(n);
  while (n--) result.v_[n] = -v_[n];
  return result;
}


template <class NT_, class AL_> 
inline Vector_<NT_,AL_> Vector_<NT_,AL_>::
compmul(const NT& x) const
{ 
  int n = d_;
  Vector_<NT_,AL_> result(n);
  while (n--) result.v_[n] = v_[n] * x;
  return result;
}


template <class NT_, class AL_> 
inline NT_ Vector_<NT_,AL_>::
operator*(const Vector_<NT_,AL_>& vec) const
{ 
  check_dimensions(vec);
  NT_ result=0;
  int n = d_;
  while (n--) result = result+v_[n]*vec.v_[n];
  return result;
}

template <class NT_, class AL_> 
inline bool Vector_<NT_,AL_>::
operator==(const Vector_<NT_,AL_>& vec)  const
{ if (vec.d_ != d_) return false;
  int i = 0;
  while ((i<d_) && (v_[i]==vec.v_[i])) i++;
  return (i==d_);
}

template <class NT_, class AL_> 
int Vector_<NT_,AL_>::
compare(const Vector_<NT_,AL_>& v1, const Vector_<NT_,AL_>& v2)
{ int i;
  v1.check_dimensions(v2);
  for(i=0; i < v1.dimension() && v1[i]==v2[i]; i++) {}
  if (i == v1.dimension()) return 0;
  return (v1[i] < v2[i]) ?  -1 : 1;
}

template <class NT_, class AL_> 
std::ostream& operator<<(std::ostream& os, const Vector_<NT_,AL_>& v)
/*{\Xbinopfunc  writes |\Mvar| componentwise to the output stream $O$.}*/
{ /* syntax: d x_0 x_1 ... x_d-1 */
    int d = v.dimension();
    switch (os.iword(CGAL::IO::mode)) {
    case CGAL::IO::BINARY:
        CGAL::write( os, d);
        for ( int i = 0; i < d; ++i)
            CGAL::write( os, v[i]);
        break;
    case CGAL::IO::ASCII:
        os << d;
        for ( int i = 0; i < d; ++i)
            os << ' ' << v[i];
        break;
    case CGAL::IO::PRETTY:
        os << "LA::Vector(" << d << " [";
        for ( int i = 0; i < d; ++i) {
            if ( i > 0)
                os << ',' << ' ';
            os << v[i];
        }
        os << "])";
        break;
    }
    return os;
}

template <class NT_, class AL_> 
std::istream& operator>>(std::istream& is, Vector_<NT_,AL_>& v)
/*{\Xbinopfunc  reads |\Mvar| componentwise from the input stream $I$.}*/
{ /* syntax: d x_0 x_1 ... x_d-1 */
  int d;
  switch (is.iword(CGAL::IO::mode)) {
    case CGAL::IO::ASCII :
    case CGAL::IO::BINARY :
      is >> d; 
      v = Vector_<NT_,AL_>(d);
      for ( int i = 0; i < d; ++i) {
          is >> v[i];
      }
      break;
    default:
      std::cerr<<"\nStream must be in ascii or binary mode"<<std::endl;
      break;
  }
  return is;
}


template <class NT_, class AL_>
typename Vector_<NT_,AL_>::allocator_type Vector_<NT_,AL_>::MM;

/*{\Ximplementation Vectors are implemented by arrays of type
|NT|. All operations on a vector |v| take time $O(|v.dimension()|)$,
except for |dimension()| and $[\ ]$ which take constant time. The space
requirement is $O(|v.dimension()|)$. }*/


} // Linear_Algebra
} // CGAL

#endif // CGAL__VECTOR___H
