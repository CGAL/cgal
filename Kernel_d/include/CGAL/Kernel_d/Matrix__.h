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

#ifndef CGAL_MATRIX___H
#define CGAL_MATRIX___H

#include <CGAL/Kernel_d/Vector__.h>
#include <new>
#include <cstddef>                 // for std::size_t, std::ptrdiff_t

namespace CGAL {
namespace Linear_Algebra {

template <typename ROW_, typename V_, typename R_, typename P_> 
class column_iterator_ {
  ROW_ row_; unsigned i_;
public:
  typedef column_iterator_ Self;
  typedef V_  value_type;
  typedef R_  reference;
  typedef P_  pointer;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::random_access_iterator_tag iterator_category;

  column_iterator_() : row_(),i_() {}
  column_iterator_(ROW_ row, unsigned i) : row_(row),i_(i) {}

  bool  operator==( const Self& x) const 
  { return row_ == x.row_ && i_ == x.i_; }
  bool  operator!=( const Self& x) const 
  { return !(*this == x); }

  R_ operator*()  const { return (**row_)[i_]; }
  P_ operator->() const { return (**row_)+i_; }

  Self& operator++() { ++row_; return *this; }
  Self  operator++(int) 
  { Self tmp = *this; ++*this; return tmp; }

  Self& operator--() { --row_; return *this; }
  Self  operator--(int) 
  { Self tmp = *this; --*this; return tmp; }

  Self operator+(difference_type i) const
  { return Self(row_+i,i_); }
  Self operator-(difference_type i) const
  { return Self(row_-i,i_); }
  difference_type operator-(const Self& x) const
  { return (row_ - x.row_); }
      
};

template <typename ROW_, typename V_, typename R_, typename P_> 
class component_iterator_ {
  ROW_ row_;  // pointer to row
  int i_, n_; // offset and limit
public:
  typedef component_iterator_ Self;
  typedef V_  value_type;
  typedef R_  reference;
  typedef P_  pointer;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag iterator_category;

  component_iterator_() : row_(),i_(),n_() {}
  component_iterator_(ROW_ row, int i, int n) 
    : row_(row),i_(i),n_(n) {}

  bool  operator==( const Self& x) const 
  { return row_==x.row_ && i_==x.i_; }
  bool  operator!=( const Self& x) const   
  { return !(*this == x); }

  R_    operator*()  const { return (**row_)[i_]; }
  P_    operator->() const { return (**row_)+i_; }

  Self& operator++() { ++i_; if (i_==n_) { ++row_; i_=0; } return *this; }
  Self  operator++(int) { Self tmp = *this; ++*this; return tmp; }
  Self& operator--() { --i_; if (i_<0) { --row_; i_=n_-1; } return *this; }
  Self  operator--(int) { Self tmp = *this; --*this; return tmp; }
     
};


/*{\Msubst
<NT_,AL_>#
<NT,AL>#
Vector_#Vector
Matrix_#Matrix
}*/
/*{\Moptions print_title=yes}*/
/*{\Moptions outfile=Matrix.man}*/
/*{\Manpage {Matrix}{}{Matrices with NT Entries}{M}}*/

template <class NT_, class AL_>
class Matrix_ 
{ 
/*{\Mdefinition An instance of data type |\Mname| is a matrix of
variables of number type |NT|. The types |\Mname| and |Vector_|
together realize many functions of basic linear algebra.}*/

public:

/*{\Mtypes 6}*/

typedef Vector_<NT_,AL_>* vector_pointer;
typedef const Vector_<NT_,AL_>* const_vector_pointer;

typedef NT_ NT;
/*{\Mtypemember the ring type of the components.}*/ 

typedef component_iterator_<vector_pointer*,NT,NT&,NT*> iterator;
/*{\Mtypemember bidirectional iterator for accessing all components
row-wise.}*/
typedef component_iterator_<vector_pointer*,NT,const NT&,const NT*> 
  const_iterator;

typedef NT* row_iterator;
/*{\Mtypemember random access iterator for accessing row
  entries.}*/
typedef const NT* row_const_iterator;

typedef column_iterator_<vector_pointer*,NT,NT&,NT*> column_iterator;
/*{\Mtypemember random access iterator for accessing column entries.}*/ 
typedef column_iterator_<vector_pointer*,NT,const NT&, const NT*> 
  column_const_iterator;

/*{\Mtext There also constant versions of the above iterators:
|const_iterator|, |row_const_iterator|, and |column_const_iterator|.}*/

class Identity {};
/*{\Mtypemember a tag class for identity initialization}*/

typedef Vector_<NT_,AL_> Vector;
/*{\Mtypemember the vector type used.}*/ 

protected:
vector_pointer* v_; int dm_,dn_; 

NT& elem(int i, int j) const { return v_[i]->v_[j]; }

typedef typename AL_::template rebind<vector_pointer>::other 
        allocator_type;
static allocator_type MM;

inline void allocate_mat_space(vector_pointer*& vi, int d)
{
  /* We use this procedure to allocate memory. We use our allocator
     memory allocation scheme. There we first get an appropriate piece
     of memory and then initialize each cell by an inplace new. */

  vi = MM.allocate(d); 
  vector_pointer* p = vi + d - 1; 
  while (p >= vi) { 
    new (p) vector_pointer*(0); p--;
  }
}

inline void deallocate_mat_space(vector_pointer*& vi, int d)
{
  /* deallocate memory via our AL_ object. */

  MM.deallocate(vi,d);
  vi = (vector_pointer*)0;
}

inline void check_dimensions(const Matrix_<NT_,AL_>& mat) const
{ 
  CGAL_assertion_msg((dm_ == mat.dm_ && dn_ == mat.dn_), 
    "Matrix::check_dimensions: incompatible matrix dimensions.") ;
}

public:

/*{\Mcreation 5}*/

Matrix_() : dm_(0),dn_(0) { v_ = (Vector**)0; }
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|.}*/
Matrix_(int n); 
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| of
dimension $n \times n$ initialized to the zero matrix.}*/
Matrix_(int m, int n); 
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| of 
dimension $m \times n$ initialized to the zero matrix.}*/
Matrix_(std::pair<int,int> p);
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| of dimension
|p.first|$\times$|p.second| initialized to the zero matrix.}*/
Matrix_(int n , const Identity&, const NT& x = NT(1) ); 
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| of
dimension $n \times n$ initialized to the identity matrix
(times |x|).}*/
Matrix_(int m, int n, const NT& x);
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| of 
dimension $m \times n$ initialized to the matrix with |x|
entries.}*/

template <class RAIterator>
void range_initialize(RAIterator first, RAIterator last,  
                      std::random_access_iterator_tag) 
{ typedef typename std::iterator_traits<RAIterator>::value_type value_type;
  typedef typename value_type::const_iterator const_iterator;
  dn_ = static_cast<int>(last-first);
  if (dn_ == 0) { dm_=0; v_=0; return; }
  dm_ = first->dimension(); 
  if (dm_ > 0) { 
    int i,j;
    allocate_mat_space(v_,dm_);
    for (i=0; i<dm_; i++) {
      v_[i] = new Vector(dn_);
      // for (int j = 0; j < dn_; j++) elem(i,j) = (*(first+j))[i];
    }
    const_iterator it;
    for (j=0; first != last; ++j, ++first) // column wise
      for (i=0, it=first->begin(); it != first->end(); ++i, ++it) // row wise
        elem(i,j) = *it;
  } else 
    v_ = (Vector**)0; 
}

template <class InputIterator>
void range_initialize(InputIterator first, InputIterator last, 
                 std::forward_iterator_tag) 
{ typedef typename std::iterator_traits<InputIterator>::value_type 
    value_type;
  std::vector<value_type> V(first,last);
  range_initialize(V.begin(),V.end(),std::random_access_iterator_tag());
}

template <class Forward_iterator>
Matrix_(Forward_iterator first, Forward_iterator last)
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|. Let $S$ be
the ordered set of $n$ column-vectors of common dimension $m$ as given
by the iterator range |[first,last)|.  |\Mvar| is initialized to an $m
\times n$ matrix with the columns as specified by $S$.  \precond
|Forward_iterator| has a value type |V| from which we require to
provide a iterator type |V::const_iterator|, to have |V::value_type ==
NT|.\\ Note that |Vector_| or |std::vector<NT>| fulfill these
requirements.}*/
{ typedef typename std::iterator_traits<Forward_iterator>::iterator_category 
    iterator_category;
  range_initialize(first,last,iterator_category()); }

Matrix_(const std::vector< Vector >& A) 
/*{\Mcreate creates an instance |\Mvar| of type |\Mname|. Let $A$ be
an array of $n$ column-vectors of common dimension $m$.  |\Mvar| is
initialized to an $m \times n$ matrix with the columns as specified by
$A$. }*/
{ range_initialize(A.begin(),A.end(),
    std::random_access_iterator_tag()); }

Matrix_(const Matrix_<NT_,AL_>&); 

Matrix_(const Vector&); 
/* creates a $d \times 1$ matrix */

Matrix_(int, int, NT**); 

Matrix_<NT_,AL_>& operator=(const Matrix_<NT_,AL_>&);

~Matrix_(); 

/*{\Moperations 3 4}*/

int row_dimension()  const {  return dm_; }
/*{\Mop returns $n$, the number of rows of |\Mvar|.}*/

int column_dimension()  const { return dn_; }
/*{\Mop returns $m$, the number of columns of |\Mvar|.}*/

std::pair<int,int> dimension() const 
/*{\Mop returns $(m,n)$, the dimension pair of |\Mvar|.}*/
{ return std::pair<int,int>(dm_,dn_); }

Vector& row(int i) const
/*{\Mop returns the $i$-th row of |\Mvar| (an $m$ - vector).\\
\precond  $0 \le i \le m - 1$. }*/
{ CGAL_assertion_msg((0<=i && i<dm_),"Matrix_: row index out of range.");
  return *v_[i]; 
}

Vector column(int i) const 
/*{\Mop returns the $i$-th column of |\Mvar| (an $n$ - vector).\\
\precond  $0 \le i \le n - 1$. }*/
{ return Vector(column_begin(i),column_end(i)); }

Vector to_vector() const 
{ 
  CGAL_assertion_msg((dn_==1), 
    "Matrix_::to_vector: cannot make vector from matrix.");
  return column(0); 
}

Vector_<NT_,AL_>& operator[](int i) const  
{ 
  CGAL_assertion_msg((0<=i && i<dm_), 
    "Matrix_::operator[]: index out of range.");
  return row(i); 
}

NT& operator()(int i, int j)
/*{\Mfunop returns $M_{ i,j }$. \\
\precond $0\le i\le m-1$ and $0\le j\le n-1$. }*/
{ CGAL_assertion_msg((0<=i && i<dm_), 
    "Matrix_::operator(): row index out of range.");
  CGAL_assertion_msg((0<=j && j<dn_), 
    "Matrix_::operator(): column index out of range.");
  return elem(i,j); 
}

NT operator()(int i, int j) const
{ 
  CGAL_assertion_msg((0<=i && i<dm_), 
    "Matrix_::operator(): row index out of range.");
  CGAL_assertion_msg((0<=j && j<dn_), 
    "Matrix_::operator(): column index out of range.");
  return elem(i,j); 
}

void swap_rows(int i, int j)
/*{\Mop swaps rows $i$ and $j$.
\precond $0\le i\le m-1$ and $0\le j\le m-1$.}*/
{ CGAL_assertion(0<=i && i<dm_ && 0<=j && j<dm_);
  std::swap(v_[i],v_[j]); 
}

void swap_columns(int i, int j) 
/*{\Mop swaps columns $i$ and $j$.
\precond $0\le i\le n-1$ and $0\le j\le n-1$.}*/
{ CGAL_assertion(0<=i && i<dn_ && 0<=j && j<dn_);
  for(int l = 0; l < dm_; l++) std::swap(elem(l,i),elem(l,j)); 
}

row_iterator row_begin(int i) 
/*{\Mop an iterator pointing to the first entry of the $i$th row.
\precond $0\le i\le m-1$.}*/
{ CGAL_assertion_msg((0<=i&&i<dm_),"Matrix: row index out of range.");
  return v_[i]->begin(); }

row_iterator row_end(int i)   
/*{\Mop an iterator pointing beyond the last entry of the $i$th row.
\precond $0\le i\le m-1$.}*/
{ CGAL_assertion_msg((0<=i&&i<dm_),"Matrix: row index out of range.");
  return v_[i]->end(); }

row_const_iterator row_begin(int i) const 
{ CGAL_assertion_msg(0<=i&&i<dm_,"Matrix: row index out of range.");
  return v_[i]->begin(); }
row_const_iterator row_end(int i) const  
{ CGAL_assertion_msg(0<=i&&i<dm_,"Matrix: row index out of range.");
  return v_[i]->end(); }

column_iterator column_begin(int i) 
/*{\Mop an iterator pointing to the first entry of the $i$th column.
\precond $0\le i\le n-1$.}*/
{ CGAL_assertion_msg(0<=i&&i<dn_,"Matrix: column index out of range.");
  return column_iterator(v_,i); }

column_iterator column_end(int i)   
/*{\Mop an iterator pointing beyond the last entry of the $i$th column.
\precond $0\le i\le n-1$.}*/
{ return column_begin(i)+dm_; }

column_const_iterator column_begin(int i) const 
{ CGAL_assertion_msg(0<=i&&i<dn_,"Matrix: column index out of range.");
  return column_const_iterator(v_,i); }
column_const_iterator column_end(int i)   const 
{ return column_begin(i)+dm_; }

iterator begin() { return iterator(v_,0,dn_); }
/*{\Mop an iterator pointing to the first entry of |\Mvar|.}*/
iterator end() { return iterator(v_+dm_,0,dn_); }
/*{\Mop an iterator pointing beyond the last entry of |\Mvar|.}*/

const_iterator begin() const { return const_iterator(v_,0,dn_); }
const_iterator end() const { return const_iterator(v_+dm_,0,dn_); }

/*{\Mtext The same operations exist for |row_const_iterator| and
|column_const_iterator|.}*/

bool  operator==(const Matrix_<NT_,AL_>& M1)  const; 
/*{\Mbinop Test for equality. }*/

bool  operator!=(const Matrix_<NT_,AL_>& M1)  const 
/*{\Mbinop Test for inequality. }*/
{ return !(*this == M1); }

/*{\Mtext \headerline{Arithmetic Operators}}*/
/*{\Mtext
\settowidth{\typewidth}{|Matrix_<NT,LA>m|}
\addtolength{\typewidth}{\colsep}
\callwidth2cm
\computewidths
\newcommand{\dimeq}[2]{ 
\\|M.row_dimension() == M1.row_dimension()| and
\\|M.column_dimension() == M1.column_dimension()|
}
}*/

Matrix_<NT_,AL_> operator+ (const Matrix_<NT_,AL_>& M1); 
/*{\Mbinop Addition. \precond \dimeq.}*/

Matrix_<NT_,AL_> operator- (const Matrix_<NT_,AL_>& M1); 
/*{\Mbinop Subtraction. \precond \dimeq.}*/

Matrix_<NT_,AL_> operator-(); // unary
/*{\Munop Negation.}*/

Matrix_<NT_,AL_>& operator-=(const Matrix_<NT_,AL_>&); 

Matrix_<NT_,AL_>& operator+=(const Matrix_<NT_,AL_>&); 

Matrix_<NT_,AL_> operator*(const Matrix_<NT_,AL_>& M1) const; 
/*{\Mbinop Multiplication. \precond \\ |\Mvar.column_dimension() = M1.row_dimension()|. }*/

Vector_<NT_,AL_> 
operator*(const Vector_<NT_,AL_>& vec) const
{  return ((*this) * Matrix_<NT_,AL_>(vec)).to_vector(); }
/*{\Mbinop  Multiplication with vector. \precond \\
|\Mvar.column_dimension() = vec.dimension()|.}*/

Matrix_<NT_,AL_> compmul(const NT& x) const; 

static int compare(const Matrix_<NT_,AL_>& M1, 
                   const Matrix_<NT_,AL_>& M2);

}; // end of class

/*{\Xtext \headerline{Input and Output}}*/

template <class NT_, class AL_> 
std::ostream&  operator<<(std::ostream& os, const Matrix_<NT_,AL_>& M);
/*{\Xbinopfunc writes matrix |\Mvar| row by row to the output stream |os|.}*/

template <class NT_, class AL_> 
std::istream&  operator>>(std::istream& is, Matrix_<NT_,AL_>& M);
/*{\Xbinopfunc reads matrix |\Mvar| row by row from the input stream |is|.}*/


template <class NT_, class AL_>
Matrix_<NT_,AL_>::
Matrix_(int dim) : dm_(dim),dn_(dim)
{ 
  CGAL_assertion_msg((dim >= 0), 
    "Matrix_::constructor: negative dimension.");
  if (dm_ > 0) { 
    allocate_mat_space(v_,dm_);
    for (int i=0; i<dm_; i++) 
      v_[i] = new Vector(dn_);
  } else 
    v_ = (Vector**)0; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_>::
Matrix_(int dim1, int dim2) : dm_(dim1),dn_(dim2)
{ 
  CGAL_assertion_msg((dim1>=0 && dim2>=0), 
    "Matrix_::constructor: negative dimension.");

  if (dm_ > 0) { 
    allocate_mat_space(v_,dm_);
    for (int i=0; i<dm_; i++) 
      v_[i] = new Vector(dn_); 
  } else 
    v_ = (Vector**)0; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_>::
Matrix_(std::pair<int,int> p) : dm_(p.first),dn_(p.second)
{ 
  CGAL_assertion_msg((dm_>=0 && dn_>=0), 
    "Matrix_::constructor: negative dimension.");
  if (dm_ > 0) { 
    allocate_mat_space(v_,dm_);
    for (int i=0; i<dm_; i++) 
      v_[i] = new Vector(dn_); 
  } else 
    v_ = (Vector**)0; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_>::
Matrix_(int dim, const Identity&, const NT& x) : dm_(dim),dn_(dim)
{ CGAL_assertion_msg((dim >= 0),
    "matrix::constructor: negative dimension.");
  if (dm_ > 0) { 
    allocate_mat_space(v_,dm_);
    for (int i=0; i<dm_; i++) 
      v_[i] = new Vector(dn_); 
    if (x!=NT(0)) for (int i=0; i<dm_; ++i) elem(i,i)=x;
  } else 
    v_ = (Vector**)0; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_>::
Matrix_(int dim1, int dim2, const NT& x) : dm_(dim1),dn_(dim2)
{ CGAL_assertion_msg((dim1>=0 && dim2>=0), 
    "Matrix_::constructor: negative dimension.");
  if (dm_ > 0) { 
    allocate_mat_space(v_,dm_);
    for (int i=0; i<dm_; ++i) v_[i] = new Vector(dn_,x);
  } else 
    v_ = (Vector**)0; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_>::
Matrix_(const Matrix_<NT_,AL_>& p) : dm_(p.dm_),dn_(p.dn_)
{ if (dm_ > 0) {  
    allocate_mat_space(v_,dm_);
    for (int i=0; i<dm_; i++) 
      v_[i] = new Vector(*p.v_[i]); 
  }
  else 
    v_ = (Vector_<NT_,AL_>**)0; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_>::
Matrix_(const Vector& v) : dm_(v.d_),dn_(1)
{ if (dm_>0) allocate_mat_space(v_,dm_);
  else v_ = (Vector_<NT_,AL_>**)0; 
  for(int i = 0; i < dm_; i++) { 
    v_[i] = new Vector(1); 
    elem(i,0) = v[i]; 
  }
}


template <class NT_, class AL_>
Matrix_<NT_,AL_>::
Matrix_(int dim1, int dim2, NT** p) : dm_(dim1),dn_(dim2)
{ 
  CGAL_assertion_msg((dim1 >= 0 && dim2 >= 0), 
    "Matrix_::constructor: negative dimension.");
  if (dm_ > 0) {
    allocate_mat_space(v_,dm_);
    for(int i=0; i<dm_; i++) { 
      v_[i] = new Vector_<NT_,AL_>(dn_); 
      for(int j=0; j<dn_; j++) 
        elem(i,j) = p[i][j]; 
    }
  } else 
    v_ = (Vector_<NT_,AL_>**)0; 
}


template <class NT_, class AL_>
Matrix_<NT_,AL_>& Matrix_<NT_,AL_>::
operator=(const Matrix_<NT_,AL_>& mat)
{ 
  if (&mat == this)
    return *this;

  int i,j; 
  if (dm_ != mat.dm_ || dn_ != mat.dn_) { 
    for(i=0; i<dm_; i++) delete v_[i]; 
    if (v_) deallocate_mat_space(v_,dm_);

    dm_ = mat.dm_; dn_ = mat.dn_; 
    if (dm_>0)
      allocate_mat_space(v_,dm_);
    for(i = 0; i < dm_; i++) 
      v_[i] = new Vector(dn_); 
  }

  for(i = 0; i < dm_; i++)
    for(j = 0; j < dn_; j++) 
      elem(i,j) = mat.elem(i,j); 
  return *this; 
}


template <class NT_, class AL_>
Matrix_<NT_,AL_>::
~Matrix_()  
{ 
  if (v_) {
    for (int i=0; i<dm_; i++) 
      delete v_[i];  
    deallocate_mat_space(v_,dm_);
  }
}


template <class NT_, class AL_>
inline bool Matrix_<NT_,AL_>::
operator==(const Matrix_<NT_,AL_>& x) const
{ 
  int i,j; 
  if (dm_ != x.dm_ || dn_ != x.dn_) 
    return false; 

  for(i = 0; i < dm_; i++)
    for(j = 0; j < dn_; j++)
      if (elem(i,j) != x.elem(i,j)) 
        return false; 
  return true; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_> Matrix_<NT_,AL_>::
operator+ (const Matrix_<NT_,AL_>& mat)
{ 
  int i,j; 
  check_dimensions(mat); 
  Matrix_<NT_,AL_> result(dm_,dn_); 
  for(i=0; i<dm_; i++)
    for(j=0; j<dn_; j++)
      result.elem(i,j) = elem(i,j) + mat.elem(i,j); 
  return result; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_> Matrix_<NT_,AL_>::
operator- (const Matrix_<NT_,AL_>& mat)
{ 
  int i,j; 
  check_dimensions(mat); 
  Matrix_<NT_,AL_> result(dm_,dn_); 
  for(i=0; i<dm_; i++)
    for(j=0; j<dn_; j++)
      result.elem(i,j) = elem(i,j) - mat.elem(i,j); 
  return result; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_> Matrix_<NT_,AL_>::
operator- ()  // unary
{ 
  int i,j; 
  Matrix_<NT_,AL_> result(dm_,dn_); 
  for(i=0; i<dm_; i++)
    for(j=0; j<dn_; j++)
      result.elem(i,j) = -elem(i,j); 
  return result; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_>& Matrix_<NT_,AL_>::
operator-= (const Matrix_<NT_,AL_>& mat) 
{ 
  int i,j; 
  check_dimensions(mat); 
  for(i=0; i<dm_; i++)
    for(j=0; j<dn_; j++)
      elem(i,j) -= mat.elem(i,j); 
  return *this; 
}

template <class NT_, class AL_>
Matrix_<NT_,AL_>& Matrix_<NT_,AL_>::
operator+= (const Matrix_<NT_,AL_>& mat) 
{ 
  int i,j; 
  check_dimensions(mat); 
  for(i=0; i<dm_; i++)
    for(j=0; j<dn_; j++)
      elem(i,j) += mat.elem(i,j); 
  return *this; 
}

template <class NT_, class AL_>
inline Matrix_<NT_,AL_> Matrix_<NT_,AL_>::
operator*(const Matrix_<NT_,AL_>& M1) const
{ CGAL_assertion_msg((dn_==M1.dm_), 
    "Matrix_::operator*: incompatible matrix types."); 
  Matrix_<NT_,AL_> result(dm_,M1.dn_); 
  int i,j,l;
  for (i=0; i<dm_; ++i)
    for (j=0; j<M1.dn_; ++j)
      for (l=0; l<dn_; ++l)
        result.elem(i,j) += elem(i,l)*M1.elem(l,j);
  return result;
}

template <class NT_, class AL_>
inline Matrix_<NT_,AL_> Matrix_<NT_,AL_>::
compmul(const NT& f) const
{ 
  int i,j; 
  Matrix_<NT_,AL_> result(dm_,dn_); 
  for(i=0; i<dm_; i++)
    for(j=0; j<dn_; j++)
      result.elem(i,j) = elem(i,j) *f; 
  return result; 
}


template <class NT, class AL>

Matrix_<NT,AL>  operator*(const NT& x, const Matrix_<NT,AL>& M)
/*{\Mbinopfunc Multiplication of every entry with |x|. }*/
{ return M.compmul(x); }

template <class NT, class AL>

Matrix_<NT,AL>  operator*(const Matrix_<NT,AL>& M, const NT& x)
/*{\Mbinopfunc Multiplication of every entry with |x|. }*/
{ return M.compmul(x); }



template <class NT_, class AL_> 
int Matrix_<NT_,AL_>::
compare(const Matrix_<NT_,AL_>& M1, const Matrix_<NT_,AL_>& M2) 
{ int i; int res;
  M1.check_dimensions(M2);
  for(i=0; i < M1.row_dimension() && 
      (res = compare(M1.row(i),M2.row(i))) != 0; i++) {}
  return res;
}


template <class NT_, class AL_> 
std::ostream&  operator<<(std::ostream& os, const Matrix_<NT_,AL_>& M)
{ 
  /* syntax: d1 d2 
             x_0,0    ... x_0,d1-1
                  d2-times
             x_d2-1,0 ... x_d2-1,d1-1 */

    int d = M.row_dimension();
    int k = M.column_dimension();
    switch (os.iword(CGAL::IO::mode)) {
    case CGAL::IO::BINARY:
        CGAL::write( os, d);
        CGAL::write( os, k);
        for ( int i = 0; i < d; ++i) {
            for ( int j = 0; j < k; ++j) {
                CGAL::write( os, M[i][j]);
            }
        }
        break;
    case CGAL::IO::ASCII:
        os << d << ' ' << k;
        for ( int i = 0; i < d; ++i) {
            for ( int j = 0; j < k; ++j) {
                os << ' ' << M[i][j];
            }
        }
        break;
    case CGAL::IO::PRETTY:
        os << "LA::Matrix((" << d << ", " << k << " [";
        for ( int i = 0; i < d; ++i) {
            for ( int j = 0; j < k; ++j) {
                if ( j != 0)
                    os << ',' << ' ';
                os << M[i][j];
            }
            if ( i != d)
                os << ",\n";
        }
        os << "])";
        break;
    }
    return os;
}

template <class NT_, class AL_> 
std::istream&  operator>>(std::istream& is, Matrix_<NT_,AL_>& M) 
{ 
  /* syntax: d1 d2 
             x_0,0  ... x_0,d1-1
                  d2-times
             x_d2,0 ... x_d2,d1-1 */

  int cdim, rdim, i;
  switch(is.iword(CGAL::IO::mode)) {
    case CGAL::IO::BINARY : 
      CGAL::read(is,rdim);
      CGAL::read(is,cdim);
      for (i=0; i<rdim*cdim; ++i)
        CGAL::read(is,M(i/rdim,i%cdim));
      break;
    case CGAL::IO::ASCII :
      is >> rdim >> cdim;
      M = Matrix_<NT_,AL_>(rdim,cdim);
      for (i=0; i<rdim*cdim; ++i)
        is >> M(i/rdim,i%cdim);
      break; 
    default:
      std::cerr<<"\nStream must be in ascii or binary mode"<<std::endl;
      break;
  }
  return is;
}

template <class NT_, class AL_>
typename Matrix_<NT_,AL_>::allocator_type Matrix_<NT_,AL_>::MM;


/*{\Ximplementation 
The data type |\Mname| is implemented by two-dimensional arrays of
variables of type |NT|. The memory layout is row oriented. Operation
|column| takes time $O(n)$, |row|, |dim1|, |dim2| take constant time,
and all other operations take time $O(nm)$.  The space requirement is
$O(nm)$.}*/


} // Linear_Algebra
} // CGAL

#endif // CGAL_MATRIX___H
