// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/Pluecker_line_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// 
// ============================================================================
#ifndef CGAL_PLUECKER_LINE_3_H
#define CGAL_PLUECKER_LINE_3_H

#include <CGAL/basic.h>
#include <algorithm>
#undef _DEBUG
#define _DEBUG 61
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

/*{\Manpage{Pluecker_line_3}{R}{Straight lines in 3-space}{pl}}*/  

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
directed straight line in the three-dimensional plane. Its
representation is based on Pluecker coordinates. (For a treatment see
the book on projected oriented geometry of Stolfi.)}*/

template <typename R_>
class Pluecker_line_3 {

/*{\Mtypes 6}*/
typedef R_ R;
/*{\Mtypemember the standard kernel type.}*/
typedef typename R::RT RT;
/*{\Mtypemember the ring type.}*/
typedef typename R::Point_3 Point_3;
/*{\Mtypemember the point type of the standard kernel.}*/
typedef typename R::Line_3  Line_3;
/*{\Mtypemember the line type of the standard kernel.}*/
typedef const RT* const_iterator;
/*{\Mtypemember iterator over Pluecker coefficients.}*/
typedef Pluecker_line_3<R> Self;
RT c_[6];

public:
/*{\Mcreation 3}*/     
Pluecker_line_3() {}
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| and
initializes it to some line.}*/

Pluecker_line_3(const Line_3& l)
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| and
initializes it to |l|.}*/
{
  Point_3 p(l.point(0)), q(l.point(1));
  c_[0] = p.hx()*q.hy() - p.hy()*q.hx();
  c_[1] = p.hx()*q.hz() - p.hz()*q.hx();
  c_[2] = p.hy()*q.hz() - p.hz()*q.hy();
  c_[3] = p.hx()*q.hw() - p.hw()*q.hx();
  c_[4] = p.hy()*q.hw() - p.hw()*q.hy();
  c_[5] = p.hz()*q.hw() - p.hw()*q.hz();
}

Pluecker_line_3(const Point_3& p, const Point_3& q)
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| and
initializes it to the oriented line through |p| and |q|.}*/
{
  c_[0] = p.hx()*q.hy() - p.hy()*q.hx();
  c_[1] = p.hx()*q.hz() - p.hz()*q.hx();
  c_[2] = p.hy()*q.hz() - p.hz()*q.hy();
  c_[3] = p.hx()*q.hw() - p.hw()*q.hx();
  c_[4] = p.hy()*q.hw() - p.hw()*q.hy();
  c_[5] = p.hz()*q.hw() - p.hw()*q.hz();

}

template <typename Forward_iterator>
Pluecker_line_3(Forward_iterator s, Forward_iterator e) 
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| and
initializes it to line with the parameters |tuple [s,e)|.}*/
{ int i=0;
  while (s!=e && i<6) c_[i++] = *s;
  CGAL_assertion(i==6);
}

/*{\Moperations 4 2 }*/

const_iterator begin() const { return c_; }
/*{\Mop returns an iterator pointing to the first
Pluecker coefficient.}*/                                                   
const_iterator end() const { return c_+6; }
/*{\Mop returns an iterator pointing beyond the last
Pluecker coefficient.}*/                                                   

Pluecker_line_3(const Self& l) 
{ std::copy(l.begin(),l.end(),c_); }

Self& operator=(const Self& l) 
{ if (&l!=this) std::copy(l.begin(),l.end(),c_); 
  return *this; }

const RT& operator[](unsigned i) const
/*{\Marrop returns constant access to the $i$th Pluecker coefficient.}*/    
{ CGAL_assertion(i<6); return c_[i]; }

int sign() const
/*{\Mop returns the sign of the first nonzero Pluecker coefficient
within the ordered tuple of coefficients.}*/
{ for (unsigned i=0; i<6; ++i) 
    if (c_[i]!=0) return CGAL_NTS sign(c_[i]); 
  return CGAL_NTS sign(c_[5]); 
}

void normalize()
/*{\Mop reduces the Pluecker coefficients to a minimal 
representation. This is done by dividing all Pluecker 
coefficients by their common gcd.}*/
{ 
  int i=0;
  while(i<6 && c_[i]==0)
    i++;
  
  if(i>5)
    return;

  RT D = c_[i];
  CGAL_assertion(D!=0);

  for(++i; i<6; ++i) 
    D = (c_[i]==0 ? D : CGAL_NTS gcd(D, c_[i]));
  if (D==0) return;
  for(int i=0; i<6; ++i) c_[i]/=D;
}

void negate()
/*{\Mop negates all Pluecker coefficients.}*/
{ for(int i=0; i<6; ++i) c_[i] = -c_[i]; }

Pluecker_line_3<R> opposite() const
/*{\Mop returns the line opposite to |\Mvar|. }*/
{ Pluecker_line_3<R> res;
  std::negate<RT> N;
  std::transform(begin(), end(), res.c_, N);
  return res;
}

static int cmp(const Pluecker_line_3<R>& l1,
               const Pluecker_line_3<R>& l2)
/*{\Mstatic returns the lexicographic order on lines defined
on their Pluecker coefficient tuples.}*/
{ for (unsigned i=0; i<5; ++i) {
    typename R::RT diff = l1[i]-l2[i];
    if ( diff != typename R::RT(0) ) return CGAL_NTS sign(diff);
  }
  return CGAL_NTS sign(l1[5]-l2[5]);
}

}; // Pluecker_line_3

/*{\Mimplementation The Pluecker coefficients of a line are stored
in a six-tuple of |RT| coefficients.}*/
                

template <typename R>
std::ostream& operator<<(std::ostream& os, const Pluecker_line_3<R>& l)
{ 
  switch( os.iword(CGAL::IO::mode) ) {
    case CGAL::IO::ASCII :
      for (unsigned i=0; i<6; ++i) os << l[i] << " "; 
      return os; 
    case CGAL::IO::BINARY :
      for (unsigned i=0; i<6; ++i) CGAL::write(os, l[i]);
      return os;
    default:
      os << "Pluecker_line_3(";
      for (unsigned i=0; i<5; ++i) os << l[i] << ", "; 
      os << l[5] << ")";
  } return os;
}

template <typename R>
Pluecker_line_3<R> categorize(const Pluecker_line_3<R>& l, int& inverted)
{ Pluecker_line_3<R> res(l);
  res.normalize();
  if ( res.sign()<0 ) { res.negate(); inverted=1; }
  else inverted=-1;
  return res;
}

struct Pluecker_line_lt {
  template <typename R>
  bool operator()(const CGAL::Pluecker_line_3<R>& l1,
                  const CGAL::Pluecker_line_3<R>& l2) const
  { return CGAL::Pluecker_line_3<R>::cmp(l1,l2)<0; }
};

CGAL_END_NAMESPACE
#endif //CGAL_PLUECKER_LINE_3_H

