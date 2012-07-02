// Copyright (c) 2011   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Triangulation_2/include/CGAL/Triangulation_hyperbolic_traits_2.h $
// $Id: Hyperboloic_isometry_2.h 57323 2010-07-05 10:07:39Z sloriot $
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_ISOMETRY_2_H
#define CGAL_HYPERBOLIC_ISOMETRY_2_H

#include <complex>

#include <CGAL/number_utils.h> 

namespace CGAL {

/*
 z -> (m*z + n)/(conj(n)*z + conj(m))
 */
template<class Gt>
class Hyperbolic_isometry_2
{
public:
  typedef Gt Geom_traits;
  typedef typename Gt::FT FT;
  typedef typename Gt::Point Point;
  typedef std::complex<FT> complex;
  
  Hyperbolic_isometry_2(const complex& m = complex(1, 0), const complex& n = complex(0, 0))
    : m_(m), n_(n)
  {
  }
  
  ~Hyperbolic_isometry_2()
  {
  }
  
  complex DoAction(const complex& z) const
  {
    return (m_*z + n_)/(std::conj(n_)*z + std::conj(m_));
  }
  
  Point DoAction(const Point& p) const
  {
    complex z = complex(p.x(), p.y());
    complex result = DoAction(z);
    
    return Point(result.real(), result.imag());
  }
  
  Hyperbolic_isometry_2 inverse() const
  {
    return Hyperbolic_isometry_2(-std::conj(m_), n_);
  }
  
  Hyperbolic_isometry_2 operator * (const Hyperbolic_isometry_2& other) const
  {
    return Hyperbolic_isometry_2(m_*other.m_ + n_*std::conj(other.n_), m_*other.n_ + n_*std::conj(other.m_));
  }
  
  void file_output(std::ostream& os) const
  {
    os << m_ << " " << n_ << std::endl;
  }
  
  complex m() const
  {
    return m_;
  }
  
  complex n() const
  {
    return n_;
  }
  
private:
  complex m_;
  complex n_;
};
  
template <class Gt>
std::ostream&
operator<<(std::ostream& os, const Hyperbolic_isometry_2<Gt> &isometry)
{
  isometry.file_output(os);
  return os;
}  
  
} // namespace CGAL

#endif // CGAL_HYPERBOLIC_ISOMETRY_2_H
