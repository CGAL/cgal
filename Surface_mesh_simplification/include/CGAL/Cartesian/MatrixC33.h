// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_CARTESIAN_MATRIXC33_H
#define CGAL_CARTESIAN_MATRIXC33_H

#include <CGAL/license/Surface_mesh_simplification.h>


#include <CGAL/Vector_3.h>
#include <CGAL/determinant.h>
#include <CGAL/number_utils.h>
#include <CGAL/Null_matrix.h>

#include <boost/optional/optional.hpp>

namespace CGAL {

template <class R_>
class MatrixC33 
{
public:
  
  typedef R_ R ;
  
  typedef typename R::FT       RT ;
  typedef typename R::Vector_3 Vector_3;

  MatrixC33 ( Null_matrix )
   :
    mR0(NULL_VECTOR)
   ,mR1(NULL_VECTOR)
   ,mR2(NULL_VECTOR)
  {}
  
  MatrixC33 ( RT const& r0x, RT const& r0y, RT const& r0z
            , RT const& r1x, RT const& r1y, RT const& r1z
            , RT const& r2x, RT const& r2y, RT const& r2z
            )
    :
     mR0(r0x,r0y,r0z)
    ,mR1(r1x,r1y,r1z)
    ,mR2(r2x,r2y,r2z)
  {}
  
  MatrixC33 ( Vector_3 const& r0, Vector_3 const& r1, Vector_3 const& r2 )
    :
     mR0(r0)
    ,mR1(r1)
    ,mR2(r2)
  {}
  
  Vector_3 const& r0() const { return mR0; }
  Vector_3 const& r1() const { return mR1; }
  Vector_3 const& r2() const { return mR2; }
  
  Vector_3& r0() { return mR0; }
  Vector_3& r1() { return mR1; }
  Vector_3& r2() { return mR2; }

  Vector_3 const& operator[] ( int row ) const { return row == 0 ? mR0 : ( row == 1 ? mR1 : mR2 ) ; }  
  Vector_3&       operator[] ( int row )       { return row == 0 ? mR0 : ( row == 1 ? mR1 : mR2 ) ; }  
  
  MatrixC33& operator+= ( MatrixC33 const& m )
  {
    mR0 = mR0 + m.r0() ;
    mR1 = mR1 + m.r1() ;
    mR2 = mR2 + m.r2() ;
    return *this ;
  }
  
  MatrixC33& operator-= ( MatrixC33 const& m )
  {
    mR0 = mR0 - m.r0() ;
    mR1 = mR1 - m.r1() ;
    mR2 = mR2 - m.r2() ;
    return *this ;
  }
    
  MatrixC33& operator*= ( RT const& c )
  {
    mR0 = mR0 * c ;
    mR1 = mR1 * c ;
    mR2 = mR2 * c ;
    return *this ;
  }
  
  MatrixC33& operator/= ( RT const& c )
  {
    mR0 = mR0 / c ;
    mR1 = mR1 / c ;
    mR2 = mR2 / c ;
    return *this ;
  }
  
  friend MatrixC33 operator+ ( MatrixC33 const& a, MatrixC33 const& b )
  {
    return MatrixC33(a.r0()+b.r0()
                    ,a.r1()+b.r1()
                    ,a.r2()+b.r2()
                    );
  }
  
  friend MatrixC33 operator- ( MatrixC33 const& a, MatrixC33 const& b )
  {
    return MatrixC33(a.r0()-b.r0()
                    ,a.r1()-b.r1()
                    ,a.r2()-b.r2()
                    );
  }
  
  friend MatrixC33 operator* ( MatrixC33 const& m, RT const& c )
  {
    return MatrixC33(m.r0()*c,m.r1()*c,m.r2()*c);
  }
  friend MatrixC33 operator* ( RT const& c, MatrixC33 const& m )
  {
    return MatrixC33(m.r0()*c,m.r1()*c,m.r2()*c);
  }
  
  friend MatrixC33 operator/ ( MatrixC33 const& m, RT const& c )
  {
    return MatrixC33(m.r0()/c,m.r1()/c,m.r2()/c);
  }
  
  friend Vector_3 operator* ( MatrixC33 const& m, Vector_3 const& v )
  {
    return Vector_3(m.r0()*v,m.r1()*v,m.r2()*v);
  }
  friend Vector_3 operator* ( Vector_3 const& v, MatrixC33 const& m )
  {
    return Vector_3(v*m.r0(),v*m.r1(),v*m.r2());
  }
 
  RT determinant() const
  {
    return CGAL::determinant(r0().x(),r0().y(),r0().z()
                            ,r1().x(),r1().y(),r1().z()
                            ,r2().x(),r2().y(),r2().z()
                            );
  }
  
  MatrixC33& transpose()
  {
    mR0 = Vector_3(r0().x(),r1().x(),r2().x());
    mR1 = Vector_3(r0().y(),r1().y(),r2().y()); 
    mR2 = Vector_3(r0().z(),r1().z(),r2().z());
    return *this ;
  }
    
private:

  Vector_3 mR0 ;
  Vector_3 mR1 ;
  Vector_3 mR2 ;
} ;

template<class R>
inline
MatrixC33<R> direct_product ( Vector_3<R> const& u, Vector_3<R> const& v )
{
  return MatrixC33<R>( v * u.x()
                     , v * u.y()
                     , v * u.z()
                     ) ;
}

template<class R>
MatrixC33<R> transposed_matrix ( MatrixC33<R> const& m )
{
  MatrixC33<R> copy = m ;
  copy.Transpose();
  return copy ;
}

template<class R>
MatrixC33<R> cofactors_matrix ( MatrixC33<R> const& m )
{
  typedef typename R::RT RT ;
  
  RT c00 =  determinant(m.r1().y(),m.r1().z(),m.r2().y(),m.r2().z());
  RT c01 = -determinant(m.r1().x(),m.r1().z(),m.r2().x(),m.r2().z());
  RT c02 =  determinant(m.r1().x(),m.r1().y(),m.r2().x(),m.r2().y());
  
  RT c10 = -determinant(m.r0().y(),m.r0().z(),m.r2().y(),m.r2().z());
  RT c11 =  determinant(m.r0().x(),m.r0().z(),m.r2().x(),m.r2().z());
  RT c12 = -determinant(m.r0().x(),m.r0().y(),m.r2().x(),m.r2().y());
  
  RT c20 =  determinant(m.r0().y(),m.r0().z(),m.r1().y(),m.r1().z());
  RT c21 = -determinant(m.r0().x(),m.r0().z(),m.r1().x(),m.r1().z());
  RT c22 =  determinant(m.r0().x(),m.r0().y(),m.r1().x(),m.r1().y());
  
  return MatrixC33<R>(c00,c01,c02
                     ,c10,c11,c12
                     ,c20,c21,c22
                     );
}

template<class R>
MatrixC33<R> adjoint_matrix ( MatrixC33<R> const& m )
{
  return cofactors_matrix(m).transpose()  ;
}

template<class R>
boost::optional< MatrixC33<R> > inverse_matrix ( MatrixC33<R> const& m )
{
  typedef typename R::RT RT ;
  
  typedef MatrixC33<R> Matrix ;
  
  typedef boost::optional<Matrix> result_type ;
   
  result_type rInverse ;
  
  RT det = m.determinant();
  
  if ( ! CGAL_NTS is_zero(det) )
  {
    RT c00 = (m.r1().y()*m.r2().z() - m.r1().z()*m.r2().y()) / det; 
    RT c01 = (m.r2().y()*m.r0().z() - m.r0().y()*m.r2().z()) / det;
    RT c02 = (m.r0().y()*m.r1().z() - m.r1().y()*m.r0().z()) / det; 
  
    RT c10 = (m.r1().z()*m.r2().x() - m.r1().x()*m.r2().z()) / det; 
    RT c11 = (m.r0().x()*m.r2().z() - m.r2().x()*m.r0().z()) / det; 
    RT c12 = (m.r1().x()*m.r0().z() - m.r0().x()*m.r1().z()) / det; 
  
    RT c20 = (m.r1().x()*m.r2().y() - m.r2().x()*m.r1().y()) / det; 
    RT c21 = (m.r2().x()*m.r0().y() - m.r0().x()*m.r2().y()) / det; 
    RT c22 = (m.r0().x()*m.r1().y() - m.r0().y()*m.r1().x()) / det; 
    
    rInverse = result_type( Matrix(c00,c01,c02
                                  ,c10,c11,c12
                                  ,c20,c21,c22
                                  )
                          ) ;
  }
  
  return rInverse ;
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_MATRIXC33_H //
// EOF //
 
 
