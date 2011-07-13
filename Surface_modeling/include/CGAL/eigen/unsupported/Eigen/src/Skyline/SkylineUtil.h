// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Guillaume Saupin <guillaume.saupin@cea.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef EIGEN_SKYLINEUTIL_H
#define EIGEN_SKYLINEUTIL_H

#ifdef NDEBUG
#define EIGEN_DBG_SKYLINE(X)
#else
#define EIGEN_DBG_SKYLINE(X) X
#endif

const unsigned int SkylineBit = 0x1200;
template<typename Lhs, typename Rhs, int ProductMode> class SkylineProduct;
enum AdditionalProductEvaluationMode {SkylineTimeDenseProduct, SkylineTimeSkylineProduct, DenseTimeSkylineProduct};
enum {IsSkyline = SkylineBit};


#define EIGEN_SKYLINE_INHERIT_ASSIGNMENT_OPERATOR(Derived, Op) \
template<typename OtherDerived> \
EIGEN_STRONG_INLINE Derived& operator Op(const Eigen::SkylineMatrixBase<OtherDerived>& other) \
{ \
  return Base::operator Op(other.derived()); \
} \
EIGEN_STRONG_INLINE Derived& operator Op(const Derived& other) \
{ \
  return Base::operator Op(other); \
}

#define EIGEN_SKYLINE_INHERIT_SCALAR_ASSIGNMENT_OPERATOR(Derived, Op) \
template<typename Other> \
EIGEN_STRONG_INLINE Derived& operator Op(const Other& scalar) \
{ \
  return Base::operator Op(scalar); \
}

#define EIGEN_SKYLINE_INHERIT_ASSIGNMENT_OPERATORS(Derived) \
  EIGEN_SKYLINE_INHERIT_ASSIGNMENT_OPERATOR(Derived, =) \
  EIGEN_SKYLINE_INHERIT_ASSIGNMENT_OPERATOR(Derived, +=) \
  EIGEN_SKYLINE_INHERIT_ASSIGNMENT_OPERATOR(Derived, -=) \
  EIGEN_SKYLINE_INHERIT_SCALAR_ASSIGNMENT_OPERATOR(Derived, *=) \
  EIGEN_SKYLINE_INHERIT_SCALAR_ASSIGNMENT_OPERATOR(Derived, /=)

#define _EIGEN_SKYLINE_GENERIC_PUBLIC_INTERFACE(Derived, BaseClass) \
  typedef BaseClass Base; \
  typedef typename Eigen::internal::traits<Derived>::Scalar Scalar; \
  typedef typename Eigen::NumTraits<Scalar>::Real RealScalar; \
  typedef typename Eigen::internal::traits<Derived>::StorageKind StorageKind; \
  typedef typename Eigen::internal::index<StorageKind>::type Index; \
  enum {  Flags = Eigen::internal::traits<Derived>::Flags, };

#define EIGEN_SKYLINE_GENERIC_PUBLIC_INTERFACE(Derived) \
  _EIGEN_SKYLINE_GENERIC_PUBLIC_INTERFACE(Derived, Eigen::SkylineMatrixBase<Derived>)

template<typename Derived> class SkylineMatrixBase;
template<typename _Scalar, int _Flags = 0> class SkylineMatrix;
template<typename _Scalar, int _Flags = 0> class DynamicSkylineMatrix;
template<typename _Scalar, int _Flags = 0> class SkylineVector;
template<typename _Scalar, int _Flags = 0> class MappedSkylineMatrix;

namespace internal {

template<typename Lhs, typename Rhs> struct skyline_product_mode;
template<typename Lhs, typename Rhs, int ProductMode = skyline_product_mode<Lhs,Rhs>::value> struct SkylineProductReturnType;

template<typename T> class eval<T,IsSkyline>
{
    typedef typename traits<T>::Scalar _Scalar;
    enum {
          _Flags = traits<T>::Flags
    };

  public:
    typedef SkylineMatrix<_Scalar, _Flags> type;
};

} // end namespace internal


#endif // EIGEN_SKYLINEUTIL_H
