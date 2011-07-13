// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Gael Guennebaud <gael.guennebaud@inria.fr>
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

#ifndef EIGEN_AUTODIFF_JACOBIAN_H
#define EIGEN_AUTODIFF_JACOBIAN_H

namespace Eigen
{

template<typename Functor> class AutoDiffJacobian : public Functor
{
public:
  AutoDiffJacobian() : Functor() {}
  AutoDiffJacobian(const Functor& f) : Functor(f) {}

  // forward constructors
  template<typename T0>
  AutoDiffJacobian(const T0& a0) : Functor(a0) {}
  template<typename T0, typename T1>
  AutoDiffJacobian(const T0& a0, const T1& a1) : Functor(a0, a1) {}
  template<typename T0, typename T1, typename T2>
  AutoDiffJacobian(const T0& a0, const T1& a1, const T1& a2) : Functor(a0, a1, a2) {}

  enum {
    InputsAtCompileTime = Functor::InputsAtCompileTime,
    ValuesAtCompileTime = Functor::ValuesAtCompileTime
  };

  typedef typename Functor::InputType InputType;
  typedef typename Functor::ValueType ValueType;
  typedef typename Functor::JacobianType JacobianType;
  typedef typename JacobianType::Scalar Scalar;
  typedef typename JacobianType::Index Index;

  typedef Matrix<Scalar,InputsAtCompileTime,1> DerivativeType;
  typedef AutoDiffScalar<DerivativeType> ActiveScalar;


  typedef Matrix<ActiveScalar, InputsAtCompileTime, 1> ActiveInput;
  typedef Matrix<ActiveScalar, ValuesAtCompileTime, 1> ActiveValue;

  void operator() (const InputType& x, ValueType* v, JacobianType* _jac=0) const
  {
    eigen_assert(v!=0);
    if (!_jac)
    {
      Functor::operator()(x, v);
      return;
    }

    JacobianType& jac = *_jac;

    ActiveInput ax = x.template cast<ActiveScalar>();
    ActiveValue av(jac.rows());

    if(InputsAtCompileTime==Dynamic)
      for (Index j=0; j<jac.rows(); j++)
        av[j].derivatives().resize(this->inputs());

    for (Index i=0; i<jac.cols(); i++)
      ax[i].derivatives() = DerivativeType::Unit(this->inputs(),i);

    Functor::operator()(ax, &av);

    for (Index i=0; i<jac.rows(); i++)
    {
      (*v)[i] = av[i].value();
      jac.row(i) = av[i].derivatives();
    }
  }
protected:

};

}

#endif // EIGEN_AUTODIFF_JACOBIAN_H
