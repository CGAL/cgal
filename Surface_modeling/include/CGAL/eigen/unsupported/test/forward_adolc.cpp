// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <g.gael@free.fr>
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

#include "main.h"
#define NUMBER_DIRECTIONS 16
#include <unsupported/Eigen/AdolcForward>

int adtl::ADOLC_numDir;

template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct TestFunc1
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  TestFunc1() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  TestFunc1(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  template<typename T>
  void operator() (const Matrix<T,InputsAtCompileTime,1>& x, Matrix<T,ValuesAtCompileTime,1>* _v) const
  {
    Matrix<T,ValuesAtCompileTime,1>& v = *_v;

    v[0] = 2 * x[0] * x[0] + x[0] * x[1];
    v[1] = 3 * x[1] * x[0] + 0.5 * x[1] * x[1];
    if(inputs()>2)
    {
      v[0] += 0.5 * x[2];
      v[1] += x[2];
    }
    if(values()>2)
    {
      v[2] = 3 * x[1] * x[0] * x[0];
    }
    if (inputs()>2 && values()>2)
      v[2] *= x[2];
  }

  void operator() (const InputType& x, ValueType* v, JacobianType* _j) const
  {
    (*this)(x, v);

    if(_j)
    {
      JacobianType& j = *_j;

      j(0,0) = 4 * x[0] + x[1];
      j(1,0) = 3 * x[1];

      j(0,1) = x[0];
      j(1,1) = 3 * x[0] + 2 * 0.5 * x[1];

      if (inputs()>2)
      {
        j(0,2) = 0.5;
        j(1,2) = 1;
      }
      if(values()>2)
      {
        j(2,0) = 3 * x[1] * 2 * x[0];
        j(2,1) = 3 * x[0] * x[0];
      }
      if (inputs()>2 && values()>2)
      {
        j(2,0) *= x[2];
        j(2,1) *= x[2];

        j(2,2) = 3 * x[1] * x[0] * x[0];
        j(2,2) = 3 * x[1] * x[0] * x[0];
      }
    }
  }
};

template<typename Func> void adolc_forward_jacobian(const Func& f)
{
    typename Func::InputType x = Func::InputType::Random(f.inputs());
    typename Func::ValueType y(f.values()), yref(f.values());
    typename Func::JacobianType j(f.values(),f.inputs()), jref(f.values(),f.inputs());

    jref.setZero();
    yref.setZero();
    f(x,&yref,&jref);
//     std::cerr << y.transpose() << "\n\n";;
//     std::cerr << j << "\n\n";;

    j.setZero();
    y.setZero();
    AdolcForwardJacobian<Func> autoj(f);
    autoj(x, &y, &j);
//     std::cerr << y.transpose() << "\n\n";;
//     std::cerr << j << "\n\n";;

    VERIFY_IS_APPROX(y, yref);
    VERIFY_IS_APPROX(j, jref);
}

void test_forward_adolc()
{
  adtl::ADOLC_numDir = NUMBER_DIRECTIONS;

  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST(( adolc_forward_jacobian(TestFunc1<double,2,2>()) ));
    CALL_SUBTEST(( adolc_forward_jacobian(TestFunc1<double,2,3>()) ));
    CALL_SUBTEST(( adolc_forward_jacobian(TestFunc1<double,3,2>()) ));
    CALL_SUBTEST(( adolc_forward_jacobian(TestFunc1<double,3,3>()) ));
    CALL_SUBTEST(( adolc_forward_jacobian(TestFunc1<double>(3,3)) ));
  }
}
