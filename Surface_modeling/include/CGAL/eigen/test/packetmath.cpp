// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2006-2008 Benoit Jacob <jacob.benoit.1@gmail.com>
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

// using namespace Eigen;

namespace Eigen {
namespace internal {
template<typename T> T negate(const T& x) { return -x; }
}
}

template<typename Scalar> bool isApproxAbs(const Scalar& a, const Scalar& b, const typename NumTraits<Scalar>::Real& refvalue)
{
  return internal::isMuchSmallerThan(a-b, refvalue);
}

template<typename Scalar> bool areApproxAbs(const Scalar* a, const Scalar* b, int size, const typename NumTraits<Scalar>::Real& refvalue)
{
  for (int i=0; i<size; ++i)
  {
    if (!isApproxAbs(a[i],b[i],refvalue))
    {
      std::cout << "a[" << i << "]: " << a[i] << " != b[" << i << "]: " << b[i] << std::endl;
      return false;
    }
  }
  return true;
}

template<typename Scalar> bool areApprox(const Scalar* a, const Scalar* b, int size)
{
  for (int i=0; i<size; ++i)
  {
    if (!internal::isApprox(a[i],b[i]))
    {
      std::cout << "a[" << i << "]: " << a[i] << " != b[" << i << "]: " << b[i] << std::endl;
      return false;
    }
  }
  return true;
}


#define CHECK_CWISE2(REFOP, POP) { \
  for (int i=0; i<PacketSize; ++i) \
    ref[i] = REFOP(data1[i], data1[i+PacketSize]); \
  internal::pstore(data2, POP(internal::pload<Packet>(data1), internal::pload<Packet>(data1+PacketSize))); \
  VERIFY(areApprox(ref, data2, PacketSize) && #POP); \
}

#define CHECK_CWISE1(REFOP, POP) { \
  for (int i=0; i<PacketSize; ++i) \
    ref[i] = REFOP(data1[i]); \
  internal::pstore(data2, POP(internal::pload<Packet>(data1))); \
  VERIFY(areApprox(ref, data2, PacketSize) && #POP); \
}

template<bool Cond,typename Packet>
struct packet_helper
{
  template<typename T>
  inline Packet load(const T* from) const { return internal::pload<Packet>(from); }

  template<typename T>
  inline void store(T* to, const Packet& x) const { internal::pstore(to,x); }
};

template<typename Packet>
struct packet_helper<false,Packet>
{
  template<typename T>
  inline T load(const T* from) const { return *from; }

  template<typename T>
  inline void store(T* to, const T& x) const { *to = x; }
};

#define CHECK_CWISE1_IF(COND, REFOP, POP) if(COND) { \
  packet_helper<COND,Packet> h; \
  for (int i=0; i<PacketSize; ++i) \
    ref[i] = REFOP(data1[i]); \
  h.store(data2, POP(h.load(data1))); \
  VERIFY(areApprox(ref, data2, PacketSize) && #POP); \
}

#define REF_ADD(a,b) ((a)+(b))
#define REF_SUB(a,b) ((a)-(b))
#define REF_MUL(a,b) ((a)*(b))
#define REF_DIV(a,b) ((a)/(b))

template<typename Scalar> void packetmath()
{
  typedef typename internal::packet_traits<Scalar>::type Packet;
  const int PacketSize = internal::packet_traits<Scalar>::size;
  typedef typename NumTraits<Scalar>::Real RealScalar;

  const int size = PacketSize*4;
  EIGEN_ALIGN16 Scalar data1[internal::packet_traits<Scalar>::size*4];
  EIGEN_ALIGN16 Scalar data2[internal::packet_traits<Scalar>::size*4];
  EIGEN_ALIGN16 Packet packets[PacketSize*2];
  EIGEN_ALIGN16 Scalar ref[internal::packet_traits<Scalar>::size*4];
  RealScalar refvalue = 0;
  for (int i=0; i<size; ++i)
  {
    data1[i] = internal::random<Scalar>()/RealScalar(PacketSize);
    data2[i] = internal::random<Scalar>()/RealScalar(PacketSize);
    refvalue = std::max(refvalue,internal::abs(data1[i]));
  }

  internal::pstore(data2, internal::pload<Packet>(data1));
  VERIFY(areApprox(data1, data2, PacketSize) && "aligned load/store");

  for (int offset=0; offset<PacketSize; ++offset)
  {
    internal::pstore(data2, internal::ploadu<Packet>(data1+offset));
    VERIFY(areApprox(data1+offset, data2, PacketSize) && "internal::ploadu");
  }

  for (int offset=0; offset<PacketSize; ++offset)
  {
    internal::pstoreu(data2+offset, internal::pload<Packet>(data1));
    VERIFY(areApprox(data1, data2+offset, PacketSize) && "internal::pstoreu");
  }

  for (int offset=0; offset<PacketSize; ++offset)
  {
    packets[0] = internal::pload<Packet>(data1);
    packets[1] = internal::pload<Packet>(data1+PacketSize);
         if (offset==0) internal::palign<0>(packets[0], packets[1]);
    else if (offset==1) internal::palign<1>(packets[0], packets[1]);
    else if (offset==2) internal::palign<2>(packets[0], packets[1]);
    else if (offset==3) internal::palign<3>(packets[0], packets[1]);
    internal::pstore(data2, packets[0]);

    for (int i=0; i<PacketSize; ++i)
      ref[i] = data1[i+offset];

    typedef Matrix<Scalar, PacketSize, 1> Vector;
    VERIFY(areApprox(ref, data2, PacketSize) && "internal::palign");
  }

  CHECK_CWISE2(REF_ADD,  internal::padd);
  CHECK_CWISE2(REF_SUB,  internal::psub);
  CHECK_CWISE2(REF_MUL,  internal::pmul);
  #ifndef EIGEN_VECTORIZE_ALTIVEC
  if (!internal::is_same<Scalar,int>::value)
    CHECK_CWISE2(REF_DIV,  internal::pdiv);
  #endif
  CHECK_CWISE1(internal::negate, internal::pnegate);
  CHECK_CWISE1(internal::conj, internal::pconj);

  for(int offset=0;offset<3;++offset)
  {
    for (int i=0; i<PacketSize; ++i)
      ref[i] = data1[offset];
    internal::pstore(data2, internal::pset1<Packet>(data1[offset]));
    VERIFY(areApprox(ref, data2, PacketSize) && "internal::pset1");
  }

  VERIFY(internal::isApprox(data1[0], internal::pfirst(internal::pload<Packet>(data1))) && "internal::pfirst");
  
  if(PacketSize>1)
  {
    for(int offset=0;offset<4;++offset)
    {
      for(int i=0;i<PacketSize/2;++i)
        ref[2*i+0] = ref[2*i+1] = data1[offset+i];
      internal::pstore(data2,internal::ploaddup<Packet>(data1+offset));
      VERIFY(areApprox(ref, data2, PacketSize) && "ploaddup");
    }
  }

  ref[0] = 0;
  for (int i=0; i<PacketSize; ++i)
    ref[0] += data1[i];
  VERIFY(isApproxAbs(ref[0], internal::predux(internal::pload<Packet>(data1)), refvalue) && "internal::predux");

  ref[0] = 1;
  for (int i=0; i<PacketSize; ++i)
    ref[0] *= data1[i];
  VERIFY(internal::isApprox(ref[0], internal::predux_mul(internal::pload<Packet>(data1))) && "internal::predux_mul");

  for (int j=0; j<PacketSize; ++j)
  {
    ref[j] = 0;
    for (int i=0; i<PacketSize; ++i)
      ref[j] += data1[i+j*PacketSize];
    packets[j] = internal::pload<Packet>(data1+j*PacketSize);
  }
  internal::pstore(data2, internal::preduxp(packets));
  VERIFY(areApproxAbs(ref, data2, PacketSize, refvalue) && "internal::preduxp");

  for (int i=0; i<PacketSize; ++i)
    ref[i] = data1[PacketSize-i-1];
  internal::pstore(data2, internal::preverse(internal::pload<Packet>(data1)));
  VERIFY(areApprox(ref, data2, PacketSize) && "internal::preverse");
}

template<typename Scalar> void packetmath_real()
{
  typedef typename internal::packet_traits<Scalar>::type Packet;
  const int PacketSize = internal::packet_traits<Scalar>::size;

  const int size = PacketSize*4;
  EIGEN_ALIGN16 Scalar data1[internal::packet_traits<Scalar>::size*4];
  EIGEN_ALIGN16 Scalar data2[internal::packet_traits<Scalar>::size*4];
  EIGEN_ALIGN16 Scalar ref[internal::packet_traits<Scalar>::size*4];

  for (int i=0; i<size; ++i)
  {
    data1[i] = internal::random<Scalar>(-1e3,1e3);
    data2[i] = internal::random<Scalar>(-1e3,1e3);
  }
  CHECK_CWISE1_IF(internal::packet_traits<Scalar>::HasSin, internal::sin, internal::psin);
  CHECK_CWISE1_IF(internal::packet_traits<Scalar>::HasCos, internal::cos, internal::pcos);
  CHECK_CWISE1_IF(internal::packet_traits<Scalar>::HasTan, internal::tan, internal::ptan);
  
  for (int i=0; i<size; ++i)
  {
    data1[i] = internal::random<Scalar>(-1,1);
    data2[i] = internal::random<Scalar>(-1,1);
  }
  CHECK_CWISE1_IF(internal::packet_traits<Scalar>::HasASin, internal::asin, internal::pasin);
  CHECK_CWISE1_IF(internal::packet_traits<Scalar>::HasACos, internal::acos, internal::pacos);

  for (int i=0; i<size; ++i)
  {
    data1[i] = internal::random<Scalar>(-87,88);
    data2[i] = internal::random<Scalar>(-87,88);
  }
  CHECK_CWISE1_IF(internal::packet_traits<Scalar>::HasExp, internal::exp, internal::pexp);

  for (int i=0; i<size; ++i)
  {
    data1[i] = internal::random<Scalar>(0,1e6);
    data2[i] = internal::random<Scalar>(0,1e6);
  }
  CHECK_CWISE1_IF(internal::packet_traits<Scalar>::HasLog, internal::log, internal::plog);
  CHECK_CWISE1_IF(internal::packet_traits<Scalar>::HasSqrt, internal::sqrt, internal::psqrt);

  ref[0] = data1[0];
  for (int i=0; i<PacketSize; ++i)
    ref[0] = std::min(ref[0],data1[i]);
  VERIFY(internal::isApprox(ref[0], internal::predux_min(internal::pload<Packet>(data1))) && "internal::predux_min");

  CHECK_CWISE2(std::min, internal::pmin);
  CHECK_CWISE2(std::max, internal::pmax);
  CHECK_CWISE1(internal::abs, internal::pabs);

  ref[0] = data1[0];
  for (int i=0; i<PacketSize; ++i)
    ref[0] = std::max(ref[0],data1[i]);
  VERIFY(internal::isApprox(ref[0], internal::predux_max(internal::pload<Packet>(data1))) && "internal::predux_max");
}

template<typename Scalar,bool ConjLhs,bool ConjRhs> void test_conj_helper(Scalar* data1, Scalar* data2, Scalar* ref, Scalar* pval)
{
  typedef typename internal::packet_traits<Scalar>::type Packet;
  const int PacketSize = internal::packet_traits<Scalar>::size;
  
  internal::conj_if<ConjLhs> cj0;
  internal::conj_if<ConjRhs> cj1;
  internal::conj_helper<Scalar,Scalar,ConjLhs,ConjRhs> cj;
  internal::conj_helper<Packet,Packet,ConjLhs,ConjRhs> pcj;
  
  for(int i=0;i<PacketSize;++i)
  {
    ref[i] = cj0(data1[i]) * cj1(data2[i]);
    VERIFY(internal::isApprox(ref[i], cj.pmul(data1[i],data2[i])) && "conj_helper pmul");
  }
  internal::pstore(pval,pcj.pmul(internal::pload<Packet>(data1),internal::pload<Packet>(data2)));
  VERIFY(areApprox(ref, pval, PacketSize) && "conj_helper pmul");
  
  for(int i=0;i<PacketSize;++i)
  {
    Scalar tmp = ref[i];
    ref[i] += cj0(data1[i]) * cj1(data2[i]);
    VERIFY(internal::isApprox(ref[i], cj.pmadd(data1[i],data2[i],tmp)) && "conj_helper pmadd");
  }
  internal::pstore(pval,pcj.pmadd(internal::pload<Packet>(data1),internal::pload<Packet>(data2),internal::pload<Packet>(pval)));
  VERIFY(areApprox(ref, pval, PacketSize) && "conj_helper pmadd");
}

template<typename Scalar> void packetmath_complex()
{
  typedef typename internal::packet_traits<Scalar>::type Packet;
  const int PacketSize = internal::packet_traits<Scalar>::size;

  const int size = PacketSize*4;
  EIGEN_ALIGN16 Scalar data1[PacketSize*4];
  EIGEN_ALIGN16 Scalar data2[PacketSize*4];
  EIGEN_ALIGN16 Scalar ref[PacketSize*4];
  EIGEN_ALIGN16 Scalar pval[PacketSize*4];

  for (int i=0; i<size; ++i)
  {
    data1[i] = internal::random<Scalar>() * Scalar(1e2);
    data2[i] = internal::random<Scalar>() * Scalar(1e2);
  }
  
  test_conj_helper<Scalar,false,false> (data1,data2,ref,pval);
  test_conj_helper<Scalar,false,true>  (data1,data2,ref,pval);
  test_conj_helper<Scalar,true,false>  (data1,data2,ref,pval);
  test_conj_helper<Scalar,true,true>   (data1,data2,ref,pval);
  
  {
    for(int i=0;i<PacketSize;++i)
      ref[i] = Scalar(std::imag(data1[i]),std::real(data1[i]));
    internal::pstore(pval,internal::pcplxflip(internal::pload<Packet>(data1)));
    VERIFY(areApprox(ref, pval, PacketSize) && "pcplxflip");
  }
  
  
}

void test_packetmath()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( packetmath<float>() );
    CALL_SUBTEST_2( packetmath<double>() );
    CALL_SUBTEST_3( packetmath<int>() );
    CALL_SUBTEST_1( packetmath<std::complex<float> >() );
    CALL_SUBTEST_2( packetmath<std::complex<double> >() );

    CALL_SUBTEST_1( packetmath_real<float>() );
    CALL_SUBTEST_2( packetmath_real<double>() );

    CALL_SUBTEST_1( packetmath_complex<std::complex<float> >() );
    CALL_SUBTEST_2( packetmath_complex<std::complex<double> >() );
  }
}
