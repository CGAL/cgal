// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
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

template<typename Scalar> bool areApprox(const Scalar* a, const Scalar* b, int size)
{
  for (int i=0; i<size; ++i)
    if (!ei_isApprox(a[i],b[i])) return false;
  return true;
}

#define CHECK_CWISE(REFOP, POP) { \
  for (int i=0; i<PacketSize; ++i) \
    ref[i] = REFOP(data1[i], data1[i+PacketSize]); \
  ei_pstore(data2, POP(ei_pload(data1), ei_pload(data1+PacketSize))); \
  VERIFY(areApprox(ref, data2, PacketSize) && #POP); \
}

#define REF_ADD(a,b) ((a)+(b))
#define REF_SUB(a,b) ((a)-(b))
#define REF_MUL(a,b) ((a)*(b))
#define REF_DIV(a,b) ((a)/(b))

namespace std {

template<> const complex<float>& min(const complex<float>& a, const complex<float>& b)
{ return a.real() < b.real() ? a : b; }

template<> const complex<float>& max(const complex<float>& a, const complex<float>& b)
{ return a.real() < b.real() ? b : a; }

}

template<typename Scalar> void packetmath()
{
  typedef typename ei_packet_traits<Scalar>::type Packet;
  const int PacketSize = ei_packet_traits<Scalar>::size;

  const int size = PacketSize*4;
  EIGEN_ALIGN_128 Scalar data1[ei_packet_traits<Scalar>::size*4];
  EIGEN_ALIGN_128 Scalar data2[ei_packet_traits<Scalar>::size*4];
  EIGEN_ALIGN_128 Packet packets[PacketSize*2];
  EIGEN_ALIGN_128 Scalar ref[ei_packet_traits<Scalar>::size*4];
  for (int i=0; i<size; ++i)
  {
    data1[i] = ei_random<Scalar>();
    data2[i] = ei_random<Scalar>();
  }

  ei_pstore(data2, ei_pload(data1));
  VERIFY(areApprox(data1, data2, PacketSize) && "aligned load/store");

  for (int offset=0; offset<PacketSize; ++offset)
  {
    ei_pstore(data2, ei_ploadu(data1+offset));
    VERIFY(areApprox(data1+offset, data2, PacketSize) && "ei_ploadu");
  }

  for (int offset=0; offset<PacketSize; ++offset)
  {
    ei_pstoreu(data2+offset, ei_pload(data1));
    VERIFY(areApprox(data1, data2+offset, PacketSize) && "ei_pstoreu");
  }

  for (int offset=0; offset<PacketSize; ++offset)
  {
    packets[0] = ei_pload(data1);
    packets[1] = ei_pload(data1+PacketSize);
         if (offset==0) ei_palign<0>(packets[0], packets[1]);
    else if (offset==1) ei_palign<1>(packets[0], packets[1]);
    else if (offset==2) ei_palign<2>(packets[0], packets[1]);
    else if (offset==3) ei_palign<3>(packets[0], packets[1]);
    ei_pstore(data2, packets[0]);

    for (int i=0; i<PacketSize; ++i)
      ref[i] = data1[i+offset];

    typedef Matrix<Scalar, PacketSize, 1> Vector;
    VERIFY(areApprox(ref, data2, PacketSize) && "ei_palign");
  }

  CHECK_CWISE(REF_ADD,  ei_padd);
  CHECK_CWISE(REF_SUB,  ei_psub);
  CHECK_CWISE(REF_MUL,  ei_pmul);
  #ifndef EIGEN_VECTORIZE_ALTIVEC
  if (!ei_is_same_type<Scalar,int>::ret)
    CHECK_CWISE(REF_DIV,  ei_pdiv);
  #endif
  CHECK_CWISE(std::min, ei_pmin);
  CHECK_CWISE(std::max, ei_pmax);

  for (int i=0; i<PacketSize; ++i)
    ref[i] = data1[0];
  ei_pstore(data2, ei_pset1(data1[0]));
  VERIFY(areApprox(ref, data2, PacketSize) && "ei_pset1");

  VERIFY(ei_isApprox(data1[0], ei_pfirst(ei_pload(data1))) && "ei_pfirst");

  ref[0] = 0;
  for (int i=0; i<PacketSize; ++i)
    ref[0] += data1[i];
  VERIFY(ei_isApprox(ref[0], ei_predux(ei_pload(data1))) && "ei_predux");

  for (int j=0; j<PacketSize; ++j)
  {
    ref[j] = 0;
    for (int i=0; i<PacketSize; ++i)
      ref[j] += data1[i+j*PacketSize];
    packets[j] = ei_pload(data1+j*PacketSize);
  }
  ei_pstore(data2, ei_preduxp(packets));
  VERIFY(areApprox(ref, data2, PacketSize) && "ei_preduxp");
}

void test_eigen2_packetmath()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( packetmath<float>() );
    CALL_SUBTEST_2( packetmath<double>() );
    CALL_SUBTEST_3( packetmath<int>() );
    CALL_SUBTEST_4( packetmath<std::complex<float> >() );
  }
}
