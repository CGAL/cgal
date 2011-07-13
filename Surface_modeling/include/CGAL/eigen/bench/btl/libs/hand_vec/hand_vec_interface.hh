//=====================================================
// File   :  hand_vec_interface.hh
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
//=====================================================
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
#ifndef HAND_VEC_INTERFACE_HH
#define HAND_VEC_INTERFACE_HH

#include <Eigen/Core>
#include "f77_interface.hh"

using namespace Eigen;

template<class real>
class hand_vec_interface : public f77_interface_base<real> {

public :

  typedef typename internal::packet_traits<real>::type Packet;
  static const int PacketSize = internal::packet_traits<real>::size;

  typedef typename f77_interface_base<real>::stl_matrix stl_matrix;
  typedef typename f77_interface_base<real>::stl_vector stl_vector;
  typedef typename f77_interface_base<real>::gene_matrix gene_matrix;
  typedef typename f77_interface_base<real>::gene_vector gene_vector;

  static void free_matrix(gene_matrix & A, int N){
    internal::aligned_free(A);
  }

  static void free_vector(gene_vector & B){
    internal::aligned_free(B);
  }

  static inline void matrix_from_stl(gene_matrix & A, stl_matrix & A_stl){
    int N = A_stl.size();
    A = (real*)internal::aligned_malloc(N*N*sizeof(real));
    for (int j=0;j<N;j++)
      for (int i=0;i<N;i++)
        A[i+N*j] = A_stl[j][i];
  }

  static inline void vector_from_stl(gene_vector & B, stl_vector & B_stl){
    int N = B_stl.size();
    B = (real*)internal::aligned_malloc(N*sizeof(real));
    for (int i=0;i<N;i++)
      B[i] = B_stl[i];
  }

  static inline std::string name() {
    #ifdef PEELING
    return "hand_vectorized_peeling";
    #else
    return "hand_vectorized";
    #endif
  }

  static inline void matrix_vector_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
  {
    asm("#begin matrix_vector_product");
    int AN = (N/PacketSize)*PacketSize;
    int ANP = (AN/(2*PacketSize))*2*PacketSize;
    int bound = (N/4)*4;
    for (int i=0;i<N;i++)
      X[i] = 0;

    for (int i=0;i<bound;i+=4)
    {
      register real* __restrict__ A0 = A + i*N;
      register real* __restrict__ A1 = A + (i+1)*N;
      register real* __restrict__ A2 = A + (i+2)*N;
      register real* __restrict__ A3 = A + (i+3)*N;

      Packet ptmp0 = internal::pset1(B[i]);
      Packet ptmp1 = internal::pset1(B[i+1]);
      Packet ptmp2 = internal::pset1(B[i+2]);
      Packet ptmp3 = internal::pset1(B[i+3]);
//       register Packet ptmp0, ptmp1, ptmp2, ptmp3;
//       asm(
//
//           "movss     (%[B],%[j],4), %[ptmp0]  \n\t"
//           "shufps   $0,%[ptmp0],%[ptmp0] \n\t"
//           "movss    4(%[B],%[j],4), %[ptmp1]  \n\t"
//           "shufps   $0,%[ptmp1],%[ptmp1] \n\t"
//           "movss    8(%[B],%[j],4), %[ptmp2]  \n\t"
//           "shufps   $0,%[ptmp2],%[ptmp2] \n\t"
//           "movss   12(%[B],%[j],4), %[ptmp3]  \n\t"
//           "shufps   $0,%[ptmp3],%[ptmp3] \n\t"
//           : [ptmp0] "=x" (ptmp0),
//             [ptmp1] "=x" (ptmp1),
//             [ptmp2] "=x" (ptmp2),
//             [ptmp3] "=x" (ptmp3)
//           : [B] "r" (B),
//             [j] "r" (size_t(i))
//           : );

      if (AN>0)
      {
//         for (size_t j = 0;j<ANP;j+=8)
//         {
//           asm(
//
//           "movaps     (%[A0],%[j],4), %%xmm8  \n\t"
//           "movaps   16(%[A0],%[j],4), %%xmm12 \n\t"
//           "movups     (%[A3],%[j],4), %%xmm11 \n\t"
//           "movups   16(%[A3],%[j],4), %%xmm15 \n\t"
//           "movups     (%[A2],%[j],4), %%xmm10 \n\t"
//           "movups   16(%[A2],%[j],4), %%xmm14 \n\t"
//           "movups     (%[A1],%[j],4), %%xmm9  \n\t"
//           "movups   16(%[A1],%[j],4), %%xmm13 \n\t"
//
//           "mulps %[ptmp0], %%xmm8  \n\t"
//           "addps (%[res0],%[j],4), %%xmm8  \n\t"
//           "mulps %[ptmp3], %%xmm11 \n\t"
//           "addps %%xmm11, %%xmm8  \n\t"
//           "mulps %[ptmp2], %%xmm10 \n\t"
//           "addps %%xmm10, %%xmm8  \n\t"
//           "mulps %[ptmp1], %%xmm9  \n\t"
//           "addps %%xmm9, %%xmm8   \n\t"
//           "movaps %%xmm8, (%[res0],%[j],4)  \n\t"
//
//           "mulps %[ptmp0], %%xmm12 \n\t"
//           "addps 16(%[res0],%[j],4), %%xmm12  \n\t"
//           "mulps %[ptmp3], %%xmm15 \n\t"
//           "addps %%xmm15, %%xmm12  \n\t"
//           "mulps %[ptmp2], %%xmm14 \n\t"
//           "addps %%xmm14, %%xmm12  \n\t"
//           "mulps %[ptmp1], %%xmm13 \n\t"
//           "addps %%xmm13, %%xmm12  \n\t"
//           "movaps %%xmm12, 16(%[res0],%[j],4) \n\t"
//           :
//           : [res0] "r" (X), [j] "r" (j),[A0] "r" (A0),
//             [A1] "r" (A1),
//             [A2] "r" (A2),
//             [A3] "r" (A3),
//             [ptmp0] "x" (ptmp0),
//             [ptmp1] "x" (ptmp1),
//             [ptmp2] "x" (ptmp2),
//             [ptmp3] "x" (ptmp3)
//           : "%xmm8", "%xmm9", "%xmm10", "%xmm11", "%xmm12", "%xmm13", "%xmm14", "%xmm15", "%r14");
//         }
          register Packet A00;
          register Packet A01;
          register Packet A02;
          register Packet A03;
          register Packet A10;
          register Packet A11;
          register Packet A12;
          register Packet A13;
          for (int j = 0;j<ANP;j+=2*PacketSize)
          {
//             A00 = internal::pload(&A0[j]);
//             A01 = internal::ploadu(&A1[j]);
//             A02 = internal::ploadu(&A2[j]);
//             A03 = internal::ploadu(&A3[j]);
//             A10 = internal::pload(&A0[j+PacketSize]);
//             A11 = internal::ploadu(&A1[j+PacketSize]);
//             A12 = internal::ploadu(&A2[j+PacketSize]);
//             A13 = internal::ploadu(&A3[j+PacketSize]);
//
//             A00 = internal::pmul(ptmp0, A00);
//             A01 = internal::pmul(ptmp1, A01);
//             A02 = internal::pmul(ptmp2, A02);
//             A03 = internal::pmul(ptmp3, A03);
//             A10 = internal::pmul(ptmp0, A10);
//             A11 = internal::pmul(ptmp1, A11);
//             A12 = internal::pmul(ptmp2, A12);
//             A13 = internal::pmul(ptmp3, A13);
//
//             A00 = internal::padd(A00,A01);
//             A02 = internal::padd(A02,A03);
//             A00 = internal::padd(A00,internal::pload(&X[j]));
//             A00 = internal::padd(A00,A02);
//             internal::pstore(&X[j],A00);
//
//             A10 = internal::padd(A10,A11);
//             A12 = internal::padd(A12,A13);
//             A10 = internal::padd(A10,internal::pload(&X[j+PacketSize]));
//             A10 = internal::padd(A10,A12);
//             internal::pstore(&X[j+PacketSize],A10);

            internal::pstore(&X[j],
              internal::padd(internal::pload(&X[j]),
                internal::padd(
                  internal::padd(internal::pmul(ptmp0,internal::pload(&A0[j])),internal::pmul(ptmp1,internal::ploadu(&A1[j]))),
                  internal::padd(internal::pmul(ptmp2,internal::ploadu(&A2[j])),internal::pmul(ptmp3,internal::ploadu(&A3[j]))) )));

            internal::pstore(&X[j+PacketSize],
              internal::padd(internal::pload(&X[j+PacketSize]),
                internal::padd(
                  internal::padd(internal::pmul(ptmp0,internal::pload(&A0[j+PacketSize])),internal::pmul(ptmp1,internal::ploadu(&A1[j+PacketSize]))),
                  internal::padd(internal::pmul(ptmp2,internal::ploadu(&A2[j+PacketSize])),internal::pmul(ptmp3,internal::ploadu(&A3[j+PacketSize]))) )));
          }
          for (int j = ANP;j<AN;j+=PacketSize)
            internal::pstore(&X[j],
              internal::padd(internal::pload(&X[j]),
                internal::padd(
                  internal::padd(internal::pmul(ptmp0,internal::pload(&A0[j])),internal::pmul(ptmp1,internal::ploadu(&A1[j]))),
                  internal::padd(internal::pmul(ptmp2,internal::ploadu(&A2[j])),internal::pmul(ptmp3,internal::ploadu(&A3[j]))) )));
      }
      // process remaining scalars
      for (int j=AN;j<N;j++)
        X[j] += internal::pfirst(ptmp0) * A0[j] + internal::pfirst(ptmp1) * A1[j] + internal::pfirst(ptmp2) * A2[j] + internal::pfirst(ptmp3) * A3[j];
    }
    for (int i=bound;i<N;i++)
    {
      real tmp0 = B[i];
      Packet ptmp0 = internal::pset1(tmp0);
      int iN0 = i*N;
      if (AN>0)
      {
        bool aligned0 = (iN0 % PacketSize) == 0;
        if (aligned0)
          for (int j = 0;j<AN;j+=PacketSize)
            internal::pstore(&X[j], internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+iN0])),internal::pload(&X[j])));
        else
          for (int j = 0;j<AN;j+=PacketSize)
            internal::pstore(&X[j], internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+iN0])),internal::pload(&X[j])));
      }
      // process remaining scalars
      for (int j=AN;j<N;j++)
        X[j] += tmp0 * A[j+iN0];
    }
    asm("#end matrix_vector_product");
  }
  
  static inline void symv(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
  {
    
//     int AN = (N/PacketSize)*PacketSize;
//     int ANP = (AN/(2*PacketSize))*2*PacketSize;
//     int bound = (N/4)*4;
    for (int i=0;i<N;i++)
      X[i] = 0;
    
    int bound = std::max(0,N-8) & 0xfffffffE;

    for (int j=0;j<bound;j+=2)
    {
      register real* __restrict__ A0 = A + j*N;
      register real* __restrict__ A1 = A + (j+1)*N;
      
      real t0 = B[j];
      Packet ptmp0 = internal::pset1(t0);
      real t1 = B[j+1];
      Packet ptmp1 = internal::pset1(t1);
      
      real t2 = 0;
      Packet ptmp2 = internal::pset1(t2);
      real t3 = 0;
      Packet ptmp3 = internal::pset1(t3);
      
      int starti = j+2;
      int alignedEnd = starti;
      int alignedStart = (starti) + internal::first_aligned(&X[starti], N-starti);
      alignedEnd = alignedStart + ((N-alignedStart)/(PacketSize))*(PacketSize);

      X[j]   += t0 * A0[j];
      X[j+1] += t1 * A1[j];
      
      X[j+1] += t0 * A0[j+1];
      t2 += A0[j+1] * B[j+1];
      
//       alignedStart = alignedEnd;
      for (int i=starti; i<alignedStart; ++i) {
        X[i] += t0 * A0[i] + t1 * A1[i];
        t2 += A0[i] * B[i];
        t3 += A1[i] * B[i];
      }
      asm("#begin symv");
      for (size_t i=alignedStart; i<alignedEnd; i+=PacketSize) {
        Packet A0i = internal::ploadu(&A0[i]);
        Packet A1i = internal::ploadu(&A1[i]);
//         Packet A0i1 = internal::ploadu(&A0[i+PacketSize]);
        Packet Xi = internal::pload(&X[i]);
        Packet Bi = internal::pload/*u*/(&B[i]);
//         Packet Xi1 = internal::pload(&X[i+PacketSize]);
//         Packet Bi1 = internal::pload/*u*/(&B[i+PacketSize]);
        Xi = internal::padd(internal::padd(Xi, internal::pmul(ptmp0, A0i)), internal::pmul(ptmp1, A1i));
        ptmp2 = internal::padd(ptmp2, internal::pmul(A0i, Bi));
        ptmp3 = internal::padd(ptmp3, internal::pmul(A1i, Bi));
//         Xi1 = internal::padd(Xi1, internal::pmul(ptmp1, A0i1));
//         ptmp2 = internal::padd(ptmp2, internal::pmul(A0i1, Bi1));
//         
        internal::pstore(&X[i],Xi);
//         internal::pstore(&X[i+PacketSize],Xi1);
//         asm(
//           "prefetchnta   64(%[A0],%[i],4)   \n\t"
//           //"movups     (%[A0],%[i],4), %%xmm8  \n\t"
//           "movsd       (%[A0],%[i],4), %%xmm8  \n\t"
//           "movhps     8(%[A0],%[i],4), %%xmm8  \n\t"
// //           "movups   16(%[A0],%[i],4), %%xmm9  \n\t"
// //           "movups   64(%[A0],%[i],4), %%xmm15  \n\t"
//           "movaps     (%[B], %[i],4), %%xmm12 \n\t"
// //           "movaps   16(%[B], %[i],4), %%xmm13 \n\t"
//           "movaps     (%[X], %[i],4), %%xmm10 \n\t"
// //           "movaps   16(%[X], %[i],4), %%xmm11 \n\t"
//           
//           "mulps %%xmm8, %%xmm12  \n\t"
// //           "mulps %%xmm9, %%xmm13  \n\t"
//           
//           "mulps %[ptmp1], %%xmm8  \n\t"
//           "addps %%xmm12, %[ptmp2]  \n\t"
//           "addps %%xmm8, %%xmm10  \n\t"
//           
//           
//           
//           
// //           "mulps %[ptmp1], %%xmm9  \n\t"
//           
// //           "addps %%xmm9, %%xmm11  \n\t"
// //           "addps %%xmm13, %[ptmp2]  \n\t"
//           
//           "movaps %%xmm10,   (%[X],%[i],4) \n\t"
// //           "movaps %%xmm11, 16(%[X],%[i],4) \n\t"
//           : 
//           : [X] "r" (X), [i] "r" (i), [A0] "r" (A0),
//             [B] "r" (B),
//             [ptmp1] "x" (ptmp1),
//             [ptmp2] "x" (ptmp2)
//           : "%xmm8", "%xmm9", "%xmm10", "%xmm11", "%xmm12", "%xmm13", "%xmm15");
      }
      asm("#end symv");
      for (int i=alignedEnd; i<N; i++) {
        X[i] += t0 * A0[i] + t1 * A1[i];
        t2 += A0[i] * B[i];
        t3 += A1[i] * B[i];
      }
      
      
      X[j]   += t2 + internal::predux(ptmp2);
      X[j+1] += t3 + internal::predux(ptmp3);
    }
    for (int j=bound;j<N;j++)
    {
      register real* __restrict__ A0 = A + j*N;
      
      real t1 = B[j];
      real t2 = 0;
      X[j] += t1 * A0[j];
      for (int i=j+1; i<N; i+=PacketSize) {
        X[i] += t1 * A0[i];
        t2 += A0[i] * B[i];
      }
      X[j] += t2;
    }
    
  }

//   static inline void matrix_vector_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
//   {
//     asm("#begin matrix_vector_product");
//     int AN = (N/PacketSize)*PacketSize;
//     int ANP = (AN/(2*PacketSize))*2*PacketSize;
//     int bound = (N/4)*4;
//     for (int i=0;i<N;i++)
//       X[i] = 0;
//
//     for (int i=0;i<bound;i+=4)
//     {
//       real tmp0 = B[i];
//       Packet ptmp0 = internal::pset1(tmp0);
//       real tmp1 = B[i+1];
//       Packet ptmp1 = internal::pset1(tmp1);
//       real tmp2 = B[i+2];
//       Packet ptmp2 = internal::pset1(tmp2);
//       real tmp3 = B[i+3];
//       Packet ptmp3 = internal::pset1(tmp3);
//       int iN0 = i*N;
//       int iN1 = (i+1)*N;
//       int iN2 = (i+2)*N;
//       int iN3 = (i+3)*N;
//       if (AN>0)
//       {
// //         int aligned0 = (iN0 % PacketSize);
//         int aligned1 = (iN1 % PacketSize);
//
//         if (aligned1==0)
//         {
//           for (int j = 0;j<AN;j+=PacketSize)
//           {
//             internal::pstore(&X[j],
//               internal::padd(internal::pload(&X[j]),
//                 internal::padd(
//                   internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+iN0])),internal::pmul(ptmp1,internal::pload(&A[j+iN1]))),
//                   internal::padd(internal::pmul(ptmp2,internal::pload(&A[j+iN2])),internal::pmul(ptmp3,internal::pload(&A[j+iN3]))) )));
//           }
//         }
//         else if (aligned1==2)
//         {
//           for (int j = 0;j<AN;j+=PacketSize)
//           {
//             internal::pstore(&X[j],
//               internal::padd(internal::pload(&X[j]),
//                 internal::padd(
//                   internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+iN0])),internal::pmul(ptmp1,internal::ploadu(&A[j+iN1]))),
//                   internal::padd(internal::pmul(ptmp2,internal::pload(&A[j+iN2])),internal::pmul(ptmp3,internal::ploadu(&A[j+iN3]))) )));
//           }
//         }
//         else
//         {
//           for (int j = 0;j<ANP;j+=2*PacketSize)
//           {
//             internal::pstore(&X[j],
//               internal::padd(internal::pload(&X[j]),
//                 internal::padd(
//                   internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+iN0])),internal::pmul(ptmp1,internal::ploadu(&A[j+iN1]))),
//                   internal::padd(internal::pmul(ptmp2,internal::ploadu(&A[j+iN2])),internal::pmul(ptmp3,internal::ploadu(&A[j+iN3]))) )));
//
//             internal::pstore(&X[j+PacketSize],
//               internal::padd(internal::pload(&X[j+PacketSize]),
//                 internal::padd(
//                   internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+PacketSize+iN0])),internal::pmul(ptmp1,internal::ploadu(&A[j+PacketSize+iN1]))),
//                   internal::padd(internal::pmul(ptmp2,internal::ploadu(&A[j+PacketSize+iN2])),internal::pmul(ptmp3,internal::ploadu(&A[j+PacketSize+iN3]))) )));
//
// //             internal::pstore(&X[j+2*PacketSize],
// //               internal::padd(internal::pload(&X[j+2*PacketSize]),
// //                 internal::padd(
// //                   internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+2*PacketSize+iN0])),internal::pmul(ptmp1,internal::ploadu(&A[j+2*PacketSize+iN1]))),
// //                   internal::padd(internal::pmul(ptmp2,internal::ploadu(&A[j+2*PacketSize+iN2])),internal::pmul(ptmp3,internal::ploadu(&A[j+2*PacketSize+iN3]))) )));
// //
// //             internal::pstore(&X[j+3*PacketSize],
// //               internal::padd(internal::pload(&X[j+3*PacketSize]),
// //                 internal::padd(
// //                   internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+3*PacketSize+iN0])),internal::pmul(ptmp1,internal::ploadu(&A[j+3*PacketSize+iN1]))),
// //                   internal::padd(internal::pmul(ptmp2,internal::ploadu(&A[j+3*PacketSize+iN2])),internal::pmul(ptmp3,internal::ploadu(&A[j+3*PacketSize+iN3]))) )));
//
//           }
//           for (int j = ANP;j<AN;j+=PacketSize)
//             internal::pstore(&X[j],
//               internal::padd(internal::pload(&X[j]),
//                 internal::padd(
//                   internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+iN0])),internal::pmul(ptmp1,internal::ploadu(&A[j+iN1]))),
//                   internal::padd(internal::pmul(ptmp2,internal::ploadu(&A[j+iN2])),internal::pmul(ptmp3,internal::ploadu(&A[j+iN3]))) )));
//         }
//       }
//       // process remaining scalars
//       for (int j=AN;j<N;j++)
//         X[j] += tmp0 * A[j+iN0] + tmp1 * A[j+iN1] + tmp2 * A[j+iN2] + tmp3 * A[j+iN3];
//     }
//     for (int i=bound;i<N;i++)
//     {
//       real tmp0 = B[i];
//       Packet ptmp0 = internal::pset1(tmp0);
//       int iN0 = i*N;
//       if (AN>0)
//       {
//         bool aligned0 = (iN0 % PacketSize) == 0;
//         if (aligned0)
//           for (int j = 0;j<AN;j+=PacketSize)
//             internal::pstore(&X[j], internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+iN0])),internal::pload(&X[j])));
//         else
//           for (int j = 0;j<AN;j+=PacketSize)
//             internal::pstore(&X[j], internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+iN0])),internal::pload(&X[j])));
//       }
//       // process remaining scalars
//       for (int j=AN;j<N;j++)
//         X[j] += tmp0 * A[j+iN0];
//     }
//     asm("#end matrix_vector_product");
//   }

//   static inline void matrix_vector_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
//   {
//     asm("#begin matrix_vector_product");
//     int AN = (N/PacketSize)*PacketSize;
//     for (int i=0;i<N;i++)
//       X[i] = 0;
//
//     for (int i=0;i<N;i+=2)
//     {
//       real tmp0 = B[i];
//       Packet ptmp0 = internal::pset1(tmp0);
//       real tmp1 = B[i+1];
//       Packet ptmp1 = internal::pset1(tmp1);
//       int iN0 = i*N;
//       int iN1 = (i+1)*N;
//       if (AN>0)
//       {
//         bool aligned0 = (iN0 % PacketSize) == 0;
//         bool aligned1 = (iN1 % PacketSize) == 0;
//
//         if (aligned0 && aligned1)
//         {
//           for (int j = 0;j<AN;j+=PacketSize)
//           {
//             internal::pstore(&X[j],
//               internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+iN0])),
//               internal::padd(internal::pmul(ptmp1,internal::pload(&A[j+iN1])),internal::pload(&X[j]))));
//           }
//         }
//         else if (aligned0)
//         {
//           for (int j = 0;j<AN;j+=PacketSize)
//           {
//             internal::pstore(&X[j],
//               internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+iN0])),
//               internal::padd(internal::pmul(ptmp1,internal::ploadu(&A[j+iN1])),internal::pload(&X[j]))));
//           }
//         }
//         else if (aligned1)
//         {
//           for (int j = 0;j<AN;j+=PacketSize)
//           {
//             internal::pstore(&X[j],
//               internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+iN0])),
//               internal::padd(internal::pmul(ptmp1,internal::pload(&A[j+iN1])),internal::pload(&X[j]))));
//           }
//         }
//         else
//         {
//           int ANP = (AN/(4*PacketSize))*4*PacketSize;
//           for (int j = 0;j<ANP;j+=4*PacketSize)
//           {
//             internal::pstore(&X[j],
//               internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+iN0])),
//               internal::padd(internal::pmul(ptmp1,internal::ploadu(&A[j+iN1])),internal::pload(&X[j]))));
//
//             internal::pstore(&X[j+PacketSize],
//               internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+PacketSize+iN0])),
//               internal::padd(internal::pmul(ptmp1,internal::ploadu(&A[j+PacketSize+iN1])),internal::pload(&X[j+PacketSize]))));
//
//             internal::pstore(&X[j+2*PacketSize],
//               internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+2*PacketSize+iN0])),
//               internal::padd(internal::pmul(ptmp1,internal::ploadu(&A[j+2*PacketSize+iN1])),internal::pload(&X[j+2*PacketSize]))));
//
//             internal::pstore(&X[j+3*PacketSize],
//               internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+3*PacketSize+iN0])),
//               internal::padd(internal::pmul(ptmp1,internal::ploadu(&A[j+3*PacketSize+iN1])),internal::pload(&X[j+3*PacketSize]))));
//           }
//           for (int j = ANP;j<AN;j+=PacketSize)
//             internal::pstore(&X[j],
//               internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+iN0])),
//               internal::padd(internal::pmul(ptmp1,internal::ploadu(&A[j+iN1])),internal::pload(&X[j]))));
//         }
//       }
//       // process remaining scalars
//       for (int j=AN;j<N;j++)
//         X[j] += tmp0 * A[j+iN0] + tmp1 * A[j+iN1];
//     }
//     int remaining = (N/2)*2;
//     for (int i=remaining;i<N;i++)
//     {
//       real tmp0 = B[i];
//       Packet ptmp0 = internal::pset1(tmp0);
//       int iN0 = i*N;
//       if (AN>0)
//       {
//         bool aligned0 = (iN0 % PacketSize) == 0;
//         if (aligned0)
//           for (int j = 0;j<AN;j+=PacketSize)
//             internal::pstore(&X[j], internal::padd(internal::pmul(ptmp0,internal::pload(&A[j+iN0])),internal::pload(&X[j])));
//         else
//           for (int j = 0;j<AN;j+=PacketSize)
//             internal::pstore(&X[j], internal::padd(internal::pmul(ptmp0,internal::ploadu(&A[j+iN0])),internal::pload(&X[j])));
//       }
//       // process remaining scalars
//       for (int j=AN;j<N;j++)
//         X[j] += tmp0 * A[j+iN0];
//     }
//     asm("#end matrix_vector_product");
//   }

//   static inline void matrix_vector_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
//   {
//     asm("#begin matrix_vector_product");
//     int AN = (N/PacketSize)*PacketSize;
//     for (int i=0;i<N;i++)
//       X[i] = 0;
//     for (int i=0;i<N;i++)
//     {
//       real tmp = B[i];
//       Packet ptmp = internal::pset1(tmp);
//       int iN = i*N;
//       if (AN>0)
//       {
//         bool aligned = (iN % PacketSize) == 0;
//         if (aligned)
//         {
//           #ifdef PEELING
//           Packet A0, A1, A2, X0, X1, X2;
//           int ANP = (AN/(8*PacketSize))*8*PacketSize;
//           for (int j = 0;j<ANP;j+=PacketSize*8)
//           {
//             A0 = internal::pload(&A[j+iN]);
//             X0 = internal::pload(&X[j]);
//             A1 = internal::pload(&A[j+PacketSize+iN]);
//             X1 = internal::pload(&X[j+PacketSize]);
//             A2 = internal::pload(&A[j+2*PacketSize+iN]);
//             X2 = internal::pload(&X[j+2*PacketSize]);
//             internal::pstore(&X[j], internal::padd(X0, internal::pmul(ptmp,A0)));
//             A0 = internal::pload(&A[j+3*PacketSize+iN]);
//             X0 = internal::pload(&X[j+3*PacketSize]);
//             internal::pstore(&X[j+PacketSize], internal::padd(internal::pload(&X1), internal::pmul(ptmp,A1)));
//             A1 = internal::pload(&A[j+4*PacketSize+iN]);
//             X1 = internal::pload(&X[j+4*PacketSize]);
//             internal::pstore(&X[j+2*PacketSize], internal::padd(internal::pload(&X2), internal::pmul(ptmp,A2)));
//             A2 = internal::pload(&A[j+5*PacketSize+iN]);
//             X2 = internal::pload(&X[j+5*PacketSize]);
//             internal::pstore(&X[j+3*PacketSize], internal::padd(internal::pload(&X0), internal::pmul(ptmp,A0)));
//             A0 = internal::pload(&A[j+6*PacketSize+iN]);
//             X0 = internal::pload(&X[j+6*PacketSize]);
//             internal::pstore(&X[j+4*PacketSize], internal::padd(internal::pload(&X1), internal::pmul(ptmp,A1)));
//             A1 = internal::pload(&A[j+7*PacketSize+iN]);
//             X1 = internal::pload(&X[j+7*PacketSize]);
//             internal::pstore(&X[j+5*PacketSize], internal::padd(internal::pload(&X2), internal::pmul(ptmp,A2)));
//             internal::pstore(&X[j+6*PacketSize], internal::padd(internal::pload(&X0), internal::pmul(ptmp,A0)));
//             internal::pstore(&X[j+7*PacketSize], internal::padd(internal::pload(&X1), internal::pmul(ptmp,A1)));
// //
// //             internal::pstore(&X[j], internal::padd(internal::pload(&X[j]), internal::pmul(ptmp,internal::pload(&A[j+iN]))));
// //             internal::pstore(&X[j+PacketSize], internal::padd(internal::pload(&X[j+PacketSize]), internal::pmul(ptmp,internal::pload(&A[j+PacketSize+iN]))));
// //             internal::pstore(&X[j+2*PacketSize], internal::padd(internal::pload(&X[j+2*PacketSize]), internal::pmul(ptmp,internal::pload(&A[j+2*PacketSize+iN]))));
// //             internal::pstore(&X[j+3*PacketSize], internal::padd(internal::pload(&X[j+3*PacketSize]), internal::pmul(ptmp,internal::pload(&A[j+3*PacketSize+iN]))));
// //             internal::pstore(&X[j+4*PacketSize], internal::padd(internal::pload(&X[j+4*PacketSize]), internal::pmul(ptmp,internal::pload(&A[j+4*PacketSize+iN]))));
// //             internal::pstore(&X[j+5*PacketSize], internal::padd(internal::pload(&X[j+5*PacketSize]), internal::pmul(ptmp,internal::pload(&A[j+5*PacketSize+iN]))));
// //             internal::pstore(&X[j+6*PacketSize], internal::padd(internal::pload(&X[j+6*PacketSize]), internal::pmul(ptmp,internal::pload(&A[j+6*PacketSize+iN]))));
// //             internal::pstore(&X[j+7*PacketSize], internal::padd(internal::pload(&X[j+7*PacketSize]), internal::pmul(ptmp,internal::pload(&A[j+7*PacketSize+iN]))));
//           }
//           for (int j = ANP;j<AN;j+=PacketSize)
//             internal::pstore(&X[j], internal::padd(internal::pload(&X[j]), internal::pmul(ptmp,internal::pload(&A[j+iN]))));
//           #else
//           for (int j = 0;j<AN;j+=PacketSize)
//             internal::pstore(&X[j], internal::padd(internal::pload(&X[j]), internal::pmul(ptmp,internal::pload(&A[j+iN]))));
//           #endif
//         }
//         else
//         {
//           #ifdef PEELING
//           int ANP = (AN/(8*PacketSize))*8*PacketSize;
//           for (int j = 0;j<ANP;j+=PacketSize*8)
//           {
//             internal::pstore(&X[j], internal::padd(internal::pload(&X[j]), internal::pmul(ptmp,internal::ploadu(&A[j+iN]))));
//             internal::pstore(&X[j+PacketSize], internal::padd(internal::pload(&X[j+PacketSize]), internal::pmul(ptmp,internal::ploadu(&A[j+PacketSize+iN]))));
//             internal::pstore(&X[j+2*PacketSize], internal::padd(internal::pload(&X[j+2*PacketSize]), internal::pmul(ptmp,internal::ploadu(&A[j+2*PacketSize+iN]))));
//             internal::pstore(&X[j+3*PacketSize], internal::padd(internal::pload(&X[j+3*PacketSize]), internal::pmul(ptmp,internal::ploadu(&A[j+3*PacketSize+iN]))));
//             internal::pstore(&X[j+4*PacketSize], internal::padd(internal::pload(&X[j+4*PacketSize]), internal::pmul(ptmp,internal::ploadu(&A[j+4*PacketSize+iN]))));
//             internal::pstore(&X[j+5*PacketSize], internal::padd(internal::pload(&X[j+5*PacketSize]), internal::pmul(ptmp,internal::ploadu(&A[j+5*PacketSize+iN]))));
//             internal::pstore(&X[j+6*PacketSize], internal::padd(internal::pload(&X[j+6*PacketSize]), internal::pmul(ptmp,internal::ploadu(&A[j+6*PacketSize+iN]))));
//             internal::pstore(&X[j+7*PacketSize], internal::padd(internal::pload(&X[j+7*PacketSize]), internal::pmul(ptmp,internal::ploadu(&A[j+7*PacketSize+iN]))));
//           }
//           for (int j = ANP;j<AN;j+=PacketSize)
//             internal::pstore(&X[j], internal::padd(internal::pload(&X[j]), internal::pmul(ptmp,internal::ploadu(&A[j+iN]))));
//           #else
//           for (int j = 0;j<AN;j+=PacketSize)
//             internal::pstore(&X[j], internal::padd(internal::pload(&X[j]), internal::pmul(ptmp,internal::ploadu(&A[j+iN]))));
//           #endif
//         }
//       }
//       // process remaining scalars
//       for (int j=AN;j<N;j++)
//         X[j] += tmp * A[j+iN];
//     }
//     asm("#end matrix_vector_product");
//   }

    static inline void atv_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
  {
    int AN = (N/PacketSize)*PacketSize;
    int bound = (N/4)*4;
    for (int i=0;i<bound;i+=4)
    {
      real tmp0 = 0;
      Packet ptmp0 = internal::pset1(real(0));
      real tmp1 = 0;
      Packet ptmp1 = internal::pset1(real(0));
      real tmp2 = 0;
      Packet ptmp2 = internal::pset1(real(0));
      real tmp3 = 0;
      Packet ptmp3 = internal::pset1(real(0));
      int iN0 = i*N;
      int iN1 = (i+1)*N;
      int iN2 = (i+2)*N;
      int iN3 = (i+3)*N;
      if (AN>0)
      {
        int align1 = (iN1 % PacketSize);
        if (align1==0)
        {
          for (int j = 0;j<AN;j+=PacketSize)
          {
            Packet b = internal::pload(&B[j]);
            ptmp0 = internal::padd(ptmp0, internal::pmul(b, internal::pload(&A[j+iN0])));
            ptmp1 = internal::padd(ptmp1, internal::pmul(b, internal::pload(&A[j+iN1])));
            ptmp2 = internal::padd(ptmp2, internal::pmul(b, internal::pload(&A[j+iN2])));
            ptmp3 = internal::padd(ptmp3, internal::pmul(b, internal::pload(&A[j+iN3])));
          }
        }
        else if (align1==2)
        {
          for (int j = 0;j<AN;j+=PacketSize)
          {
            Packet b = internal::pload(&B[j]);
            ptmp0 = internal::padd(ptmp0, internal::pmul(b, internal::pload(&A[j+iN0])));
            ptmp1 = internal::padd(ptmp1, internal::pmul(b, internal::ploadu(&A[j+iN1])));
            ptmp2 = internal::padd(ptmp2, internal::pmul(b, internal::pload(&A[j+iN2])));
            ptmp3 = internal::padd(ptmp3, internal::pmul(b, internal::ploadu(&A[j+iN3])));
          }
        }
        else
        {
          for (int j = 0;j<AN;j+=PacketSize)
          {
            Packet b = internal::pload(&B[j]);
            ptmp0 = internal::padd(ptmp0, internal::pmul(b, internal::pload(&A[j+iN0])));
            ptmp1 = internal::padd(ptmp1, internal::pmul(b, internal::ploadu(&A[j+iN1])));
            ptmp2 = internal::padd(ptmp2, internal::pmul(b, internal::ploadu(&A[j+iN2])));
            ptmp3 = internal::padd(ptmp3, internal::pmul(b, internal::ploadu(&A[j+iN3])));
          }
        }
        tmp0 = internal::predux(ptmp0);
        tmp1 = internal::predux(ptmp1);
        tmp2 = internal::predux(ptmp2);
        tmp3 = internal::predux(ptmp3);
      }
      // process remaining scalars
      for (int j=AN;j<N;j++)
      {
        tmp0 += B[j] * A[j+iN0];
        tmp1 += B[j] * A[j+iN1];
        tmp2 += B[j] * A[j+iN2];
        tmp3 += B[j] * A[j+iN3];
      }
      X[i+0] = tmp0;
      X[i+1] = tmp1;
      X[i+2] = tmp2;
      X[i+3] = tmp3;
    }

    for (int i=bound;i<N;i++)
    {
      real tmp0 = 0;
      Packet ptmp0 = internal::pset1(real(0));
      int iN0 = i*N;
      if (AN>0)
      {
        if (iN0 % PacketSize==0)
          for (int j = 0;j<AN;j+=PacketSize)
            ptmp0 = internal::padd(ptmp0, internal::pmul(internal::pload(&B[j]), internal::pload(&A[j+iN0])));
        else
          for (int j = 0;j<AN;j+=PacketSize)
            ptmp0 = internal::padd(ptmp0, internal::pmul(internal::pload(&B[j]), internal::ploadu(&A[j+iN0])));
        tmp0 = internal::predux(ptmp0);
      }
      // process remaining scalars
      for (int j=AN;j<N;j++)
        tmp0 += B[j] * A[j+iN0];
      X[i+0] = tmp0;
    }
  }

//   static inline void atv_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
//   {
//     int AN = (N/PacketSize)*PacketSize;
//     for (int i=0;i<N;i++)
//       X[i] = 0;
//     for (int i=0;i<N;i++)
//     {
//       real tmp = 0;
//       Packet ptmp = internal::pset1(real(0));
//       int iN = i*N;
//       if (AN>0)
//       {
//         bool aligned = (iN % PacketSize) == 0;
//         if (aligned)
//         {
//           #ifdef PEELING
//           int ANP = (AN/(8*PacketSize))*8*PacketSize;
//           for (int j = 0;j<ANP;j+=PacketSize*8)
//           {
//             ptmp =
//               internal::padd(internal::pmul(internal::pload(&B[j]), internal::pload(&A[j+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+PacketSize]), internal::pload(&A[j+PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+2*PacketSize]), internal::pload(&A[j+2*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+3*PacketSize]), internal::pload(&A[j+3*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+4*PacketSize]), internal::pload(&A[j+4*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+5*PacketSize]), internal::pload(&A[j+5*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+6*PacketSize]), internal::pload(&A[j+6*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+7*PacketSize]), internal::pload(&A[j+7*PacketSize+iN])),
//               ptmp))))))));
//           }
//           for (int j = ANP;j<AN;j+=PacketSize)
//             ptmp = internal::padd(ptmp, internal::pmul(internal::pload(&B[j]), internal::pload(&A[j+iN])));
//           #else
//           for (int j = 0;j<AN;j+=PacketSize)
//             ptmp = internal::padd(ptmp, internal::pmul(internal::pload(&B[j]), internal::pload(&A[j+iN])));
//           #endif
//         }
//         else
//         {
//           #ifdef PEELING
//           int ANP = (AN/(8*PacketSize))*8*PacketSize;
//           for (int j = 0;j<ANP;j+=PacketSize*8)
//           {
//             ptmp =
//               internal::padd(internal::pmul(internal::pload(&B[j]), internal::ploadu(&A[j+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+PacketSize]), internal::ploadu(&A[j+PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+2*PacketSize]), internal::ploadu(&A[j+2*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+3*PacketSize]), internal::ploadu(&A[j+3*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+4*PacketSize]), internal::ploadu(&A[j+4*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+5*PacketSize]), internal::ploadu(&A[j+5*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+6*PacketSize]), internal::ploadu(&A[j+6*PacketSize+iN])),
//               internal::padd(internal::pmul(internal::pload(&B[j+7*PacketSize]), internal::ploadu(&A[j+7*PacketSize+iN])),
//               ptmp))))))));
//           }
//           for (int j = ANP;j<AN;j+=PacketSize)
//             ptmp = internal::padd(ptmp, internal::pmul(internal::pload(&B[j]), internal::ploadu(&A[j+iN])));
//           #else
//           for (int j = 0;j<AN;j+=PacketSize)
//             ptmp = internal::padd(ptmp, internal::pmul(internal::pload(&B[j]), internal::ploadu(&A[j+iN])));
//           #endif
//         }
//         tmp = internal::predux(ptmp);
//       }
//       // process remaining scalars
//       for (int j=AN;j<N;j++)
//         tmp += B[j] * A[j+iN];
//       X[i] = tmp;
//     }
//   }

  static inline void axpy(real coef, const gene_vector & X, gene_vector & Y, int N){
    int AN = (N/PacketSize)*PacketSize;
    if (AN>0)
    {
      Packet pcoef = internal::pset1(coef);
      #ifdef PEELING
      const int peelSize = 3;
      int ANP = (AN/(peelSize*PacketSize))*peelSize*PacketSize;
      float* X1 = X + PacketSize;
      float* Y1 = Y + PacketSize;
      float* X2 = X + 2*PacketSize;
      float* Y2 = Y + 2*PacketSize;
      Packet x0,x1,x2,y0,y1,y2;
      for (int j = 0;j<ANP;j+=PacketSize*peelSize)
      {
        x0 = internal::pload(X+j);
        x1 = internal::pload(X1+j);
        x2 = internal::pload(X2+j);

        y0 = internal::pload(Y+j);
        y1 = internal::pload(Y1+j);
        y2 = internal::pload(Y2+j);

        y0 = internal::pmadd(pcoef, x0, y0);
        y1 = internal::pmadd(pcoef, x1, y1);
        y2 = internal::pmadd(pcoef, x2, y2);

        internal::pstore(Y+j,  y0);
        internal::pstore(Y1+j, y1);
        internal::pstore(Y2+j, y2);
//         internal::pstore(&Y[j+2*PacketSize], internal::padd(internal::pload(&Y[j+2*PacketSize]), internal::pmul(pcoef,internal::pload(&X[j+2*PacketSize]))));
//         internal::pstore(&Y[j+3*PacketSize], internal::padd(internal::pload(&Y[j+3*PacketSize]), internal::pmul(pcoef,internal::pload(&X[j+3*PacketSize]))));
//         internal::pstore(&Y[j+4*PacketSize], internal::padd(internal::pload(&Y[j+4*PacketSize]), internal::pmul(pcoef,internal::pload(&X[j+4*PacketSize]))));
//         internal::pstore(&Y[j+5*PacketSize], internal::padd(internal::pload(&Y[j+5*PacketSize]), internal::pmul(pcoef,internal::pload(&X[j+5*PacketSize]))));
//         internal::pstore(&Y[j+6*PacketSize], internal::padd(internal::pload(&Y[j+6*PacketSize]), internal::pmul(pcoef,internal::pload(&X[j+6*PacketSize]))));
//         internal::pstore(&Y[j+7*PacketSize], internal::padd(internal::pload(&Y[j+7*PacketSize]), internal::pmul(pcoef,internal::pload(&X[j+7*PacketSize]))));
      }
      for (int j = ANP;j<AN;j+=PacketSize)
        internal::pstore(&Y[j], internal::padd(internal::pload(&Y[j]), internal::pmul(pcoef,internal::pload(&X[j]))));
      #else
      for (int j = 0;j<AN;j+=PacketSize)
        internal::pstore(&Y[j], internal::padd(internal::pload(&Y[j]), internal::pmul(pcoef,internal::pload(&X[j]))));
      #endif
    }
    // process remaining scalars
    for (int i=AN;i<N;i++)
      Y[i] += coef * X[i];
  }


};

#endif
