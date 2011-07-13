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

static int nb_load;
static int nb_loadu;
static int nb_store;
static int nb_storeu;

#define EIGEN_DEBUG_ALIGNED_LOAD    { nb_load++;    }
#define EIGEN_DEBUG_UNALIGNED_LOAD  { nb_loadu++;   }
#define EIGEN_DEBUG_ALIGNED_STORE   { nb_store++;   }
#define EIGEN_DEBUG_UNALIGNED_STORE { nb_storeu++;  }

#define VERIFY_ALIGNED_UNALIGNED_COUNT(XPR,AL,UL,AS,US) {\
    nb_load = nb_loadu = nb_store = nb_storeu = 0; \
    XPR; \
    if(!(nb_load==AL && nb_loadu==UL && nb_store==AS && nb_storeu==US)) \
      std::cerr << " >> " << nb_load << ", " << nb_loadu << ", " << nb_store << ", " << nb_storeu << "\n"; \
    VERIFY( (#XPR) && nb_load==AL && nb_loadu==UL && nb_store==AS && nb_storeu==US ); \
  }


#include "main.h"

void test_unalignedcount()
{
  #ifdef EIGEN_VECTORIZE_SSE
  VectorXf a(40), b(40);
  VERIFY_ALIGNED_UNALIGNED_COUNT(a += b, 20, 0, 10, 0);
  VERIFY_ALIGNED_UNALIGNED_COUNT(a.segment(0,40) += b.segment(0,40), 10, 10, 10, 0);
  VERIFY_ALIGNED_UNALIGNED_COUNT(a.segment(0,40) -= b.segment(0,40), 10, 10, 10, 0);
  VERIFY_ALIGNED_UNALIGNED_COUNT(a.segment(0,40) *= 3.5, 10, 0, 10, 0);
  VERIFY_ALIGNED_UNALIGNED_COUNT(a.segment(0,40) /= 3.5, 10, 0, 10, 0);
  #else
  // The following line is to eliminate "variable not used" warnings
  nb_load = nb_loadu = nb_store = nb_storeu = 0;
  #endif
}
