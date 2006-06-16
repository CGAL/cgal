// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#include "CGAL/benchmark_basic.hpp"
#include "CGAL/Bench.hpp"

CGAL_BENCHMARK_BEGIN_NAMESPACE

bool Bench_base::m_got_signal = false;
int Bench_base::m_name_length = 32;

CGAL_BENCHMARK_END_NAMESPACE
