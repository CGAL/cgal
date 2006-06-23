// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
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
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_BENCHMARK_CONFIG_HPP
#define CGAL_BENCHMARK_CONFIG_HPP

#ifdef CGAL_BENCHMARK_CFG_NO_NAMESPACE
#  define CGAL_BENCHMARK_USING_NAMESPACE_STD
#  define CGAL_BENCHMARK_STD
#  define CGAL_BENCHMARK std
#else
#  define CGAL_BENCHMARK_USING_NAMESPACE_STD using namespace std;
#  define CGAL_BENCHMARK_STD std
#  ifndef CGAL_BENCHMARK_USE_NAMESPACE
#    define CGAL_BENCHMARK_USE_NAMESPACE 1
#  endif
#endif

#if CGAL_BENCHMARK_USE_NAMESPACE
#  define CGAL_BENCHMARK_BEGIN_NAMESPACE namespace CGAL { namespace benchmark {
#  define CGAL_BENCHMARK_END_NAMESPACE }}
#else
#  define CGAL_BENCHMARK_BEGIN_NAMESPACE
#  define CGAL_BENCHMARK_END_NAMESPACE
#endif

#endif
