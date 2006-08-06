// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_PROFILE_COUNTER_H
#define CGAL_PROFILE_COUNTER_H

// This file contains the class Profile_counter which is able to keep track
// of a number, and prints a message in the destructor.
// Typically, it can be used as a profile counter in a static variable.

#include <CGAL/basic.h>
#include <string>

CGAL_BEGIN_NAMESPACE

struct Profile_counter
{
    Profile_counter(const std::string & ss)
      : i(0), s(ss) {}

    void operator++() { ++i; }

    ~Profile_counter()
    {
        std::cerr << "[CGAL::Profile_counter] ";
        std::cerr.width(10);
        std::cerr << i << " " << s << std::endl;
    }

private:
    unsigned int i;
    const std::string s;
};

#ifdef CGAL_PROFILE
#  define CGAL_PROFILER(Y) { static CGAL::Profile_counter tmp(Y); ++tmp; }
#else
#  define CGAL_PROFILER(Y)
#endif

CGAL_END_NAMESPACE

#endif // CGAL_PROFILE_COUNTER_H
