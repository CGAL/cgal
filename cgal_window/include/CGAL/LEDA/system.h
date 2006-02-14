// Copyright (c) 1999  Martin-Luther-University Halle-Wittenberg (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Matthias Baesken, Algorithmic Solutions

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#if defined(unix) || defined(__unix__) || defined(__unix) || defined(_AIX)
#include <CGAL/LEDA/sys/unix.h>

#elif defined(__CYGWIN32__)
#include <CGAL/LEDA/sys/cygwin32.h>

#elif defined(__WIN32__) || defined(_WIN32) || defined(__NT__)
#include <CGAL/LEDA/sys/win32.h>



#else
// unknown system

#endif

