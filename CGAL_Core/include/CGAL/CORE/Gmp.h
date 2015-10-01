/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * $URL$
 * $Id$
 ***************************************************************************/

// CORE LIBRARY FILE
#ifndef _CORE_GMP_H_
#define _CORE_GMP_H_

#include <CGAL/CORE/Impl.h>
#include <gmp.h>

namespace CORE { 

CGAL_CORE_EXPORT std::ostream& io_write (std::ostream &, mpz_srcptr);
CGAL_CORE_EXPORT std::ostream& io_write (std::ostream &, mpq_srcptr);
CGAL_CORE_EXPORT std::istream& io_read (std::istream &, mpz_ptr);
CGAL_CORE_EXPORT std::istream& io_read (std::istream &, mpq_ptr);
//std::ostream& operator<< (std::ostream &, mpz_srcptr);
//std::ostream& operator<< (std::ostream &, mpq_srcptr);
//std::istream& operator>> (std::istream &, mpz_ptr);
//std::istream& operator>> (std::istream &, mpq_ptr);

} //namespace CORE

#ifdef CGAL_HEADER_ONLY
#include <CGAL/CORE/Gmp_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // _CORE_GMP_H_
