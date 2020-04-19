/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
