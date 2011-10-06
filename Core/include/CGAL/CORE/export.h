/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/); you may
 * redistribute it under the terms of the Q Public License version 1.0.
 * See the file LICENSE.QPL distributed with CORE.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: Expr.h
 * Synopsis: a class of Expression in Level 3
 * 
 * Written by 
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Igor Pechtchanski <pechtcha@cs.nyu.edu>
 *       Vijay Karamcheti <vijayk@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu> 
 *       Vikram Sharma<sharma@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/

// Author(s)     : Andreas Fabri

#ifndef CGAL_CORE_EXPORT_H
#define CGAL_CORE_EXPORT_H

#include <boost/config.hpp>

#if defined(BOOST_MSVC) && ( ! defined(CGAL_EXPORTS) ) && defined(CGAL_BUILD_SHARED_LIB)

#if defined(CGAL_Core_EXPORTS) // add by CMake or in cpp files of the dll
#define	CGAL_CORE_EXPORT __declspec (dllexport)
#define CGAL_CORE_EXPIMP_TEMPLATE
#else
#define CGAL_CORE_EXPORT __declspec (dllimport)
#define CGAL_CORE_EXPIMP_TEMPLATE extern
#endif
 
#else 

#define  CGAL_CORE_EXPORT
#define CGAL_CORE_EXPIMP_TEMPLATE 
#endif

#endif //  CGAL_CORE_EXPORT_H


