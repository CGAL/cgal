// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_DEBUG_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_DEBUG_H

#include <CGAL/Straight_skeleton_3/internal/StackTrace.h>

#include <iostream>

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// #define CGAL_SS3_DUMP_FILES

#define CGAL_SS3_PROFILE_FILTERING_MECHANISMS
#define CGAL_SS3_RUN_TIMERS
#define CGAL_SS3_EXIT_ASAP

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

#define CGAL_SS3_TRACE_STREAM std::cout

/* just a helper for code location */
#define CGAL_SS3_TRACE_LOC CGAL_SS3_TRACE_STREAM << "[DEBUG] " << __FILE__ << ":" << __LINE__ << ": ";

#ifndef CGAL_SS3_TRACE_VERBOSITY
#  define CGAL_SS3_TRACE_VERBOSITY 4
#endif

#define CGAL_SS3_ENABLE_TRACE // generic
// #define CGAL_SS3_TRAITS_ENABLE_TRACE // traits & kernel
#define CGAL_SS3_HDS_ENABLE_TRACE // Polyhedron and related classes
// #define CGAL_SS3_SKEL_DS_ENABLE_TRACE // Skeleton and related classes
#define CGAL_SS3_IO_ENABLE_TRACE // DB, IO, ...
#define CGAL_SS3_TRANSF_ENABLE_TRACE // Polyhedron transformation (facet merging, perturbations, etc.)
#define CGAL_SS3_SPLITTER_ENABLE_TRACE // vertex splitters
#define CGAL_SS3_ALGO_ENABLE_TRACE // supporting algorithms
#define CGAL_SS3_CORE_ENABLE_TRACE // main algo

// -------------------------------------------------------------------------------------------------

// Generic
#ifdef CGAL_SS3_ENABLE_TRACE
# define CGAL_SS3_TRACE_CODE(code) code
# define CGAL_SS3_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------

#ifdef CGAL_SS3_TRAITS_ENABLE_TRACE
# define CGAL_SS3_TRAITS_TRACE_CODE(code) code
# define CGAL_SS3_TRAITS_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_TRAITS_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_TRAITS_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------

// Polyhedron and related classes
#ifdef CGAL_SS3_HDS_ENABLE_TRACE
# define CGAL_SS3_HDS_TRACE_CODE(code) code
# define CGAL_SS3_HDS_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_HDS_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_HDS_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------

// Skeleton and related classes
#ifdef CGAL_SS3_SKEL_DS_ENABLE_TRACE
# define CGAL_SS3_SKEL_DS_TRACE_CODE(code) code
# define CGAL_SS3_SKEL_DS_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_SKEL_DS_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_SKEL_DS_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------

// Reading and writing data
#ifdef CGAL_SS3_IO_ENABLE_TRACE
# define CGAL_SS3_IO_TRACE_CODE(code) code
# define CGAL_SS3_IO_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_IO_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_IO_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------

// Transformation of polyhedra
#ifdef CGAL_SS3_TRANSF_ENABLE_TRACE
# define CGAL_SS3_TRANSF_TRACE_CODE(code) code
# define CGAL_SS3_TRANSF_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_TRANSF_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_TRANSF_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------

// Vertex splitters
#ifdef CGAL_SS3_SPLITTER_ENABLE_TRACE
# define CGAL_SS3_SPLITTER_TRACE_CODE(code) code
# define CGAL_SS3_SPLITTER_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_SPLITTER_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_SPLITTER_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------

// Supporting algorithms such as self-intersection detection
#ifdef CGAL_SS3_ALGO_ENABLE_TRACE
# define CGAL_SS3_ALGO_TRACE_CODE(code) code
# define CGAL_SS3_ALGO_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_ALGO_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_ALGO_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------

// Main construction algorithms
#ifdef CGAL_SS3_CORE_ENABLE_TRACE
# define CGAL_SS3_CORE_TRACE_CODE(code) code
# define CGAL_SS3_CORE_TRACE(m) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_CORE_TRACE_V(l,m) if (l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
# define CGAL_SS3_CORE_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_SS3_TRACE_VERBOSITY) /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << m << std::endl;
#endif

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

#ifndef CGAL_SS3_TRACE_CODE
# define CGAL_SS3_TRACE_CODE(code)
#endif

#ifndef CGAL_SS3_ENABLE_TRACE
# define CGAL_SS3_TRACE_CODE(code) code
# define CGAL_SS3_TRACE(m)
# define CGAL_SS3_TRACE_V(l,m)
# define CGAL_SS3_TRACE_IF(c,l,m)
#endif

#ifndef CGAL_SS3_TRAITS_ENABLE_TRACE
# define CGAL_SS3_TRAITS_TRACE_CODE(code) code
# define CGAL_SS3_TRAITS_TRACE(m)
# define CGAL_SS3_TRAITS_TRACE_V(l,m)
# define CGAL_SS3_TRAITS_TRACE_IF(c,l,m)
#endif

#ifndef CGAL_SS3_HDS_ENABLE_TRACE
# define CGAL_SS3_HDS_TRACE_CODE(code)
# define CGAL_SS3_HDS_TRACE(m)
# define CGAL_SS3_HDS_TRACE_V(l,m)
# define CGAL_SS3_HDS_TRACE_IF(c,l,m)
#endif

#ifndef CGAL_SS3_SKEL_DS_ENABLE_TRACE
# define CGAL_SS3_SKEL_DS_TRACE_CODE(code)
# define CGAL_SS3_SKEL_DS_TRACE(m)
# define CGAL_SS3_SKEL_DS_TRACE_V(l,m)
# define CGAL_SS3_SKEL_DS_TRACE_IF(c,l,m)
#endif

#ifndef CGAL_SS3_IO_ENABLE_TRACE
# define CGAL_SS3_IO_TRACE_CODE(code)
# define CGAL_SS3_IO_TRACE(m)
# define CGAL_SS3_IO_TRACE_V(l,m)
# define CGAL_SS3_IO_TRACE_IF(c,l,m)
#endif

#ifndef CGAL_SS3_TRANSF_ENABLE_TRACE
# define CGAL_SS3_TRANSF_TRACE_CODE(code)
# define CGAL_SS3_TRANSF_TRACE(m)
# define CGAL_SS3_TRANSF_TRACE_V(l,m)
# define CGAL_SS3_TRANSF_TRACE_IF(c,l,m)
#endif

#ifndef CGAL_SS3_SPLITTER_ENABLE_TRACE
# define CGAL_SS3_SPLITTER_TRACE_CODE(code)
# define CGAL_SS3_SPLITTER_TRACE(m)
# define CGAL_SS3_SPLITTER_TRACE_V(l,m)
# define CGAL_SS3_SPLITTER_TRACE_IF(c,l,m)
#endif

#ifndef CGAL_SS3_ALGO_ENABLE_TRACE
# define CGAL_SS3_ALGO_TRACE_CODE(code)
# define CGAL_SS3_ALGO_TRACE(m)
# define CGAL_SS3_ALGO_TRACE_V(l,m)
# define CGAL_SS3_ALGO_TRACE_IF(c,l,m)
#endif

#ifndef CGAL_SS3_CORE_ENABLE_TRACE
# define CGAL_SS3_CORE_TRACE_CODE(code)
# define CGAL_SS3_CORE_TRACE(m)
# define CGAL_SS3_CORE_TRACE_V(l,m)
# define CGAL_SS3_CORE_TRACE_IF(c,l,m)
#endif


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

#ifdef DEBUG

/* macros that warn if a smart pointer is invalid */
#define CGAL_SS3_DEBUG_SPTR(sptr) if (!sptr) { /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << "shared pointer is invalid: " << (#sptr) << std::endl; util::StackTrace::print(CGAL_SS3_TRACE_STREAM); }
#define CGAL_SS3_DEBUG_WPTR(wptr) if (wptr.expired()) { /*CGAL_SS3_TRACE_LOC*/ CGAL_SS3_TRACE_STREAM << "weak pointer is expired: " << (#wptr) << std::endl; util::StackTrace::print(CGAL_SS3_TRACE_STREAM); }

#else

#define CGAL_SS3_DEBUG_SPTR(sptr)
#define CGAL_SS3_DEBUG_WPTR(wptr)

#endif

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_DEBUG_H */
