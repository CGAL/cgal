// Copyright (c) 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CERES_SUPPORT_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CERES_SUPPORT_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <boost/config.hpp>

// Ignore warnings on Windows of type
// " class XXXX needs to have dll-interface to be used by clients of struct YYYY" in an interaction
// between Ceres and Eigen, which doesn't make much sense as Eigen is header-only.

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4251)
#endif


// If a user is using minilog instead of glog, the verbosity cannot be controlled via ceres options...
// Since we really don't want to hear about it, override that value with '0'
#ifdef MINIGLOG_MAX_LOG_LEVEL
#  define MINIGLOG_MAX_LOG_LEVEL_WAS_DEFINED
#  pragma push_macro("MINIGLOG_MAX_LOG_LEVEL")
#  undef MINIGLOG_MAX_LOG_LEVEL
#endif

#define MAX_LOG_LEVEL 0

#include "ceres/ceres.h"

#ifdef MINIGLOG_MAX_LOG_LEVEL_WAS_DEFINED
#  pragma pop_macro("MINIGLOG_MAX_LOG_LEVEL")
#  undef MINIGLOG_MAX_LOG_LEVEL_WAS_DEFINED
#else
#  undef MINIGLOG_MAX_LOG_LEVEL
#endif


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CERES_SUPPORT_H
