// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Clement Jamin
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_3_PROFILING_TOOLS_H
#define CGAL_MESH_3_PROFILING_TOOLS_H

#include <CGAL/license/Mesh_3.h>


// TBB timers
#ifdef CGAL_LINKED_WITH_TBB
  #include <tbb/tick_count.h>
  struct WallClockTimer
  {
    tbb::tick_count t;
    WallClockTimer()
    {
      t = tbb::tick_count::now();
    }
    void reset()
    {
      t = tbb::tick_count::now();
    }
    double elapsed() const
    {
      return (tbb::tick_count::now() - t).seconds();
    }
  };

#else
  #include <CGAL/Real_timer.h>

  struct WallClockTimer
  {
    CGAL::Real_timer t;
    WallClockTimer()
    {
      t.start();
    }
    void reset()
    {
      t.reset();
    }
    double elapsed() const
    {
      return t.time();
    }
  };
#endif

#endif // CGAL_MESH_3_PROFILING_TOOLS_H
