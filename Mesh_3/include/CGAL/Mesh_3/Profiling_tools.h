// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Clement Jamin
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_PROFILING_TOOLS_H
#define CGAL_MESH_3_PROFILING_TOOLS_H

// TBB timers
#ifdef CGAL_LINKED_WITH_TBB
  #include <tbb/tick_count.h>
  struct Wall_clock_timer
  {
    tbb::tick_count t;
    Wall_clock_timer()
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

  struct Wall_clock_timer
  {
    CGAL::Real_timer t;
    Wall_clock_timer()
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
