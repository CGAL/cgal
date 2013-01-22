// Copyright (c) 2013
// INRIA Sophia-Antipolis (France).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_TBB_H
#define CGAL_TBB_H

#ifdef CGAL_LINKED_WITH_TBB

#include <tbb/tbb.h>
#include <tbb/task_scheduler_init.h>

namespace CGAL {

/*!
\ingroup Installation
 
\brief The class TBB_configuration provides control over TBB.
*/
class TBB_configuration
{
public:
  /*!
  \brief    Set the maximum number of threads that TBB may create.
  \details  If max_num_threads <= 0, TBB will decide the number of threads, 
            which is typically the number of hardware threads.
  */
  static void set_max_number_of_threads(int max_num_threads = -1)
  {
    static tbb::task_scheduler_init init(tbb::task_scheduler_init::deferred);
    
    if (init.is_active())
      init.terminate();

    if (max_num_threads > 0)
      init.initialize(max_num_threads);
    else
      init.initialize();
  }

  /*!
  \brief    Returns the number of threads TBB scheduler would create if initialized by 
            default.
  \details  This number is typically the number of hardware threads (logical 
            processors).
  */
  static int get_default_number_of_threads()
  {
    return tbb::task_scheduler_init::default_num_threads();
  }

protected:
  // We don't want this class to be instantiated by an user
  TBB_configuration() {}
};

} //namespace CGAL

#endif // CGAL_LINKED_WITH_TBB
#endif // CGAL_TBB_H
