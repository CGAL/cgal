// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>

#ifndef CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H
#define CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H

#ifndef CGAL_PM_BOUNDING_BOX_BASE_H
#include <CGAL/Pm_bounding_box_base.h>
#endif

#ifndef CGAL_LEDA_WINDOW_H
#include <CGAL/leda_window.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Planar_map>
inline Window_stream& operator<<(Window_stream& os,
                          const Pm_bounding_box_base<Planar_map> &b)
{}  

CGAL_END_NAMESPACE

#endif // CGAL_IO_PM_BOUNDING_BOX_BASE_WINDOW_STREAM_H


