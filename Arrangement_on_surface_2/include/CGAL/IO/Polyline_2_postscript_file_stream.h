// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
//
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_POLYLINE_POSTSCRIPT_FILE_STREAM_H
#define CGAL_POLYLINE_POSTSCRIPT_FILE_STREAM_H

/*! \file
 * Postscript output stream for the Polyline_2<SegmentTraits> class.
 */

#include <CGAL/IO/Postscript_file_stream.h>
#include <CGAL/Arr_traits_2/Polyline_2.h>

CGAL_BEGIN_NAMESPACE

/*! Export a Polyline_2 object to a Postscript stream
 * \param ps_stream the Postscript stream
 * \param polyline the polyline curve
 * \return the Postscript stream
 */
template<typename SegmentTraits> Postscript_file_stream &
operator<<(Postscript_file_stream & ps_stream,
           const _Polyline_2<SegmentTraits> & polyline)
{
  unsigned int i;
  for (i = 0; i < polyline.size(); ++i) ps_stream << polyline[i];
  return ps_stream;
}

CGAL_END_NAMESPACE

#endif
