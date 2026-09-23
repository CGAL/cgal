// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//

#ifndef CGAL_POLYHEDRON_3_TO_LCC_H
#define CGAL_POLYHEDRON_3_TO_LCC_H

#include <CGAL/license/Polyhedron.h>


#define CGAL_DEPRECATED_HEADER "<CGAL/Polyhedron_3_to_lcc.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/import_face_graph_to_lcc.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/assertions.h>
#include <iostream>
#include <map>
#include <CGAL/config.h>

namespace CGAL {

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \deprecated This function is deprecated since CGAL 6.2. Use `import_face_graph_to_lcc()` instead.
*/
template< class LCC, class Polyhedron >
CGAL_DEPRECATED
typename LCC::Dart_descriptor
import_from_polyhedron_3(LCC& alcc, const Polyhedron &apoly)
{
  return import_face_graph_to_lcc(apoly, alcc);
}


/*!
  \deprecated This function is deprecated since CGAL 6.3. Use `import_face_graph_to_lcc()` instead.
*/
template< class LCC, class Polyhedron >
CGAL_DEPRECATED
typename LCC::Dart_descriptor
polyhedron_3_to_lcc(LCC& alcc, const Polyhedron &apoly)
{
  return import_face_graph_to_lcc(apoly, alcc);
}


/** converts a Polyhedron_3 read into a flux into 3D linear cell complex.
 * @param alcc the linear cell complex where Polyhedron_3 will be converted.
 * @param ais the istream where read the Polyhedron_3.
 * @return A dart created during the conversion.
 */
template < class LCC >
typename LCC::Dart_descriptor
polyhedron_3_flux_to_lcc(LCC& alcc, std::istream& ais)
{
  if (!ais.good())
  {
    std::cout << "Error reading flux." << std::endl;
    return LCC::null_descriptor;
  }
  CGAL::Polyhedron_3<typename LCC::Traits> P;
  ais >> P;
  return polyhedron_3_to_lcc<LCC, CGAL::Polyhedron_3
                                  <typename LCC::Traits> > (alcc, P);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_IMPORT_FROM_POLYHEDRON_3_H
