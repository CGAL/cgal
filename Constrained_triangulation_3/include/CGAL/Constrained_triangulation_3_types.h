// Copyright (c) 2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CT_3_TYPES_H
#define CGAL_CT_3_TYPES_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/config.h>

#include <functional>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace CGAL {

/**
 * @ingroup PkgConstrainedTriangulation3Classes
 * \brief Signed integral type to store the index of constraints.
 */
using CDT_3_signed_index = int; // must be signed

/**
 * @ingroup PkgConstrainedTriangulation3Classes
 * \brief Exception type thrown when a constrained Delaunay triangulation
 * cannot be restored after the insertion of constraints.
 */

 class Constrained_triangulation_insertion_exception : public std::domain_error
 {
 public:
   /** \brief Constructor.
    * \param failed_faces a vector of indices of the faces that could not be processed.
    */
   Constrained_triangulation_insertion_exception(const std::vector<CDT_3_signed_index>& failed_faces)
       : std::domain_error(std::invoke([&]() {
         std::ostringstream oss;
         oss << "Constrained Delaunay triangulation could not be restored: " << m_failed_faces.size()
             << " faces could not be processed:\n";
         for(const CDT_3_signed_index& index : m_failed_faces) {
           oss << " " << index << "\n";
         }
         return oss.str();
       }))
       , m_failed_faces(failed_faces) {}

   /** \brief Returns the vector of indices of the faces that could not be processed.
    */
   const std::vector<CDT_3_signed_index>& failed_faces() const {
     return m_failed_faces;
   }

 private:
   std::vector<CDT_3_signed_index> m_failed_faces;
 };

} // namespace CGAL

#endif // CGAL_CT_3_TYPES_H
