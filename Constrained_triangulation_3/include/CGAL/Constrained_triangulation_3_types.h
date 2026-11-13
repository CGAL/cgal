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
#include <utility>
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
   Constrained_triangulation_insertion_exception(std::vector<CDT_3_signed_index> failed_faces)
       : std::domain_error(std::invoke([&]() {
         std::ostringstream oss;
         oss << "Constrained Delaunay triangulation could not be restored: " << failed_faces.size()
             << " faces could not be processed:\n";
         for(const CDT_3_signed_index& index : failed_faces) {
           oss << "  face #" << index << "\n";
         }
         return oss.str();
       }))
       , m_failed_faces(std::move(failed_faces)) {}

   /** \brief Returns the vector of indices of the faces that could not be processed.
    */
   const std::vector<CDT_3_signed_index>& failed_faces() const {
     return m_failed_faces;
   }

 private:
   std::vector<CDT_3_signed_index> m_failed_faces;
 };

/**
 * @ingroup PkgConstrainedTriangulation3Classes
 * \brief Exception type thrown when a piecewise linear complex (PLC) facet
 * is not coplanar.
 */
class Non_planar_plc_facet_exception : public std::exception
{
public:
  /** \brief Default constructor.
   */
  Non_planar_plc_facet_exception() = default;

  /** \brief Constructor taking the index of the faulty PLC facet.
   * \param polygon_id the index of the non-planar PLC facet.
   */
  explicit Non_planar_plc_facet_exception(CDT_3_signed_index polygon_id)
      : m_polygon_id(polygon_id),
        m_msg("ERROR: PLC facet #" + std::to_string(polygon_id) + " is not coplanar")
  {}

  /** \brief Returns the error message.
   */
  const char* what() const noexcept override
  {
    return m_msg.empty() ? "ERROR: PLC facet is not coplanar" : m_msg.c_str();
  }

  /** \brief Returns the index of the faulty PLC facet if available, or `-1` otherwise.
   */
  CDT_3_signed_index polygon_id() const noexcept {
    return m_polygon_id;
  }

private:
  CDT_3_signed_index m_polygon_id = -1;
  std::string m_msg;
};

} // namespace CGAL

#endif // CGAL_CT_3_TYPES_H
