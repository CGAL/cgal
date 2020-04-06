// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_LABEL_H
#define CGAL_CLASSIFICATION_LABEL_H

#include <CGAL/license/Classification.h>
#include <CGAL/IO/Color.h>

#include <memory>

namespace CGAL {

namespace Classification {

/// \cond SKIP_IN_MANUAL
class Label_set;
/// \endcond

/*!
\ingroup PkgClassificationLabel

\brief %Classification label (for example: vegetation, ground, etc.)
defined as a set of relationships with classification features.

*/
class Label
{
private:

  std::string m_name;
  std::size_t m_index;
  std::size_t m_standard_index;
  CGAL::Color m_color;

  friend Label_set;

public:

  /// \cond SKIP_IN_MANUAL
  // Undocumented: Labels should be created by the set
  Label (std::string name, std::size_t index, std::size_t standard_index,
         const CGAL::Color& color)
    : m_name (name), m_index (index), m_standard_index (standard_index)
    , m_color (color)
  { }
  /// \endcond

  /*!
    returns the name of the classification label (e.g. vegetation).
  */
  const std::string& name() const { return m_name; }

  void set_name (const std::string& name) { m_name = name; }

  /*!
    returns the index of the classification label in the label set.

    \note This index cannot be changed by the user and is handled directly by the label set.
  */
  std::size_t index() const { return m_index; }

  /*!
    returns the standard index of the classification label (i.e. index in the ASPRS standard).

    \note This index is purely user-oriented and is not used by the classification algorithms.
  */
  std::size_t standard_index() const { return m_standard_index; }

  void set_standard_index(std::size_t idx) { m_standard_index = idx; }

  /*!
    returns the color used to represent the label.

    \note The color is purely user-oriented and is not used by the
    classification algorithms. It is not to be confused with a color
    attribute embedded in a data set which _can_ be used (see
    `Color_channel`).
  */
  const CGAL::Color& color() const { return m_color; }
  void set_color (const Color& color) { m_color = color; }
};

#ifdef DOXYGEN_RUNNING
/*!
  \ingroup PkgClassificationLabel

  \brief %Handle to a classification `Label`.

  \cgalModels Handle
*/
class Label_handle { };
#else
typedef std::shared_ptr<Label> Label_handle;
#endif

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_LABEL_H
