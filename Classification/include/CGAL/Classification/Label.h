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

#include <boost/shared_ptr.hpp>

namespace CGAL {

namespace Classification {

/*!
\ingroup PkgClassificationLabel

\brief %Classification label (for example: vegetation, ground, etc.)
defined as a set of relationships with classification features.

*/
class Label
{
private:

  std::string m_name;

public:

  /*! 
    \param name name of the classification label (e.g. vegetation).
  */ 
  Label (std::string name) : m_name (name) { }

  const std::string& name() const { return m_name; }

  /// \cond SKIP_IN_MANUAL
  void set_name (const std::string& name) { m_name = name; }
  /// \endcond
};

#ifdef DOXYGEN_RUNNING
/*!
  \ingroup PkgClassificationLabel

  \brief %Handle to a classification `Label`.

  \cgalModels Handle
*/
class Label_handle { };
#else
typedef boost::shared_ptr<Label> Label_handle;
#endif

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_LABEL_H
