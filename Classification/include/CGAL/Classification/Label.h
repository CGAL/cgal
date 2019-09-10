// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0+
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
