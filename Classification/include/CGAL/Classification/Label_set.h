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

#ifndef CGAL_CLASSIFICATION_LABEL_SET_H
#define CGAL_CLASSIFICATION_LABEL_SET_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Label.h>

#include <vector>

namespace CGAL {

namespace Classification {
  
/*!
\ingroup PkgClassificationLabel

\brief Set of `Label` used as input by classification
algorithms.

*/
class Label_set
{
  typedef std::vector<Label_handle> Base;
  Base m_labels;
  
public:
  
  Label_set() { }
  
  /// \cond SKIP_IN_MANUAL
  virtual ~Label_set() { }
  /// \endcond

  /*!
    \brief Adds a label.

    \note Names are not used for identification: two labels in the
    same set can have the same name (but not the same handle).

    \param name name of the label.

    \return a handle to the newly added label.
  */
  Label_handle add (const char* name)
  {
    Label_handle out (new Classification::Label (name));
    m_labels.push_back (out);
    return out;
  }

  /*!
    \brief Removes a label.

    \param label the handle to the label that must be removed.

    \return `true` if the label was correctly removed,
    `false` if its handle was not found.
  */ 
  bool remove (Label_handle label)
  {
    std::size_t idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i] == label)
      {
        m_labels.erase (m_labels.begin() + i);
        idx = i;
        break;
      }
    if (idx == (std::size_t)(-1))
      return false;

    return true;
  }

  /*!
    \brief Returns how many labels are defined.
  */  
  std::size_t size () const
  {
    return m_labels.size();
  }
  
  /*!
    \brief Returns the \f$i^{th}\f$ label.
  */  
  Label_handle operator[] (std::size_t i) const
  {
    return m_labels[i];
  }


  /*!
    \brief Removes all labels.
  */
  void clear ()
  {
    m_labels.clear();
  }


};



} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_LABEL_SET_H
