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

#ifndef CGAL_CLASSIFICATION_LABEL_SET_H
#define CGAL_CLASSIFICATION_LABEL_SET_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Label.h>

#include <vector>

namespace CGAL {

namespace Classification {

/*!
\ingroup PkgClassificationLabel

\brief sets of `Label` used as input by classification
algorithms.

*/
class Label_set
{
  typedef std::vector<Label_handle> Base;
  Base m_labels;

public:

  Label_set() { }

  /*!
    \brief initializes the set with the provided `labels` names.
  */
  Label_set (const std::initializer_list<const char*>& labels)
  {
    for (const char* l : labels)
      add(l);
  }

  /// \cond SKIP_IN_MANUAL
  virtual ~Label_set() { }
  /// \endcond

  /*!
    \brief adds a label.

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
    \brief removes a label.

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
    \brief returns how many labels are defined.
  */
  std::size_t size () const
  {
    return m_labels.size();
  }

  /*!
    \brief returns the \f$i^{th}\f$ label.
  */
  Label_handle operator[] (std::size_t i) const
  {
    return m_labels[i];
  }


  /*!
    \brief removes all labels.
  */
  void clear ()
  {
    m_labels.clear();
  }


};



} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_LABEL_SET_H
