// Copyright (c) 2018 GeometryFactory Sarl (France).
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

#ifndef CGAL_CLASSIFICATION_FEATURES_EIGENVALUE_H
#define CGAL_CLASSIFICATION_FEATURES_EIGENVALUE_H

#include <CGAL/license/Classification.h>

#include <vector>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Local_eigen_analysis.h>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance matrix of a
    local neighborhood. 

    Its default name is "eigenvalue_0", "eigenvalue_1" or
    "eigenvalue_2", depending on which eigenvalue is chosen in the
    constructor.
  */
class Eigenvalue : public Feature_base
{
protected:

  const Classification::Local_eigen_analysis& eigen;
  unsigned int m_idx;
  
public:

  /*!
    Constructs the feature.

    \tparam Input model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
    \param idx index of the eigen value (0, 1 or 2).
  */
  template <typename InputRange>
  Eigenvalue (const InputRange& input,
              const Classification::Local_eigen_analysis& eigen,
              unsigned int idx)
    : eigen (eigen), m_idx (idx)
  {
    CGAL_USE(input);
    std::ostringstream oss;
    oss << "eigenvalue" << idx;
    this->set_name (oss.str());
  }

  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
    return eigen.eigenvalue(pt_index)[std::size_t(m_idx)];
  }
  /// \endcond

};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURES_EIGENVALUE_H
