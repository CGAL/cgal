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

#ifndef CGAL_CLASSIFICATION_FEATURES_EIGEN_H
#define CGAL_CLASSIFICATION_FEATURES_EIGEN_H

#include <CGAL/license/Classification.h>

#include <vector>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Local_eigen_analysis.h>

namespace CGAL {

namespace Classification {

namespace Feature {

/// \cond SKIP_IN_MANUAL
class Eigen_feature : public Feature_base
{
protected:
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
  std::vector<float> attrib;
#else
  const Classification::Local_eigen_analysis& eigen;
#endif
  
public:
  template <typename InputRange>
  Eigen_feature (const InputRange&,
                 const Classification::Local_eigen_analysis& eigen)
#ifndef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    : eigen (eigen)
#endif
  {
  }

#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
  virtual void init (std::size_t size, const Classification::Local_eigen_analysis& eigen)
  {
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      attrib.push_back (get_value (eigen, i));
  }
#else
  virtual void init (std::size_t, const Classification::Local_eigen_analysis&)
  {
  }
#endif
  
  virtual float get_value (const Classification::Local_eigen_analysis& eigen, std::size_t i) = 0;
  virtual float value (std::size_t pt_index)
  {
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    return attrib[pt_index];
#else
    return get_value(eigen, pt_index);
#endif
  }

};
/// \endcond
  
  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance matrix of a
    local neighborhood. Linearity is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_1 - \lambda_2}{\lambda_1}
    \f]

    Its default name is "linearity".
  */
class Linearity
#ifdef DOXYGEN_RUNNING
  : public Feature_base
#else
  : public Eigen_feature
#endif
{
public:
  /*!
    Constructs the feature.

    \tparam Input model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \param input point range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Linearity (const InputRange& input,
             const Local_eigen_analysis& eigen) : Eigen_feature (input, eigen)
  {
    this->set_name("linearity");
    this->init(input.size(), eigen);
  }

  /// \cond SKIP_IN_MANUAL
  virtual float get_value (const Local_eigen_analysis& eigen, std::size_t i)
  {
    const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[2] < 1e-15)
      return 0.;
    else
      return ((ev[2] - ev[1]) / ev[2]);
  }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance matrix of a
    local neighborhood. Planarity is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_2 - \lambda_3}{\lambda_1}
    \f]

    Its default name is "planarity".
  */
class Planarity
#ifdef DOXYGEN_RUNNING
  : public Feature_base
#else
  : public Eigen_feature
#endif
{
public:
  /*!
    Constructs the feature.

    \param input point range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Planarity (const InputRange& input,
             const Local_eigen_analysis& eigen)
    : Eigen_feature(input, eigen)
  {
    this->set_name("planarity");
    this->init(input.size(), eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual float get_value (const Local_eigen_analysis& eigen, std::size_t i)
  {
    const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[2] < 1e-15)
      return 0.;
    else
      return ((ev[1] - ev[0]) / ev[2]);
  }
  /// \endcond
 
};

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance matrix of a
    local neighborhood. Sphericity is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_3}{\lambda_1}
    \f]

    Its default name is "sphericity".
  */
class Sphericity
#ifdef DOXYGEN_RUNNING
  : public Feature_base
#else
  : public Eigen_feature
#endif
{
public:
  /*!
    Constructs the feature.

    \param input point range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Sphericity (const InputRange& input,
              const Local_eigen_analysis& eigen)
    : Eigen_feature(input, eigen)
  {
    this->set_name("sphericity");
    this->init(input.size(), eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual float get_value (const Local_eigen_analysis& eigen, std::size_t i)
  {
    const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[2] < 1e-15)
      return 0.;
    else
      return (ev[0] / ev[2]);
  }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance matrix of a
    local neighborhood. Omnivariance is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    (\lambda_1 \times \lambda_2 \times \lambda_3)^{\frac{1}{3}}
    \f]

    Its default name is "omnivariance".
  */
class Omnivariance
#ifdef DOXYGEN_RUNNING
  : public Feature_base
#else
  : public Eigen_feature
#endif
{
public:
  /*!
    Constructs the feature.

    \param input point range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Omnivariance (const InputRange& input,
                const Local_eigen_analysis& eigen)
    : Eigen_feature(input, eigen)
  {
    this->set_name("omnivariance");
    this->init(input.size(), eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual float get_value (const Local_eigen_analysis& eigen, std::size_t i)
  {
    const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    return (std::pow (CGAL::abs(ev[0] * ev[1] * ev[2]), 0.333333333f));
  }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance matrix of a
    local neighborhood. Anisotropy is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_1 - \lambda_3}{\lambda_1}
    \f]
    
    Its default name is "anisotropy".
  */
class Anisotropy
#ifdef DOXYGEN_RUNNING
  : public Feature_base
#else
  : public Eigen_feature
#endif
{
public:
  /*!
    Constructs the feature.

    \param input point range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Anisotropy (const InputRange& input,
              const Local_eigen_analysis& eigen)
    : Eigen_feature(input, eigen)
  {
    this->set_name("anisotropy");
    this->init(input.size(), eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual float get_value (const Local_eigen_analysis& eigen, std::size_t i)
  {
    const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[2] < 1e-15)
      return 0.;
    else
      return ((ev[2] - ev[0]) / ev[2]);
  }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance matrix of a
    local neighborhood. Eigentropy is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    - \sum_{i=1}^3 \lambda_i \times \log{\lambda_i}
    \f]
    
    Its default name is "eigentropy".
  */
class Eigentropy
#ifdef DOXYGEN_RUNNING
  : public Feature_base
#else
  : public Eigen_feature
#endif
{
public:
  /*!
    Constructs the feature.

    \param input point range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Eigentropy (const InputRange& input,
              const Local_eigen_analysis& eigen)
    : Eigen_feature(input, eigen)
  {
    this->set_name("eigentropy");
    this->init(input.size(), eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual float get_value (const Local_eigen_analysis& eigen, std::size_t i)
  {
    const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[0] < 1e-15
        || ev[1] < 1e-15
        || ev[2] < 1e-15)
      return 0.;
    else
      return (- ev[0] * std::log(ev[0])
              - ev[1] * std::log(ev[1])
              - ev[2] * std::log(ev[2]));
  }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance matrix of a
    local neighborhood. The sum of the eigenvalues is defined, for the
    3 eigenvalues \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$,
    as:

    \f[
    \lambda_1 + \lambda_2 + \lambda_3
    \f]

    Its default name is "sum_eigen".
  */
class Sum_eigenvalues
#ifdef DOXYGEN_RUNNING
  : public Feature_base
#else
  : public Eigen_feature
#endif
{
public:
  /*!
    Constructs the feature.

    \param input point range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Sum_eigenvalues (const InputRange& input,
                   const Local_eigen_analysis& eigen)
    : Eigen_feature(input, eigen)
  {
    this->set_name("sum_eigen");
    this->init(input.size(), eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual float get_value (const Local_eigen_analysis& eigen, std::size_t i)
  {
    return eigen.sum_of_eigenvalues(i);
  }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on the eigenvalues of the covariance
    matrix of a local neighborhood. Surface variation is defined, for
    the 3 eigenvalues \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge
    0\f$, as:

    \f[
    \frac{\lambda_3}{\lambda_1 + \lambda_2 + \lambda_3}
    \f]

    Its default name is "surface_variation".
  */
class Surface_variation
#ifdef DOXYGEN_RUNNING
  : public Feature_base
#else
  : public Eigen_feature
#endif
{
public:
  /*!
    Constructs the feature.

    \param input point range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Surface_variation (const InputRange& input,
                     const Local_eigen_analysis& eigen)
    : Eigen_feature(input, eigen)
  {
    this->set_name("surface_variation");
    this->init(input.size(), eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual float get_value (const Local_eigen_analysis& eigen, std::size_t i)
  {
    const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[0] + ev[1] + ev[2] < 1e-15)
      return 0.;
    else
      return (ev[0] / (ev[0] + ev[1] + ev[2]));
  }
  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURES_EIGEN_H
