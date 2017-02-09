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
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_ATTRIBUTES_EIGEN_H
#define CGAL_CLASSIFICATION_ATTRIBUTES_EIGEN_H

#include <vector>

#include <CGAL/Classification/Local_eigen_analysis.h>

namespace CGAL {

namespace Classification {

namespace Attribute {

template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Eigen_attribute : public Attribute_base
{
protected:
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;

#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_ATTRIBUTES
  std::vector<double> attrib;
#else
  Local_eigen_analysis& eigen;
#endif
  
public:
  Eigen_attribute (const PointRange&,
                   Local_eigen_analysis& eigen)
#ifndef CGAL_CLASSIFICATION_PRECOMPUTE_ATTRIBUTES
    : eigen (eigen)
#endif
  {
  }

  virtual void init (const PointRange& input, Local_eigen_analysis& eigen)
  {
    this->set_weight(1.);
#ifndef CGAL_CLASSIFICATION_PRECOMPUTE_ATTRIBUTES
    std::vector<double> attrib;
#endif
    attrib.reserve (input.size());
    for (std::size_t i = 0; i < input.size(); ++ i)
      attrib.push_back (get_value (eigen, i));

    this->compute_mean_max (attrib, mean, this->max);
  }
  
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i) = 0;
  virtual double value (std::size_t pt_index)
  {
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_ATTRIBUTES
    return attrib[pt_index];
#else
    return get_value(eigen, pt_index);
#endif
  }
  virtual std::string name() { return "eigen_attribute"; }

};
  
  /*!
    \ingroup PkgClassificationAttributes

    %Attribute based on the eigenvalues of the covariance matrix of a
    local neighborhood. Linearity is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_1 - \lambda_2}{\lambda_1}
    \f]
    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Linearity : public Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits>
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  typedef Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits> Base;
public:
  /*!
    Constructs the attribute.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Linearity (const PointRange& input,
             Local_eigen_analysis& eigen) : Base (input, eigen)
  {
    this->init(input, eigen);
  }

  /// \cond SKIP_IN_MANUAL
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i)
  {
    const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[2] < 1e-15)
      return 0.;
    else
      return ((ev[2] - ev[1]) / ev[2]);
  }
  virtual std::string name() { return "linearity"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    %Attribute based on the eigenvalues of the covariance matrix of a
    local neighborhood. Planarity is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_2 - \lambda_3}{\lambda_1}
    \f]
    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Planarity : public Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits>
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  typedef Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits> Base;
public:
  /*!
    Constructs the attribute.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Planarity (const PointRange& input,
             Local_eigen_analysis& eigen)
    : Base(input, eigen)
  {
    this->init(input, eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i)
  {
    const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[2] < 1e-15)
      return 0.;
    else
      return ((ev[1] - ev[0]) / ev[2]);
  }
  virtual std::string name() { return "planarity"; }
  /// \endcond
 
};

  /*!
    \ingroup PkgClassificationAttributes

    %Attribute based on the eigenvalues of the covariance matrix of a
    local neighborhood. Sphericity is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_3}{\lambda_1}
    \f]
    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Sphericity : public Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits>
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  typedef Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits> Base;
public:
  /*!
    Constructs the attribute.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Sphericity (const PointRange& input,
              Local_eigen_analysis& eigen)
    : Base(input, eigen)
  {
    this->init(input, eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i)
  {
    const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[2] < 1e-15)
      return 0.;
    else
      return (ev[0] / ev[2]);
  }
  virtual std::string name() { return "sphericity"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    %Attribute based on the eigenvalues of the covariance matrix of a
    local neighborhood. Omnivariance is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    (\lambda_1 \times \lambda_2 \times \lambda_3)^{\frac{1}{3}}
    \f]
    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Omnivariance : public Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits>
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  typedef Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits> Base;
public:
  /*!
    Constructs the attribute.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Omnivariance (const PointRange& input,
                Local_eigen_analysis& eigen)
    : Base(input, eigen)
  {
    this->init(input, eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i)
  {
    const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    return (std::pow (std::fabs(ev[0] * ev[1] * ev[2]), 0.333333333));
  }
  virtual std::string name() { return "omnivariance"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    %Attribute based on the eigenvalues of the covariance matrix of a
    local neighborhood. Anisotropy is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_1 - \lambda_3}{\lambda_1}
    \f]
    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Anisotropy : public Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits>
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  typedef Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits> Base;
public:
  /*!
    Constructs the attribute.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Anisotropy (const PointRange& input,
              Local_eigen_analysis& eigen)
    : Base(input, eigen)
  {
    this->init(input, eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i)
  {
    const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[2] < 1e-15)
      return 0.;
    else
      return ((ev[2] - ev[0]) / ev[2]);
  }
  virtual std::string name() { return "anisotropy"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    %Attribute based on the eigenvalues of the covariance matrix of a
    local neighborhood. Eigentropy is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    - \sum_{i=1}^3 \lambda_i \times \log{\lambda_i}
    \f]
    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Eigentropy : public Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits>
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  typedef Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits> Base;
public:
  /*!
    Constructs the attribute.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Eigentropy (const PointRange& input,
              Local_eigen_analysis& eigen)
    : Base(input, eigen)
  {
    this->init(input, eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i)
  {
    const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[0] < 1e-15
        || ev[1] < 1e-15
        || ev[2] < 1e-15)
      return 0.;
    else
      return (- ev[0] * std::log(ev[0])
              - ev[1] * std::log(ev[1])
              - ev[2] * std::log(ev[2]));
  }
  virtual std::string name() { return "eigentropy"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    %Attribute based on the eigenvalues of the covariance matrix of a
    local neighborhood. The sum of the eigenvalues is defined, for the
    3 eigenvalues \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$,
    as:

    \f[
    \lambda_1 + \lambda_2 + \lambda_3
    \f]
    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Sum_eigenvalues : public Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits>
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  typedef Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits> Base;
public:
  /*!
    Constructs the attribute.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Sum_eigenvalues (const PointRange& input,
                   Local_eigen_analysis& eigen)
    : Base(input, eigen)
  {
    this->init(input, eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i)
  {
    return eigen.sum_of_eigenvalues(i);
  }
  virtual std::string name() { return "sum_eigen"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    %Attribute based on the eigenvalues of the covariance
    matrix of a local neighborhood. Surface variation is defined, for
    the 3 eigenvalues \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge
    0\f$, as:

    \f[
    \frac{\lambda_3}{\lambda_1 + \lambda_2 + \lambda_3}
    \f]
    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Surface_variation : public Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits>
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  typedef Eigen_attribute<Geom_traits, PointRange, PointMap, DiagonalizeTraits> Base;
public:
  /*!
    Constructs the attribute.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Surface_variation (const PointRange& input,
                     Local_eigen_analysis& eigen)
    : Base(input, eigen)
  {
    this->init(input, eigen);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double get_value (Local_eigen_analysis& eigen, std::size_t i)
  {
    const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
    if (ev[0] + ev[1] + ev[2] < 1e-15)
      return 0.;
    else
      return (ev[0] / (ev[0] + ev[1] + ev[2]));
  }
  virtual std::string name() { return "surface_variation"; }
  /// \endcond
};

} // namespace Attribute

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTES_EIGEN_H
