// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
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

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on the eigenvalues of the covariance matrix
    of a local neighborhood.

    Linearity is defined, for the 3 eigenvalues \f$\lambda_1 \ge
    \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_1 - \lambda_2}{\lambda_1}
    \f]
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Linearity : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Linearity (RandomAccessIterator begin,
             RandomAccessIterator end,
             Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back ((ev[2] - ev[1]) / ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "linearity"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on the eigenvalues of the covariance matrix
    of a local neighborhood.

    Planarity is defined, for the 3 eigenvalues \f$\lambda_1 \ge
    \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_2 - \lambda_3}{\lambda_1}
    \f]
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Planarity : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Planarity (RandomAccessIterator begin,
             RandomAccessIterator end,
             Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back ((ev[1] - ev[0]) / ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "planarity"; }
  /// \endcond
 
};

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on the eigenvalues of the covariance matrix
    of a local neighborhood.

    Sphericity is defined, for the 3 eigenvalues \f$\lambda_1 \ge
    \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_3}{\lambda_1}
    \f]
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Sphericity : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Sphericity (RandomAccessIterator begin,
              RandomAccessIterator end,
              Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back (ev[0] / ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "sphericity"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on the eigenvalues of the covariance matrix
    of a local neighborhood.

    Omnivariance is defined, for the 3 eigenvalues \f$\lambda_1 \ge
    \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    (\lambda_1 \times \lambda_2 \times \lambda_3)^{\frac{1}{3}}
    \f]
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Omnivariance : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Omnivariance (RandomAccessIterator begin,
                RandomAccessIterator end,
                Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back (std::pow (std::fabs(ev[0] * ev[1] * ev[2]), 0.333333333));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "omnivariance"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on the eigenvalues of the covariance matrix
    of a local neighborhood.

    Anisotropy is defined, for the 3 eigenvalues \f$\lambda_1 \ge
    \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_1 - \lambda_3}{\lambda_1}
    \f]
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Anisotropy : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Anisotropy (RandomAccessIterator begin,
              RandomAccessIterator end,
              Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back ((ev[2] - ev[0]) / ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "anisotropy"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on the eigenvalues of the covariance matrix
    of a local neighborhood.

    Eigentropy is defined, for the 3 eigenvalues \f$\lambda_1 \ge
    \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    - \sum_{i=1}^3 \lambda_i \times \log{\lambda_i}
    \f]
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Eigentropy : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Eigentropy (RandomAccessIterator begin,
              RandomAccessIterator end,
              Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[0] < 1e-15
            || ev[1] < 1e-15
            || ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back (- ev[0] * std::log(ev[0])
                            - ev[1] * std::log(ev[1])
                            - ev[2] * std::log(ev[2]));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "eigentropy"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on the eigenvalues of the covariance matrix
    of a local neighborhood.

    The sum of the eigenvalues is defined, for the 3 eigenvalues
    \f$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \lambda_1 + \lambda_2 + \lambda_3
    \f]
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Sum_eigenvalues : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Sum_eigenvalues (RandomAccessIterator begin,
                   RandomAccessIterator end,
                   Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      attrib.push_back (eigen.sum_of_eigenvalues(i));

    this->compute_mean_max (attrib, mean, this->max);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "sum_eigen"; }
  /// \endcond
};

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on the eigenvalues of the covariance matrix
    of a local neighborhood.

    Surface variation is defined, for the 3 eigenvalues \f$\lambda_1
    \ge \lambda_2 \ge \lambda_3 \ge 0\f$, as:

    \f[
    \frac{\lambda_3}{\lambda_1 + \lambda_2 + \lambda_3}
    \f]
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Surface_variation : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Surface_variation (RandomAccessIterator begin,
                     RandomAccessIterator end,
                     Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);

        if (ev[0] + ev[1] + ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back (ev[0] / (ev[0] + ev[1] + ev[2]));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "surface_variation"; }
  /// \endcond
};

} // namespace Attribute

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTES_EIGEN_H
