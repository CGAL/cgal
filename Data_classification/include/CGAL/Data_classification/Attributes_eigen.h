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

#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTES_EIGEN_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTES_EIGEN_H

#include <vector>

#include <CGAL/Data_classification/Local_eigen_analysis.h>

namespace CGAL {

namespace Data_classification {


  /*!
    \ingroup PkgDataClassification

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
class Attribute_linearity : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_linearity (RandomAccessIterator begin,
                       RandomAccessIterator end,
                       Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
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
    \ingroup PkgDataClassification

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
class Attribute_planarity : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_planarity (RandomAccessIterator begin,
                       RandomAccessIterator end,
                       Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
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
    \ingroup PkgDataClassification

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
class Attribute_sphericity : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_sphericity (RandomAccessIterator begin,
                        RandomAccessIterator end,
                        Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
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
    \ingroup PkgDataClassification

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
class Attribute_omnivariance : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_omnivariance (RandomAccessIterator begin,
                          RandomAccessIterator end,
                          Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
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
    \ingroup PkgDataClassification

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
class Attribute_anisotropy : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_anisotropy (RandomAccessIterator begin,
                        RandomAccessIterator end,
                        Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
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
    \ingroup PkgDataClassification

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
class Attribute_eigentropy : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_eigentropy (RandomAccessIterator begin,
                        RandomAccessIterator end,
                        Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
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
    \ingroup PkgDataClassification

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
class Attribute_sum_eigenvalues : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_sum_eigenvalues (RandomAccessIterator begin,
                             RandomAccessIterator end,
                             Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      attrib.push_back (eigen.sum_eigenvalues(i));

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
    \ingroup PkgDataClassification

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
class Attribute_surface_variation : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_surface_variation (RandomAccessIterator begin,
                               RandomAccessIterator end,
                               Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
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


} // namespace Data_classification

} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTES_EIGEN_H
