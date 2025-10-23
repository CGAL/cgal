// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>

#ifndef CGAL_HDVF_TRAITS_D_H
#define CGAL_HDVF_TRAITS_D_H

#include <CGAL/license/HDVF.h>
#include <CGAL/Bbox_d.h>
#include <CGAL/Dimension.h>
#include <functional>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <Eigen/SVD>

namespace CGAL {
namespace Homological_discrete_vector_field {

  /*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Hdvf_traits_d` implements the `HDVFTraits` concept for dD data, using a geometric kernel `K`.

 @tparam K a geometric traits class. Must be either `Epick_d` or `Epeck_d` with a fixed dimension.

 \cgalModels{HDVFTraits}

 */

template <typename K>
struct Hdvf_traits_d {
    using Dimension = Dimension_tag< K::Dimension::value >;
    typedef K Kernel;
    typedef typename K::Point_d Point;
    typedef typename K::FT FT;
    typedef CGAL::Bbox_d<Dimension> Bbox;
    typedef typename Exact_predicates_inexact_constructions_kernel::Point_3 Point3;
    static std::function<Point3(const Point&)> to_point3;
    
    // Set of standard projection operators
    // R^d -> R^3 projection on first coordinates (x,y,z)
    static std::function<Point3(const Point&)> default_projection ;
    
    // Projection on a given 3D affine plane (defined by three unit orthogonal vectors + point)
    static std::function<typename Exact_predicates_inexact_constructions_kernel::Point_3(const typename K::Point_d&)> plane_projection_builder(const typename K::Point_d& p0, const typename K::Vector_d& d1, const typename K::Vector_d& d2, const typename K::Vector_d& d3) {
        typedef Exact_predicates_inexact_constructions_kernel::Point_3 Point_3;
        typedef typename K::Point_d Point_d;
        typedef typename K::Vector_d Vector_d;
        
        std::function<Point_3(const Point_d&)> project = [&](const Point_d& p) {
            Vector_d tmp(p-p0);
            return Point_3(tmp*d1, tmp*d2, tmp*d3);
        };
        return project;
    }
    // Define a method computing the PCA projection frame
    static std::function<typename Exact_predicates_inexact_constructions_kernel::Point_3(const typename K::Point_d&)> pca_frame_builder(const std::vector<typename K::Point_d>& pts) {
        typedef Exact_predicates_inexact_constructions_kernel::Point_3 Point_3;
        typedef typename K::Point_d Point_d;
        typedef typename K::Vector_d Vector_d;
        typedef double FT;
        typedef Eigen_vector<FT> Vector;
        typedef Eigen_matrix<FT> Matrix;
        
        // Barycenter of points
        Vector_d bary;
        for (Point p : pts)
            bary += (p-ORIGIN) ;
        Point_d barycenter = ORIGIN + bary/pts.size() ;
        // Define the PCA matrix by SVD
        Matrix A(pts.size(), K::Dimension::value);
        // Fill the matrix (p[i]-barycenter along line i)
        for (int i=0; i<pts.size(); ++i) {
            Vector_d tmp(pts.at(i)-barycenter);
            for (int j=0; j<K::Dimension::value; ++j) {
                A.set(i,j,tmp[j]);
            }
        }
        // Compute SVD decomposition
#if EIGEN_VERSION_AT_LEAST(3,4,90)
        Eigen::JacobiSVD<Matrix::EigenType, Eigen::ComputeThinU | Eigen::ComputeThinV> jacobiSvd(A.eigen_object());
#else
        Eigen::JacobiSVD<Matrix::EigenType> jacobiSvd(A.eigen_object(), ::Eigen::ComputeThinU | ::Eigen::ComputeThinV);
#endif
        std::vector<std::vector<double>> d(3);
        std::vector<Vector_d> dirs(3);
        // Fill d[i] vector
        for (int i=0; i<3; ++i) {
            d.at(i).resize(K::Dimension::value);
            for (int j=0; j<K::Dimension::value; ++j) {
                d.at(i).at(j) = jacobiSvd.matrixV()(j,i);
            }
            dirs.at(i) = Vector_d(d.at(i).begin(), d.at(i).end());
        }
        return plane_projection_builder(barycenter, dirs.at(0), dirs.at(1), dirs.at(2));
    }
};

template <typename K>
std::function<typename Exact_predicates_inexact_constructions_kernel::Point_3(const typename K::Point_d&)> Hdvf_traits_d<K>::to_point3 = Hdvf_traits_d<K>::default_projection;

template <typename K>
std::function<typename Exact_predicates_inexact_constructions_kernel::Point_3(const typename K::Point_d&)> Hdvf_traits_d<K>::default_projection =
    [](const K::Point_d& p) { return typename Exact_predicates_inexact_constructions_kernel::Point_3(p[0], p[1], p[2]); };




} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */

#endif // CGAL_HDVF_TRAITS_D_H
