// Copyright (c) 2022 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Daniel Zint
//                 Julian Stahl

#ifndef CGAL_DUAL_CONTOURING_3_INTERNAL_DUAL_CONTOURING_3_H
#define CGAL_DUAL_CONTOURING_3_INTERNAL_DUAL_CONTOURING_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Origin.h>
#include <CGAL/centroid.h>

#include <Eigen/SVD>
#include <array>
#include <map>
#include <mutex>
#include <vector>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

namespace Positioning {
/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Computes the vertex position for a point in Dual Contouring using Quadric Error Metrics and the SVD pseudo
 * inverse.
 *
 * \details
 *
 * \tparam use_bbox clamp vertex position to the bounding box of the cell
 *
 */
template <bool use_bbox = false>
class QEM_SVD {
public:
    /**
     * \ingroup PkgIsosurfacing3Ref
     *
     * \brief Computes the vertex position for a point in Dual Contouring.
     *
     * \details
     *
     * \tparam Domain_ must be a model of `IsosurfacingDomainWithGradient`.
     *
     * \param domain the domain providing input data and its topology
     * \param iso_value value of the isosurface
     * \param cell the cell within the domain for which the vertex position ins computed
     * \param point the point position of the vertex that belongs to that cell
     *
     * \return true, if the voxel intersects the isosurface
     */
    template <class Domain_>
    bool position(const Domain_& domain, const typename Domain_::FT iso_value,
                  const typename Domain_::Cell_descriptor& cell, typename Domain_::Point& point) const {
        typedef typename Domain_::Point Point;
        typedef typename Domain_::Geom_traits::Vector_3 Vector;
        typedef typename Domain_::FT FT;

        typename Domain_::Cell_vertices vertices = domain.cell_vertices(cell);

        std::vector<Point> pos(vertices.size());
        std::transform(vertices.begin(), vertices.end(), pos.begin(),
                       [&](const auto& v) { return domain.position(v); });

        point = CGAL::centroid(pos.begin(), pos.end(), CGAL::Dimension_tag<0>());  // set point to cell center

        // compute edge intersections
        std::vector<Point> edge_intersections;
        std::vector<Vector> edge_intersection_normals;

        for (const auto& edge : domain.cell_edges(cell)) {
            const auto& edge_vertices = domain.edge_vertices(edge);
            const auto& v0 = edge_vertices[0];
            const auto& v1 = edge_vertices[1];

            const auto& val0 = domain.value(v0);
            const auto& val1 = domain.value(v1);

            const auto& p0 = domain.position(v0);
            const auto& p1 = domain.position(v1);

            if ((val0 <= iso_value) != (val1 <= iso_value)) {  // this edge is intersected by the isosurface
                const FT u = (val0 - iso_value) / (val0 - val1);
                const Point p_lerp = CGAL::ORIGIN + ((1 - u) * (p0 - CGAL::ORIGIN) + u * (p1 - CGAL::ORIGIN));
                edge_intersections.push_back(p_lerp);
                edge_intersection_normals.push_back(domain.gradient(p_lerp));
            }
        }

        if (edge_intersections.empty()) {
            return false;
        }

        // SVD QEM
        Eigen::Matrix3d A;
        A.setZero();
        Eigen::Vector3d rhs;
        rhs.setZero();
        for (std::size_t i = 0; i < edge_intersections.size(); ++i) {
            Eigen::Vector3d n_k = {edge_intersection_normals[i].x(), edge_intersection_normals[i].y(),
                                   edge_intersection_normals[i].z()};
            Eigen::Vector3d p_k = {edge_intersections[i].x(), edge_intersections[i].y(), edge_intersections[i].z()};
            double d_k = n_k.transpose() * p_k;

            Eigen::Matrix3d A_k = n_k * n_k.transpose();
            Eigen::Vector3d b_k = d_k * n_k;
            A += A_k;
            rhs += b_k;
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        // set threshold as in Peter Lindstrom's paper, "Out-of-Core
        // Simplification of Large Polygonal Models"
        svd.setThreshold(1e-3);

        // Init x hat
        Eigen::Vector3d x_hat;
        x_hat << point.x(), point.y(), point.z();

        // Lindstrom formula for QEM new position for singular matrices
        Eigen::VectorXd v_svd = x_hat + svd.solve(rhs - A * x_hat);

        point = Point(v_svd[0], v_svd[1], v_svd[2]);

        // bbox
        if (use_bbox) {
            CGAL::Bbox_3 bbox = pos[0].bbox() + pos[7].bbox();  // TODO remove[0],[7]

            FT x = (std::min<FT>)((std::max<FT>)(point.x(), bbox.xmin()), bbox.xmax());
            FT y = (std::min<FT>)((std::max<FT>)(point.y(), bbox.ymin()), bbox.ymax());
            FT z = (std::min<FT>)((std::max<FT>)(point.z(), bbox.zmin()), bbox.zmax());
            point = Point(x, y, z);
        }

        return true;
    }
};

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Returns cell center.
 */
class Cell_center {
public:
    /**
     * \ingroup PkgIsosurfacing3Ref
     *
     * \brief Computes the vertex position for a point in Dual Contouring.
     *
     * \details
     *
     * \tparam Domain_ must be a model of `IsosurfacingDomainWithGradient`.
     *
     * \param domain the domain providing input data and its topology
     * \param iso_value value of the isosurface
     * \param cell the cell within the domain for which the vertex position ins computed
     * \param point the point position of the vertex that belongs to that cell
     *
     * \return true, if the voxel intersects the isosurface
     */
    template <class Domain_>
    bool position(const Domain_& domain, const typename Domain_::FT iso_value,
                  const typename Domain_::Cell_descriptor& vh, typename Domain_::Point& point) const {
        typedef typename Domain_::Point Point;
        typedef typename Domain_::Geom_traits::Vector_3 Vector;
        typedef typename Domain_::FT FT;

        typename Domain_::Cell_vertices vertices = domain.cell_vertices(vh);

        std::vector<Point> pos(vertices.size());
        std::transform(vertices.begin(), vertices.end(), pos.begin(),
                       [&](const auto& v) { return domain.position(v); });

        point = CGAL::centroid(pos.begin(), pos.end(), CGAL::Dimension_tag<0>());  // set point to cell center

        bool allSmaller = true;
        bool allGreater = true;
        for (const auto& v : vertices) {
            const bool& b = domain.value(v) <= iso_value;
            allSmaller = allSmaller && b;
            allGreater = allGreater && !b;
        }

        if (allSmaller || allGreater) {
            return false;
        }

        return true;
    }
};

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Computes the centroid of all cell edge intersections with the isosurface.
 */
class Centroid_of_edge_intersections {
public:
    /**
     * \ingroup PkgIsosurfacing3Ref
     *
     * \brief Computes the vertex position for a point in Dual Contouring.
     *
     * \details
     *
     * \tparam Domain_ must be a model of `IsosurfacingDomainWithGradient`.
     *
     * \param domain the domain providing input data and its topology
     * \param iso_value value of the isosurface
     * \param cell the cell within the domain for which the vertex position ins computed
     * \param point the point position of the vertex that belongs to that cell
     *
     * \return true, if the voxel intersects the isosurface
     */
    template <class Domain_>
    bool position(const Domain_& domain, const typename Domain_::FT iso_value,
                  const typename Domain_::Cell_descriptor& cell, typename Domain_::Point& point) const {
        typedef typename Domain_::Point Point;
        typedef typename Domain_::Geom_traits::Vector_3 Vector;
        typedef typename Domain_::FT FT;

        typename Domain_::Cell_vertices vertices = domain.cell_vertices(cell);

        // compute edge intersections
        std::vector<Point> edge_intersections;

        for (const auto& edge : domain.cell_edges(cell)) {
            const auto& edge_vertices = domain.edge_vertices(edge);
            const auto& v0 = edge_vertices[0];
            const auto& v1 = edge_vertices[1];

            const auto& val0 = domain.value(v0);
            const auto& val1 = domain.value(v1);

            const auto& p0 = domain.position(v0);
            const auto& p1 = domain.position(v1);

            if ((val0 <= iso_value) != (val1 <= iso_value)) {  // this edge is intersected by the isosurface
                const FT u = (val0 - iso_value) / (val0 - val1);
                const Point p_lerp = CGAL::ORIGIN + ((1 - u) * (p0 - CGAL::ORIGIN) + u * (p1 - CGAL::ORIGIN));
                edge_intersections.push_back(p_lerp);
            }
        }

        if (edge_intersections.empty()) {
            return false;
        }

        point = CGAL::centroid(edge_intersections.begin(), edge_intersections.end(),
                               CGAL::Dimension_tag<0>());  // set point to center of edge intersections

        return true;
    }
};
}  // namespace Positioning

template <class Domain_, class Positioning_>
class Dual_contouring_vertex_positioning {
private:
    typedef Domain_ Domain;
    typedef Positioning_ Positioning;

    typedef typename Domain::FT FT;
    typedef typename Domain::Point Point;
    typedef typename Domain::Cell_descriptor Cell_descriptor;

public:
    Dual_contouring_vertex_positioning(const Domain& domain, FT iso_value, const Positioning& positioning)
        : domain(domain), iso_value(iso_value), positioning(positioning), points_counter(0) {}

    void operator()(const Cell_descriptor& v) {
        // compute dc-vertices
        Point p;
        if (positioning.position(domain, iso_value, v, p)) {

            std::lock_guard<std::mutex> lock(mutex);
            map_voxel_to_point[v] = p;
            map_voxel_to_point_id[v] = points_counter++;
        }
    }

    // private:
    const Domain& domain;
    FT iso_value;
    const Positioning& positioning;

    std::map<Cell_descriptor, std::size_t> map_voxel_to_point_id;
    std::map<Cell_descriptor, Point> map_voxel_to_point;
    std::size_t points_counter;

    std::mutex mutex;
};

template <class Domain_>
class Dual_contouring_face_generation {
private:
    typedef Domain_ Domain;

    typedef typename Domain_::FT FT;
    typedef typename Domain_::Edge_descriptor Edge_descriptor;
    typedef typename Domain_::Cell_descriptor Cell_descriptor;

public:
    Dual_contouring_face_generation(const Domain& domain, FT iso_value) : domain(domain), iso_value(iso_value) {}

    void operator()(const Edge_descriptor& e) {
        // save all faces
        const auto& vertices = domain.edge_vertices(e);
        const FT s0 = domain.value(vertices[0]);
        const FT s1 = domain.value(vertices[1]);

        if (s0 <= iso_value && s1 > iso_value) {
            const auto& voxels = domain.cells_incident_to_edge(e);

            std::lock_guard<std::mutex> lock(mutex);
            faces[e].insert(faces[e].begin(), voxels.begin(), voxels.end());

        } else if (s1 <= iso_value && s0 > iso_value) {
            const auto& voxels = domain.cells_incident_to_edge(e);

            std::lock_guard<std::mutex> lock(mutex);
            faces[e].insert(faces[e].begin(), voxels.rbegin(), voxels.rend());
        }
    }

    // private:
    std::map<Edge_descriptor, std::vector<Cell_descriptor>> faces;

    const Domain& domain;
    FT iso_value;

    std::mutex mutex;
};

}  // namespace internal
}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_DUAL_CONTOURING_3_INTERNAL_DUAL_CONTOURING_3_H
