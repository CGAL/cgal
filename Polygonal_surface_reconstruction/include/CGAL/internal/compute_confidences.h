// Copyright (c) 2018  Liangliang Nan. All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Liangliang Nan

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_CANDIDATE_CONFIDENCES_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_CANDIDATE_CONFIDENCES_H

#include <CGAL/license/Polygonal_surface_reconstruction.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/assertions.h>
#include <CGAL/internal/point_set_with_planes.h>
#include <CGAL/internal/alpha_shape_mesh.h>
#include <CGAL/internal/hypothesis.h>
#include <CGAL/internal/parameters.h>


// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

namespace CGAL {
namespace internal {

                /**
                *        Computes the confidences of the candidate faces.
                */

                template <typename Kernel>
                class Candidate_confidences
                {
                private:
                        typedef typename Kernel::FT                                FT;
                        typedef typename Kernel::Point_3                        Point;
                        typedef typename Kernel::Point_2                        Point2;
                        typedef typename Kernel::Vector_3                        Vector;
                        typedef typename Kernel::Line_3                                Line;
                        typedef typename Kernel::Segment_3                        Segment;
                        typedef typename Kernel::Plane_3                        Plane;
                        typedef CGAL::Polygon_2<Kernel>                                Polygon;
                        typedef internal::Planar_segment<Kernel>                Planar_segment;
                        typedef internal::Point_set_with_planes<Kernel>                Point_set;
                        typedef CGAL::Surface_mesh<Point>                        Polygon_mesh;
                        typedef typename Polygon_mesh::Face_index                Face_descriptor;
                        typedef typename Polygon_mesh::Edge_index                Edge_descriptor;
                        typedef typename Polygon_mesh::Vertex_index                Vertex_descriptor;
                        typedef typename Polygon_mesh::Halfedge_index                Halfedge_descriptor;

                public:
                        Candidate_confidences() {}
                        ~Candidate_confidences() {}

                        /// Computes the confidence values for each face
                        /// - supporting point number:        stored as property 'f:num_supporting_points'
                        /// - face area:                                stored as property 'f:face_area'
                        /// - covered area:                                stored as property 'f:covered_area'
                        void compute(const Point_set& point_set, Polygon_mesh& mesh);

                private:
                        // Returns the indices of the supporting point for 'face'
                        std::vector<std::size_t> supporting_points(Face_descriptor face, const Polygon_mesh& mesh, const Point_set& point_set);

                        inline FT triangle_area(const Point& p1, const Point& p2, const Point& p3) const {
                                const Vector& orth = CGAL::cross_product(Vector(p1, p2), Vector(p1, p3));
                                return std::sqrt(orth.squared_length()) * FT(0.5);
                        }

                        FT face_area(Face_descriptor f, const Polygon_mesh& mesh) const {
                                FT result(0);

                                const typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = mesh.points();

                                Halfedge_around_face_circulator<Polygon_mesh> cir(mesh.halfedge(f), mesh), done(cir);
                                Halfedge_descriptor p_hd = *cir;
                                Vertex_descriptor p_vd = mesh.target(p_hd);
                                const Point& p = coords[p_vd];
                                ++cir;

                                do {
                                        Halfedge_descriptor q_hd = *cir;
                                        Vertex_descriptor q_vd = mesh.target(q_hd);
                                        const Point& q = coords[q_vd];

                                        Halfedge_descriptor r_hd = mesh.next(q_hd);
                                        Vertex_descriptor r_vd = mesh.target(r_hd);
                                        const Point& r = coords[r_vd];

                                        result += triangle_area(p, q, r);

                                        ++cir;
                                } while (cir != done);

                                return result;
                        }
                };


                //////////////////////////////////////////////////////////////////////////

                // implementation

                template <typename Kernel>
                std::vector<std::size_t> Candidate_confidences<Kernel>::supporting_points(Face_descriptor face, const Polygon_mesh& mesh, const Point_set& point_set) {
                        std::vector<std::size_t> indices;

                        if (face == Polygon_mesh::null_face())
                                return indices;

                        // The supporting planar segment of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, Planar_segment*> face_supporting_segments =
                                mesh.template property_map<Face_descriptor, Planar_segment*>("f:supp_segment").first;

                        Planar_segment* segment = face_supporting_segments[face];
                        if (segment == nullptr)
                                return indices;

                        // The supporting plane of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
                                mesh.template property_map<Face_descriptor, const Plane*>("f:supp_plane").first;
                        // We do everything by projecting the point onto the face's supporting plane
                        const Plane* supporting_plane = face_supporting_planes[face];
                        CGAL_assertion(supporting_plane == segment->supporting_plane());

                        Polygon plg; // The projection of the face onto it supporting plane
                        const typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = mesh.points();
                        Halfedge_around_face_circulator<Polygon_mesh> cir(mesh.halfedge(face), mesh), done(cir);
                        do {
                                Halfedge_descriptor hd = *cir;
                                Vertex_descriptor vd = mesh.target(hd);
                                const Point& p = coords[vd];
                                const Point2& q = supporting_plane->to_2d(p);

                                // Removes duplicated vertices
                                // The last point in the polygon
                                if (!plg.is_empty()) {
                                        const Point2& r = plg[plg.size() - 1];
                                        if (CGAL::squared_distance(q, r) < CGAL::snap_squared_distance_threshold<FT>())
                                                continue;
                                }
                                plg.push_back(q);

                                ++cir;
                        } while (cir != done);

                        if (plg.size() < 3 || !plg.is_simple())
                                return indices;

                        const typename Point_set::Point_map& points = point_set.point_map();
                        for (std::size_t i = 0; i < segment->size(); ++i) {
                                std::size_t idx = segment->at(i);
                                const Point& p = points[idx];
                                if (plg.bounded_side(supporting_plane->to_2d(p)) == CGAL::ON_BOUNDED_SIDE)
                                        indices.push_back(idx);
                        }

                        return indices;
                }

                template <typename Kernel>
                void Candidate_confidences<Kernel>::compute(const Point_set& point_set, Polygon_mesh& mesh) {
                        const unsigned int K = 6;

                        const typename Point_set::Point_map& points = point_set.point_map();
                        FT avg_spacing = compute_average_spacing<Concurrency_tag>(points, K);

                        // The number of supporting points of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, std::size_t> face_num_supporting_points =
                                mesh.template add_property_map<Face_descriptor, std::size_t>("f:num_supporting_points").first;

                        // The area of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, FT> face_areas =
                                mesh.template add_property_map<Face_descriptor, FT>("f:face_area").first;

                        // The point covered area of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, FT> face_covered_areas =
                                mesh.template add_property_map<Face_descriptor, FT>("f:covered_area").first;

                        // The supporting plane of each face
                        typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
                                mesh.template property_map<Face_descriptor, const Plane*>("f:supp_plane").first;

                        FT degenerate_face_area_threshold = CGAL::snap_squared_distance_threshold<FT>() * CGAL::snap_squared_distance_threshold<FT>();

                        for(auto f : mesh.faces()) {
                                const Plane* supporting_plane = face_supporting_planes[f];
                                // Face area
                                FT area = face_area(f, mesh);
                                face_areas[f] = area;

                                if (area > degenerate_face_area_threshold) {
                                        const std::vector<std::size_t>& indices = supporting_points(f, mesh, point_set);
                                        face_num_supporting_points[f] = indices.size();

                                        std::vector<Point> pts;
                                        for (std::size_t i = 0; i < indices.size(); ++i) {
                                                std::size_t idx = indices[i];
                                                const Point& p = points[idx];
                                                pts.push_back(p);
                                        }

                                        FT covered_area(0);
                                        Alpha_shape_mesh<Kernel> alpha_mesh(pts.begin(), pts.end(), *supporting_plane);
                                        Polygon_mesh covering_mesh;
                                        FT radius = avg_spacing * FT(5.0);
                                        if (alpha_mesh.extract_mesh(radius * radius, covering_mesh)) {
                                                // We cannot use the area of the 3D faces, because the alpha shape mesh is
                                                // not perfectly planar
                                                const typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = covering_mesh.points();
                                                for(auto face : covering_mesh.faces()) {
                                                        // We have to use the projected version
                                                        Polygon plg; // the projection of the face onto it supporting plane
                                                        Halfedge_around_face_circulator<Polygon_mesh> cir(covering_mesh.halfedge(face), covering_mesh), done(cir);
                                                        do {
                                                                Halfedge_descriptor hd = *cir;
                                                                Vertex_descriptor vd = covering_mesh.target(hd);
                                                                const Point& p = coords[vd];
                                                                const Point2& q = supporting_plane->to_2d(p);
                                                                plg.push_back(q);
                                                                ++cir;
                                                        } while (cir != done);
                                                        covered_area += std::abs(plg.area());
                                                }
                                        }

                                        face_covered_areas[f] = covered_area;
                                        if (covered_area > area)
                                                face_covered_areas[f] = area;
                                }
                                else { // For tiny faces, we can simple assign zero supporting points
                                        face_num_supporting_points[f] = 0;
                                        face_covered_areas[f] = FT(0.0);
                                }
                        }
                }

        } //namespace internal

} //namespace CGAL


#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_CANDIDATE_CONFIDENCES_H
