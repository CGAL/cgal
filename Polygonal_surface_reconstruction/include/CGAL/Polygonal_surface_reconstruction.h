// Copyright (c) 2018  Liangliang Nan. All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Liangliang Nan

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Polygonal_surface_reconstruction.h>

#include <CGAL/bounding_box.h>
#include <CGAL/property_map.h>
#include <CGAL/internal/hypothesis.h>
#include <CGAL/internal/compute_confidences.h>
#include <CGAL/internal/point_set_with_planes.h>

#include <unordered_map>

/*!
\file Polygonal_surface_reconstruction.h
*/

namespace CGAL {

        /*!
        \ingroup PkgPolygonalSurfaceReconstruction

        \brief

        Implementation of the Polygonal Surface Reconstruction method.

        Given a set of 3D points with unoriented normals sampled
        from the outer boundary of a piecewise planar object, the Polygonal Surface
        Reconstruction method \cgalCite{nan2017polyfit} outputs a simplified and
        watertight surface mesh interpolating the input point set.

        The method first generates a set of face candidates by intersecting the planar
        primitives. Then an optimal subset of the candidate faces is selected through
        optimization under hard constraints that enforce the final model to be manifold
        and watertight.

        The reconstruction assumes the planar segmentation of the point cloud is
        provided in the input.

        \tparam GeomTraits a geometric traits class, model of Kernel
        */
        template <class GeomTraits>
        class Polygonal_surface_reconstruction
        {
        public:

                /// \name Types

                typedef typename GeomTraits::FT                                FT;                        ///< number type.
                typedef typename GeomTraits::Point_3                Point;                ///< point type.
                typedef typename GeomTraits::Vector_3                Vector;                ///< vector type.
                typedef typename GeomTraits::Plane_3                Plane;                ///< plane type.

        private:

                // Polygon mesh storing the candidate faces
                typedef CGAL::Surface_mesh<Point>                                Polygon_mesh;

                typedef typename Polygon_mesh::Face_index                Face_descriptor;
                typedef typename Polygon_mesh::Vertex_index                Vertex_descriptor;
                typedef typename Polygon_mesh::Edge_index                Edge_descriptor;
                typedef typename Polygon_mesh::Halfedge_index        Halfedge_descriptor;

                // Public methods
        public:

                /// \name Creation

                /*!
                Creates a Polygonal Surface Reconstruction object.
                After construction, candidate faces are generated and point/face confidence values are
                computed, allowing to reuse them in the subsequent reconstruction step with different parameters.

                \tparam PointRange is the range of input points, model of `ConstRange`.
                \tparam PointMap is a model of `ReadablePropertyMap` with value        type `GeomTraits::Point_3`.
                \tparam NormalMap is a model of `ReadablePropertyMap` with value type `GeomTraits::Vector_3`.
                \tparam IndexMap is a model of `ReadablePropertyMap` with value        type `int`.

                \param points range of input points.
                \param point_map property map: value_type of `typename PointRange::const_iterator` -> `Point_3`
                \param normal_map property map: value_type of `typename PointRange::const_iterator` -> `Vector_3`
                \param index_map property map: value_type of `typename PointRange::const_iterator` -> `int`,
                denoting the index of the plane it belongs to (-1 if the point is not assigned to a plane)
                */
                template <
                        typename PointRange,
                        typename PointMap,
                        typename NormalMap,
                        typename IndexMap
                >
                        Polygonal_surface_reconstruction(
                                const PointRange& points,
                                PointMap point_map,
                                NormalMap normal_map,
                                IndexMap index_map
                        );


                /// \name Operations

                /** Reconstructs a watertight polygonal mesh model.

                \tparam MixedIntegerProgramTraits a model of `MixedIntegerProgramTraits`

                \tparam PolygonMesh a model of `MutableFaceGraph`

                \return `true` if the reconstruction succeeded, `false` otherwise.
                */
                template <typename MixedIntegerProgramTraits, typename PolygonMesh>
                bool reconstruct(
                        PolygonMesh& output_mesh,                ///< the final reconstruction result
                        double wt_fitting = 0.43,                ///< weight for the data fitting term.
                        double wt_coverage = 0.27,                ///< weight for the point coverage term.
                        double wt_complexity = 0.30                ///< weight for the model complexity term.
                );

                /*! Gives the user the possibility to access the intermediate candidate faces
                        (i.e., the faces induced by the intersection of the supporting planes).
                \tparam PolygonMesh a model of `MutableFaceGraph`.
                */
                template <typename PolygonMesh>
                void output_candidate_faces(PolygonMesh& candidate_faces) const;

                /// Gets the error message (if reconstruction failed).
                const std::string& error_message() const { return error_message_; }

                // Data members.
        private:
                internal::Hypothesis<GeomTraits> hypothesis_;

                // The generated candidate faces stored as a polygon mesh
                Polygon_mesh        candidate_faces_;

                std::string                error_message_;

        private: // Copying is not allowed
                Polygonal_surface_reconstruction(const Polygonal_surface_reconstruction& psr);

        }; // end of Polygonal_surface_reconstruction


        //////////////////////////////////////////////////////////////////////////s

        // implementations

        template <class GeomTraits>

        template <
                typename PointRange,
                typename PointMap,
                typename NormalMap,
                typename IndexMap
        >
                Polygonal_surface_reconstruction<GeomTraits>::Polygonal_surface_reconstruction(
                        const PointRange& points,
                        PointMap point_map,
                        NormalMap normal_map,
                        IndexMap index_map
                ) : error_message_("")
        {
                if (points.empty()) {
                        error_message_ = "empty input points";
                        return;
                }

                typedef internal::Planar_segment<GeomTraits>                        Planar_segment;
                typedef internal::Point_set_with_planes<GeomTraits>                Point_set_with_planes;

                Point_set_with_planes point_set(points, point_map, normal_map, index_map);

                const std::vector< Planar_segment* >& planar_segments = point_set.planar_segments();
                if (planar_segments.size() < 4) {
                        error_message_ = "at least 4 planes required to reconstruct a closed surface mesh (only "
                                + std::to_string(planar_segments.size()) + " provided)";
                        return;
                }

                hypothesis_.generate(point_set, candidate_faces_);

                typedef internal::Candidate_confidences<GeomTraits>                Candidate_confidences;
                Candidate_confidences conf;
                conf.compute(point_set, candidate_faces_);
        }


        template <class GeomTraits>
        template <typename PolygonMesh>
        void Polygonal_surface_reconstruction<GeomTraits>::output_candidate_faces(PolygonMesh& candidate_faces) const {
                candidate_faces.clear();        // make sure it is empty.
                CGAL::copy_face_graph(candidate_faces_, candidate_faces);
        }


        template <class GeomTraits>
        template <typename MixedIntegerProgramTraits, typename PolygonMesh>
        bool Polygonal_surface_reconstruction<GeomTraits>::reconstruct(
                PolygonMesh& output_mesh,
                double wt_fitting /* = 0.43 */,
                double wt_coverage /* = 0.27 */,
                double wt_complexity /* = 0.30 */)
        {
                if (!error_message_.empty()) { // an error has occurred in the constructor
                        return false;
                }

                if (candidate_faces_.num_faces() < 4) {
                        error_message_ = "at least 4 candidate faces required to reconstruct a closed surface mesh (only "
                                + std::to_string(candidate_faces_.num_faces()) + " computed)";
                        return false;
                }

                typedef typename internal::Hypothesis<GeomTraits>::Adjacency Adjacency;
                const Adjacency& adjacency = hypothesis_.extract_adjacency(candidate_faces_);

                // Internal data structure
                Polygon_mesh target_mesh = candidate_faces_;

                // The number of supporting points of each face
                typename Polygon_mesh::template Property_map<Face_descriptor, std::size_t> face_num_supporting_points =
                        target_mesh.template add_property_map<Face_descriptor, std::size_t>("f:num_supporting_points").first;

                // The area of each face
                typename Polygon_mesh::template Property_map<Face_descriptor, FT> face_areas =
                        target_mesh.template add_property_map<Face_descriptor, FT>("f:face_area").first;

                // The point covered area of each face
                typename Polygon_mesh::template Property_map<Face_descriptor, FT> face_covered_areas =
                        target_mesh.template add_property_map<Face_descriptor, FT>("f:covered_area").first;

                // The supporting plane of each face
                typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
                        target_mesh.template add_property_map<Face_descriptor, const Plane*>("f:supp_plane").first;

                // Gives each face an index
                typename Polygon_mesh::template Property_map<Face_descriptor, std::size_t> face_indices =
                        target_mesh.template add_property_map<Face_descriptor, std::size_t>("f:index").first;

                double total_points = 0.0;
                std::size_t idx = 0;
                for(auto f : target_mesh.faces()) {
                        total_points += face_num_supporting_points[f];
                        face_indices[f] = idx;
                        ++idx;
                }


                typedef MixedIntegerProgramTraits                                                                MIP_Solver;
                typedef typename MixedIntegerProgramTraits::Variable                        Variable;
                typedef typename MixedIntegerProgramTraits::Linear_objective        Linear_objective;
                typedef typename MixedIntegerProgramTraits::Linear_constraint        Linear_constraint;

                MIP_Solver solver;

                // Adds variables

                // Binary variables:
                // x[0] ... x[num_faces - 1] : binary labels of all the input faces
                // x[num_faces] ... x[num_faces + num_edges - 1] : binary labels of all the intersecting edges (remain or not)
                // x[num_faces + num_edges] ... x[num_faces + num_edges + num_edges] : binary labels of corner edges (sharp edge of not)

                std::size_t num_faces = target_mesh.number_of_faces();
                std::size_t num_edges(0);

                typedef typename internal::Hypothesis<GeomTraits>::Intersection        Intersection;

                std::unordered_map<const Intersection*, std::size_t> edge_usage_status;        // keep or remove an intersecting edges
                for (std::size_t i = 0; i < adjacency.size(); ++i) {
                        const Intersection& fan = adjacency[i];
                        if (fan.size() == 4) {
                                std::size_t var_idx = num_faces + num_edges;
                                edge_usage_status[&fan] = var_idx;
                                ++num_edges;
                        }
                }

                std::size_t total_variables = num_faces + num_edges + num_edges;

                const std::vector<Variable*>& variables = solver.create_variables(total_variables);
                for (std::size_t i = 0; i < total_variables; ++i) {
                        Variable* v = variables[i];
                        v->set_variable_type(Variable::BINARY);
                }

                // Adds objective

                const typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = target_mesh.points();
                std::vector<Point> vertices(target_mesh.number_of_vertices());
                idx = 0;
                for(auto v : target_mesh.vertices()) {
                        vertices[idx] = coords[v];
                        ++idx;
                }

                typedef typename GeomTraits::Iso_cuboid_3                Box;

                const Box& box = CGAL::bounding_box(vertices.begin(), vertices.end());
                FT dx = box.xmax() - box.xmin();
                FT dy = box.ymax() - box.ymin();
                FT dz = box.zmax() - box.zmin();
                FT box_area = FT(2.0) * (dx * dy + dy * dz + dz * dx);

                // Chooses a better scale: all actual values multiplied by total number of points
                double coeff_data_fitting = wt_fitting;
                double coeff_coverage = total_points * wt_coverage / box_area;
                double coeff_complexity = total_points * wt_complexity / double(adjacency.size());

                Linear_objective * objective = solver.create_objective(Linear_objective::MINIMIZE);

                std::unordered_map<const Intersection*, std::size_t> edge_sharp_status;        // the edge is sharp or not
                std::size_t num_sharp_edges = 0;
                for (std::size_t i = 0; i < adjacency.size(); ++i) {
                        const Intersection& fan = adjacency[i];
                        if (fan.size() == 4) {
                                std::size_t var_idx = num_faces + num_edges + num_sharp_edges;
                                edge_sharp_status[&fan] = var_idx;

                                // Accumulates model complexity term
                                objective->add_coefficient(variables[var_idx], coeff_complexity);
                                ++num_sharp_edges;
                        }
                }
                CGAL_assertion(num_edges == num_sharp_edges);

                for(auto f : target_mesh.faces()) {
                        std::size_t var_idx = face_indices[f];

                        // Accumulates data fitting term
                        std::size_t num = face_num_supporting_points[f];
                        objective->add_coefficient(variables[var_idx], -coeff_data_fitting * num);

                        // Accumulates model coverage term
                        double uncovered_area = (face_areas[f] - face_covered_areas[f]);
                        objective->add_coefficient(variables[var_idx], coeff_coverage * uncovered_area);
                }

                // Adds constraints: the number of faces associated with an edge must be either 2 or 0
                std::size_t var_edge_used_idx = 0;
                for (std::size_t i = 0; i < adjacency.size(); ++i) {
                        Linear_constraint* c = solver.create_constraint(0.0, 0.0);
                        const Intersection& fan = adjacency[i];
                        for (std::size_t j = 0; j < fan.size(); ++j) {
                                Face_descriptor f = target_mesh.face(fan[j]);
                                std::size_t var_idx = face_indices[f];
                                c->add_coefficient(variables[var_idx], 1.0);
                        }

                        if (fan.size() == 4) {
                                std::size_t var_idx = num_faces + var_edge_used_idx;
                                c->add_coefficient(variables[var_idx], -2.0);  //
                                ++var_edge_used_idx;
                        }
                        else { // boundary edge
                                   // will be set to 0 (i.e., we don't allow open surface)
                        }
                }

                // Adds constraints: for the sharp edges. The explanation of posing this constraint can be found here:
                // https://user-images.githubusercontent.com/15526536/30185644-12085a9c-942b-11e7-831d-290dd2a4d50c.png
                double M = 1.0;
                for (std::size_t i = 0; i < adjacency.size(); ++i) {
                        const Intersection& fan = adjacency[i];
                        if (fan.size() != 4)
                                continue;

                        // If an edge is sharp, the edge must be selected first:
                        // X[var_edge_usage_idx] >= X[var_edge_sharp_idx]
                        Linear_constraint* c = solver.create_constraint(0.0);
                        std::size_t var_edge_usage_idx = edge_usage_status[&fan];
                        c->add_coefficient(variables[var_edge_usage_idx], 1.0);
                        std::size_t var_edge_sharp_idx = edge_sharp_status[&fan];
                        c->add_coefficient(variables[var_edge_sharp_idx], -1.0);

                        for (std::size_t j = 0; j < fan.size(); ++j) {
                                Face_descriptor f1 = target_mesh.face(fan[j]);
                                const Plane* plane1 = face_supporting_planes[f1];
                                std::size_t fid1 = face_indices[f1];
                                for (std::size_t k = j + 1; k < fan.size(); ++k) {
                                        Face_descriptor f2 = target_mesh.face(fan[k]);
                                        const Plane* plane2 = face_supporting_planes[f2];
                                        std::size_t fid2 = face_indices[f2];

                                        if (plane1 != plane2) {
                                                // The constraint is:
                                                //X[var_edge_sharp_idx] + M * (3 - (X[fid1] + X[fid2] + X[var_edge_usage_idx])) >= 1
                                                // which equals to
                                                //X[var_edge_sharp_idx] - M * X[fid1] - M * X[fid2] - M * X[var_edge_usage_idx] >= 1 - 3M
                                                c = solver.create_constraint(1.0 - 3.0 * M);
                                                c->add_coefficient(variables[var_edge_sharp_idx], 1.0);
                                                c->add_coefficient(variables[fid1], -M);
                                                c->add_coefficient(variables[fid2], -M);
                                                c->add_coefficient(variables[var_edge_usage_idx], -M);
                                        }
                                }
                        }
                }

                // Optimization

                if (solver.solve()) {

                        // Marks results
                        const std::vector<double>& X = solver.solution();

                        std::vector<Face_descriptor> to_delete;
                        std::size_t f_idx(0);
                        for(auto f : target_mesh.faces()) {
                                if (static_cast<int>(std::round(X[f_idx])) == 0)
                                        to_delete.push_back(f);
                                ++f_idx;
                        }

                        for (std::size_t i = 0; i < to_delete.size(); ++i) {
                                Face_descriptor f = to_delete[i];
                                Halfedge_descriptor h = target_mesh.halfedge(f);
                                Euler::remove_face(h, target_mesh);
                        }

                        // Marks the sharp edges
                        typename Polygon_mesh::template Property_map<Edge_descriptor, bool> edge_is_sharp =
                                target_mesh.template add_property_map<Edge_descriptor, bool>("e:sharp_edges").first;
                        for (auto e : target_mesh.edges())
                                edge_is_sharp[e] = false;

                        for (std::size_t i = 0; i < adjacency.size(); ++i) {
                                const Intersection& fan = adjacency[i];
                                if (fan.size() != 4)
                                        continue;

                                std::size_t idx_sharp_var = edge_sharp_status[&fan];
                                if (static_cast<int>(X[idx_sharp_var]) == 1) {
                                        for (std::size_t j = 0; j < fan.size(); ++j) {
                                                Halfedge_descriptor h = fan[j];
                                                Face_descriptor f = target_mesh.face(h);
                                                if (f != Polygon_mesh::null_face()) { // some faces may be deleted
                                                        std::size_t fid = face_indices[f];
                                                        if (static_cast<int>(std::round(X[fid])) == 1) {
                                                                Edge_descriptor e = target_mesh.edge(h);
                                                                edge_is_sharp[e] = true;
                                                                break;
                                                        }
                                                }
                                        }
                                }
                        }

                        // Converts from internal data structure to the required `PolygonMesh`.
                        output_mesh.clear();        // make sure it is empty.
                        CGAL::copy_face_graph(target_mesh, output_mesh);
                }
                else {
                        error_message_ = "solving the binary program failed";
                        return false;
                }

                return true;
        }

} //namespace CGAL

#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H
