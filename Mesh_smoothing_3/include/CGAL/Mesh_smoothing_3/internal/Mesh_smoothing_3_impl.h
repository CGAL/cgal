// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : François Protais

#ifndef CGAL_MESH_SMOOTHING_3_INTERNAL_MESH_SMOOTHING_3_IMPL_H
#define CGAL_MESH_SMOOTHING_3_INTERNAL_MESH_SMOOTHING_3_IMPL_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <CGAL/Mesh_smoothing_3/internal/utils/math_functions.h>

namespace CGAL {

namespace Mesh_smoothing_3 {

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::Mesh_smoother(TetrahedralMesh &mesh, BoundaryMesh const &boundary, EdgeNetwork const &edge_network)
: _mesh(mesh)
, _boundary(boundary)
, _edge_network(edge_network)
{}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_locked_boundary(bool locked) {
    _lock_boundary = locked;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_verbose(bool verbose) {
    _verbose = verbose;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_max_number_of_iteration(unsigned number_of_iterations) {
    _max_number_of_iteration = number_of_iterations;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag> // necessay
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_optimization_mode(Parameters::OPTIMIZATION_MODE mode) {
    _optimization_mode = mode;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_predicates_mode(Parameters::PREDICATES_MODE mode) {
    _predicates_mode = mode;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_maximum_running_time(double time_limit) {
    _time_limit = time_limit;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_maximum_number_of_metric_evaluations(unsigned max_nb_metric_evaluations) {
    _max_nb_metric_evaluations = max_nb_metric_evaluations;
}


template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_minimum_valid_edge_size(double val) {
    _min_valid_edge_size = val;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_boundary_weight(double val) {
    _boundary_weight = val;
}


template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_boundary_weight(Parameters::BOUNDARY_WEIGHTING_MODE mode) {
    switch (mode) {
        case Parameters::DEFAULT:
            _boundary_weight = 1;
            break;
        case Parameters::STRONG:
            _boundary_weight = 10;
            break;
        case Parameters::SOFT:
            _boundary_weight = 1e-3;
            break;
    }
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
unsigned Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::get_total_number_of_lbfgs_iterations() const {
    return _nb_lbfgs_iterations;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
unsigned Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::get_number_of_invalid_steps_measured_by_predicates() const {
    return _nb_predicates_invalid_steps;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_vertex_lock(Vertex_descriptor vertex_descriptor, bool locked) {
    auto res = _user_locks.emplace(vertex_descriptor, std::array<bool, 3>{false, false, false});
    (*res.first).second = {locked, locked, locked};
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_vertex_dim_lock(Vertex_descriptor vertex_descriptor, unsigned dimension, bool locked) {
    assert(dimension < 3);
    auto res = _user_locks.emplace(vertex_descriptor, std::array<bool, 3>{false, false, false});
    (*res.first).second[dimension] = locked;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_vertices_dim_locks(Vertex_descriptor_map<std::array<bool, 3>> const &vertex_dimension_locks)  {
    _user_locks = vertex_dimension_locks;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::clear_locks()  {
    _user_locks.clear();
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_boundary_query(Boundary_point_query boundary_query) {
    _lock_boundary = false;
    _boundary_query_type = POINT_QUERY;
    _boundary_point_query = boundary_query;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_boundary_query(Boundary_polygon_query boundary_query) {
    _lock_boundary = false;
    _boundary_query_type = ELEMENT_QUERY;
    _boundary_polygon_query = boundary_query;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_boundary_query(Boundary_point_batch_query boundary_query) {
    _lock_boundary = false;
    _boundary_query_type = BATCH_POINT_QUERY;
    _boundary_point_batch_query = boundary_query;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_boundary_query(Boundary_polygon_batch_query boundary_query) {
    _lock_boundary = false;
    _boundary_query_type = BATCH_ELEMENT_QUERY;
    _boundary_polygon_batch_query = boundary_query;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_curves_query(Curve_point_query curve_query) {
    _lock_boundary = false;
    _curve_query_type = POINT_QUERY;
    _curve_point_query = curve_query;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_curves_query(Curve_segment_query curve_query) {
    _lock_boundary = false;
    _curve_query_type = ELEMENT_QUERY;
    _curve_segment_query = curve_query;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_curves_query(Curve_point_batch_query curve_query) {
    _lock_boundary = false;
    _curve_query_type = BATCH_POINT_QUERY;
    _curve_point_batch_query = curve_query;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_curves_query(Curve_segment_batch_query curve_query) {
    _lock_boundary = false;
    _curve_query_type = BATCH_ELEMENT_QUERY;
    _curve_segment_batch_query = curve_query;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_vertex_target_position(Vertex_descriptor v, Point_3 const &target_position, double weight) {
    _vertex_target_positions.push_back({v, target_position, weight});
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_vertex_target_positions(std::vector<std::tuple<Vertex_descriptor, Point_3, double>> const &target_positions) {
    _vertex_target_positions = target_positions;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::clear_vertex_target_positions() {
    _vertex_target_positions.clear();
}


template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
bool Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::is_vert_locked(Vertex_descriptor Vertex_descriptor) const {
    auto iter = _user_locks.find(Vertex_descriptor);
    if (iter == _user_locks.end()) return false;
    return (*iter).second[0] && (*iter).second[1] && (*iter).second[2];
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_update_validator(Update_validator func) {
    _update_validator = func;
}


template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::check_refs() {
    // when is it that the input tetrahedra is wrong, or is it just tangled?
    // I didn't find a good way to sanitize the inputs because of the variety of input garbage :(

    // bool input_contains_invalid_references = false;
    // for (auto cell : _mesh.cell_range()) {
    //     auto ref = _mesh.cell_reference_shape(cell);
    //     std::array<Eigen::Vector3d, 4> tet;
    //     for (unsigned i = 0; i < 4; ++i) {
    //         tet[i] = convert_to_eigen(ref[i]);
    //     }
    //     if ((tet[1]-tet[0]).cross(tet[2]-tet[0]).dot(tet[3]-tet[0]) <= 0) {
    //         input_contains_invalid_references = true;
    //     }
    // }
    // assert(!input_contains_invalid_references);
    // if (input_contains_invalid_references && _verbose) {
    //     CGAL::IO::Color_stream_guard bright_red(std::cout, CGAL::IO::Ansi_color::BrightBlue);
    //     std::cout << "Mesh_smoother: some tetrahedron reference shapes are invalid (negative volume). Results may be incorrect." << std::endl;
    // }
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::clear_internal_data() {
    _compressed_coords.resize(0);
    _compressed_locks.clear();
    _vertex_original_to_compressed.clear();
    _cell_original_to_compressed.clear();
    _tetrahedra.clear();
    _tetrahedron_refs.clear();
    _vert2tet_corner.clear();

    _bnd_faces.clear();
    _vert_and_face_corners.clear();
    _face_surface_id.clear();

    _curve_edges.clear();
    _curve_ids.clear();
    _point_targets.clear();

    _callback_initialized = false;
    _callback_vertex_map_data.clear();
    _callback_cell_map_data.clear();

    _current_coords_to_check.clear();
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::create_compress_sorted_data() {
    clear_internal_data();
    if (_verbose) std::cout << "Mesh_smoother: copying and compressing data " << std::endl;
    std::size_t nb_points = 0;
    _vertex_original_to_compressed.reserve(_mesh.nb_vertices());
    std::vector<unsigned> nb_tet_on_verts;
    std::vector<double> temp_coordinates_storage;
    nb_tet_on_verts.reserve(_mesh.nb_vertices());
    _compressed_locks.reserve(3*_mesh.nb_vertices());
    temp_coordinates_storage.reserve(3*_mesh.nb_vertices());

    auto get_compressed_point_id = [&](Vertex_descriptor vertex_descriptor) {
        auto res = _vertex_original_to_compressed.emplace(vertex_descriptor, nb_points);
        if (res.second) {
            ++nb_points;
            nb_tet_on_verts.push_back(0);
            Point_3 coordinates = _mesh.vertex_coordinates(vertex_descriptor);
            auto iterator = _user_locks.find(vertex_descriptor);
            for (unsigned i = 0; i < 3; ++i) {
                temp_coordinates_storage.push_back(coordinates[i]);
                _compressed_locks.push_back(iterator == _user_locks.end() ? false : (*iterator).second[i]);
            }
        }
        return (*res.first).second;
    };

    // initialize volume

    _tetrahedra.reserve(_mesh.nb_cells());
    _tetrahedron_refs.reserve(_mesh.nb_cells());
    for (auto cell : _mesh.cell_range()) {
        bool hasUnlocked = false;
        auto vertex_indices = _mesh.cell_vertices(cell);
        for (unsigned i = 0; i < 4; ++i) {
            if (!is_vert_locked(vertex_indices[i])) {
                hasUnlocked = true;
                break;
            }
        }
        if (!hasUnlocked) continue;

        auto cell_ref = _mesh.cell_reference_shape(cell);
        std::array<Eigen::Vector3d, 4> cell_ref_eigen;
        std::array<unsigned, 4> compressed_index;
        for (unsigned i = 0; i < 4; ++i) {
            cell_ref_eigen[i] = convert_to_eigen(cell_ref[i]);
            compressed_index[i] = get_compressed_point_id(vertex_indices[i]);
            ++nb_tet_on_verts[compressed_index[i]];
        }
        _cell_original_to_compressed.emplace(cell, _tetrahedra.size());
        _tetrahedra.push_back(compressed_index);
        _tetrahedron_refs.push_back(Mesh_smoothing_3_internal::Math_functions::transform_coordinates_to_gradient_base(cell_ref_eigen));
    }

    _tetrahedra.shrink_to_fit();
    _tetrahedron_refs.shrink_to_fit();

    // computing vertex -> tet data structure
    std::vector<unsigned> curr_id(nb_points,0);
    _vert2tet_corner.reserve(nb_points);
    for (unsigned v = 0; v < nb_points; ++v) {
        _vert2tet_corner.push_back(std::vector<unsigned>(nb_tet_on_verts[v]));
    }
    for (unsigned t = 0; t < _tetrahedra.size(); ++t) {
        for (unsigned i = 0; i < 4; ++i) {
            unsigned v = _tetrahedra[t][i];
            _vert2tet_corner[v][curr_id[v]++] = 4*t+i;
        }
    }

    // initialize boundary

    std::vector<unsigned> currFace(10);
    _bnd_faces.reserve(_boundary.nb_faces());
    _face_surface_id.reserve(_boundary.nb_faces());
    std::unordered_map<unsigned, std::vector<std::array<unsigned, 2>>> vert2faces; // todo: maybe use a map from boost
    for (auto face : _boundary.face_range()) {
        bool hasUnlocked = false;
        for (auto vertex_descriptor : _boundary.face_vertices(face)) {
            if (!is_vert_locked(vertex_descriptor)) {
                hasUnlocked = true;
                break;
            }
        }
        if (!hasUnlocked) continue;
        currFace.clear();
        for (auto vertex_descriptor : _boundary.face_vertices(face)) {
            currFace.push_back(get_compressed_point_id(vertex_descriptor));
        }
        for (unsigned i = 0; i < currFace.size(); ++i) {
            if (_lock_boundary) {
                for (unsigned d = 0; d < 3; ++d) {
                    _compressed_locks[3*currFace[i]+d] = true;
                }
            }
            std::array<unsigned, 2> location_pair = {static_cast<unsigned>(_bnd_faces.size()), i};
            auto res = vert2faces.emplace(currFace[i], std::vector<std::array<unsigned, 2>>{location_pair});
            if (!res.second) (*res.first).second.push_back(location_pair);
        }
        _bnd_faces.push_back(currFace);
        _face_surface_id.push_back(_boundary.patch_id(face));
    }
    // not ordered, do a sort? then I lose memory alignment? check if becomes bottleneck
    _vert_and_face_corners.assign(vert2faces.begin(), vert2faces.end());


    // initialize curve network

    _curve_edges.reserve(_edge_network.nb_edges());
    _curve_ids.reserve(_edge_network.nb_edges());
    for (auto edge : _edge_network.edge_range()) {
        Vertex_descriptor v0 = _edge_network.edge_vertex(edge, 0);
        Vertex_descriptor v1 = _edge_network.edge_vertex(edge, 1);
        if (is_vert_locked(v0) && is_vert_locked(v1)) continue;
        _curve_edges.push_back({get_compressed_point_id(v0), get_compressed_point_id(v1)});
        _curve_ids.push_back(_edge_network.curve_id(edge));
    }



    // initialize rest of data
    _point_targets.reserve(_vertex_target_positions.size());

    _compressed_coords = Eigen::Map<Eigen::VectorXd>(temp_coordinates_storage.data(), static_cast<Eigen::Index>(temp_coordinates_storage.size()));

    if (_update_validator != nullptr) {
        _current_coords_to_check.reserve(_mesh.nb_vertices());
        for (auto [v, id] : _vertex_original_to_compressed) {
            _current_coords_to_check.emplace(v, _mesh.vertex_coordinates(v));
        }
    }
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::initialize_boundary() {
    std::vector<unsigned> currFace;
    _bnd_faces.reserve(_boundary.nb_faces());
    _face_surface_id.reserve(_boundary.nb_faces());
    std::unordered_map<unsigned, std::vector<std::array<unsigned, 2>>> vert2faces; // todo: maybe use a map from boost
    for (auto face : _boundary.face_range()) {
        currFace.clear();
        for (auto vertex_descriptor : _boundary.face_vertices(face)) {
            auto iterator = _vertex_original_to_compressed.find(static_cast<Vertex_descriptor>(vertex_descriptor));
            // it is not dramatic for Tetrahedral_conformal_smoother to not have the perfect face
            if (iterator == _vertex_original_to_compressed.end()) continue;
            currFace.push_back(iterator->second);
        }
        if (currFace.empty()) continue;
        for (unsigned i = 0; i < currFace.size(); ++i) {
            if (_lock_boundary) {
                for (unsigned d = 0; d < 3; ++d) {
                    _compressed_locks[3*currFace[i]+d] = true;
                }
            }
            std::array<unsigned, 2> location_pair = {static_cast<unsigned>(_bnd_faces.size()), i};
            auto res = vert2faces.emplace(currFace[i], std::vector<std::array<unsigned, 2>>{location_pair});
            if (!res.second) (*res.first).second.push_back(location_pair);
        }
        _bnd_faces.push_back(currFace);
        _face_surface_id.push_back(_boundary.patch_id(face));
    }
    // not ordered, do a sort? then I lose memory alignment? check if becomes bottleneck
    _vert_and_face_corners.assign(vert2faces.begin(), vert2faces.end());
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::initialize_curve_network() {
    std::vector<unsigned> currFace;
    _curve_edges.reserve(_edge_network.nb_edges());
    _curve_ids.reserve(_edge_network.nb_edges());
    for (auto edge : _edge_network.edge_range()) {
        std::array<unsigned, 2> edge_vertices;
        bool not_in = false;
        for (unsigned i = 0; i < 2; ++i) {
            Vertex_descriptor v = _edge_network.edge_vertex(edge, i);
            auto iterator = _vertex_original_to_compressed.find(static_cast<Vertex_descriptor>(v));
            // if we do not find the vertex, we just ignore the edge
            if (iterator == _vertex_original_to_compressed.end()) {
                not_in = true;
                break;
            }
            edge_vertices[i] = iterator->second;
        }
        if (not_in) continue;
        _curve_edges.push_back(edge_vertices);
        _curve_ids.push_back(_edge_network.curve_id(edge));
    }

    _point_targets.reserve(_vertex_target_positions.size());
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::initialise_point_targets() {
    for (auto [v, target, weight] : _vertex_target_positions) {
        auto iterator = _vertex_original_to_compressed.find(static_cast<Vertex_descriptor>(v));
        if (iterator == _vertex_original_to_compressed.end()) continue;
        _point_targets.push_back({iterator->second, convert_to_inner(target), weight});
    }
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::rescale_geometry() {
    if (_compressed_coords.size() == 0) return;
    std::size_t nb_vertices = static_cast<std::size_t>(_compressed_coords.size()) / 3;
    Eigen::Vector3d mini(_compressed_coords[0], _compressed_coords[1], _compressed_coords[2]);
    Eigen::Vector3d maxi(_compressed_coords[0], _compressed_coords[1], _compressed_coords[2]);
    Eigen::Vector3d center(0.,0.,0.);
    for (unsigned i = 0; i < nb_vertices; ++i) {
        for (unsigned d = 0; d < 3;++d) {
            mini[d] = (std::min)(mini[d], _compressed_coords[3*i+d]);
            maxi[d] = (std::max)(maxi[d], _compressed_coords[3*i+d]);
            center[d] += _compressed_coords[3*i+d];
        }
    }

    center /= static_cast<double>(nb_vertices);
    _scale = (std::max)((std::max)(maxi[0] - mini[0], maxi[1] - mini[1]), maxi[2]-mini[2]) / _rescaled_bbox_scale;
    if (_scale < 1e-14) { // todo figure out a better thresh old here
        _scale = 1.;
        return;
    }
    _shift = center;

    for (int i = 0; i < static_cast<int>(nb_vertices); ++i) {
        for (int d = 0; d < 3;++d) {
            _compressed_coords[3*i+d] = (_compressed_coords[3*i+d] - _shift[d]) / _scale;
        }
    }
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::update_mesh_coordinates() {
    for (auto [vertex_descriptor, compressed_id] : _vertex_original_to_compressed) {
        if (_compressed_locks[3*compressed_id] && _compressed_locks[3*compressed_id+1] && _compressed_locks[3*compressed_id+2]) continue;
        Eigen::Vector3d pt = Mesh_smoothing_3_internal::Math_functions::sub_col_vector(_compressed_coords, compressed_id);
        _mesh. set_vertex_coordinates(vertex_descriptor, convert_to_user(pt));
    }
}


template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_callback_function(Callback_function callback_function, Callback_setting setting) {
    _callback_function = callback_function;
    _callback_setting = setting;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::initialize_callback() {
    if (_callback_initialized) return;
    _callback_vertex_map_data.reserve(_vertex_original_to_compressed.size());
    _callback_cell_map_data.reserve(_cell_original_to_compressed.size());
    Vertex_status default_vertex_status;
    for (auto [vertex_descriptor, compressed_id] : _vertex_original_to_compressed) {
        _callback_vertex_map_data.emplace(vertex_descriptor, default_vertex_status);
    }
    Cell_status default_cell_status;
    for (auto [cell_descriptor, compressed_id] : _cell_original_to_compressed) {
        _callback_cell_map_data.emplace(cell_descriptor, default_cell_status);
    }
    _callback_initialized = true;
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
bool Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::run_callback(Iteration_status const &status, std::vector<Vertex_status> const &vertex_status, std::vector<Cell_status> const &cell_status) {
    initialize_callback();
    update_mesh_coordinates();
    _callback_status = status;
    _callback_status.scaling_factor = _scale;
    for (auto [vertex_descriptor, compressed_id] : _vertex_original_to_compressed) {
        _callback_vertex_map_data.at(vertex_descriptor) = vertex_status[compressed_id];
    }
    for (auto [cell_descriptor, compressed_id] : _cell_original_to_compressed) {
        _callback_cell_map_data.at(cell_descriptor) = cell_status[compressed_id];
    }
    return _callback_function(_callback_status, _callback_vertex_map_data, _callback_cell_map_data);
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::initialise_boundary_query(Tetrahedral_mesh_smoother &smoother) {
    switch (_boundary_query_type) {
    case NONE:
        smoother.set_boundary_with_singular_query(
            _bnd_faces,
            &_vert_and_face_corners,
            nullptr,
            _face_surface_id
        );
        break;
    case POINT_QUERY:
        {
            typename Tetrahedral_mesh_smoother::Boundary_query boundary_query = [&](std::vector<Eigen::Vector3d> const &coords, Surface_patch_index surface_id) {
                Eigen::Vector3d center = Eigen::Vector3d::Zero();
                for (auto const &coord : coords) {
                    center += coord;
                }
                center /= static_cast<double>(coords.size());
                double radius = 0;
                for (auto const &coord : coords) {
                    radius += (coord-center).norm();
                }
                radius /= static_cast<double>(coords.size());
                radius /= _scale;
                auto [user_point, user_normal, weight] =  _boundary_point_query(convert_to_user(center), surface_id, radius);
                Eigen::Vector3d proj = convert_to_inner(user_point);
                Eigen::Vector3d normal = convert_to_eigen(user_normal);
                normal.normalize();
                return typename Tetrahedral_mesh_smoother::Plane{ proj, normal, weight };
            };
            smoother.set_boundary_with_singular_query(
                _bnd_faces,
                &_vert_and_face_corners,
                boundary_query,
                _face_surface_id
            );
        }
        break;
    case ELEMENT_QUERY:
        {
            typename Tetrahedral_mesh_smoother::Boundary_query boundary_query = [&](std::vector<Eigen::Vector3d> const &coords, Surface_patch_index surface_id) {
                std::vector<Point_3> polygon(coords.size());
                for (unsigned i = 0; i < coords.size(); ++i) {
                    polygon[i] = convert_to_user(coords[i]);
                }
                auto [user_point, user_normal, weight] =  _boundary_polygon_query(polygon, surface_id);
                Eigen::Vector3d proj = convert_to_inner(user_point);
                Eigen::Vector3d normal = convert_to_eigen(user_normal);
                normal.normalize();
                return typename Tetrahedral_mesh_smoother::Plane{ proj, normal, weight };
            };
            smoother.set_boundary_with_singular_query(
                _bnd_faces,
                &_vert_and_face_corners,
                boundary_query,
                _face_surface_id
            );
        }
        break;
    case BATCH_POINT_QUERY:
        {
            typename Tetrahedral_mesh_smoother::Boundary_batch_query boundary_query = [&](std::vector<std::vector<Eigen::Vector3d>> const &coords, std::vector<Surface_patch_index> surface_ids, std::vector<typename Tetrahedral_mesh_smoother::Plane> &results) {
                if (coords.size() != _boundary_batch_info_points.size()) {
                    _boundary_batch_info_points.resize(coords.size());
                    _boundary_batch_info_radii.resize(coords.size());
                    _boundary_batch_planes.resize(coords.size());
                }
                Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, coords.size(), [&](std::size_t uiter) {
                    Eigen::Vector3d center = Eigen::Vector3d::Zero();
                    for (auto const &coord : coords[uiter]) {
                        center += coord;
                    }
                    center /= static_cast<double>(coords[uiter].size());
                    double radius = 0;
                    for (auto const &coord : coords[uiter]) {
                        radius += (coord-center).norm();
                    }
                    radius /= static_cast<double>(coords[uiter].size());
                    radius /= _scale;
                    _boundary_batch_info_points[uiter] = convert_to_user(center);
                    _boundary_batch_info_radii[uiter] = radius;
                });

                _boundary_point_batch_query(_boundary_batch_info_points, surface_ids, _boundary_batch_info_radii, _boundary_batch_planes);

                Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, coords.size(), [&](std::size_t t) {
                    auto [user_point, user_normal, weight] =  _boundary_batch_planes[t];
                    Eigen::Vector3d proj = convert_to_inner(user_point);
                    Eigen::Vector3d normal = convert_to_eigen(user_normal);
                    normal.normalize();
                    results[t] = { proj, normal, weight };
                });
            };
            smoother.set_boundary_with_batch_query(
                _bnd_faces,
                &_vert_and_face_corners,
                boundary_query,
                _face_surface_id
            );
        }
        break;
    case BATCH_ELEMENT_QUERY:
        {
            typename Tetrahedral_mesh_smoother::Boundary_batch_query boundary_query = [&](std::vector<std::vector<Eigen::Vector3d>> const &coords, std::vector<Surface_patch_index> surface_ids, std::vector<typename Tetrahedral_mesh_smoother::Plane> &results) {
                if (coords.size() != _boundary_batch_info_polygons.size()) {
                    _boundary_batch_info_polygons.reserve(coords.size());
                    for (unsigned i = 0; i < coords.size(); ++i) {
                        _boundary_batch_info_polygons.push_back(std::vector<Point_3>(coords[i].size()));
                    }
                    _boundary_batch_planes.resize(coords.size());
                }
                Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, coords.size(), [&](std::size_t t) {
                    for (unsigned i = 0; i < coords[t].size(); ++i) {
                        _boundary_batch_info_polygons[t][i] = convert_to_user(coords[t][i]);
                    }
                });

                _boundary_polygon_batch_query(_boundary_batch_info_polygons, surface_ids, _boundary_batch_planes);

                Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, coords.size(), [&](std::size_t t) {
                    auto [user_point, user_normal, weight] =  _boundary_batch_planes[t];
                    Eigen::Vector3d proj = convert_to_inner(user_point);
                    Eigen::Vector3d normal = convert_to_eigen(user_normal);
                    normal.normalize();
                    results[t] = { proj, normal, weight };
                });

            };
            smoother.set_boundary_with_batch_query(
                _bnd_faces,
                &_vert_and_face_corners,
                boundary_query,
                _face_surface_id
            );
        }
        break;
    }
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::initialise_curve_queries(Tetrahedral_mesh_smoother &smoother) {
    switch (_curve_query_type) {
    case NONE:
        break;
    case POINT_QUERY:
        {
            typename Tetrahedral_mesh_smoother::Curve_query curve_query = [&](std::array<Eigen::Vector3d, 2> const &coords, Curve_index curve_id) {
                Eigen::Vector3d center = 0.5*(coords[0] + coords[1]);
                double radius = (coords[1] - coords[0]).norm() / _scale;
                auto [user_point, user_tangent, weight] =  _curve_point_query(convert_to_user(center), curve_id, radius);
                Eigen::Vector3d proj = convert_to_inner(user_point);
                Eigen::Vector3d tangent = convert_to_eigen(user_tangent);
                tangent.normalize();
                return typename Tetrahedral_mesh_smoother::Curve_tangent{ proj, tangent, weight };
            };
            smoother.set_curve_network_with_singular_query(
                _curve_edges,
                _curve_ids,
                curve_query
            );
        }
        break;
    case ELEMENT_QUERY:
        {
            typename Tetrahedral_mesh_smoother::Curve_query curve_query = [&](std::array<Eigen::Vector3d, 2> const &coords, Curve_index curve_id) {
                std::array<Point_3, 2> edge = { convert_to_user(coords[0]),  convert_to_user(coords[1])};
                auto [user_point, user_tangent, weight] =  _curve_segment_query(edge, curve_id);
                Eigen::Vector3d proj = convert_to_inner(user_point);
                Eigen::Vector3d tangent = convert_to_eigen(user_tangent);
                tangent.normalize();
                return typename Tetrahedral_mesh_smoother::Curve_tangent{ proj, tangent, weight };
            };
            smoother.set_curve_network_with_singular_query(
                _curve_edges,
                _curve_ids,
                curve_query
            );
        }
        break;
    case BATCH_POINT_QUERY:
        {
            typename Tetrahedral_mesh_smoother::Curve_batch_query curve_query = [&](std::vector<std::array<Eigen::Vector3d, 2>> const &edges, std::vector<Curve_index> const &, std::vector<typename Tetrahedral_mesh_smoother::Curve_tangent> &results) {
                if (edges.size() != _curve_batch_info_points.size()) {
                    _curve_batch_info_points.resize(edges.size());
                    _curve_batch_info_radii.resize(edges.size());
                    _curve_batch_tangents.resize(edges.size());
                }
                for (unsigned e = 0; e < edges.size(); ++e) {
                    _curve_batch_info_points[e] = convert_to_user(0.5*(edges[e][0] + edges[e][1]));
                    _curve_batch_info_radii[e] = (edges[e][1] - edges[e][0]).norm() / _scale;
                }

                _curve_point_batch_query(_curve_batch_info_points, _curve_ids, _curve_batch_info_radii, _curve_batch_tangents);

                for (unsigned e = 0; e < edges.size(); ++e) {
                    auto [user_point, user_tangent, weight] =  _curve_batch_tangents[e];
                    Eigen::Vector3d proj = convert_to_inner(user_point);
                    Eigen::Vector3d tangent = convert_to_eigen(user_tangent);
                    tangent.normalize();
                    results[e] = typename Tetrahedral_mesh_smoother::Curve_tangent{ proj, tangent, weight };
                }
            };
            smoother.set_curve_network_with_batch_query(
                _curve_edges,
                _curve_ids,
                curve_query
            );
        }
        break;
    case BATCH_ELEMENT_QUERY:
        {
            typename Tetrahedral_mesh_smoother::Curve_batch_query curve_query = [&](std::vector<std::array<Eigen::Vector3d, 2>> const &edges, std::vector<Curve_index> const &, std::vector<typename Tetrahedral_mesh_smoother::Curve_tangent> &results) {
                if (edges.size() != _curve_batch_info_edges.size()) {
                    _curve_batch_info_edges.resize(edges.size());
                    _curve_batch_tangents.resize(edges.size());
                }
                for (unsigned e = 0; e < edges.size(); ++e) {
                    _curve_batch_info_edges[e] = {convert_to_user(edges[e][0]), convert_to_user(edges[e][1])};
                }

                _curve_segment_batch_query(_curve_batch_info_edges, _curve_ids, _curve_batch_tangents);

                for (unsigned e = 0; e < edges.size(); ++e) {
                    auto [user_point, user_tangent, weight] =  _curve_batch_tangents[e];
                    Eigen::Vector3d proj = convert_to_inner(user_point);
                    Eigen::Vector3d tangent = convert_to_eigen(user_tangent);
                    tangent.normalize();
                    results[e] = { proj, tangent, weight };
                }
            };
            smoother.set_curve_network_with_batch_query(
                _curve_edges,
                _curve_ids,
                curve_query
            );
        }
        break;
    }
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::initialise_smoother(Tetrahedral_mesh_smoother &smoother) {
    initialise_boundary_query(smoother);
    initialise_curve_queries(smoother);

    initialise_point_targets();
    smoother.set_quadratic_target_positions(_point_targets);

    if (_callback_function != nullptr) {
        smoother.callback_function = [&](Iteration_status const &status, std::vector<Vertex_status> const &vertex_data, std::vector<Cell_status> const &cell_data)
            {
                return run_callback(status, vertex_data, cell_data);
            };
    }

    if (_update_validator != nullptr) {
        smoother.set_validation_query(
            [&](Eigen::VectorXd const &coords) {
                for (auto &[v, pt]: _current_coords_to_check) {
                    Eigen::Vector3d live_pt = Mesh_smoothing_3_internal::Math_functions::sub_col_vector(coords, _vertex_original_to_compressed.at(v));
                    pt = convert_to_user(live_pt);
                }
                return _update_validator(_current_coords_to_check);
            }
        );
    }

    smoother.exact_predicate_status = _predicates_mode > Parameters::NO_CHECK;
    smoother.exact_predicate_optimization_check = _predicates_mode > Parameters::STATUS_UPDATE;
    smoother.exact_predicate_linesearch_enforcement = _predicates_mode > Parameters::CHECK_IN_LBFGS;

    smoother.callback_setting = _callback_setting;
    smoother.verbose = _verbose;
    smoother.min_valid_edge_size = _min_valid_edge_size / _scale;
    smoother.boundary_weight = _boundary_weight;
    smoother.fine_time_logging = _verbose && _tetrahedra.size() > 15'000'000; // todo: add a number of core specific threshold here
    smoother.time_limit = _time_limit;
    smoother.max_nb_metric_evaluations = _max_nb_metric_evaluations;

    smoother.smoother_status = &_smoothing_status;
}

namespace internal {
    inline std::pair<unsigned, unsigned> get_nb_inverted_cells(std::vector<double> const &determinants, std::vector<bool> const &exact_checks) {
        unsigned nb_det_inverted = 0;
        unsigned nb_exactly_inverted = 0;
        for (unsigned i = 0; i < determinants.size(); ++i) {
            if (determinants[i] <= 0) ++nb_det_inverted;
            if (!exact_checks[i]) ++nb_exactly_inverted;
        }

        return {nb_det_inverted, nb_exactly_inverted};
    }
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
std::pair<unsigned, unsigned> Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::get_nb_of_inverted_cells(Cell_descriptor_map<std::pair<double, bool>> *cell_determinants) {
    check_refs();
    create_compress_sorted_data();
    rescale_geometry();
    Tetrahedral_mesh_smoother smoother(_compressed_coords, _compressed_locks, _tetrahedra, _tetrahedron_refs, _vert2tet_corner);

    auto determinants = smoother.get_determinants();
    auto exact_checks = smoother.get_predicate_is_positive();

    if (cell_determinants != nullptr) {
        for (auto [cell_descriptor, compressed_id] : _cell_original_to_compressed) {
            (*cell_determinants)[cell_descriptor] = {determinants[compressed_id], exact_checks[compressed_id]};
        }
    }

    return get_nb_inverted_cells(determinants, exact_checks);
}

template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
void Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::set_input_smoothing_status(Smoothing_status const &status) {
    _smoothing_status = status;
}


template<typename TetrahedralMesh, typename BoundaryMesh, typename EdgeNetwork, typename ConcurrencyTag>
Smoothing_status Mesh_smoother<TetrahedralMesh, BoundaryMesh, EdgeNetwork, ConcurrencyTag>::run() {
    if (!_smoothing_status.in_progress()) _smoothing_status = Smoothing_status();

    check_refs();
    create_compress_sorted_data();
    rescale_geometry();
    Tetrahedral_mesh_smoother smoother(_compressed_coords, _compressed_locks, _tetrahedra, _tetrahedron_refs, _vert2tet_corner);

    auto [nb_det_inverted, nb_exactly_inverted] = internal::get_nb_inverted_cells(smoother.get_determinants(), smoother.get_predicate_is_positive());
    _smoothing_status.nb_initial_invalid_elements = _predicates_mode > Parameters::NO_CHECK ? nb_exactly_inverted : nb_det_inverted;

    if constexpr(!std::is_convertible_v<ConcurrencyTag, CGAL::Parallel_tag>) {
        Eigen::setNbThreads(1);
        assert(Eigen::nbThreads() == 1);
    }

    initialise_smoother(smoother);
    _smoothing_status.add_time(true);
    switch (_optimization_mode)
    {
    case Parameters::CONFORMAL:
        smoother.run_untangling(_max_number_of_iteration);
        break;
    case Parameters::MIXED_WITH_SIZE:
        assert(false /*Not implemented yet*/);
        break;
    case Parameters::LAPLACIAN:
        smoother.run_laplacian_gradient_descent(_max_number_of_iteration);
        break;
    case Parameters::MIN_MAXING_CONFORMAL:
        smoother.run_quality_maximization(_max_number_of_iteration);
        break;

    default:
        assert(false /*Unreachable*/);
        break;
    }

    update_mesh_coordinates();
    _nb_lbfgs_iterations = smoother.number_of_lbfgs_iter;
    _nb_predicates_invalid_steps = smoother.get_number_of_invalid_steps_with_predicates();

    std::tie(nb_det_inverted, nb_exactly_inverted) = internal::get_nb_inverted_cells(smoother.get_determinants(), smoother.get_predicate_is_positive());
    _smoothing_status.nb_invalid_elements = _predicates_mode > Parameters::NO_CHECK ? nb_exactly_inverted : nb_det_inverted;
    _smoothing_status.add_time(true);
    _smoothing_status.conclude();
    if constexpr(!std::is_convertible_v<ConcurrencyTag, CGAL::Parallel_tag>) {
        Eigen::setNbThreads(0);
    }
    return _smoothing_status;
}

} } // end of CGAL::Mesh_smoothing_3 namespace

#endif //CGAL_MESH_SMOOTHING_3_INTERNAL_MESH_SMOOTHING_3_IMPL_H
