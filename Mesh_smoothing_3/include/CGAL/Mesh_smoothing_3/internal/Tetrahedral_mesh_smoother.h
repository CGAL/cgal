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

#ifndef CGAL_MESH_SMOOTHING_3_INTERNAL_TETRAHEDRAL_MESH_SMOOTHER_H
#define CGAL_MESH_SMOOTHING_3_INTERNAL_TETRAHEDRAL_MESH_SMOOTHER_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <CGAL/Mesh_smoothing_3/internal/utils/math_functions.h>
#include <CGAL/Mesh_smoothing_3/internal/utils/Function_minimizer.h>
#include <CGAL/Mesh_smoothing_3/internal/utils/log_time.h>
#include <CGAL/Mesh_smoothing_3/internal/utils/loops.h>


#include <CGAL/IO/Color_ostream.h>

#include <Eigen/Eigen>

#include <vector>
#include <array>
#include <set>
#include <functional>
#include <random>
#include <fstream>


#include <CGAL/Mesh_smoothing_3/Smoothing_status.h>


namespace CGAL {

namespace Mesh_smoothing_3_internal {

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
class Tetrahedral_mesh_smoother {
public:
    Tetrahedral_mesh_smoother(
        Eigen::VectorXd &coords,
        std::vector<bool> const &locks,
        std::vector<std::array<unsigned, 4>> const &tetrahedra,
        std::vector<std::array<Eigen::Vector3d, 4>> const &tet_inv_grad,
        std::vector<std::vector<unsigned>> const &vert2tet_corner
    );

    using Plane = std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>; // point, normal, weight
    using Boundary_query = std::function<Plane (std::vector<Eigen::Vector3d> const &poly, Surface_patch_index surface_id)>;
    using Boundary_batch_query = std::function<void (std::vector<std::vector<Eigen::Vector3d>> const &polys, std::vector<Surface_patch_index> const &surface_id, std::vector<Plane> &results)>;


    void set_boundary_without_query(
        std::vector<std::vector<unsigned>> const &bnd_faces,
        std::vector<std::pair<unsigned, std::vector<std::array<unsigned, 2>>>> const *vert_and_face_corners, // pointer to highlight that the code keeps but not copy
        std::vector<Surface_patch_index> const &face_ids
    );


    void set_boundary_with_singular_query(
        std::vector<std::vector<unsigned>> const &bnd_faces,
        std::vector<std::pair<unsigned, std::vector<std::array<unsigned, 2>>>> const *vert_and_face_corners, // pointer to highlight that the code keeps but not copy
        Boundary_query boundary_query,
        std::vector<Surface_patch_index> const &face_ids
    );

    void set_boundary_with_batch_query(
        std::vector<std::vector<unsigned>> const &bnd_faces,
        std::vector<std::pair<unsigned, std::vector<std::array<unsigned, 2>>>> const *vert_and_face_corners, // pointer to highlight that the code keeps but not copy
        Boundary_batch_query boundary_query,
        std::vector<Surface_patch_index> const &face_ids
    );


    using Curve_tangent = std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>; // point, normal, weight

    using Curve_query = std::function<Curve_tangent (std::array<Eigen::Vector3d, 2> const &edge, Curve_index curve_id)>;
    using Curve_batch_query = std::function<void (std::vector<std::array<Eigen::Vector3d, 2>> const &edges, std::vector<Curve_index> const &curve_ids, std::vector<Curve_tangent> &results)>;
    void set_curve_network_with_singular_query(
        std::vector<std::array<unsigned, 2>> const &edges,
        std::vector<Curve_index> const &edge_ids,
        Curve_query query
    );

    void set_curve_network_with_batch_query(
        std::vector<std::array<unsigned, 2>> const &edges,
        std::vector<Curve_index> const &edge_ids,
        Curve_batch_query query
    );

    void set_quadratic_target_positions(std::vector<std::tuple<unsigned, Eigen::Vector3d, double>> const &targets);


    using Size_query = std::function<double (unsigned vertex_id, Eigen::Vector3d coord)>;
    using Size_batch_query = std::function<void (std::vector<unsigned> const &vertices, std::vector<Eigen::Vector3d> const &coords, std::vector<double> &sizes)>;

    void set_target_sizing(Size_query query); // todo: implement and use.
    void set_target_sizing(Size_batch_query query);

    using Validation_query = std::function<double (Eigen::VectorXd const &coords)>;

    void set_validation_query(Validation_query query);


    unsigned max_lbfgs_iter = 500;

    void set_starting_untangling_epsilon(double eps) { start_untangle_eps = eps; };




    double min_valid_edge_size = 1e-6;

    double boundary_weight = 1.;

    bool verbose = false;
    bool fine_time_logging = false;

    bool laplacian_precond = false;

    void run_laplacian_gradient_descent(unsigned max_number_iter = 100);

    bool run_untangling(unsigned max_number_iter = 1000);

    bool run_quality_maximization(unsigned max_number_iter = 100);

    std::vector<double> const &get_determinants() const {return _determinants; }
    double min_det() const {return _det_min; }
    double num_inverted() const {return _nb_inverted; }

    std::vector<double> const &get_conformal_energies() const {return _conformal_energies; }
    double get_max_conformal_energy() const {return _conformal_energy_max; }

    double start_untangle_eps = -1;

    unsigned number_of_outer_iter = 0;
    unsigned number_of_lbfgs_iter = 0;

    double time_limit = 0.;
    unsigned max_nb_metric_evaluations = 0;
    CGAL::Mesh_smoothing_3::Smoothing_status * smoother_status = nullptr;

public:
    bool exact_predicate_status = false;
    bool exact_predicate_optimization_check = false;
    bool exact_predicate_linesearch_enforcement = false;

    std::vector<bool> const &get_predicate_is_positive() {
        evaluate_exact_predicates();
        return _predicate_is_positive;
    }

    unsigned get_number_of_invalid_steps_with_predicates() const { return _predicates_nb_invalid_steps; }

public:

    enum DEBUG_CALLBACK_SETTING {NOTHING, OUTER_ITER, LBFGS_ITER} callback_setting = NOTHING;
    enum OPTIMIZATION_TYPE {UNDEFINED, UNTANGLING, STIFFENING, LAPLACIAN, INFLATION};

    struct LBFGS_status {
        unsigned iter = (std::numeric_limits<unsigned>::max)();
        double step = 0.;
        unsigned nbEval = 0;
        bool enabled() const { return iter != (std::numeric_limits<unsigned>::max)(); }
    };


    struct Iteration_status {
        OPTIMIZATION_TYPE opt = UNDEFINED;
        bool boundary_enabled = false;
        double min_edge_size = 0.;
        double boundary_weight = 1.;
        double scaling_factor = 1.;

        unsigned outer_iter_nb = 0;
        LBFGS_status lbfgs_status;
        bool is_in_lbfgs() const { return lbfgs_status.enabled(); }

        double smoothing_energy = 0.;
        double boundary_energy = 0.;
        double min_det = 0.;
        unsigned nb_inverted = 0;
        double opt_parameter = 0.;
    };

    struct Vertex_status {
        std::array<bool, 3> lock = {false, false, false};
        bool is_locked() const { return lock[0] && lock[1] && lock[2]; }
        double local_edge_size = 0.;
        Eigen::Vector3d smoothing_gradient = Eigen::Vector3d::Zero();
        Eigen::Vector3d boundary_gradient = Eigen::Vector3d::Zero();
        Eigen::Vector3d lbfgs_gradient = Eigen::Vector3d::Zero();
    };

    struct Tetrahedron_status {
        double energy_value = 0.;
        double weight = 0.;
        double det = 0.;
    };

    using Callback_function = std::function<bool (Iteration_status const &status, std::vector<Vertex_status> const&vertex_data, std::vector<Tetrahedron_status> const&tet_data)>;

    Callback_function callback_function = nullptr;

private:
    Iteration_status _callback_status;
    std::vector<Vertex_status> _callback_vert_storage;
    std::vector<Tetrahedron_status> _callback_tet_storage;
    Eigen::VectorXd _callback_smoothing_gradient;
    Eigen::VectorXd _callback_boundary_gradient;
    bool run_callback(OPTIMIZATION_TYPE opt_type, unsigned iter, LBFGS_status lbfgs_status = {(std::numeric_limits<unsigned>::max)(), 0., 0}, Eigen::VectorXd const *g = nullptr);


private:
    // mesh data
    struct Tet_storage {
        // input
        std::array<unsigned, 4> verts;
        std::array<Eigen::Vector3d, 4> ig;
        Eigen::Matrix3d compute_jacobian(Eigen::VectorXd const &coords) const;

        // internal usage
        bool skip = false;
        double local_edge_size = 1.;
        double det_estimation = 1.;
        double fval;
        std::array<Eigen::Vector3d, 4> vert_grad;
    };

    Eigen::VectorXd &_coords;
    std::vector<bool> _locks;
    std::vector<Tet_storage> _tet_storage;
    std::vector<std::vector<unsigned>> const &_vert2tet_corner;

    std::size_t nb_vertices() const { return _locks.size() / 3; }
    bool vertex_is_locked(unsigned v) const { return _locks[3*v+0]&&_locks[3*v+1]&&_locks[3*v+2]; }


    unsigned _nb_inverted;
    double _det_min;
    std::vector<double> _determinants;
    double _conformal_energy_max;
    std::vector<double> _conformal_energies;
    bool _collapsed_area_detected = false;
    void compute_determinants();

    void update_local_size();
    std::vector<double> _local_size;

    std::vector<bool> _predicate_is_positive;
    unsigned _predicates_nb_invalid_steps = 0;
    unsigned evaluate_exact_predicates(Eigen::VectorXd const *x = nullptr);
    double _predicate_infinite_energy = 1e100;

    void gather_energy_gradient(Eigen::VectorXd &g) const;

    Validation_query _validation_query = nullptr;
private:
    // untangling / smoothing
    void update_untangling_eps(double decrease_rate);
    double _untangling_eps;
    double _untangling_ref_eps = 1e-3;
    double _NO_UNTANGLING_EPS = 1e-30;
    double _untangling_max_energy = 1.;
    double regularized_untangling_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr);
    double power_mips_untangling_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr);
    double mips_untangling_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr);
    double sym_dirichlet_untangling_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr);

    inline double untangling_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr) {
        if (exact_predicate_linesearch_enforcement) {
            unsigned nb_invalid = evaluate_exact_predicates(&x);
            if (nb_invalid != 0) return _predicate_infinite_energy;
        }
        if (_validation_query != nullptr) {
            if (!_validation_query(x)) return _predicate_infinite_energy;
        }
        return power_mips_untangling_energy(x, g);
    }

private:
    double laplacian_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr);

    Eigen::SparseMatrix<double> compute_diffusion_matrix(double w, bool graph = true);

private:
    // quality maximization
    void update_qis_tau(double decrease_rate);
    double _qis_tau;
    double const _QIS_CONVERGENCE_THRESHOLD = 1e-3;
    double _QIS_INFINITE_ENERGY = 1e60;
    double qis_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr);

private:
    // boundary
    struct Boundary_poly {
        Boundary_poly(unsigned size, Surface_patch_index id)
        : verts(size)
        , surface_id(id)
        , vert_grad(size, Eigen::Vector3d::Zero())
        {}

        // input
        std::vector<unsigned> verts;
        Surface_patch_index surface_id;

        Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
        Eigen::Vector3d pt = Eigen::Vector3d::Zero();
        Eigen::Vector3d prev_center = Eigen::Vector3d::Zero();
        double max_drift = 0;
        double avg_edge_size = 1.;
        double weight = 1.;

        std::vector<Eigen::Vector3d> vert_grad;
    };

    double ACCEPTED_LOCAL_VARIATION_FROM_BOUNDARY = 0.02; // used as weight
    double MAX_DRIFT_BEFORE_UPDATE = 0.3;
    double MAXIMUM_RELATIVE_DISTANCE_RE_WEIGHTING = 1./3.;

    std::vector<Boundary_poly> _bnd_poly;
    std::vector<std::pair<unsigned, std::vector<std::array<unsigned, 2>>>> const *_vert_and_face_corners;
    Boundary_query _boundary_query = nullptr;

    Boundary_batch_query _boundary_batch_query = nullptr;
    bool _boundary_batch_mode = false;
    std::vector<std::vector<Eigen::Vector3d>> _boundary_live_coords;
    std::vector<Surface_patch_index> _boundary_batch_ids;
    std::vector<Plane> _boundary_batch_query_results;

    bool has_bnd_terms() const;

    void update_boundary_info(Eigen::VectorXd const &x, bool reset = false);
    double boundary_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr);

private:
    // curve and target positions

    struct Edge_data {
        Edge_data(std::array<unsigned, 2> verts_, Curve_index curve_id_)
        : verts(verts_)
        , curve_id(curve_id_)
        {}
        std::array<unsigned, 2> verts;
        Curve_index curve_id;

        Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
        Eigen::Vector3d pt = Eigen::Vector3d::Zero();
        Eigen::Vector3d prev_center = Eigen::Vector3d::Zero();
        double max_drift = 0;
        double weight = 1.;
    };

    std::vector<Edge_data> _edge_data;
    Curve_query _curve_query = nullptr;

    Curve_batch_query _curve_batch_query = nullptr;
    bool _curve_batch_mode = false;
    std::vector<std::array<Eigen::Vector3d, 2>> _edge_live_coords;
    std::vector<Curve_index> _edge_batch_ids;
    std::vector<Curve_tangent> _edge_batch_query_results;

    std::vector<double> _point_prev_target_weight;
    std::vector<std::tuple<unsigned, Eigen::Vector3d, double>> _point_target_position;

    bool has_curves_and_points_terms() const;
    void update_curves_and_points_info(Eigen::VectorXd const &x, bool reset = false);
    double curves_and_points_energies(Eigen::VectorXd const &x, Eigen::VectorXd *g = nullptr);

};

// ==============================================================================================//
//                                           Implementation                                      //
// ==============================================================================================//

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::Tetrahedral_mesh_smoother(
    Eigen::VectorXd &coords,
    std::vector<bool> const &locks,
    std::vector<std::array<unsigned, 4>> const &tetrahedra,
    std::vector<std::array<Eigen::Vector3d, 4>> const &tet_inv_grad,
    std::vector<std::vector<unsigned>> const &vert2tet_corner
)
: _coords(coords)
, _locks(locks)
, _tet_storage(tetrahedra.size())
, _vert2tet_corner(vert2tet_corner)
, _determinants(tetrahedra.size(), 0)
, _conformal_energies(tetrahedra.size(), 1)
, _local_size(nb_vertices(), min_valid_edge_size)
, _predicate_is_positive(tetrahedra.size(), true)
{
    for (unsigned t = 0; t < tetrahedra.size(); ++t) {
        _tet_storage[t].verts = tetrahedra[t];
        _tet_storage[t].ig = tet_inv_grad[t];
    }
    compute_determinants();
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::set_boundary_without_query(
    std::vector<std::vector<unsigned>> const &bnd_faces,
    std::vector<std::pair<unsigned, std::vector<std::array<unsigned, 2>>>> const *vert_and_face_corners,
    std::vector<Surface_patch_index> const &ids
)
{
    _vert_and_face_corners = vert_and_face_corners;
    _bnd_poly.clear();
    _bnd_poly.reserve(bnd_faces.size());
    _boundary_live_coords.clear();
    _boundary_live_coords.reserve(bnd_faces.size());
    _boundary_batch_ids.clear();
    _boundary_batch_ids.reserve(bnd_faces.size());
    _boundary_batch_query_results.clear();
    _boundary_batch_query_results.reserve(bnd_faces.size());
    assert(ids.size() == bnd_faces.size());
    for (unsigned f = 0; f < bnd_faces.size(); ++f) {
        if (bnd_faces[f].empty()) continue;
        _bnd_poly.push_back(Boundary_poly(static_cast<unsigned>(bnd_faces[f].size()), ids[f]));
        _bnd_poly.back().verts = bnd_faces[f];
        _boundary_live_coords.push_back(std::vector<Eigen::Vector3d>(bnd_faces[f].size()));
        _boundary_batch_ids.push_back(ids[f]);
    }
}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::set_boundary_with_singular_query(
    std::vector<std::vector<unsigned>> const &bnd_faces,
    std::vector<std::pair<unsigned, std::vector<std::array<unsigned, 2>>>> const *vert_and_face_corners,
    Boundary_query boundary_query,
    std::vector<Surface_patch_index> const &ids
)
{
    set_boundary_without_query(bnd_faces, vert_and_face_corners, ids);
    _boundary_batch_mode = false;
    _boundary_query = boundary_query;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::set_boundary_with_batch_query(
    std::vector<std::vector<unsigned>> const &bnd_faces,
    std::vector<std::pair<unsigned, std::vector<std::array<unsigned, 2>>>> const *vert_and_face_corners,
    Boundary_batch_query boundary_query,
    std::vector<Surface_patch_index> const &ids
)
{
    set_boundary_without_query(bnd_faces, vert_and_face_corners, ids);
    _boundary_batch_mode = true;
    _boundary_batch_query = boundary_query;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::set_curve_network_with_singular_query(
        std::vector<std::array<unsigned, 2>> const &edges,
        std::vector<Curve_index> const &edge_ids,
        Curve_query query
)
{
    assert(edges.size() == edge_ids.size());
    _edge_data.reserve(edges.size());
    for (unsigned i = 0; i < edges.size(); ++i) {
        _edge_data.push_back(Edge_data(edges[i], edge_ids[i]));
    }
    _curve_batch_mode = false;
    _curve_query = query;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::set_curve_network_with_batch_query(
        std::vector<std::array<unsigned, 2>> const &edges,
        std::vector<Curve_index> const &edge_ids,
        Curve_batch_query query
)
{
    assert(edges.size() == edge_ids.size());
    _edge_data.reserve(edges.size());
    for (unsigned i = 0; i < edges.size(); ++i) {
        _edge_data.push_back(Edge_data(edges[i], edge_ids[i]));
    }
    _curve_batch_mode = true;
    _curve_batch_query = query;
    _edge_live_coords.resize(edges.size());
    _edge_batch_ids = edge_ids;
    _edge_batch_query_results.resize(edges.size());
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::set_quadratic_target_positions(
        std::vector<std::tuple<unsigned, Eigen::Vector3d, double>> const &targets
)
{
    _point_prev_target_weight.resize(targets.size());
    _point_target_position = targets;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::set_validation_query(Validation_query query) {
    _validation_query = query;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline Eigen::Matrix3d Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::Tet_storage::compute_jacobian(Eigen::VectorXd const &coords) const {
        return    ig[0] * Math_functions::sub_line_vector(coords,verts[0])
                + ig[1] * Math_functions::sub_line_vector(coords,verts[1])
                + ig[2] * Math_functions::sub_line_vector(coords,verts[2])
                + ig[3] * Math_functions::sub_line_vector(coords,verts[3]);
}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::compute_determinants()
{
    double min_trace = min_valid_edge_size * min_valid_edge_size;

    struct Reduction {
        double det_min = (std::numeric_limits<double>::max)();
        unsigned nb_inverted = 0;
        double conformal_energy_max = 0;
        bool collapsed_area_detected = false;

        void join(Reduction const &other) {
            det_min = (std::min)(det_min, other.det_min);
            nb_inverted += other.nb_inverted;
            conformal_energy_max = (std::max)(conformal_energy_max, other.conformal_energy_max);
            collapsed_area_detected = collapsed_area_detected || other.collapsed_area_detected;
        }
    };

    Reduction reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _tet_storage.size(), [&](std::size_t t, Reduction& r) {
        Tet_storage const &tet = _tet_storage[t];

        Eigen::Matrix3d J = tet.compute_jacobian(_coords);
        _determinants[t] = J.determinant();
        r.nb_inverted += _determinants[t] <= 0;
        r.det_min = (std::min)(_determinants[t], r.det_min);
        double trace = J.squaredNorm();
        r.collapsed_area_detected = r.collapsed_area_detected || (trace < min_trace);

        if (_determinants[t] <= 0) {
            _conformal_energies[t] = (std::numeric_limits<double>::max)();
            r.conformal_energy_max = _conformal_energies[t];
            return;
        }
        double power_trace = std::pow(trace, 3./2.);
        _conformal_energies[t] = power_trace / _determinants[t];
        r.conformal_energy_max = (std::max)(_conformal_energies[t], r.conformal_energy_max);
    });

    _det_min = reduction.det_min;
    _nb_inverted = reduction.nb_inverted;
    _conformal_energy_max = reduction.conformal_energy_max;
    _collapsed_area_detected = reduction.collapsed_area_detected;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
Eigen::SparseMatrix<double> Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::compute_diffusion_matrix(double w, bool graph) {
    unsigned n = 3 * nb_vertices();
    double l = w / (1+w);

    Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, _tet_storage.size(), [&](std::size_t t) {
        Tet_storage &tet = _tet_storage[t];
        Eigen::Matrix3d J = tet.compute_jacobian(_coords);
        Eigen::Matrix3d dfdJ = 2*J.transpose();
        for (unsigned v = 0; v < 4; ++v) {
            tet.vert_grad[v] = dfdJ * tet.ig[v];
        }
    });

    Eigen::SparseMatrix<double> H(n,n);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(4*4*3*_tet_storage.size());

    // identity term
    for (unsigned v = 0; v < nb_vertices(); ++v) {
        for (unsigned d = 0; d < 3; ++d) {
            unsigned v_d = 3*v+d;
            triplets.push_back({static_cast<int>(v_d), static_cast<int>(v_d), 1-l});
        }
    }
    if (graph) {
        std::set<unsigned> neigh;
        for (unsigned v = 0; v < nb_vertices(); ++v) {
            neigh.clear();
            for (unsigned tc : _vert2tet_corner[v]) {
                unsigned tet_id = tc / 4;
                Tet_storage const &tet = _tet_storage[tet_id];
                unsigned loc_id = tc % 4;
                for (unsigned i = 0; i < 4;++i) {
                    if (loc_id == i) continue;
                    neigh.insert(tet.verts[i]);
                }
            }

            for (unsigned d = 0; d < 3; ++d) {
                unsigned v_d = 3*v+d;
                triplets.push_back({static_cast<int>(v_d), static_cast<int>(v_d), l*neigh.size()});
                for (unsigned other_v : neigh) {
                    int other_v_d = 3*other_v+d;
                    triplets.push_back({static_cast<int>(v_d), static_cast<int>(other_v_d), -l});
                }
            }
        }
    }
    else {
        for (unsigned v = 0; v < nb_vertices(); ++v) {
            for (unsigned tc : _vert2tet_corner[v]) {
                unsigned tet_id = tc / 4;
                Tet_storage const &tet = _tet_storage[tet_id];
                unsigned loc_id = tc % 4;
                for (unsigned d = 0; d < 3; ++d) {
                    unsigned v_d = 3*v+d;
                    for (unsigned i = 0; i < 4;++i) {
                        int other_v_d = 3*tet.verts[i]+d;
                        double lapl_w = l*2*tet.ig[i].dot(tet.ig[loc_id]);
                        triplets.push_back({static_cast<int>(v_d), static_cast<int>(other_v_d), lapl_w});
                    }
                }
            }
        }
    }
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::update_local_size() {
    Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, nb_vertices(), [&](std::size_t v) {
        // does geometric average on all edge of the vertex
        // pass to log to avoid double overflow with high number of multiplications
        double log_avg = 0.;
        Eigen::Vector3d x = Math_functions::sub_col_vector(_coords, v);
        for (unsigned tc : _vert2tet_corner[v]) {
            for (unsigned j = 1; j < 4; ++j){
                Eigen::Vector3d y = Math_functions::sub_col_vector(_coords, _tet_storage[tc/4].verts[(tc%4+j)%4]);
                log_avg += std::log10((std::max)((x-y).norm(), min_valid_edge_size));
            }
        }
        double nbEdges = static_cast<double>(_vert2tet_corner[v].size()) * 3;
        log_avg /= nbEdges;
        _local_size[v] = std::pow(10., log_avg);
    });

    double min_allowed_det = min_valid_edge_size * min_valid_edge_size * min_valid_edge_size;

    Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, _tet_storage.size(), [&](std::size_t t) {
        Tet_storage &tet = _tet_storage[t];
        double avg_edge_size = 1;
        for(unsigned i = 0; i < 4; ++i) {
            avg_edge_size *= _local_size[tet.verts[i]];
        }
        avg_edge_size = std::pow(avg_edge_size, 1./4.);
        tet.local_edge_size = avg_edge_size;
        tet.det_estimation = (std::max)(avg_edge_size * avg_edge_size * avg_edge_size, min_allowed_det);
    });
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
unsigned Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::evaluate_exact_predicates(Eigen::VectorXd const *x) {

    Eigen::VectorXd const & coords = x == nullptr ? _coords : *x;

    struct Reduction {
        unsigned invalid_tet = 0;
        void join(Reduction const &other) { invalid_tet += other.invalid_tet; }
    };

    Reduction reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _tet_storage.size(), [&](unsigned t, Reduction &reduction) {
        Tet_storage const &tet = _tet_storage[t];
        bool exact_check = Math_functions::strictly_positive_tetrahedra({
            Math_functions::sub_line_vector(coords,tet.verts[0]),
            Math_functions::sub_line_vector(coords,tet.verts[1]),
            Math_functions::sub_line_vector(coords,tet.verts[2]),
            Math_functions::sub_line_vector(coords,tet.verts[3])
        });
        if (x == nullptr) _predicate_is_positive[t] = exact_check;
        if (!exact_check) ++reduction.invalid_tet;
    });

    return reduction.invalid_tet;
}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
bool Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::run_callback(OPTIMIZATION_TYPE opt_type, unsigned iter, LBFGS_status lbfgs_status, Eigen::VectorXd const *g) {
    if (callback_function == nullptr) return false;
    if (callback_setting == DEBUG_CALLBACK_SETTING::NOTHING) return false;
    if (callback_setting == DEBUG_CALLBACK_SETTING::OUTER_ITER && lbfgs_status.enabled()) return false;

    _callback_vert_storage.resize(nb_vertices());
    _callback_tet_storage.resize(_tet_storage.size());
    _callback_smoothing_gradient.resize(_coords.size());
    _callback_boundary_gradient.resize(_coords.size());

    _callback_status.opt = opt_type;
    _callback_status.boundary_enabled = has_bnd_terms();
    _callback_status.min_edge_size = min_valid_edge_size;
    _callback_status.boundary_weight = boundary_weight;

    _callback_status.outer_iter_nb = iter;
    _callback_status.lbfgs_status = lbfgs_status;

    _callback_status.smoothing_energy =  opt_type == UNTANGLING ? untangling_energy(_coords, &_callback_smoothing_gradient)
                                      : (opt_type == STIFFENING ? qis_energy(_coords, &_callback_smoothing_gradient)
                                      :                           laplacian_energy(_coords, &_callback_smoothing_gradient));

    _callback_status.boundary_energy = has_bnd_terms() ? boundary_energy(_coords, &_callback_boundary_gradient): 0.;
    compute_determinants();

    _callback_status.min_det = _det_min;
    _callback_status.nb_inverted = _nb_inverted;
    _callback_status.opt_parameter = opt_type == UNTANGLING ? _untangling_eps : (opt_type == STIFFENING ? _qis_tau : 0.);

    Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, nb_vertices(), [&](std::size_t v) {
        _callback_vert_storage[v].lock[0] = _locks[3*v+0];
        _callback_vert_storage[v].lock[1] = _locks[3*v+1];
        _callback_vert_storage[v].lock[2] = _locks[3*v+2];
        _callback_vert_storage[v].local_edge_size = _local_size[v];
        _callback_vert_storage[v].smoothing_gradient = Math_functions::sub_line_vector(_callback_smoothing_gradient, v);
        _callback_vert_storage[v].boundary_gradient = Math_functions::sub_line_vector(_callback_boundary_gradient, v);
        if (g != nullptr) _callback_vert_storage[v].lbfgs_gradient = Math_functions::sub_line_vector(*g, v);
    });

    Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, _tet_storage.size(), [&](std::size_t t) {
        _callback_tet_storage[t].energy_value = _tet_storage[t].fval;
        _callback_tet_storage[t].weight = _tet_storage[t].local_edge_size;
        _callback_tet_storage[t].det = _determinants[t];
    });

    return callback_function(_callback_status, _callback_vert_storage, _callback_tet_storage);
}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::gather_energy_gradient(Eigen::VectorXd &g) const {
    Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, nb_vertices(), [&](std::size_t v) {
        for (unsigned d = 0; d < 3; ++d) {
            if (_locks[3*v+d]) continue;
            for (unsigned tc : _vert2tet_corner[v]) {
                g(3*v+d) += _tet_storage[tc/4].vert_grad[tc%4](d);
            }
        }
    });

}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::update_untangling_eps(double decrease_rate) {
    // mix of foldover and 1999 epsilons.
    if (_det_min > 0) {
        _untangling_eps = _NO_UNTANGLING_EPS;
        return;
    }

    struct Reduction {
        double weighted_det_min = (std::numeric_limits<double>::max)();
        void join(const Reduction& other) {
            weighted_det_min = (std::min)(weighted_det_min, other.weighted_det_min);
        }
    };

    auto reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _tet_storage.size(), [&](std::size_t t, Reduction &r) {
        r.weighted_det_min = (std::min)(r.weighted_det_min, _determinants[t] / _tet_storage[t].det_estimation);
        return r;
    });
    double weighted_det_min = reduction.weighted_det_min;

    double _1999_eps = std::sqrt(1e-18 + 4*1e-2* weighted_det_min*weighted_det_min);

    constexpr double forced_decrease_rate = 0.1;
    double sigma = (std::max)(forced_decrease_rate, 1 - decrease_rate);
    double mu = (1-sigma) * Math_functions::chi(_untangling_eps, weighted_det_min);
    double foldover_eps = 2 * std::sqrt(mu*(mu - weighted_det_min));

    _untangling_eps = (std::min)(foldover_eps, _1999_eps);

    // todo: I don't remember how I obtained this formula, but it was working well for hard untangling. To study again.
    // double c = sqrt(_untangling_eps*_untangling_eps+weighted_det_min*weighted_det_min);
    // double decrease = (1-sigma*c/(std::abs(weighted_det_min) + c));
    // double custom_eps = decrease * _untangling_eps;
    // _untangling_eps = (std::min)(_untangling_eps, custom_eps);
    return;

}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline double Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::mips_untangling_energy(Eigen::VectorXd const &x, Eigen::VectorXd *grad) {

    struct Reduction {
        double untangling_max_energy = 1.;
        double F = 0;
        void join(const Reduction& other) {
            untangling_max_energy = (std::max)(untangling_max_energy, other.untangling_max_energy);
            F += other.F;
        }
    };

    Reduction reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _tet_storage.size(), [&](std::size_t t, Reduction& r) {
        Tet_storage &tet = _tet_storage[t];
        if (tet.skip) return;
        double scaled_epsilon = tet.det_estimation*_untangling_eps;
        double weight = tet.local_edge_size;

        Eigen::Matrix3d J = tet.compute_jacobian(x);
        double det = J.determinant();

        double c1 = Math_functions::chi(scaled_epsilon, det);
        double inv_c1 = 1./c1;
        double cbrt_c1 = std::cbrt(c1);
        double trace = J.squaredNorm();

        double f =  trace * cbrt_c1 * inv_c1; //f / std::pow(c1, 2./3.);

        tet.fval = f;
        r.untangling_max_energy = (std::max)(r.untangling_max_energy, tet.fval);

        r.F +=  weight * f;

        if (grad == nullptr) return;

        double c3 = Math_functions::chi_deriv(scaled_epsilon, det);

        Eigen::Matrix3d K = Math_functions::dual_basis(J);

        Eigen::Matrix3d dfdJ = J * (weight * 2. * cbrt_c1 * inv_c1)
                             - K * (weight * (2./3.) * f * c3 * inv_c1 );

        dfdJ.transposeInPlace();
        for (unsigned v = 0; v < 4; ++v) {
            tet.vert_grad[v] = dfdJ * tet.ig[v];
        }
    });


    if (grad != nullptr)
        gather_energy_gradient(*grad);

    _untangling_max_energy = reduction.untangling_max_energy;

    return reduction.F;
}



template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline double Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::power_mips_untangling_energy(Eigen::VectorXd const &x, Eigen::VectorXd *grad) {
    struct Reduction {
        double untangling_max_energy = 1.;
        double F = 0;
        void join(const Reduction& other) {
            untangling_max_energy = (std::max)(untangling_max_energy, other.untangling_max_energy);
            F += other.F;
        }
    };

    Reduction reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _tet_storage.size(), [&](std::size_t t, Reduction& r) {
        Tet_storage &tet = _tet_storage[t];
        if (tet.skip) return;
        double scaled_epsilon = tet.det_estimation*_untangling_eps;
        double weight = tet.local_edge_size;

        Eigen::Matrix3d J = tet.compute_jacobian(x);
        double det = J.determinant();

        double c1 = Math_functions::chi(scaled_epsilon, det);
        double trace = J.squaredNorm();
        double root_trace = std::sqrt(trace);
        double pow_trace = trace*root_trace; //std::pow(trace, 3. / 2.);

        double f = pow_trace / c1;

        tet.fval = f;
        r.untangling_max_energy = (std::max)(r.untangling_max_energy, tet.fval);

        r.F +=  weight * f;

        if (grad == nullptr) return;

        double c2 = Math_functions::chi_deriv(scaled_epsilon, det);

        Eigen::Matrix3d K = Math_functions::dual_basis(J);

        Eigen::Matrix3d dfdJ = J * (weight * 3. * root_trace / c1)
                             - K * (weight * (f * c2/c1));

        dfdJ.transposeInPlace();
        for (unsigned v = 0; v < 4; ++v) {
            tet.vert_grad[v] = dfdJ * tet.ig[v];
        }
    });

    if (grad != nullptr)
        gather_energy_gradient(*grad);

    _untangling_max_energy = reduction.untangling_max_energy;

    return reduction.F;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline double Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::regularized_untangling_energy(Eigen::VectorXd const &x, Eigen::VectorXd *grad) {
    double det_regularization_weight = 1e-3;

    struct Reduction {
        double untangling_max_energy = 1.;
        double F = 0;
        void join(const Reduction& other) {
            untangling_max_energy = (std::max)(untangling_max_energy, other.untangling_max_energy);
            F += other.F;
        }
    };

    Reduction reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _tet_storage.size(), [&](std::size_t t, Reduction& r) {
        Tet_storage &tet = _tet_storage[t];
        if (tet.skip) return;
        double delta = 1./tet.det_estimation;
        double scaled_epsilon = tet.det_estimation*_untangling_eps;
        double weight = tet.local_edge_size;

        Eigen::Matrix3d J = tet.compute_jacobian(x);
        double det = J.determinant();

        double c1 = Math_functions::chi(scaled_epsilon, det);
        double trace = J.squaredNorm();
        double root_trace = std::sqrt(trace);
        double pow_trace = trace*root_trace; //std::pow(trace, 3. / 2.);

        double f = pow_trace / c1;

        tet.fval = f;
        r.untangling_max_energy = (std::max)(r.untangling_max_energy, tet.fval);

        double d_g = delta * det;
        double c1_g = Math_functions::chi(_untangling_eps, d_g);


        double g = (1+d_g*d_g) / c1_g;

        r.F +=  weight * (f + det_regularization_weight * g);

        if (grad == nullptr) return;

        double c2 = Math_functions::chi_deriv(scaled_epsilon, det);
        double c2_g = Math_functions::chi_deriv(_untangling_eps, d_g);

        Eigen::Matrix3d K = Math_functions::dual_basis(J);

        Eigen::Matrix3d dfdJ = J * (weight * 3. * root_trace / c1)
                             - K * (weight * (f * c2/c1
                                             + det_regularization_weight * delta * (g * c2_g - 2*d_g)/c1_g));

        dfdJ.transposeInPlace();
        for (unsigned v = 0; v < 4; ++v) {
            tet.vert_grad[v] = dfdJ * tet.ig[v];
        }
    });
    if (grad != nullptr)
        gather_energy_gradient(*grad);

    _untangling_max_energy = reduction.untangling_max_energy;

    return reduction.F;
}



template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline double Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::laplacian_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g) {
    struct Reduction {
        double F = 0;
        void join(const Reduction& other) {
            F += other.F;
        }
    };
    Reduction reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _tet_storage.size(), [&](std::size_t t, Reduction& r) {
        Tet_storage &tet = _tet_storage[t];
        Eigen::Matrix3d J = tet.compute_jacobian(x);
        double f = J.squaredNorm();
        r.F += f;
        tet.fval = f;

        if (g == nullptr) return;

        Eigen::Matrix3d dfdJ = 2*J.transpose();
        for (unsigned v = 0; v < 4; ++v) {
            tet.vert_grad[v] = dfdJ * tet.ig[v];
        }
    });
    if (g != nullptr)
        gather_energy_gradient(*g);

    return reduction.F;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::update_qis_tau(double decrease_rate) {
    double sigma = (std::max)(1.-decrease_rate, 1e-1);
    _qis_tau = _qis_tau + sigma*(1.-_qis_tau*_conformal_energy_max)/_conformal_energy_max;
}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline double Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::qis_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g) {
    struct Reduction {
        double F = 0;
        unsigned above_max_qual = 0;
        void join(const Reduction& other) {
            F += other.F;
            above_max_qual += other.above_max_qual;
        }
    };
    Reduction reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _tet_storage.size(), [&](std::size_t t, Reduction& r) {
        Tet_storage &tet = _tet_storage[t];

        Eigen::Matrix3d J = tet.compute_jacobian(x);
        double det = J.determinant();
        if (det <= 0) {
            ++r.above_max_qual;
            return;
        }
        double trace = J.squaredNorm();
        double root_trace = std::sqrt(trace);
        double pow_trace = trace*root_trace; //std::pow(trace, 3. / 2.);

        double f = pow_trace / det;
        r.above_max_qual += _qis_tau * f >= 1.;
        if (r.above_max_qual) return;

        double denom = (1-_qis_tau*f);
        double qis = f / denom;
        tet.fval = qis;
        r.F +=  tet.local_edge_size * qis;
        if (g == nullptr) return;
        Eigen::Matrix3d K = Math_functions::dual_basis(J);
        double dqis_dx = 1./(denom*denom);
        Eigen::Matrix3d dqis_dJ = J * (tet.local_edge_size * dqis_dx * 3. * root_trace / det) - K * (tet.local_edge_size * dqis_dx * f / det);
        dqis_dJ.transposeInPlace();
        for (unsigned v = 0; v < 4; ++v) {
            tet.vert_grad[v] = dqis_dJ * tet.ig[v];
        }
    });
    if (reduction.above_max_qual > 0 || reduction.F > _QIS_INFINITE_ENERGY) {
        if (g != nullptr) g->setZero();
        return _QIS_INFINITE_ENERGY;
    }

    if (g != nullptr)
        gather_energy_gradient(*g);

    return reduction.F;
}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline bool Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::has_bnd_terms() const {
    bool has_terms = !_bnd_poly.empty();
    bool has_query = _boundary_batch_mode ? (_boundary_batch_query != nullptr) : (_boundary_query != nullptr);
    return has_terms && has_query;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::update_boundary_info(Eigen::VectorXd const &x, bool reset) {
    if (!has_bnd_terms()) return;

    auto update_poly_coord = [&](unsigned t) {
        Boundary_poly const &poly = _bnd_poly[t];
        for (unsigned i=0; i<poly.verts.size(); ++i) {
            _boundary_live_coords[t][i] = Math_functions::sub_col_vector(x, poly.verts[i]);
        }
    };

    if (_boundary_batch_mode) {
        Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, _bnd_poly.size(), [&](std::size_t t) {
            update_poly_coord(t);
        });
        _boundary_batch_query(_boundary_live_coords, _boundary_batch_ids, _boundary_batch_query_results);
    }

    Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, _bnd_poly.size(), [&](std::size_t t) {
        Boundary_poly &poly = _bnd_poly[t];
        Eigen::Vector3d center = Eigen::Vector3d::Zero();
        for (unsigned v : poly.verts) {
            center += Math_functions::sub_col_vector(x, v);
        }
        center /= static_cast<double>(poly.verts.size());
        if (!reset && (center - poly.prev_center).norm()<poly.max_drift) return;

        poly.prev_center = center;
        poly.avg_edge_size = 0.;
        for (unsigned i = 0; i < poly.verts.size(); ++i) {
            poly.avg_edge_size += log10(_local_size[poly.verts[i]]);
        }
        poly.avg_edge_size /= poly.verts.size();
        poly.avg_edge_size = std::pow(10., poly.avg_edge_size);

        if (!_boundary_batch_mode) update_poly_coord(t);
        auto [pt, n, weight] = _boundary_batch_mode ? _boundary_batch_query_results[t] : _boundary_query(_boundary_live_coords[t], poly.surface_id);

        poly.max_drift = MAX_DRIFT_BEFORE_UPDATE*poly.avg_edge_size;
        poly.A = weight* n * n.transpose();
        poly.pt = pt;
        if (reset) {
            double max_poly_relative_dist = 0.;
            for (unsigned v : poly.verts) {
                Eigen::Vector3d dir = poly.pt - Math_functions::sub_col_vector(x, v);
                double dist = dir.transpose() * poly.A * dir;
                dist /= _local_size[v]*_local_size[v];
                max_poly_relative_dist = (std::max)(max_poly_relative_dist, dist);
            }
            poly.weight = 1./(std::max)(MAXIMUM_RELATIVE_DISTANCE_RE_WEIGHTING,max_poly_relative_dist);
        }
    });

}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline double Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::boundary_energy(Eigen::VectorXd const &x, Eigen::VectorXd *g) {
    if (!has_bnd_terms()) return 0.;
    double drift = boundary_weight/(ACCEPTED_LOCAL_VARIATION_FROM_BOUNDARY*ACCEPTED_LOCAL_VARIATION_FROM_BOUNDARY);
    struct Reduction {
        double F = 0.;
        void join(const Reduction& other) {
            F += other.F;
        }
    };
    Reduction reduction = Mesh_smoothing_3_internal::reduce<ConcurrencyTag, Reduction>(0, _bnd_poly.size(), [&](std::size_t t, Reduction& r) {
        Boundary_poly &poly = _bnd_poly[t];
        for (unsigned tc = 0; tc < poly.verts.size(); ++tc) {
            double regul = poly.weight * drift/_local_size[poly.verts[tc]];
            Eigen::Vector3d pt = Math_functions::sub_col_vector(x, poly.verts[tc]);
            Eigen::Vector3d dn = pt-poly.pt;
            double dist = (dn.transpose() * poly.A * dn);
            r.F += regul * dist;
            if (g == nullptr) continue;
            poly.vert_grad[tc] = regul * (2 * poly.A *dn);
        }
    });
    if (g != nullptr) {
        Mesh_smoothing_3_internal::for_each<ConcurrencyTag>(0, _vert_and_face_corners->size(), [&](std::size_t i) {
            unsigned v = (*_vert_and_face_corners)[i].first;
            for (unsigned d = 0; d < 3; ++d) {
                if (_locks[3*v+d]) continue;
                for (auto [t, tc] : (*_vert_and_face_corners)[i].second) {
                    (*g)(3*v+d) += _bnd_poly[t].vert_grad[tc](d);
                }
            }
        });
    }
    return reduction.F;
}


template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline bool Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::has_curves_and_points_terms() const {
    return !_edge_data.empty() || !_point_prev_target_weight.empty();
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::update_curves_and_points_info(Eigen::VectorXd const &x, bool reset) {
   if (!has_curves_and_points_terms()) return;


    if (_curve_batch_mode) {
        for (unsigned e = 0; e < _edge_data.size(); ++e) {
            _edge_live_coords[e][0] = Math_functions::sub_col_vector(x, _edge_data[e].verts[0]);
            _edge_live_coords[e][1] = Math_functions::sub_col_vector(x, _edge_data[e].verts[1]);
        }
        _curve_batch_query(_edge_live_coords, _edge_batch_ids, _edge_batch_query_results);
    }

    for (unsigned e = 0; e < _edge_data.size(); ++e) {
        Edge_data &edge = _edge_data[e];
        Eigen::Vector3d pt0 = Math_functions::sub_col_vector(x, edge.verts[0]);
        Eigen::Vector3d pt1 = Math_functions::sub_col_vector(x, edge.verts[1]);
        Eigen::Vector3d center = 0.5*(pt0+pt1);

        if (!reset && (center - edge.prev_center).norm()<edge.max_drift) continue;

        edge.prev_center = center;
        double avg_edge_size = std::sqrt(_local_size[edge.verts[0]]*_local_size[edge.verts[1]]);

        auto [pt, n, weight] = _curve_batch_mode ? _edge_batch_query_results[e] : _curve_query({pt0, pt1}, edge.curve_id);

        edge.max_drift = MAX_DRIFT_BEFORE_UPDATE*avg_edge_size;
        edge.A = weight * (Eigen::Matrix3d::Identity() - n * n.transpose());
        edge.pt = pt;

        if (reset)  {
            double max_edge_relative_dist = 0.;
            for (unsigned v : edge.verts) {
                Eigen::Vector3d dir = edge.pt - Math_functions::sub_col_vector(x, v);
                double dist = dir.transpose() * edge.A * dir;
                dist /= _local_size[v]*_local_size[v];
                max_edge_relative_dist = (std::max)(max_edge_relative_dist, dist);
            }
            edge.weight = 1./(std::max)(MAXIMUM_RELATIVE_DISTANCE_RE_WEIGHTING,max_edge_relative_dist);
        }

    }

    for (unsigned i = 0; i < _point_target_position.size(); ++i) {
        double curr_dist = (Math_functions::sub_col_vector(x, std::get<0>(_point_target_position[i])) - std::get<1>(_point_target_position[i])).squaredNorm();
         if (reset) _point_prev_target_weight[i] = 1./(std::max)(MAXIMUM_RELATIVE_DISTANCE_RE_WEIGHTING,curr_dist);
    }
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline double Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::curves_and_points_energies(Eigen::VectorXd const &x, Eigen::VectorXd *g) {
    if (!has_curves_and_points_terms()) return 0.;
    double F = 0;
    double drift = boundary_weight/(ACCEPTED_LOCAL_VARIATION_FROM_BOUNDARY*ACCEPTED_LOCAL_VARIATION_FROM_BOUNDARY);
    for (unsigned e = 0; e < _edge_data.size(); ++e) {
        Edge_data &edge = _edge_data[e];
        for (unsigned ev = 0; ev < 2; ++ev) {
            unsigned v = edge.verts[ev];
            double regul = edge.weight * drift/_local_size[v];
            Eigen::Vector3d pt = Math_functions::sub_col_vector(x, v);
            Eigen::Vector3d dn = pt-edge.pt;
            double dist = (dn.transpose() * edge.A * dn);
            F += regul * dist;
            if (g == nullptr) continue;
            Eigen::Vector3d grad = regul * (2 * edge.A *dn);
            for (unsigned d = 0; d < 3; ++d) {
                if (_locks[3*v+d]) continue;
                (*g)(3*v+d) += grad[d];
            }
        }
    }

    for (unsigned i = 0; i < _point_target_position.size(); ++i) {
        unsigned v = std::get<0>(_point_target_position[i]);
        Eigen::Vector3d target = std::get<1>(_point_target_position[i]);
        double weight = std::get<2>(_point_target_position[i]) * _point_prev_target_weight[i] * drift / _local_size[v];

        Eigen::Vector3d dn = Math_functions::sub_col_vector(x, v) - target;
        F += weight * dn.squaredNorm();

        if (g == nullptr) continue;
        Eigen::Vector3d grad = weight * (2*dn);
        for (unsigned d = 0; d < 3; ++d) {
            if (_locks[3*v+d]) continue;
            (*g)(3*v+d) += grad[d];
        }

    }

    return F;
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
void Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::run_laplacian_gradient_descent(unsigned max_number_iter) {
    if (verbose) std::cout << "==== Tetrahedral_mesh_smoother  ====" << "\n";
    if (verbose) std::cout << "----   running laplacian smoothing    ----" << "\n";
    if (verbose) std::cout << "Nb of optimization variables (vertices x3): " << _coords.size()  << std::endl;
    if (verbose) std::cout << "Nb of energy terms (tetrahedra): " << _tet_storage.size() << std::endl;
    if (verbose) std::cout << "Nb of used OpenMP core: " << Eigen::nbThreads() << std::endl;
    if (verbose) std::cout << "Nb iterations: " << max_number_iter << std::endl;
    Time_log logging("Tetrahedral_mesh_smoother:: laplacian");
    if (verbose && has_bnd_terms()) std::cout << "WARNING: Laplacian run ignores boundaries"  << std::endl;

    number_of_outer_iter = 1;
    number_of_lbfgs_iter = 0;

    run_callback(LAPLACIAN, 0);
    Function_minimizer opt ([&](Eigen::VectorXd const &x, Eigen::VectorXd &g) -> double {
        g.setZero();
        return laplacian_energy(x, &g);
    });
    opt._call_back = [&](Eigen::VectorXd const &, Eigen::VectorXd const & g,double,double step,unsigned iter,unsigned nbEval) {
        bool stop_required = run_callback(LAPLACIAN, 0, {iter, step, nbEval}, &g);
        ++number_of_lbfgs_iter;
        if (verbose) std::cout << "." << std::flush;
        if (verbose && stop_required) std::cout << "--- CALLBACK REQUIRED STOPPING ---" << std::flush;
        return stop_required;
    };
    opt.set_max_iter(max_number_iter);
    opt._parameters.min_step = 1e-40;
    opt._parameters.past = 8;
    opt._parameters.mem_size = 8;
    if (verbose) std::cout << " Running laplacian smoothing: " << std::endl;
    opt.lbfgs_optimize(_coords);
    if (verbose) std::cout << " done." << std::endl;
    compute_determinants();
    if (verbose) logging.log_total_time();
}

template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline bool Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::run_untangling(unsigned max_number_iter) {
    if (verbose) std::cout << "==== Tetrahedral_mesh_smoother untangling ====" << "\n";
    if (_coords.size() == 0 || _tet_storage.empty()) {
        if (verbose) std::cout << "No variables to optimize." << std::endl;
        return true;
    }

    number_of_outer_iter = 0;
    number_of_lbfgs_iter = 0;


    Time_log logging("Tetrahedral_mesh_smoother");
    compute_determinants();
    update_local_size();
    update_boundary_info(_coords, true);
    update_curves_and_points_info(_coords, true);
    if (verbose) std::cout << "Nb of optimization variables (vertices x3): " << _coords.size()  << std::endl;
    if (verbose) std::cout << "Nb of energy terms (tetrahedra): " << _tet_storage.size() << std::endl;
    if (verbose) std::cout << "Nb of used OpenMP core: " << Eigen::nbThreads( ) << std::endl;

    if (exact_predicate_status) {
        unsigned curr_invalid = evaluate_exact_predicates();
        if (verbose) std::cout << "Exact predicate check: " << curr_invalid << " invalid tetrahedra" << std::endl;
        if (curr_invalid != 0) {
            exact_predicate_linesearch_enforcement = false; // would not move otherwise
            exact_predicate_optimization_check = false; // would be way too talkative
        }
        _predicates_nb_invalid_steps = 0;
    }

    _NO_UNTANGLING_EPS = 1e-12;
    _untangling_ref_eps = 1e-1;

    if (start_untangle_eps > 0) {
        _untangling_ref_eps = start_untangle_eps;
    }
    _untangling_eps = _untangling_ref_eps;
    update_untangling_eps(1.);
    bool res = _det_min > 0;
    if (verbose) std::cout << "Initial energy: " << untangling_energy(_coords) << " | eps: " << _untangling_eps << " detmin: " << _det_min << " ninv: " << _nb_inverted << std::endl;
    if (verbose && has_bnd_terms()) std::cout << "Boundary energy: " << boundary_energy(_coords)  << std::endl;
    if (verbose && has_curves_and_points_terms()) std::cout << "Curves and point energy: " << curves_and_points_energies(_coords)  << std::endl;




    if (run_callback(UNTANGLING, 0)) {
        if (verbose) std::cout << "Early callback stop. Returning false." << std::endl;
        return false;
    }

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> linear_solver;
    if (laplacian_precond) {
        if (verbose) std::cout << "Computing diffusion pre-conditioning"<< std::endl;
        double diffusion_weight = 10;
        unsigned diffusion_power = 2;
        if (verbose) std::cout << "    diffusion weight = " << diffusion_weight << std::endl;
        if (verbose) std::cout << "    diffusion power = " << diffusion_power << std::endl;
        Time_log linear_timing_logger("Diffusion_inversion");
        Eigen::SparseMatrix<double> H = compute_diffusion_matrix(diffusion_weight);
        if (diffusion_power != 1) {
            Eigen::SparseMatrix<double> res = H * H;
            for (unsigned i = 0; i < diffusion_power - 2; ++i) {
            res = res * H;
            }
            H = res;
        }
        linear_solver.setTolerance(1e-8);
        linear_solver.compute(H);
        if (verbose) logging.log_sub_step("decomposition done" );
    }


    auto precond = [&](Eigen::VectorXd const &x, Eigen::VectorXd const &d) {
        if (!laplacian_precond) return d;
        if (verbose) logging.restart();
        Eigen::VectorXd md = linear_solver.solve(d);
        for (unsigned i = 0; i < static_cast<unsigned>(x.size()); ++i) {
            if (_locks[i]) md[i] = 0;
        }
        double scale_back = 1.;
        for (unsigned v = 0; v < nb_vertices(); ++v) {
            Eigen::Vector3d loc_d = Math_functions::sub_line_vector(md, v);
            double move = loc_d.norm();
            if (move > 10*_local_size[v]) {
                scale_back = (std::min)(scale_back, _local_size[v]/move);
            }
        }
        if (scale_back < 0.1) md *= 10.*scale_back;
        if (verbose) logging.log_sub_step("Inverse applied");
        return md;
    };

    smoother_status->return_code = CGAL::Mesh_smoothing_3::Smoothing_return_code::MAX_ITERATIONS_REACHED; // default if no other stopping criteria is met
    smoother_status->add_time(true);

    bool prev_res = false;
    for (unsigned iter = 0; iter < max_number_iter; ++iter) {
        if (verbose) std::cout << "Optimization iteration #" << iter << "\n";
        if (verbose) std::cout << "    curr eps: " << _untangling_eps << std::endl;
        ++number_of_outer_iter;
        ++smoother_status->nb_iterations;

        double e_prev = untangling_energy(_coords);
        double b_prev = boundary_energy(_coords);
        double e_curve_points_prev = curves_and_points_energies(_coords);
        if (!fine_time_logging && verbose) std::cout << "    optimizing"<<std::flush;
        Time_log fineLogging("optimization iter " + std::to_string(iter));
        Function_minimizer opt ([&](Eigen::VectorXd const &x, Eigen::VectorXd &g) -> double {
            g.setZero();
            double f = untangling_energy(x, &g);
            f += boundary_energy(x, &g);
            f += curves_and_points_energies(x, &g);
            return f;
        });
        opt._call_back = [&](Eigen::VectorXd const & x, Eigen::VectorXd const & g,double f,double step,unsigned bfgs_iter,unsigned nbEval) {
            bool stop_required = run_callback(UNTANGLING, iter, {bfgs_iter, step, nbEval}, &g);

            smoother_status->add_time();
            ++smoother_status->nb_vertex_updates;
            smoother_status->nb_metric_evaluations += nbEval;

            if (time_limit > 0 && smoother_status->total_time > time_limit) stop_required = true;
            if (max_nb_metric_evaluations > 0 && smoother_status->nb_metric_evaluations > static_cast<unsigned>(max_nb_metric_evaluations)) stop_required = true;

            bool significant_step = false;
            if (exact_predicate_optimization_check) {
                unsigned curr_invalid = evaluate_exact_predicates();
                if (curr_invalid != 0) {
                    ++_predicates_nb_invalid_steps;
                    significant_step = true;
                }
            }

            update_boundary_info(x, false);
            update_curves_and_points_info(x, false);
            ++number_of_lbfgs_iter;
            if (fine_time_logging)  fineLogging.log_sub_step("LBFGS iter " + std::to_string(bfgs_iter), "Energy: " + std::to_string(f));
            else if (verbose) std::cout << (significant_step ? "!" : ".") << std::flush;
            if (verbose && stop_required) std::cout << "--- CALLBACK REQUIRED STOPPING ---" << std::flush;
            return stop_required;
        };
        opt._precond = precond;

        opt.set_max_iter(max_lbfgs_iter);
        opt._parameters.min_step = 1e-40;
        opt._parameters.past = 8;
        opt._parameters.mem_size = 8;
        bool opt_res = opt.lbfgs_optimize(_coords);
        if (!fine_time_logging && verbose) std::cout << " done." << std::endl;
        if (fine_time_logging) fineLogging.log_total_time();
        double e = untangling_energy(_coords);
        double b = boundary_energy(_coords);
        double e_curve_points = curves_and_points_energies(_coords);
        compute_determinants();
        if (verbose) std::cout << "    E: " << e_prev << " -> " << e  << " | "  << "fmax: " << _untangling_max_energy << ", detmin: " << _det_min << " (eps=" << _untangling_eps << ") ninv: " << _nb_inverted << "\n";
        if (verbose && has_bnd_terms()) std::cout << "    B: " << b_prev << " -> " << b << std::endl;;
        if (verbose && has_curves_and_points_terms()) std::cout << "    Curves&Points: " << e_curve_points_prev << " -> " << e_curve_points << std::endl;;
        if (verbose) std::cout << "    Status: " << opt.get_message() << std::endl;

        if (time_limit > 0 && smoother_status->total_time > time_limit) {
            if (verbose) std::cout << "Time limit reached, stopping." << " ( " << smoother_status->total_time << "s / " << time_limit << "s )" << std::endl;
            smoother_status->return_code = CGAL::Mesh_smoothing_3::Smoothing_return_code::TIME_LIMIT_REACHED;
            break;
        }
        if (max_nb_metric_evaluations > 0 && smoother_status->nb_metric_evaluations > static_cast<unsigned>(max_nb_metric_evaluations)) {
            if (verbose) std::cout << "Max number of metric evaluations reached, stopping." << " ( " << smoother_status->nb_metric_evaluations << " / " << max_nb_metric_evaluations << " )" << std::endl;
            smoother_status->return_code = CGAL::Mesh_smoothing_3::Smoothing_return_code::MAX_NUMBER_OF_METRIC_EVALUATIONS_REACHED;
            break;
        }

        double improvement_ratio = e/e_prev;
        if (iter == 0) {
            // reseting the eps after the first iter to account for degenerated configurations
            improvement_ratio = 1;
            _untangling_eps = _untangling_ref_eps;
        }
        update_untangling_eps(improvement_ratio);
        res = _det_min > 0;
        if (verbose) logging.log_sub_step("optimization iter " + std::to_string(iter+1) );
        update_local_size();
        if (has_bnd_terms()) {
            if (verbose) std::cout << "    re assessing boundary to improve matching" << std::endl;
            update_boundary_info(_coords, true);
        }
        if (has_curves_and_points_terms()) {
            if (verbose) std::cout << "    updating curves and target point data" << std::endl;
            update_curves_and_points_info(_coords, true);
        }

        if (exact_predicate_status) {
            unsigned curr_invalid = evaluate_exact_predicates();
            if (verbose) std::cout << "Exact predicate check: " << curr_invalid << " invalid tetrahedra" << std::endl;
        }

        if (run_callback(UNTANGLING, iter+1)) {
            smoother_status->return_code = CGAL::Mesh_smoothing_3::Smoothing_return_code::USER_ABORT;
            if (verbose) std::cout << "Callback required stop, breaking." << std::endl;
            break;
        }

        if (prev_res && res && opt_res) {
            smoother_status->return_code = CGAL::Mesh_smoothing_3::Smoothing_return_code::CONVERGENCE_REACHED;
            if (verbose) std::cout << "Optimization converged and no tangled elements." << std::endl;
            break;
        };
        prev_res = res;
    }
    if (verbose) logging.log_total_time();

    smoother_status->add_time();
    return res;
}


// TODO need update for curve and point targets
template<typename Surface_patch_index, typename Curve_index, typename ConcurrencyTag>
inline bool Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index, ConcurrencyTag>::run_quality_maximization(unsigned max_number_iter) {
    if (verbose) std::cout << "==== Tetrahedral_mesh_smoother quality maximization ====" << "\n";
    if (_coords.size() == 0 || _tet_storage.empty()) {
        if (verbose) std::cout << "No variables to optimize." << std::endl;
        return true;
    }
    number_of_outer_iter = 0;
    number_of_lbfgs_iter = 0;

    Time_log logging("Tetrahedral_mesh_smoother");
    if (_det_min <= 0) {
        {
            CGAL::IO::Color_stream_guard bright_red(std::cout, CGAL::IO::Ansi_color::BrightRed);
            std::cout << "Inverted elements detected, Tetrahedral_mesh_smoother will first try to untangle them." << std::endl;
        }
        run_untangling(1500);
        if (_det_min <= 0) {
            CGAL::IO::Color_stream_guard bright_red(std::cout, CGAL::IO::Ansi_color::BrightRed);
            std::cout << "Input mesh still contains inverted elements, Tetrahedral_mesh_smoother cannot run its quality maximization routine." << std::endl;
            return false;
        }
        {
            CGAL::IO::Color_stream_guard bright_red(std::cout, CGAL::IO::Ansi_color::BrightBlue);
            std::cout << "Untangling was success. Now maximizing worst quality." << std::endl;
        }
    }
    else {
        CGAL::IO::Color_stream_guard bright_red(std::cout, CGAL::IO::Ansi_color::BrightRed);
        std::cout << "Running standard elliptic energy before quality maximization." << std::endl;
        run_untangling(2);
    }
    update_local_size();
    update_boundary_info(_coords, true);
    if (verbose) std::cout << "Nb of optimization variables (vertices x3): " << _coords.size()  << std::endl;
    if (verbose) std::cout << "Nb of energy terms (tetrahedra): " << _tet_storage.size() << std::endl;

    _qis_tau = 1./(1.01*_conformal_energy_max);
    if (verbose) std::cout << "Initial energy: " << qis_energy(_coords) << " | worst quality: " << _conformal_energy_max << " detmin: " << _det_min << std::endl;
    if (verbose && has_bnd_terms()) std::cout << "Boundary energy: " << boundary_energy(_coords)  << std::endl;

    bool res = false;

    if (run_callback(STIFFENING, 0)) {
        if (verbose) std::cout << "Early callback stop. Returning false." << std::endl;
        return false;
    }
    double prev_max_f = _conformal_energy_max;
    double init_max_f = _conformal_energy_max;

    for (unsigned iter = 0; iter < max_number_iter; ++iter) {
        if (verbose) std::cout << "Optimization iteration #" << iter << "\n";
        if (verbose) std::cout << "    curr quality bound: " << 1./_qis_tau << std::endl;;
        ++number_of_outer_iter;

        double e_prev = qis_energy(_coords); // both with same eps
        double b_prev = boundary_energy(_coords);
        if (!fine_time_logging && verbose) std::cout << "    optimizing"<<std::flush;
        Time_log fineLogging("optimization iter " + std::to_string(iter));
        Function_minimizer opt ([&](Eigen::VectorXd const &x, Eigen::VectorXd &g){
            g.setZero();
            return qis_energy(x, &g) + boundary_energy(x, &g);
        });
        opt._call_back = [&](Eigen::VectorXd const & x, Eigen::VectorXd const & g,double f,double step,unsigned bfgs_iter,unsigned nbEval) {
            bool stop_required = run_callback(STIFFENING, iter, {bfgs_iter, step, nbEval}, &g);
            update_boundary_info(x, false);
            ++number_of_lbfgs_iter;
            if (fine_time_logging)  fineLogging.log_sub_step("LBFGS iter " + std::to_string(bfgs_iter), "Energy: " + std::to_string(f));
            else if (verbose) std::cout << "." << std::flush;
            if (verbose && stop_required) std::cout << "--- CALLBACK REQUIRED STOPPING ---" << std::flush;
            return stop_required;
        };

        opt.set_max_iter(max_lbfgs_iter);
        opt._parameters.min_step = 1e-40;
        opt._parameters.past = 8;
        opt._parameters.mem_size = 8;
        bool opt_res = opt.lbfgs_optimize(_coords);
        if (!fine_time_logging && verbose) std::cout << " done." << std::endl;
        if (fine_time_logging) fineLogging.log_total_time();
        double e = qis_energy(_coords);
        double b = boundary_energy(_coords);
        compute_determinants();
        if (verbose) std::cout << "    E: " << e_prev << " -> " << e  << " | "  << "fmax: " << prev_max_f << " -> " << _conformal_energy_max << " | quality bound: " <<  1./_qis_tau << ", detmin: " << _det_min << "\n";
        if (verbose && has_bnd_terms()) std::cout << "    B: " << b_prev << " -> " << b << std::endl;;
        if (verbose) std::cout << "    Status: " << opt.get_message() << std::endl;

        update_qis_tau(e/e_prev);
        if (verbose) logging.log_sub_step("optimization iter " + std::to_string(iter+1) );
        update_local_size();
        if (has_bnd_terms()) {
            if (verbose) std::cout << "    re assessing boundary to improve matching" << std::endl;
            update_boundary_info(_coords, true);
        }
        if (run_callback(STIFFENING, iter+1)) {
            std::cout << "Callback required stop, breaking." << std::endl;
            break;
        }
        res = opt_res && ((prev_max_f-_conformal_energy_max)/prev_max_f < _QIS_CONVERGENCE_THRESHOLD);
        if (res) {
            break;
        };
        prev_max_f = _conformal_energy_max;
    }
    if (res) {
        if (verbose) std::cout << "Optimization converged." << std::endl;
        if (verbose) std::cout << "Static threshold used: " << _QIS_CONVERGENCE_THRESHOLD << std::endl;
    }
    else {
        if (verbose) std::cout << "Optimization failed to fully converge in " << max_number_iter << " iterations. " << std::endl;
    }
    if (verbose) std::cout << "Worst quality improvement: " << init_max_f << " -> " << _conformal_energy_max << std::endl;

    if (verbose) logging.log_total_time();
    return res;
}

} } // end of CGAL::Mesh_smoothing_3_internal namespace

#endif // CGAL_MESH_SMOOTHING_3_INTERNAL_TETRAHEDRAL_MESH_SMOOTHER_H
