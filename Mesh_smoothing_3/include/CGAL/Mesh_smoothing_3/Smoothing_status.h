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

#ifndef CGAL_MESH_SMOOTHING_3_SMOOTHING_STATUS_H
#define CGAL_MESH_SMOOTHING_3_SMOOTHING_STATUS_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <chrono>

namespace CGAL {

namespace Mesh_smoothing_3 {

// Re-defining the Parallel_if_available_tag specifically for Mesh_smoothing_3
// as it can automatically link with OpenMP if available, which is not the case for the rest of CGAL.
#if defined(CGAL_LINKED_WITH_TBB) || (defined(_OPENMP) && _OPENMP >= 201307)
using Parallel_if_available_tag = CGAL::Parallel_tag;
#else
using Parallel_if_available_tag = CGAL::Sequential_tag;
#endif

/*!
 * \ingroup pkgMeshSmoothing3Functions
 *
 * \brief Return code of the mesh smoothing.
 */
enum class Smoothing_return_code
{
    ALL_VERTICES_FROZEN, ///< All vertices are frozen.
    CONVERGENCE_REACHED, ///< Convergence reached.
    MAX_ITERATIONS_REACHED, ///< Maximum number of iterations reached.
    MAX_NUMBER_OF_METRIC_EVALUATIONS_REACHED, ///< Maximum number of metric evaluations reached.
    TIME_LIMIT_REACHED, ///< Time limit reached.
    USER_ABORT, ///< Smoothing was aborted by the user.
};


/*!
 * \ingroup pkgMeshSmoothing3Functions
 *
 * \brief Status of the mesh smoothing process.
 */
struct Smoothing_status {
    Smoothing_return_code return_code = Smoothing_return_code::ALL_VERTICES_FROZEN; ///< Return code of the smoothing algorithm.
    unsigned nb_iterations = 0; ///< Number of smoothing/untangling iterations performed.
    unsigned nb_vertex_updates = 0; ///< Number of times vertices were updated.
    unsigned nb_metric_evaluations = 0; ///< Number of times the quality metric was evaluated for optimization.
    unsigned nb_initial_invalid_elements = 0; ///< Number of negatively oriented elements in the mesh at the beginning of the smoothing process.
    unsigned nb_invalid_elements = 0; ///< Number of negatively oriented elements in the mesh at the end of the smoothing process.
    double total_time = 0.; ///< Total time spent in the smoothing algorithm, in seconds.
    double pre_processing_time = 0.; ///< Time spent in the preprocessing step, in seconds.
    double optimization_time = 0.; ///< Time spent in the optimization step, in seconds.

    bool valid_mesh() const { return nb_invalid_elements == 0; } ///< returns true if the mesh is valid (no negatively oriented elements).


// internal usage

    void add_time(bool pre_processing = false) {
        _running = true;
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(now - previous).count()) * 1e-6;
        if (pre_processing) pre_processing_time += time;
        else optimization_time += time;
        total_time += time;
        previous = now;
    }

    Smoothing_status() {
        previous = std::chrono::steady_clock::now();
    }

    bool in_progress() const { return _running; }
    void conclude() { _running = false; }

private:
    std::chrono::steady_clock::time_point previous;
    bool _running = false;
};

} } // end of CGAL::Mesh_smoothing_3 namespace

#endif // CGAL_MESH_SMOOTHING_3_SMOOTHING_STATUS_H
