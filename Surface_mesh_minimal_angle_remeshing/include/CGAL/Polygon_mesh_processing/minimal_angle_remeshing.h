// Copyright (c) 2019  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Kaimo Hu

#ifndef CGAL_POLYGON_MESH_PROCESSING_MINIMAL_ANGLE_REMESHING_H
#define CGAL_POLYGON_MESH_PROCESSING_MINIMAL_ANGLE_REMESHING_H

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

// TODO(kaimo Hu): replace line 20 with line 19 when published
// #include <CGAL/Polygon_mesh_processing/internal/minangle_remesh_impl.h>
#include "internal/minangle_remesh_impl.h"
// TODO(kaimo Hu): unccomment the following to use the parameters interface
// #include "../../../BGL/include/CGAL/boost/graph/parameters_interface.h"

namespace CGAL{
namespace Polygon_mesh_processing {

/*!
* \ingroup PMP_meshing_grp
* @brief remeshes a triangule mesh using a minimal angle optimization approach.
* This operation applies carefully prioritized local operations to greedily 
* search for the coarsest mesh with minimal interior angle threshold and the 
* approximation error threshold defined by the user. The local operations 
* include edge split, edge collapse, edge flip (optional), vertex relocation, 
* and projection to the initial triangular mesh.
*
* @tparam TriangleMesh model of `MutableFaceGraph`.
*         The descriptor types `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         and `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` must be
*         models of `Hashable`.
*         If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
*         and no `face_index_map` is given
*         as a named parameter, then the internal one must be initialized
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param tm the triangle mesh to be remeshed
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" 
* among the ones listed below
*
*
* \cgalNamedParamsBegin
* 
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, 
*                       using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the
*                     vertex point type.}
*     \cgalParamExtra{Exact constructions kernels are not supported by this function.}
*   \cgalParamNEnd
* 
*   \cgalParamNBegin{max_error_threshold}
*     \cgalParamDescription{The approximation error, expressed as the percentage of
*                           diagonal length of the input surface mesh.}
*     \cgalParamType{double}
*     \cgalParamDefault{0.2}
*   \cgalParamNEnd
*   \cgalParamNBegin{min_angle_threshold}
*     \cgalParamDescription{The minimal angle that the remesh surface mesh should achieve.}
*     \cgalParamType{double}
*     \cgalParamDefault{30}
*   \cgalParamNEnd
*   \cgalParamNBegin{max_mesh_complexity}
*     \cgalParamDescription{The maximal mesh complexity that the remesh surface mesh should 
*                           maintain, expressed as the number of vertices of the remesh 
*                           surface mesh.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{100000000}
*   \cgalParamNEnd
*   \cgalParamNBegin{smooth_angle_delta}
*     \cgalParamDescription{The minimal step for angle improvement. If a local operator 
*                           improves the minimal angle less than this value, then it is 
*                           considered the minimal angle has not been improved.}
*     \cgalParamType{double}
*     \cgalParamDefault{0.1}
*   \cgalParamNEnd
* 
*   \cgalParamNBegin{apply_edge_flip}
*     \cgalParamDescription{Indicates whether we apply the local operator `Edge_flip` when
*                           improving the minimal angle. Experiments show that enable the edge
*                           flipping achieves much better results}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*   \cgalParamNBegin{edge_flip_strategy}
*     \cgalParamDescription{Indicates whether we want to improve the minimal angle, or the vertex 
*                           valences when flipping an edge. Experiments show that 
*                           `Improve_the_minimal_angle` is better.}
*     \cgalParamType{EdgeFlipStrategy}
*     \cgalParamDefault{EdgeFlipStrategy::k_improve_angle}
*   \cgalParamNEnd
*   \cgalParamNBegin{flip_after_split_and_collapse}
*     \cgalParamDescription{Indicates whether we perform `Edge_flip` after applying edge split or
*                           edge collapse, such that the local valance of the new generated vertex
*                           can be improved. Experiments show that enabling this is better.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*   \cgalParamNBegin{relocate_after_local_operations}
*     \cgalParamDescription{Indicates whether we perform `Vertex_relocate` after applying other
*                           local operators, such as edge split, edge collapse. Note we do not
*                           perform `Vertex_relocate` after edge flip in order to perverse sharp
*                           features. Experiments show that enabling this is better than not.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
* 
*   \cgalParamNBegin{relocate_strategy}
*     \cgalParamDescription{The options for `Vertex_relocation`, we now have barycenter or CVT 
*                           barycenter. Barycenter is simply the average of the one-ring vertices,
*                           which is also called Laplacian-operation. Instead, the CVT barycenter 
*                           is the Centroid Voronoi Tessellation center of the one-ring vertices. 
*                           Experiments show that the CVT barycenter is better.}
*     \cgalParamType{RelocateStrategy}
*     \cgalParamDefault{RelocateStrategy::k_cvt_barycenter}
*   \cgalParamNEnd
*   \cgalParamNBegin{keep_vertex_in_one_ring}
*     \cgalParamDescription{Indicate whether we want to keep the vertex in one-ring region when
*                           performing `Edge_collapse` or `Vertex_relocate`. This was indicated
*                           in Fig. 4 of our paper. Experiments show that disabling this is much
*                           better. However, we may have some folder-overs in some results.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*   \cgalParamNBegin{use_local_aabb_tree}
*     \cgalParamDescription{Indicates whether we construct the local AABB tree or not when 
*                           simulating the edge collapse operator. If it is true, we use the 
*                           constructed local AABB tree to update the links from input surface 
*                           mesh to remesh surface mesh; Otherwise, we update the links directly. 
*                           Experiments show that the efficiencies for both cases are almost the 
*                           same.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*   \cgalParamNBegin{collapsed_list_size}
*     \cgalParamDescription{In some rare cases, the software traps into local dead loops (e.g., 
*                           edge split -> edge collapse -> edge split ->edge collapse...). Hence, 
*                           we maintain a global collapse list, such that if an edge is collapsed
*                           not long ago (recorded in the list), we deny the edge collapse when it
*                           can be in later steps. This parameter indicates the maximum size of
*                           the collapse list. The larger this value is, the less possibility that 
*                           the software traps into local dead loops, but higher computational and 
*                           memory cost. Experiments show that this strategy is really effective.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{10}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{decrease_max_errors}
*     \cgalParamDescription{Due to the randomness of sampling, sometimes the max error threshold 
*                           cannot be strictly bounded to the threshold. If this parameter is 
*                           enabled, we explicitly decrease the max error as long as we find the 
*                           Hausdorff distance between the input and the remesh exceeds the 
*                           threshold. This will make sure the approximation error is strictly 
*                           bounded regardless of the sampling randomness.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*   \cgalParamNBegin{verbose_progress}
*     \cgalParamDescription{If this parameter is enabled, we output the detailed information when 
*                           each local operator is applied. The detailed information may include 
*                           current maximal error, minimal angle, or the size of the dynamic 
*                           priority queue.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*   \cgalParamNBegin{apply_initial_mesh_simplification}
*     \cgalParamDescription{If this parameter is enabled, we apply the initial mesh simplification
*                           before improving the minimal angles. This reduces the mesh complexity 
*                           around 20-30% on average.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*   \cgalParamNBegin{apply_final_vertex_relocation}
*     \cgalParamDescription{If this parameter is enabled, we apply the final vertex relocation 
*                           after improving the minimal angles to the min angle threshold. This 
*                           will improve the average quality of the triangles with respect to the 
*                           parameter `Smooth angle delta`. The smaller `Smooth angle delta` is 
*                           set, the better results we get along with the higher computational 
*                           cost.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{samples_per_face_in}
*     \cgalParamDescription{The number of samples per facet on input surface mesh.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{10}
*   \cgalParamNEnd
*   \cgalParamNBegin{samples_per_face_out}
*     \cgalParamDescription{The number of samples per facet on remesh surface mesh.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{10}
*   \cgalParamNEnd
*   \cgalParamNBegin{max_samples_per_area}
*     \cgalParamDescription{The maximal number of samplers per unit area. This parameter is 
*                           used to avoid too dense sampling.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{10000}
*   \cgalParamNEnd
*   \cgalParamNBegin{min_samples_per_triangle}
*     \cgalParamDescription{The minimal number of samples per triangle. This parameter is used 
*                           to guarantee that each triangle, no matter how small it is, has 
*                           certain number of samples on it.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{1}
*   \cgalParamNEnd
* 
*   \cgalParamNBegin{bvd_iteration_count}
*     \cgalParamDescription{The large this parameter is, the more uniform the samples in the 
*                           facets are. However, the computational cost will be dramatically 
*                           higher as well.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{1}
*   \cgalParamNEnd
*   \cgalParamNBegin{sample_number_strategy}
*     \cgalParamDescription{The option of sample number strategy on remesh surface mesh. If it 
*                           is `Fixed`, then the number of samples per facet is roughly the same 
*                           as the one on the input surface mesh; If it is `Variable`, then the 
*                           number of samples per facet is variable with respect to the size of 
*                           facets of the remesh surface mesh: the more facets the remesh surface 
*                           mesh has, the smaller number of samples we generate on each facet of 
*                           the remesh surface mesh. The `Variable` option makes the total 
*                           samples on the input surface mesh and remesh surface mesh roughly the 
*                           same.}
*     \cgalParamType{SampleNumberStrategy}
*     \cgalParamDefault{SampleNumberStrategy::k_fixed}
*   \cgalParamNEnd
*   \cgalParamNBegin{sample_strategy}
*     \cgalParamDescription{The option of sample strategy. If it is `uniform`, then the number 
*                           of samples per facet is proportional to its area; If it is `adaptive`, 
*                           then the number of samples per facet is roughly the same.}
*     \cgalParamType{SampleStrategy}
*     \cgalParamDefault{SampleStrategy::k_adaptive}
*   \cgalParamNEnd
*   \cgalParamNBegin{use_stratified_sampling}
*     \cgalParamDescription{If this parameter is enabled, the vertex samples, edge samples and 
*                           the facet samples partition the area of the surface mesh, 
*                           respectively (we call the partition area as the capacity of each 
*                           sample); If it is disabled, all the three types of samples will 
*                           partition the area of the surface mesh together. In the former case, 
*                           the feature samples (include the vertex samples and the edge samples) 
*                           own higher area weights; while in the latter case, all the samples 
*                           have the same area weights. Enabling this parameter preserves 
*                           features better, but sacrifices the overall approximation error, 
*                           because the facet samples dominates all the samples.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{sum_theta}
*     \cgalParamDescription{The maximal value of Gaussian curvature, expressed as the times of PI. 
*                           If a Gaussian curvature exceeds some value, it will be clamped to this 
*                           value.}
*     \cgalParamType{double}
*     \cgalParamDefault{1.0}
*   \cgalParamNEnd
*   \cgalParamNBegin{sum_delta}
*     \cgalParamDescription{The scale of the Gaussian curvature. If dividing a Gaussian curvature 
*                           with this scale and the result exceeds Sum theta (PI), then it is 
*                           clamped to Sum theta (PI).}
*     \cgalParamType{double}
*     \cgalParamDefault{0.5}
*   \cgalParamNEnd
*   \cgalParamNBegin{dihedral_theta}
*     \cgalParamDescription{The maximal value of large dihedral angle value, expressed as the 
*                           times of PI. If the large dihedral angle exceeds some value, it will 
*                           be clamped to this value.}
*     \cgalParamType{double}
*     \cgalParamDefault{1.0}
*   \cgalParamNEnd
*   \cgalParamNBegin{dihedral_delta}
*     \cgalParamDescription{The scale of the large dihedral angle. If dividing the large dihedral 
*                           angle with this scale and the result exceeds Dihedral theta (PI), 
*                           then it is clamped to Dihedral theta (PI).}
*     \cgalParamType{double}
*     \cgalParamDefault{0.5}
*   \cgalParamNEnd
* 
*   \cgalParamNBegin{feature_difference_delta}
*     \cgalParamDescription{This parameter is only used when getting the initial position of the 
*                           vertex after applying edge collapse. If the feature intensity 
*                           difference between the two end points is smaller than their maximal 
*                           value multiplied with this parameter, then the midpoint will be 
*                           selected as the initial position; otherwise, the end point with 
*                           higher feature intensity will be selected.}
*     \cgalParamType{double}
*     \cgalParamDefault{0.15}
*   \cgalParamNEnd
*   \cgalParamNBegin{feature_control_delta}
*     \cgalParamDescription{This parameter is used for vertex classification. Please refer to 
*                           Fig. 10 in the paper for more detailed explanation.}
*     \cgalParamType{double}
*     \cgalParamDefault{0.5}
*   \cgalParamNEnd
*   \cgalParamNBegin{inherit_element_types}
*     \cgalParamDescription{If this parameter is enabled, we calculate the edge types (crease 
*                           edge or non-crease edge) of the remesh surface mesh in advance, 
*                           and then maintain the edge types explicitly during the whole 
*                           remeshing process. Otherwise, we calculate the edge types according 
*                           to the given parameters each time when a local operator is applied.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*   \cgalParamNBegin{use_feature_intensity_weights}
*     \cgalParamDescription{If this parameter is enabled, the weight for each sample is its 
*                           according area multiplied by its feature intensity; Otherwise, 
*                           the weight is its according area. Enabling this parameter keeps the 
*                           features better.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
* 
*   \cgalParamNBegin{vertex_optimize_count}
*     \cgalParamDescription{The number of iterations we perform when optimizing a vertex position. 
*                           Please refer to the first row of Fig. 15 in the paper for more 
*                           details.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{2}
*   \cgalParamNEnd
*   \cgalParamNBegin{vertex_optimize_ratio}
*     \cgalParamDescription{The ratio to get to the optimal position when optimizing a vertex 
*                           position. Refer to the second row of Fig. 15 in the paper for more 
*                           details.}
*     \cgalParamType{double}
*     \cgalParamDefault{0.9}
*   \cgalParamNEnd
*   \cgalParamNBegin{stencil_ring_size}
*     \cgalParamDescription{The stencil ring size when collecting the samples from input surface 
*                           mesh to remesh surface mesh. Please refer to the gray region of 
*                           Fig. 3 as well as the discussion in Sec. 6.1 in the paper for more 
*                           details.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{1}
*   \cgalParamNEnd
*   \cgalParamNBegin{optimize_strategy}
*     \cgalParamDescription{Options for vertex position optimization. If it is `Approximation`, 
*                           then the vertex position is just the optimal position we calculated; 
*                           otherwise, the vertex position is the projection of the optimal 
*                           position to the input surface mesh. Please refer to Sec. 6.2 in the 
*                           paper for more details.}
*     \cgalParamType{OptimizeStrategy}
*     \cgalParamDefault{OptimizeStrategy::k_approximation}
*   \cgalParamNEnd
* 
*   \cgalParamNBegin{face_optimize_type}
*     \cgalParamDescription{The facet sample types used for optimizing a vertex position.}
*     \cgalParamType{OptimizeType}
*     \cgalParamDefault{OptimizeType::k_both}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{edge_optimize_type}
*     \cgalParamDescription{The edge sample types used for optimizing a vertex position.}
*     \cgalParamType{OptimizeType}
*     \cgalParamDefault{OptimizeType::k_both}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_optimize_type}
*     \cgalParamDescription{The vertex sample types used for optimizing a vertex position.
*                           `Optimize after local operations`: If this parameter is enabled, 
*                           we perform the vertex position optimization after each local 
*                           operator is applied.}
*     \cgalParamType{OptimizeType}
*     \cgalParamDefault{OptimizeType::k_both}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{optimize_after_local_operations}
*     \cgalParamDescription{If this parameter is enabled, we perform the vertex position 
*                           optimization after each local operator is applied.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*/
template <class TriangleMesh, class NamedParameters>
void minimal_angle_remeshing(TriangleMesh& tm, const NamedParameters& np)
{
  // step 1: construct the internal remesher
  using parameters::get_parameter;
  using parameters::choose_parameter;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  GT gt = choose_parameter(get_parameter(np, internal_np::geom_traits), GT());
  typedef internal::Minangle_remesher<GT> Minangle_remesher;
  Minangle_remesher remesher;

  // step 2: setup the parameters
  // 2.1: general parameters
  //double max_error_threshold = choose_parameter(get_parameter(np, 
  //  internal_np::max_error_threshold), 0.2);
  //remesher.set_max_error_threshold(max_error_threshold);
  remesher.set_max_error_threshold(0.2);
  //double min_angle_threshold = choose_parameter(get_parameter(np,
  //  internal_np::min_angle_threshold), 30.0);
  //remesher.set_min_angle_threshold(min_angle_threshold);
  remesher.set_min_angle_threshold(30);
  //int max_mesh_complexity = choose_parameter(get_parameter(np,
  //  internal_np::max_mesh_complexity), 100000000);
  //remesher.set_max_mesh_complexity(max_mesh_complexity);
  remesher.set_max_mesh_complexity(100000000);
  //double smooth_angle_delta = choose_parameter(get_parameter(np,
  //  internal_np::smooth_angle_delta), 0.1);
  //remesher.set_smooth_angle_delta(smooth_angle_delta);
  remesher.set_smooth_angle_delta(0.1);

  //bool apply_edge_flip = choose_parameter(get_parameter(np,
  //  internal_np::apply_edge_flip), true);
  //remesher.set_apply_edge_flip(apply_edge_flip);
  remesher.set_apply_edge_flip(true);
  //EdgeFlipStrategy edge_flip_strategy = choose_parameter(get_parameter(np,
  //  internal_np::edge_flip_strategy), EdgeFlipStrategy::k_improve_angle);
  //remesher.set_edge_flip_strategy(edge_flip_strategy);
  remesher.set_edge_flip_strategy(EdgeFlipStrategy::k_improve_angle);
  //bool flip_after_split_and_collapse = choose_parameter(get_parameter(np,
  //  internal_np::flip_after_split_and_collapse), true);
  //remesher.set_flip_after_split_and_collapse(flip_after_split_and_collapse);
  remesher.set_flip_after_split_and_collapse(true);
  //bool relocate_after_local_operations = choose_parameter(get_parameter(np,
  //  internal_np::relocate_after_local_operations), true);
  //remesher.set_relocate_after_local_operations(relocate_after_local_operations);
  remesher.set_relocate_after_local_operations(true);

  //RelocateStrategy relocate_strategy = choose_parameter(get_parameter(np,
  //  internal_np::relocate_strategy), RelocateStrategy::k_cvt_barycenter);
  //remesher.set_relocate_strategy(relocate_strategy);
  remesher.set_relocate_strategy(RelocateStrategy::k_cvt_barycenter);
  //bool keep_vertex_in_one_ring = choose_parameter(get_parameter(np,
  //  internal_np::keep_vertex_in_one_ring), false);
  //remesher.set_keep_vertex_in_one_ring(keep_vertex_in_one_ring);
  remesher.set_keep_vertex_in_one_ring(false);
  //bool use_local_aabb_tree = choose_parameter(get_parameter(np,
  //  internal_np::use_local_aabb_tree), true);
  //remesher.set_use_local_aabb_tree(use_local_aabb_tree);
  remesher.set_use_local_aabb_tree(true);
  //int collapsed_list_size = choose_parameter(get_parameter(np,
  //  internal_np::collapsed_list_size), 10);
  //remesher.set_collapsed_list_size(collapsed_list_size);
  remesher.set_collapsed_list_size(10);

  //bool decrease_max_errors = choose_parameter(get_parameter(np,
  //  internal_np::decrease_max_errors), true);
  //remesher.set_decrease_max_errors(decrease_max_errors);
  remesher.set_decrease_max_errors(true);
  //bool verbose_progress = choose_parameter(get_parameter(np,
  //  internal_np::verbose_progress), true);
  //remesher.set_verbose_progress(verbose_progress);
  remesher.set_verbose_progress(true);
  //bool apply_initial_mesh_simplification = choose_parameter(get_parameter(np,
  //  internal_np::apply_initial_mesh_simplification), true);
  //remesher.set_apply_initial_mesh_simplification(apply_initial_mesh_simplification);
  remesher.set_apply_initial_mesh_simplification(true);
  //bool apply_final_vertex_relocation = choose_parameter(get_parameter(np,
  //  internal_np::apply_final_vertex_relocation), true);
  //remesher.set_apply_final_vertex_relocation(apply_final_vertex_relocation);
  remesher.set_apply_final_vertex_relocation(true);

  //// 2.2: sample parameters
  //int samples_per_face_in = choose_parameter(get_parameter(np,
  //  internal_np::samples_per_face_in), 10);
  //remesher.set_samples_per_face_in(samples_per_face_in);
  remesher.set_samples_per_face_in(10);
  //int samples_per_face_out = choose_parameter(get_parameter(np,
  //  internal_np::samples_per_face_out), 10);
  //remesher.set_samples_per_face_out(samples_per_face_out);
  remesher.set_samples_per_face_out(10);
  //int max_samples_per_area = choose_parameter(get_parameter(np,
  //  internal_np::max_samples_per_area), 10000);
  //remesher.set_max_samples_per_area(max_samples_per_area);
  remesher.set_max_samples_per_area(10000);
  //int min_samples_per_triangle = choose_parameter(get_parameter(np,
  //  internal_np::min_samples_per_triangle), 1);
  //remesher.set_min_samples_per_triangle(min_samples_per_triangle);
  remesher.set_min_samples_per_triangle(1);

  //int bvd_iteration_count = choose_parameter(get_parameter(np,
  //  internal_np::bvd_iteration_count), 1);
  //remesher.set_bvd_iteration_count(bvd_iteration_count);
  remesher.set_bvd_iteration_count(1);
  //SampleNumberStrategy sample_number_strategy = choose_parameter(get_parameter(np,
  //  internal_np::sample_number_strategy), SampleNumberStrategy::k_fixed);
  //remesher.set_sample_number_strategy(sample_number_strategy);
  remesher.set_sample_number_strategy(SampleNumberStrategy::k_fixed);
  //SampleStrategy sample_strategy = choose_parameter(get_parameter(np,
  //  internal_np::sample_strategy), SampleStrategy::k_adaptive);
  //remesher.set_sample_strategy(sample_strategy);
  remesher.set_sample_strategy(SampleStrategy::k_adaptive);
  //bool use_stratified_sampling = choose_parameter(get_parameter(np,
  //  internal_np::use_stratified_sampling), false);
  //remesher.set_use_stratified_sampling(use_stratified_sampling);
  remesher.set_use_stratified_sampling(false);

  //// 2.3: feature function parameters
  //double sum_theta = choose_parameter(get_parameter(np,
  //  internal_np::sum_theta), 1.0);
  //remesher.set_sum_theta(sum_theta);
  remesher.set_sum_theta(1.0);
  //double sum_delta = choose_parameter(get_parameter(np,
  //  internal_np::sum_delta), 0.5);
  //remesher.set_sum_delta(sum_delta);
  remesher.set_sum_delta(0.5);
  //double dihedral_theta = choose_parameter(get_parameter(np,
  //  internal_np::dihedral_theta), 1.0);
  //remesher.set_dihedral_theta(dihedral_theta);
  remesher.set_dihedral_delta(1.0);
  //double dihedral_delta = choose_parameter(get_parameter(np,
  //  internal_np::dihedral_delta), 0.5);
  //remesher.set_dihedral_delta(dihedral_delta);
  remesher.set_dihedral_delta(0.5);

  //double feature_difference_delta = choose_parameter(get_parameter(np,
  //  internal_np::feature_difference_delta), 0.15);
  //remesher.set_feature_difference_delta(feature_difference_delta);
  remesher.set_feature_difference_delta(0.15);
  //double feature_control_delta = choose_parameter(get_parameter(np,
  //  internal_np::feature_control_delta), 0.5);
  //remesher.set_feature_control_delta(feature_control_delta);
  remesher.set_feature_control_delta(0.5);
  //bool inherit_element_types = choose_parameter(get_parameter(np,
  //  internal_np::inherit_element_types), false);
  //remesher.set_inherit_element_types(inherit_element_types);
  remesher.set_inherit_element_types(false);
  //bool use_feature_intensity_weights = choose_parameter(get_parameter(np,
  //  internal_np::use_feature_intensity_weights), false);
  //remesher.set_use_feature_intensity_weights(use_feature_intensity_weights);
  remesher.set_use_feature_intensity_weights(false);

  //// 2.4 vertex relocate parameters
  //int vertex_optimize_count = choose_parameter(get_parameter(np,
  //  internal_np::vertex_optimize_count), 2);
  //remesher.set_vertex_optimize_count(vertex_optimize_count);
  remesher.set_vertex_optimize_count(2);
  //double vertex_optimize_ratio = choose_parameter(get_parameter(np,
  //  internal_np::vertex_optimize_ratio), 0.9);
  //remesher.set_vertex_optimize_ratio(vertex_optimize_ratio);
  remesher.set_vertex_optimize_ratio(0.9);
  //int stencil_ring_size = choose_parameter(get_parameter(np,
  //  internal_np::stencil_ring_size), 1);
  //remesher.set_stencil_ring_size(stencil_ring_size);
  remesher.set_stencil_ring_size(1);
  //OptimizeStrategy optimize_strategy = choose_parameter(get_parameter(np,
  //  internal_np::optimize_strategy), OptimizeStrategy::k_approximation);
  //remesher.set_optimize_strategy(optimize_strategy);
  remesher.set_optimize_strategy(OptimizeStrategy::k_approximation);

  //OptimizeType face_optimize_type = choose_parameter(get_parameter(np,
  //  internal_np::face_optimize_type), OptimizeType::k_both);
  //remesher.set_face_optimize_type(face_optimize_type);
  remesher.set_face_optimize_type(OptimizeType::k_both);
  //OptimizeType edge_optimize_type = choose_parameter(get_parameter(np,
  //  internal_np::edge_optimize_type), OptimizeType::k_both);
  //remesher.set_edge_optimize_type(edge_optimize_type);
  remesher.set_edge_optimize_type(OptimizeType::k_both);
  //OptimizeType vertex_optimize_type = choose_parameter(get_parameter(np,
  //  internal_np::vertex_optimize_type), OptimizeType::k_both);
  //remesher.set_vertex_optimize_type(vertex_optimize_type);
  remesher.set_vertex_optimize_type(OptimizeType::k_both);
  //bool optimize_after_local_operations = choose_parameter(get_parameter(np,
  //  internal_np::optimize_after_local_operations), true);
  //remesher.set_optimize_after_local_operations(optimize_after_local_operations);
  remesher.set_optimize_after_local_operations(true);

  // step 3: set the input/output, and perform remeshing
  TriangleMesh tm_original(tm);
  remesher.set_input(&tm_original, remesher.get_verbose_progress());
  remesher.set_remesh(&tm, remesher.get_verbose_progress());
  remesher.maximize_minimal_angle();
};

template <class TriangleMesh>
void minimal_angle_remeshing(TriangleMesh& tm)
{
  minimal_angle_remeshing(tm, parameters::all_default());
};

} } // end of CGAL::Polygon_mesh_processing

#endif
