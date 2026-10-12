// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_FEATURE_GRAPH_INTERNAL_DGTAL_INTERFACE_H
#define CGAL_FEATURE_GRAPH_INTERNAL_DGTAL_INTERFACE_H

#include <CGAL/license/Feature_graph.h>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/helpers/StdDefs.h>
// This file must be included before CGAL files that uses Eigen. ??
#include <DGtal/dec/ATSolver2D.h>
#include <DGtal/dec/DiscreteExteriorCalculusFactory.h>

#include <CGAL/Feature_graph/internal/convert_image_to_surface_mesh.h>

namespace CGAL
{
namespace Feature_graph
{
namespace internal
{

// Option 1 : Translate Mesh : face_descriptor -> graph<DigitalSurface>::vertex_descriptor, halfedge_descriptor -> graph<DigitalSurface>::edge_descriptor, ...
// template < class TDigitalSurfaceContainer >
// struct DGtal_mesh_face_graph : public boost::graph_traits< DGtal::DigitalSurface< TDigitalSurfaceContainer > >
// {
//     typedef boost::graph_traits< DGtal::DigitalSurface< TDigitalSurfaceContainer > > Base;
//     typedef typename Base::vertex_descriptor face_descriptor;
// }

// Option 2 : Translate Image and use DGtal's surface creation
// Image -> Point predicate with `bool operator()( PointVector< dim, Integer > )`

// Option 3 : Translate CGAL mesh and use AT on it

// Option 4 : Construct DigitalSurface from CGAL Mesh

// Option 5 : Recode AT

template <typename Image>
struct LabelledImageToBinaryWrapper
{
  typedef std::array<int, 3> Point;
  const Image* image;
  const float* label;
  LabelledImageToBinaryWrapper(const Image& image, const float& label)
    : image(&image), label(&label)
  { }
  template <typename IntegerVector_3>
  bool operator()(const IntegerVector_3& coords) const
  {
    return image->value(coords[0], coords[1], coords[2]) == (*label);
  }
  // LabelledImageToBinaryWrapper& operator=(const LabelledImageToBinaryWrapper& other)
  // {
  //   image = other.image;
  //   label = other.label;
  //   return *this;
  // }
};

template <typename KSpace, typename Image>
void make_DGtal_surface_from_image(const KSpace& kSpace, const Image& image,
  std::vector<DGtal::CountedPtr<typename DGtal::Shortcuts<KSpace>::DigitalSurface>>& surfaces,
  std::vector<float>& surfaces_domain)
{
  typedef DGtal::Shortcuts<KSpace> SH;
  typedef DGtal::ShortcutsGeometry<KSpace> SHG;
  DGtal::Parameters params = SH::defaultParameters() | SHG::defaultParameters();
  // find all image labels
  std::set<float> labels;
  for (std::size_t z = 0; z < image.zdim(); z++)
  for (std::size_t y = 0; y < image.ydim(); y++)
  for (std::size_t x = 0; x < image.xdim(); x++)
  {
    labels.insert(image.value(x,y,z));
  }
  // create a surface per label
  for (const float& label : labels)
  {
    DGtal::CountedPtr<LabelledImageToBinaryWrapper<Image>> binary_image_ptr(new LabelledImageToBinaryWrapper<Image>(image, label));
    surfaces.push_back(SH::makeDigitalSurface(binary_image_ptr, kSpace, params));
    surfaces_domain.push_back(label);
  }
}

template <typename KSpace>
struct Similar_normal_predicate
{
  typedef DGtal::Shortcuts<KSpace> SH;
  typedef typename SH::RealVector RealVector;
  typedef typename SH::RealVectors RealVectors;

  typedef typename SH::SCell Vertex;

  RealVector base_normal;
  RealVectors normals;
  const std::map<Vertex, std::size_t>* surfel_scell_index;
  double cos_angle_threshold;

  Similar_normal_predicate(const double &cos_angle_threshold,
    const RealVector &base_normal,
    const RealVectors &normals,
    const std::map<Vertex, std::size_t> &surfel_scell_index)
    : cos_angle_threshold(cos_angle_threshold), base_normal(base_normal), normals(normals), surfel_scell_index(&surfel_scell_index)
  {}

  bool operator()(const Vertex &surfel) const
  {
    const std::size_t &surfel_id = surfel_scell_index->at(surfel);
    const RealVector &surfel_normal = normals[surfel_id];
    const double cos_normal_angle = base_normal.dot(surfel_normal);
    return cos_normal_angle > cos_angle_threshold;
  }
};

// regularizes the normal estimation.
//  For each surfel, a patch is created.
//  it constains surfels that have normal close to the initial surfel, constructed by breadth first search.
//  If the patch has a few number of surfel, then invert the normal of the initial surfel.
template <typename KSpace>
void regularize_inversion_on_normals(
  const typename DGtal::Shortcuts<KSpace>::DigitalSurface &surface,
  const typename DGtal::Shortcuts<KSpace>::SurfelRange &surfels,
  typename DGtal::Shortcuts<KSpace>::RealVectors &normals,
  const double angle_threshold = boost::math::constants::pi<double>()*.25,
  const std::size_t number_of_surfel_threshold = 8)
{
  typedef DGtal::Shortcuts<KSpace> SH;
  typedef typename SH::SCell SCell;
  typedef typename SH::RealVector RealVector;
  typedef typename SH::RealVectors RealVectors;
  const double &cos_angle_threshold = std::cos(angle_threshold);

  std::map<SCell, std::size_t> surfel_scell_index;
  for (std::size_t surfel_id = 0, nb_surfel = surfels.size(); surfel_id < nb_surfel; surfel_id++)
  {
    surfel_scell_index.emplace(surfels[surfel_id], surfel_id);
  }
  for (std::size_t surfel_id = 0, nb_surfel = surfels.size(); surfel_id < nb_surfel; surfel_id++)
  {
    const SCell& surfel_scell = surfels[surfel_id];
    const RealVector &surfel_normal = normals[surfel_id];
    Similar_normal_predicate<KSpace> normal_predicate(cos_angle_threshold, surfel_normal, normals, surfel_scell_index);
    DGtal::BreadthFirstVisitor<typename SH::DigitalSurface> visitor( surface, surfel_scell );
    std::size_t visited_surfel_number = 0;
    while (!visitor.finished() && visited_surfel_number < number_of_surfel_threshold)
    {
      visited_surfel_number++;
      visitor.expand(normal_predicate);
    }
    if (visited_surfel_number < number_of_surfel_threshold)
    {
      for (const SCell &surfel_inverted : visitor.markedVertices())
      {
        normals[surfel_scell_index.at(surfel_inverted)] *= -1.0;
      }
    }
  }
}

template <typename Integer, typename Point>
std::size_t get_Khalimsky_coordinate_from_index_coords(const Integer& max_x, const Integer& max_y, const Point& coords)
{
  return coords[2]*max_x*max_y*4 + coords[1]*max_x*2 + coords[0];
}

// normals are for surfaces
// values are for all elements
template <typename KSpace, typename SharpnessPropertyMap, typename Sharpness_type, typename SurfelNormalPropertyMap, typename NormalType>
void compute_AT_normals_and_values_from_surface(
  const KSpace& kSpace,
  const typename DGtal::Shortcuts<KSpace>::DigitalSurface &surface,
  const typename DGtal::Shortcuts<KSpace>::CellRange &pointels,
  const typename DGtal::Shortcuts<KSpace>::CellRange &linels,
  const typename DGtal::Shortcuts<KSpace>::CellRange &surfels,
  const typename DGtal::Shortcuts<KSpace>::SurfelRange &signed_surfels,
  const typename DGtal::Shortcuts<KSpace>::RealVectors &normals,
  SharpnessPropertyMap sharpness_property_map,
  SurfelNormalPropertyMap surfel_normal_property_map,
  const double alpha = 0.05,
  const double lambda = 0.02,
  const double width_start = 2.0,
  const double width_end = 0.25,
  const double width_divide = 2.0)
{
  typedef typename DGtal::Shortcuts<KSpace> SH;

  // DGtal::CanonicSCellEmbedder<KSpace> embedder = SH::getSCellEmbedder(kSpace);
  // DGtal::CanonicCellEmbedder<KSpace> uembedder = SH::getCellEmbedder(kSpace);

  // create solver
  typedef DGtal::DiscreteExteriorCalculusFactory<DGtal::EigenLinearAlgebraBackend> CalculusFactory;
  const auto calculus = CalculusFactory::createFromNSCells<2>(signed_surfels.begin(), signed_surfels.end());
  DGtal::ATSolver2D<KSpace> at_solver(calculus, 1);
  at_solver.initInputVectorFieldU2(normals, signed_surfels.cbegin(), signed_surfels.cend());
  at_solver.setUp(alpha, lambda);
  at_solver.solveGammaConvergence(width_start, width_end, width_divide);

  // Vector_field_output_translator normal_translator
      // <SurfelNormalPropertyMap, KSpace, NormalType>
      // (surfel_normal_property_map, kSpace, surfels, embedder);
  // at_solver.getOutputVectorFieldU2(normal_translator, surfels.cbegin(), surfels.cend());
  // out_surfaces.resize(surfels.size());
  // std::transform(surfels.cbegin(), surfels.cend(), out_surfaces.begin(),
  // [embedder](const DGtal_types::SH::SCell &c)
  // { return embedder(c); });

  // recover normals on surfels
  std::vector<typename SH::RealVector> out_surfel_normals;
  out_surfel_normals.resize(normals.size());
  at_solver.getOutputVectorFieldU2(out_surfel_normals, signed_surfels.cbegin(), signed_surfels.cend());

  int max_x = kSpace.max(0)+1;
  int max_y = kSpace.max(1)+1;

  for (std::size_t surfel_id = 0, number_of_signed_surfels = signed_surfels.size(); surfel_id < number_of_signed_surfels; surfel_id++)
  {
    const typename SH::SCell& surfel_cell = signed_surfels[surfel_id];
    const typename SH::RealVector& normal = out_surfel_normals[surfel_id];
    const std::size_t surfel_coordinate_index = get_Khalimsky_coordinate_from_index_coords(max_x, max_y, kSpace.sKCoords(surfel_cell));
    put(surfel_normal_property_map, surfel_coordinate_index, NormalType(normal[0], normal[1], normal[2]));
  }

  // recover sharpness on pointels
  std::vector<double> out_pointel_values;
  out_pointel_values.resize(pointels.size());
  at_solver.getOutputScalarFieldV0(out_pointel_values, pointels.cbegin(), pointels.cend(), at_solver.Average);
  for (std::size_t pointel_id = 0, number_of_pointels = pointels.size(); pointel_id < number_of_pointels; pointel_id++)
  {
    const typename SH::Cell& pointel_cell = pointels[pointel_id];
    const double& value = 1.0 - out_pointel_values[pointel_id]; // AT values: 0 is sharp, 1 is smooth. We want: 1 is sharp, 0 is smooth
    const std::size_t pointel_coordinate_index = get_Khalimsky_coordinate_from_index_coords(max_x, max_y, kSpace.uKCoords(pointel_cell));
    put(sharpness_property_map, pointel_coordinate_index, Sharpness_type(value));
  }

  // recover sharpness on linel
  std::vector<double> out_linel_values;
  out_linel_values.resize(linels.size());
  at_solver.getOutputScalarFieldV0(out_linel_values, linels.cbegin(), linels.cend(), at_solver.Average);
  for (std::size_t linel_id = 0, number_of_linels = linels.size(); linel_id < number_of_linels; linel_id++)
  {
    const typename SH::Cell& linel_cell = linels[linel_id];
    const double& value = 1.0 - out_linel_values[linel_id]; // AT values: 0 is sharp, 1 is smooth. We want: 1 is sharp, 0 is smooth
    const std::size_t linel_coordinate_index = get_Khalimsky_coordinate_from_index_coords(max_x, max_y, kSpace.uKCoords(linel_cell));
    put(sharpness_property_map, linel_coordinate_index, Sharpness_type(value));
  }

  // recover sharpness on surfel
  std::vector<double> out_surfel_values;
  out_surfel_values.resize(surfels.size());
  at_solver.getOutputScalarFieldV0(out_surfel_values, surfels.cbegin(), surfels.cend(), at_solver.Average);
  for (std::size_t surfel_id = 0, number_of_surfels = surfels.size(); surfel_id < number_of_surfels; surfel_id++)
  {
    const typename SH::Cell& surfel_cell = surfels[surfel_id];
    const double& value = 1.0 - out_surfel_values[surfel_id]; // AT values: 0 is sharp, 1 is smooth. We want: 1 is sharp, 0 is smooth
    const std::size_t surfel_coordinate_index = get_Khalimsky_coordinate_from_index_coords(max_x, max_y, kSpace.uKCoords(surfel_cell));
    put(sharpness_property_map, surfel_coordinate_index, Sharpness_type(value));
  }



  // DGtal_types::SH::Scalars features(pointels.size());
  // // put sharpness in property map
  // //    first determine values for linels and surfel (strategy : average of adjacent pointels)
  // std::unordered_map<std::size_t, double> linel_values; // holds sharpnes * 2
  // std::unordered_map<std::size_t, double> surfel_values; // holds sharpnes * 4
  // for (std::size_t pointel_id = 0, number_of_pointels = pointels.size(); pointel_id < number_of_pointels; pointel_id++)
  // {
  //   const SH::SCell& pointel_cell = pointels[pointel_id];
  //   const double& sharpness = out_values[pointel_id];
  //   const SH::RealPoint& point = embedder(surfel_cell);
  //   // for each linel adjacent to the pointel (may be outside of the surface)
  //   for (int axis = 0; axis < 3; axis++)
  //     for (int change = -1; change < 2; change += 2) // change is -1 then +1
  //     {
  //       SH::RealPoint linel_point = point;
  //       linel_point[axis] += change;
  //       if (linel_point[axis] < 0 || linel_point[axis] > kSpace.max(axis)) // Khalimsky coordinate indices are not defined outside of the space
  //         continue;
  //       const std::size_t linel_coordinate_index = get_Khalimsky_coordinate_index_from_DGtal_point(kSpace, linel_point);
  //       linel_values[linel_coordinate_index] += sharpness;
  //     }
  //   // for each surfel adjacent to the pointel (may be outside of the surface)
  //   for (int axis = 0; axis < 3; axis++)
  //   {
  //     int axis_1 = (axis+1)%3;
  //     int axis_2 = (axis+2)%3;
  //     // change on the two other axis
  //     for (int change_axis_1 = -1; change_axis_1 < 2; change_axis_1 += 2) // change_axis_1 is -1 then +1
  //       for (int change_axis_2 = -1; change_axis_2 < 2; change_axis_2 += 2) // change_axis_2 is -1 then +1
  //       {
  //         SH::RealPoint surfel_point = point;
  //         surfel_point[axis_1] += change_axis_1;
  //         if (surfel_point[axis_1] < 0 || surfel_point[axis_1] > kSpace.max(axis_1)) // Khalimsky coordinate indices are not defined outside of the space
  //           continue;
  //         surfel_point[axis_2] += change_axis_2;
  //         if (surfel_point[axis_2] < 0 || surfel_point[axis_2] > kSpace.max(axis_2)) // Khalimsky coordinate indices are not defined outside of the space
  //           continue;
  //         const std::size_t surfel_coordinate_index = get_Khalimsky_coordinate_index_from_DGtal_point(kSpace, surfel_point);
  //         surfel_values[surfel_coordinate_index] += sharpness;
  //       }
  //   }
  //   // add pointel to map
  //   const std::size_t pointel_coordinate_index = get_Khalimsky_coordinate_index_from_DGtal_point(kSpace, point);
  //   put(pointel_normal_property_map, pointel_coordinate_index, sharpness);
  // }
  // //   add linel and surfels to map
  // for (const auto& [linel_index, value] : linel_values)
  //   put(pointel_normal_property_map, linel_index, value * .5); // divide by 2
  // for (const auto& [surfel_index, value] : surfel_values)
  //   put(pointel_normal_property_map, surfel_index, value * .25); // divide by 4

  // out_valuePoints.resize(pointels.size());
  // std::transform(pointels.cbegin(), pointels.cend(), out_valuePoints.begin(),
  // [surface_info](const DGtal_types::SH::Cell &v)
  // { return uembedder(v); });
  // out_values.resize(features.size());
  // std::transform(features.cbegin(), features.cend(), out_values.begin(),
  // [](const double &v) { return v; });
}

// struct DGtal_mesh_face_graph
// {
//     struct DGtal_mesh_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
//                                                  public virtual boost::vertex_list_graph_tag,
//                                                  public virtual boost::edge_list_graph_tag,
//                                                  public virtual boost::adjacency_graph_tag
//     {};

//     typedef std::size_t vertex_descriptor;
//     typedef std::size_t halfedge_descriptor;
//     typedef std::size_t edge_descriptor;
//     typedef std::size_t face_descriptor;
//     typedef boost::undirected_tag directed_category;
//     typedef boost::disallow_parallel_edge_tag edge_parallel_category;

// };

} /* namespace internal */

// SharpnessPropertyMap assigns a sharpness value for each element of the image
//      It must define `put(SharpnessPropertyMap map, std::size_t element_index, SharpnesValue sharpness)`
// SurfelNormalPropertyMap assigns a normal for each surfel of the image
//      It must define `put(SurfelNormalPropertyMap map, std::size_t surfel_index, NormalVector normal)`
template <typename K, typename Image, typename SharpnessPropertyMap, typename Sharpness_type, typename SurfelNormalPropertyMap, typename NormalVector>
void compute_AmbrosioTortorelli_on_labelled_image(const Image& image,
  SharpnessPropertyMap sharpness_property_map,
  SurfelNormalPropertyMap surfel_normal_property_map,
  const double alpha = 0.05,
  const double lambda = 0.02,
  const double width_start = 2.0,
  const double width_end = 0.25,
  const double width_divide = 2.0,
  const float& value_outside = 0)
{
  // make a DGtal surface for each volume of the labeled image
  typedef DGtal::Z3i::KSpace KSpace;
  typedef typename DGtal::Shortcuts<KSpace> SH;
  typedef typename DGtal::ShortcutsGeometry<KSpace> SHG;
  typedef typename SH::DigitalSurface DigitalSurface;
  typedef typename KSpace::Point Integer_point;
  DGtal::Parameters params = SH::defaultParameters() | SHG::defaultParameters();
  Integer_point lastPoint;
  lastPoint[0] = image.xdim() - 1;
  lastPoint[1] = image.ydim() - 1;
  lastPoint[2] = image.zdim() - 1;
  KSpace kSpace = SH::getKSpace(Integer_point::zero, lastPoint, params);
  std::vector<DGtal::CountedPtr<DigitalSurface>> surfaces;
  std::vector<float> surfaces_domain;
  internal::make_DGtal_surface_from_image<KSpace, Image>(kSpace, image, surfaces, surfaces_domain);
  // compute sharpness and normals
  for (std::size_t surface_index = 0, number_of_surfaces = surfaces.size(); surface_index < number_of_surfaces; surface_index++)
  {
    DGtal::CountedPtr<DigitalSurface> surface_pointer = surfaces[surface_index];
    const float& domain = surfaces_domain[surface_index];
    internal::LabelledImageToBinaryWrapper binary_image(image, domain);
    // compute initial normals
    typename SH::SurfelRange signed_surfels = SH::getSurfelRange(surface_pointer, params);
    typename SH::RealVectors normals = SHG::getIINormalVectors(binary_image, kSpace, signed_surfels, params);
    internal::regularize_inversion_on_normals<KSpace>(*surface_pointer, signed_surfels, normals);
    // compute AmbrosioTortorelli
    typename SH::CellRange pointels = SH::getCellRange(surface_pointer, 0);
    typename SH::CellRange linels = SH::getCellRange(surface_pointer, 1);
    // typename SH::CellRange surfels = SH::getCellRange(surface_pointer, 2);
    typename SH::CellRange surfels;
    surfels.resize(signed_surfels.size());
    std::transform( signed_surfels.cbegin(), signed_surfels.cend(), surfels.begin(),
                  [&] (const typename KSpace::SCell& c) { return kSpace.uCell( kSpace.sKCoords(c) ); } );
    internal::compute_AT_normals_and_values_from_surface<KSpace, SharpnessPropertyMap, Sharpness_type, SurfelNormalPropertyMap, NormalVector>(
      /*In */ kSpace, *surface_pointer,
      /*Out*/ pointels, linels, surfels, signed_surfels, normals, sharpness_property_map, surfel_normal_property_map,
      /*Optional In*/ alpha, lambda, width_start, width_end, width_divide);
    }
    // // add values to property map
    // internal::convert_dgtal_coordinates(pointels);
    // internal::convert_dgtal_coordinates(surfels);
    // internal::make_pointel_value_map(pointels, values, point_values);
    // internal::add_values_of_pointels_to_graph(pointels, values, surface_graph);
    // internal::add_normals_of_surfel_to_graph(surfels, normals, surface_graph);
  }
} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_INTERNAL_DGTAL_INTERFACE_H