// Copyright (c) 2026 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sven Oesau
//
#ifndef CGAL_IO_WRITE_VG_H
#define CGAL_IO_WRITE_VG_H

#include <CGAL/IO/Color.h>
#include <CGAL/property_map.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <iostream>
#include <fstream>

namespace CGAL {
namespace IO {

/*!
 * \ingroup PkgStreamSupportIoFuncsVG
 *
 * \brief writes the content of `points` and `regions` to the output stream `os` in the VG (Vertex Group) format.
 *
 * \attention The format does not support binary streams.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam RegionRange
 *    a model of `ForwardRange` whose value type is `std::pair<Primitive, std::vector<item>>`,
 *    where item is the key type of `Point_map`. `Primitive` is the geometric primitive and needs to implement
 *    `operator <<` for the output stream `os`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param os the output stream
 * \param points point set
 * \param regions a range of regions providing their algebraic primitive and elements of the point cloud.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{serializer}
 *     \cgalParamDescription{serializes `Primitive` into a string to be written to `os` and provides the type of the primitive (unsigned int) and the number of written parameters.}
 *     \cgalParamType{a function with the signature `const char* (Primitive &p, unsigned int &type, std::size_t& number_of_parameters)`}
 *     \cgalParamDefault{The `Primitive` is not exported into `os`. Default handler are provided for `geom_traits::Plane_3` and `geom_traits::Sphere_3`}
 *     \cgalParamExtra{The type of primitive is defined as 0: Plane, 1: Cylinder, 2: Sphere, 3: Cone, 4: Torus, 5: General}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the point set `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
 *                    of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
 *     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{normal_map}
 *     \cgalParamDescription{a property map associating normals to the elements of the point range}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
 *                    of the iterator of `PointRange` and whose value type is `geom_traits::Vector_3`}
 *     \cgalParamDefault{If this parameter is omitted, no normals are written to `os`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{color_map}
 *     \cgalParamDescription{a property map associating colors to the elements of the point range}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
 *                    of the iterator of `PointRange` and whose value type is `CGAL::Color`}
 *     \cgalParamDefault{If this parameter is omitted, no colors are written to `os`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{labels}
 *     \cgalParamDescription{a range of labels providing one label for each region.}
 *     \cgalParamType{a model of the concept `RandomAccessContainer` whose value type is `std::string` or `const char *`}
 *     \cgalParamDefault{No group labels are written to `os`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{colors}
 *     \cgalParamDescription{a range of colors providing one color for each region.}
 *     \cgalParamType{a model of the concept `RandomAccessContainer` whose value type is `CGAL::Color`}
 *     \cgalParamDefault{No group colors are written to `os`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{the precision of the stream `os`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the writing was successful, `false` otherwise.
 */

template <typename PointRange, typename RegionRange, typename NamedParameters>
bool write_VG(std::ostream& os,
              const PointRange& points,
              const RegionRange& regions,
              const NamedParameters& np = parameters::default_values()) {
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, NamedParameters>;
  using PointMap = typename NP_helper::Const_point_map;
  using Point = typename boost::property_traits<PointMap>::value_type;
  using GeomTraits = typename NP_helper::Geom_traits;
  using Region = typename std::iterator_traits<typename RegionRange::iterator>::value_type;
  using Primitive = typename Region::first_type;

  PointMap point_map = NP_helper::get_const_point_map(points, np);

  if (!os.good()) {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  set_stream_precision_from_NP(os, np);

  os << "num_points: " << points.size() << "\n";

  for (const auto& p : points)
    os << get(point_map, p) << "\n";

  os << "\n";

  if constexpr (!(is_default_parameter<NamedParameters, internal_np::color_map_t>::value)) {
    using DummyColorMap = Constant_property_map<Point, Color>;
    using Color_map = typename internal_np::Lookup_named_param_def<internal_np::color_map_t, NamedParameters, DummyColorMap> ::type;
    Color_map color_map = parameters::choose_parameter<Color_map>(parameters::get_parameter(np, internal_np::color_map));
    os << "num_colors: " << points.size() << "\n";
    for (const auto& p : points)
      os << get(color_map, p) << "\n";
    os << "\n";
  }

  if constexpr (!(is_default_parameter<NamedParameters, internal_np::normal_t>::value)) {
    typename NP_helper::Normal_map normal_map = NP_helper::get_normal_map(points, np);
    os << "num_normals: " << points.size() << "\n";
    for (const auto& p : points)
      os << get(normal_map, p) << "\n";
    os << "\n";
  }

  if (!regions.empty()) {
    os << "num_groups: " << regions.size() << "\n";
    std::size_t idx = 0;
    for (const Region& r : regions) {
      /*
      * group_type: type     # integer denoting the of the segment (0: PLANE, 1: CYLINDER, 2: SPHERE, 3: CONE, 4: TORUS, 5: GENERAL)
      * num_group_parameters: NUM_GROUP_PARAMETERS    # integer number denoting the number of floating point values representing the segment (e.g., 4 for planes)
      * group_parameters: float[NUM_GROUP_PARAMETERS] # a sequence of NUM_GROUP_PARAMETERS floating point numbers (e.g., a, b, c, and d for a plane)
      * group_label: label   # the label (a string) of the segment
      * group_color: r g b   # 3 floating point numbers denoting the color of this segment
      * group_num_points: N  # N is an integer denoting the number of points in this segment (can be 0)
      * id1 ... idN          # N integer numbers denoting the indices of the points in this segment
      * num_children: num    # a segment/primitive/object may contain subsegment (that has the same representation as this segment)
      */
      if constexpr (!(is_default_parameter<NamedParameters, internal_np::serializer_t>::value)) {
        using Serializer = typename internal_np::Lookup_named_param_def<internal_np::serializer_t, NamedParameters, std::string>::type;
        Serializer serializer = parameters::choose_parameter<Serializer>(parameters::get_parameter(np, internal_np::serializer));
        unsigned int type = 5;
        std::size_t number_of_parameters = 0;
        std::string str = serializer(r.first, type, number_of_parameters);
        os << "group_type: " << type << "\n";
        os << "num_group_parameters: " << number_of_parameters << "\n";
        os << "group_parameters: " << str << "\n";
      }
      else if constexpr (std::is_same<typename GeomTraits::Plane_3, Primitive>::value) {
        os << "group_type: 0\n";
        os << "num_group_parameters: 4\n";
        os << "group_parameters: " << r.first << "\n";
      }
      else if constexpr (std::is_same<typename GeomTraits::Sphere_3, Primitive>::value) {
        os << "group_type: 2\n";
        os << "num_group_parameters: 4\n";
        os << "group_parameters: " << r.first << "\n";
      }
      else {
        os << "group_type: 5\n";
        os << "num_group_parameters: 0\n";
        os << "group_parameters: \n";
      }

      if constexpr (!(is_default_parameter<NamedParameters, internal_np::labels_t>::value)) {
        using Labels = typename internal_np::Lookup_named_param_def<internal_np::labels_t, NamedParameters, std::vector<std::string>>::type;
        Labels labels = parameters::choose_parameter(parameters::get_parameter(np, internal_np::labels), Labels());
        if (idx < labels.size())
          os << "group_label: " << labels[idx] << "\n";
      }

      if constexpr (!(is_default_parameter<NamedParameters, internal_np::colors_t>::value)) {
        using Colors = typename internal_np::Lookup_named_param_def<internal_np::colors_t, NamedParameters, std::vector<Color>>::type;
        Colors colors = parameters::choose_parameter(parameters::get_parameter(np, internal_np::colors), Colors());
        if (idx < colors.size())
          os << "group_color: " << colors[idx] << "\n";
      }

      os << "group_num_points: " << r.second.size() << "\n";
      for (const auto& id : r.second) {
        os << id << " ";
      }

      os << "\nnum_children: 0\n";
      idx++;
    }
  }

  return os.good();
}

template <typename PointRange, typename Regions, typename NamedParameters>
bool write_VG(const std::string& filename,
              const PointRange& points,
              const Regions& regions,
              const NamedParameters& np = parameters::default_values()) {
  std::ofstream os(filename);
  return write_VG(os, points, regions, np);
}

} // namespace IO
} // namespace CGAL

#endif // CGAL_IO_WRITE_VG_H
