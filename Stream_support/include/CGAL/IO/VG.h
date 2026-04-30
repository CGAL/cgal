// Copyright (c) 2026 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sven Oesau

#ifndef CGAL_IO_VG_H
#define CGAL_IO_VG_H


#include <CGAL/IO/Color.h>
#include <CGAL/property_map.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <boost/tti/has_member_function.hpp>
#include <iostream>
#include <fstream>
#include <string>

namespace CGAL {
namespace IO {
namespace internal {

inline
void trim(std::string& s) {
  auto l = [](int c) {return !std::isspace(c); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), l));
  s.erase(std::find_if(s.rbegin(), s.rend(), l).base(), s.end());
}

inline
bool next_line(std::istream& is, std::string& line) {
  do {
    if (!std::getline(is, line))
      return false;
    trim(line);
  } while (line.empty() || line[0] == '#');
  return true;
}

BOOST_TTI_HAS_MEMBER_FUNCTION(reserve)
template <typename Range, bool has_reserve = has_member_function_reserve<Range, void, boost::mpl::vector<std::size_t>>::value>
struct Reserver {
};

template <typename Range>
struct Reserver<Range, false> {
  void operator()(Range&, std::size_t) {}
};

template <typename Range>
struct Reserver<Range, true> {
  void operator()(Range& r, std::size_t n) { r.reserve(n); }
};

template <typename Map>
bool read_elements(std::istream& is, Map cmap, std::size_t num, std::string& line) {
  using Type = typename boost::property_traits<Map>::value_type;
  std::istringstream ss;
  for (std::size_t i = 0; i < num; ++i) {
    if (!internal::next_line(is, line))
      return false;

    if (!std::isdigit(line[0]))
      return false;

    ss.clear();
    ss.str(line);

    double a, b, c;
    if (ss >> iformat(a) >> iformat(b) >> iformat(c))
      put(cmap, i, Type(a, b, c));
  }
  next_line(is, line);

  return true;
}

template <typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_colors(std::ifstream& is, std::size_t num_colors, std::string& line, const CGAL_NP_CLASS& np) {
  using parameters::is_default_parameter;
  if constexpr (!(is_default_parameter<CGAL_NP_CLASS, internal_np::color_map_t>::value)) {
    using Color_map = typename internal_np::Lookup_named_param_def<internal_np::color_map_t, CGAL_NP_CLASS, Constant_property_map<std::size_t, Color>>::type;
    Color_map color_map = parameters::choose_parameter<Color_map>(parameters::get_parameter(np, internal_np::color_map));
    return internal::read_elements(is, color_map, num_colors, line);
  }
  else while (internal::next_line(is, line) && std::isdigit(line[0]));

  return true;
}

template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_normals(std::ifstream& is, std::size_t num_normals, std::string& line, PointRange points, const CGAL_NP_CLASS& np) {
  using parameters::is_default_parameter;
  using NP_helper = Point_set_processing_3_np_helper<PointRange, CGAL_NP_CLASS>;
  if constexpr (!(is_default_parameter<CGAL_NP_CLASS, internal_np::normal_t>::value)) {
    typename NP_helper::Normal_map normal_map = NP_helper::get_normal_map(points, np);
    return internal::read_elements(is, normal_map, num_normals, line);
  }
  else while (internal::next_line(is, line) && (std::isdigit(line[0]) || line[0] == '-'));

  return true;
}

} // namespace internal

/*!
 * \ingroup PkgStreamSupportIoFuncsVG
 *
 * \brief reads the content of the input stream `is` into `points` and `regions` according to the VG (Vertex Group) format.
 *
 * \attention The format does not support binary streams.
 *
 * \tparam PointRange
 *         a model of the concepts `RandomAccessContainer` and `BackInsertionSequence` whose value type is the point type
 * \tparam RegionOutputIterator type of the output iterator who must accept `std::pair<Primitive, std::vector<item>>`,
 *         where item is the key type of `Point_map`. `Primitive` needs either to implement `operator >>` to be read the input stream
 *         `is` or to provide the named parameter `constructor` described below to be created from `std::vector<FT>`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param is input stream
 * \param points point range to which the read points are written
 * \param rout output iterator to which the read regions are written
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{normal_map}
 *     \cgalParamDescription{a property map associating normals to the elements inserted into `points`}
 *     \cgalParamType{a model of `WritablePropertyMap` whose key type is `std::size_t` and whose value type is `geom_traits::Vector_3`}
 *     \cgalParamDefault{If this parameter is omitted, no normals are read.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{color_map}
 *     \cgalParamDescription{a property map associating colors to the elements inserted into `points`}
 *     \cgalParamType{a model of `WritablePropertyMap` whose key type is `std::size_t` and whose value type is `CGAL::Color`}
 *     \cgalParamDefault{If this parameter is omitted, no colors are read.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{constructor}
 *     \cgalParamDescription{a function taking the group type and `std::vector<FT>` as input to create an instance of `Primitive`}
 *     \cgalParamType{a function with the signature `Primitive (int type, std::vector<FT>&)` or `Primitive (int type, std::vector<FT>&&)`}
 *     \cgalParamDefault{Not provided, the `operator >>` of `Primitive` is used.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{labels}
 *     \cgalParamDescription{a range of labels providing one label for each region.}
 *     \cgalParamType{a model of the concept `RandomAccessContainer` whose value type is `std::string` or `const char *`}
 *     \cgalParamDefault{No group labels are read.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{colors}
 *     \cgalParamDescription{a range of colors providing one color for each region.}
 *     \cgalParamType{a model of the concept `RandomAccessContainer` whose value type is `CGAL::Color`}
 *     \cgalParamDefault{No group colors are read.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the writing was successful, `false` otherwise.
 */

template <typename PointRange, typename RegionOutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_VG(std::ifstream& is,
  PointRange& points,
  RegionOutputIterator rout,
  const CGAL_NP_CLASS& np = parameters::default_values()) {
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  using Point = typename std::iterator_traits<typename PointRange::iterator>::value_type;
  using GeomTraits = typename CGAL::Kernel_traits<Point>::type;
  using Region = typename std::iterator_traits<typename RegionOutputIterator::container_type::iterator>::value_type;
  using FT = typename GeomTraits::FT;

  if (!is.good()) {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  std::string line;
  if (!internal::next_line(is, line))
    return false;

  std::istringstream ss(line);
  std::string token;

  ss >> token;

  if (token != "num_points:" && token != "num_points")
    return false;

  std::size_t num_points;
  ss >> num_points;

  if (num_points == 0)
    return false;

  internal::Reserver<PointRange> r;
  r(points, num_points);

  std::size_t added = 0;
  while (internal::next_line(is, line)) {
    FT x, y, z;
    ss.clear();
    ss.str(line);
    if (ss >> iformat(x) >> iformat(y) >> iformat(z)) {
      points.push_back(Point(x, y, z));
      added++;
    }
    else {
      break;
    }
  }

  if (added != num_points) {
    std::cerr << "Error: could only read " << added << " points from file instead of specified " << num_points << std::endl;
    return false;
  }

  ss.clear();
  ss.str(line);
  ss >> token;

  while (token != "num_groups" && token != "num_groups:") {
    if (token == "num_colors:" || token == "num_colors") {
      std::size_t num_colors;
      ss >> num_colors;
      internal::read_colors(is, num_colors, line, np);
    }
    else
      if (token == "num_normals" || token == "num_normals:") {
        std::size_t num_normals;
        ss >> num_normals;
        internal::read_normals(is, num_normals, line, points, np);
      }
      else {
        std::cerr << "Error: Expected \"num_normals:\", \"num_colors:\" or \"num_groups:\", but found " << line << std::endl;
        return false;
      }

    ss.clear();
    ss.str(line);
    ss >> token;
  }

  std::size_t num_groups;
  ss >> num_groups;

  if (token != "num_groups" && token != "num_groups:") {
    std::cerr << "Error: expected 'num_groups' token, but got '" << token << "' instead" << std::endl;
    return false;
  }

  if (num_groups == 0)
    return true;

  for (std::size_t i = 0; i < num_groups; ++i) {
    int num_group_parameters = -1;
    unsigned int type = -1;
    Region region;
    while (internal::next_line(is, line)) {
      ss.clear();
      ss.str(line);

      if (!(ss >> token)) {
        std::cerr << "Error: expected group information" << std::endl;
        return false;
      }

      if (std::isdigit(token[0])) {
        std::cerr << "Error: expected group attribute, but found digit " << token << std::endl;
        return false;
      }

      if (token == "num_children:" || token == "num_children")
        continue;

      if (token == "group_type:" || token == "group_type") {
        ss >> type;
        continue;
      }

      if (token == "group_label:" || token == "group_label") {
        if constexpr (!(is_default_parameter<CGAL_NP_CLASS, internal_np::labels_t>::value)) {
          using Labels = typename internal_np::Lookup_named_param_def<internal_np::labels_t, CGAL_NP_CLASS, std::vector<std::string>>::type;
          Labels labels = parameters::choose_parameter<Labels>(parameters::get_parameter(np, internal_np::labels));
          labels.resize(i + 1); // fill missing labels with empty strings if there are not enough labels provided
          ss >> labels[i];
        }
        continue;
      }

      if (token == "group_colors:" || token == "group_colors") {
        if constexpr (!(is_default_parameter<CGAL_NP_CLASS, internal_np::colors_t>::value)) {
          using Colors = typename internal_np::Lookup_named_param_def<internal_np::colors_t, CGAL_NP_CLASS, std::vector<Color>>::type;
          Colors colors = parameters::choose_parameter<Colors>(parameters::get_parameter(np, internal_np::colors));
          colors.resize(i + 1); // fill missing colors with empty strings if there are not enough colors provided
          float r, g, b;
          ss >> r >> g >> b;
          r = std::clamp(r, 0.f, 1.f);
          g = std::clamp(g, 0.f, 1.f);
          b = std::clamp(b, 0.f, 1.f);
          colors[i] = Color(int(r * 255 + 0.5), int(g * 255 + 0.5), int(b * 255 + 0.5));
        }
        continue;
      }

      if (token == "num_group_parameters:" || token == "num_group_parameters") {
        if (!(ss >> num_group_parameters)) {
          std::cerr << "Error: expected number of group parameters after 'num_group_parameters' token, but got '" << line << "' instead" << std::endl;
          return false;
        }
        continue;
      }

      if (token == "group_parameters:" || token == "group_parameters") {
        if (num_group_parameters < 0) {
          std::cerr << "Error: expected 'num_group_parameters' token before 'group_parameters' token" << std::endl;
          return false;
        }

        if (num_group_parameters > 0) {
          if constexpr (is_default_parameter<CGAL_NP_CLASS, internal_np::constructor_t>::value) {
            //typename GeomTraits::Plane_3 pl;
            ss >> region.first;
          }
          else {
            using Constructor = typename internal_np::Lookup_named_param_def<internal_np::constructor_t, CGAL_NP_CLASS, int>::type;
            Constructor constructor = parameters::choose_parameter<Constructor>(parameters::get_parameter(np, internal_np::constructor));
            std::string params;
            std::getline(ss, params);
            region.first = std::move(constructor(type, params));
          }
        }
        continue;
      }

      if (token == "group_num_points:" || token == "group_num_points") {
        std::size_t num_group_points;
        if (!(ss >> num_group_points)) {
          std::cerr << "Error: expected number of group points after 'group_num_points' token, but got '" << line << "' instead" << std::endl;
          return false;
        }

        std::size_t id;
        for (std::size_t j = 0; j < num_group_points; ++j) {
          if (!(ss >> id))
            break;
          region.second.push_back(id);
        }

        if (region.second.size() < num_group_points) {
          // if not all point indices were on the same line, read the remaining ones from the following line
          internal::next_line(is, line);
          ss.clear();
          ss.str(line);

          for (std::size_t j = region.second.size(); j < num_group_points; ++j) {
            std::size_t id;
            if (!(ss >> id)) {
              std::cerr << "Error: expected point index after 'group_num_points' token, but got '" << line << "' instead" << std::endl;
              return false;
            }
            region.second.push_back(id);
          }
        }

        *rout++ = std::move(region);

        break;
      }
    }
  }

  return true;
}
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

template <typename PointRange, typename RegionOutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_VG(const std::string& filename,
  PointRange& points,
  RegionOutputIterator rout,
  const CGAL_NP_CLASS& np = parameters::default_values()) {
  std::ifstream is(filename);
  return read_VG(is, points, rout, np);
}

} // namespace IO
} // namespace CGAL

#endif
