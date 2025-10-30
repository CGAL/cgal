// Copyright (c) 2017  GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_IO_LAS_READ_LAS_POINTS_H
#define CGAL_IO_LAS_READ_LAS_POINTS_H

#include <CGAL/config.h>

#include <CGAL/IO/LAS/Las_property.h>
#include <CGAL/IO/LAS.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/type_traits/is_iterator.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>

#ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable:4251) // DLL warning from LASlib
#endif

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#define USE_AS_DLL 1
#include <lasreader_las.hpp>
#undef USE_AS_DLL

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif

#ifdef BOOST_MSVC
#  pragma warning(pop)
#endif

namespace CGAL {

namespace IO {



namespace internal {

namespace LAS {

inline void get_value(const LASpoint& r, double& v, LAS_property::X&)
{ v = r.get_x(); }
inline void get_value(const LASpoint& r, double& v, LAS_property::Y&)
{ v = r.get_y(); }
inline void get_value(const LASpoint& r, double& v, LAS_property::Z&)
{ v = r.get_z(); }
inline void get_value(const LASpoint& r, unsigned short& v, LAS_property::Intensity&)
{ v = r.get_intensity(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::Return_number&)
{ v = r.get_return_number(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::Number_of_returns&)
{ v = r.get_number_of_returns(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::Scan_direction_flag&)
{ v = r.get_scan_direction_flag(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::Edge_of_flight_line&)
{ v = r.get_edge_of_flight_line(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::Classification&)
{ v = r.get_classification(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::Synthetic_flag&)
{ v = r.get_synthetic_flag(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::Keypoint_flag&)
{ v = r.get_keypoint_flag(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::Withheld_flag&)
{ v = r.get_withheld_flag(); }
inline void get_value(const LASpoint& r, float& v, LAS_property::Scan_angle&)
{ v = r.get_scan_angle(); }
inline void get_value(const LASpoint& r, unsigned char& v, LAS_property::User_data&)
{ v = r.get_user_data(); }
inline void get_value(const LASpoint& r, unsigned short& v, LAS_property::Point_source_ID&)
{ v = r.get_point_source_ID(); }
inline void get_value(const LASpoint& r, unsigned int& v, LAS_property::Deleted_flag&)
{ v = r.get_deleted_flag(); }
inline void get_value(const LASpoint& r, double& v, LAS_property::GPS_time&)
{ v = r.get_gps_time(); }
inline void get_value(const LASpoint& r, unsigned short& v, LAS_property::R&)
{ v = r.get_R(); }
inline void get_value(const LASpoint& r, unsigned short& v, LAS_property::G&)
{ v = r.get_G(); }
inline void get_value(const LASpoint& r, unsigned short& v, LAS_property::B&)
{ v = r.get_B(); }
inline void get_value(const LASpoint& r, unsigned short& v, LAS_property::I&)
{ v = r.get_I(); }


template <std::size_t N>
struct Filler
{
  template <class Value_tuple, class LAS_property_tuple>
  static void fill(const LASpoint& r, Value_tuple& values, LAS_property_tuple wrappers)
  {
    get_value(r, std::get<N>(values), std::get<N+2>(wrappers));
    Filler<N-1>::fill(r, values, wrappers);
  }
};

template<int ...>
struct seq { };

template<int N, int ...S>
struct gens : gens<N-1, N-1, S...> { };

template<int ...S>
struct gens<0, S...> {
  typedef seq<S...> type;
};

template<class ValueType, class Functor, class Tuple, int ...S>
ValueType call_functor(Functor f, Tuple t, seq<S...>) {
  return f(std::get<S>(t) ...);
}

template <class ValueType, class Functor, typename ... T>
ValueType call_functor(Functor f, std::tuple<T...>& t)
{
  return call_functor<ValueType>(f, t, typename gens<sizeof...(T)>::type());
}

template<>
struct Filler<0>
{
  template <class Value_tuple, class LAS_property_tuple>
  static void fill(const LASpoint& r, Value_tuple& values, LAS_property_tuple wrappers)
  {
    get_value(r, std::get<0>(values), std::get<2>(wrappers));
  }
};

template <typename OutputValueType, typename PropertyMap, typename T, LAS_property::Id::Id id>
void process_properties (const LASpoint& reader, OutputValueType& new_element,
                         std::pair<PropertyMap, LAS_property::Base<T,id> >&& current);

template <typename OutputValueType, typename PropertyMap, typename T, LAS_property::Id::Id id,
          typename NextPropertyBinder, typename ... PropertyMapBinders>
void process_properties (const LASpoint& reader, OutputValueType& new_element,
                         std::pair<PropertyMap, LAS_property::Base<T,id> >&& current,
                         NextPropertyBinder&& next,
                         PropertyMapBinders&& ... properties);

template <typename OutputValueType,
          typename PropertyMap,
          typename Constructor,
          typename ... T,
          LAS_property::Id::Id ... id>
void process_properties (const LASpoint& reader, OutputValueType& new_element,
                         std::tuple<PropertyMap, Constructor, LAS_property::Base<T,id>...>&& current)
{
  typedef typename PropertyMap::value_type PmapValueType;
  std::tuple<T...> values;
  Filler<sizeof...(T)-1>::fill(reader, values, current);
  PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
  put (std::get<0>(current), new_element, new_value);
}

template <typename OutputValueType,
          typename PropertyMap,
          typename Constructor,
          typename ... T,
          LAS_property::Id::Id ... id,
          typename NextPropertyBinder,
          typename ... PropertyMapBinders>
void process_properties (const LASpoint& reader, OutputValueType& new_element,
                         std::tuple<PropertyMap, Constructor, LAS_property::Base<T,id>...>&& current,
                         NextPropertyBinder&& next,
                         PropertyMapBinders&& ... properties)
{
  typedef typename PropertyMap::value_type PmapValueType;
  std::tuple<T...> values;
  Filler<sizeof...(T)-1>::fill(reader, values, current);
  PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
  put (std::get<0>(current), new_element, new_value);

  process_properties (reader, new_element, std::forward<NextPropertyBinder>(next),
                      std::forward<PropertyMapBinders>(properties)...);
}

template <typename OutputValueType, typename PropertyMap, typename T, LAS_property::Id::Id id>
void process_properties (const LASpoint& reader, OutputValueType& new_element,
                         std::pair<PropertyMap, LAS_property::Base<T,id> >&& current)
{
  T new_value = T();
  get_value (reader, new_value, current.second);
  put (current.first, new_element, new_value);
}

template <typename OutputValueType, typename PropertyMap, typename T, LAS_property::Id::Id id,
          typename NextPropertyBinder, typename ... PropertyMapBinders>
void process_properties (const LASpoint& reader, OutputValueType& new_element,
                         std::pair<PropertyMap, LAS_property::Base<T,id> >&& current,
                         NextPropertyBinder&& next,
                         PropertyMapBinders&& ... properties)
{
  T new_value = T();
  get_value (reader, new_value, current.second);
  put (current.first, new_element, new_value);
  process_properties (reader, new_element, std::forward<NextPropertyBinder>(next),
                      std::forward<PropertyMapBinders>(properties)...);
}

} // namespace LAS
} // namespace internal


// documenation in ../LAS.h
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename ... PropertyHandler>
bool read_LAS_with_properties(std::istream& is,
                              PointOutputIterator output,
                              PropertyHandler&& ... properties)
{
  typedef OutputIteratorValueType Enriched_point;

  if(!is)
    return false;

#if LAS_TOOLS_VERSION < 240319
  LASreaderLAS lasreader;
#else
  LASreaderLAS lasreader(nullptr);
#endif
  lasreader.open(is);

  while(lasreader.read_point())
  {
    const LASpoint& laspoint = lasreader.point;
    Enriched_point new_point;

    internal::LAS::process_properties (laspoint, new_point, std::forward<PropertyHandler>(properties)...);

    *(output ++) = new_point;
  }

  lasreader.close();

  return true;

}

/// \cond SKIP_IN_MANUAL

template <typename OutputIterator,
          typename ... PropertyHandler>
bool read_LAS_with_properties(std::istream& is,
                              OutputIterator output,
                              PropertyHandler&& ... properties)
{
  return read_LAS_with_properties<typename value_type_traits<OutputIterator>::type>(is, output, std::forward<PropertyHandler>(properties)...);
}

/// \endcond

// documenation in ../LAS.h
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool read_LAS(std::istream& is,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::type PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  return read_LAS_with_properties(is, output, make_las_point_reader(point_map));
}

/// \cond SKIP_IN_MANUAL

template <typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_LAS(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  return read_LAS<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

/// \endcond

// documentation in ../LAS.h
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool read_LAS(const std::string& filename,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np)
{
  std::ifstream is(filename, std::ios::binary);
  CGAL::IO::set_mode(is, CGAL::IO::BINARY);
  return read_LAS<OutputIteratorValueType>(is, output, np);
}

/// \cond SKIP_IN_MANUAL

template <typename OutputIterator,typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_LAS(const std::string& fname, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values())
{
  std::ifstream is(fname, std::ios::binary);
  CGAL::IO::set_mode(is, CGAL::IO::BINARY);
  return read_LAS<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

/// \endcond

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_LAS_READ_LAS_POINTS_H
