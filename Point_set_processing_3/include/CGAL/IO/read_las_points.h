// Copyright (c) 2017  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_READ_LAS_POINTS_H
#define CGAL_READ_LAS_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/config.h>

#include <tuple>

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <sstream>
#include <string>

#ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable:4251) // DLL warning from LASlib
#endif

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#define USE_AS_DLL
#include <lasreader_las.hpp>
#undef USE_AS_DLL

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif

#ifdef BOOST_MSVC
#  pragma warning(pop)
#endif


namespace CGAL {


  /// \cond SKIP_IN_MANUAL
  namespace LAS_property
  {
    namespace Id
    {
      enum Id
      {
        X,
        Y,
        Z,
        Intensity,
        Return_number,
        Number_of_returns,
        Scan_direction_flag,
        Edge_of_flight_line,
        Classification,
        Synthetic_flag,
        Keypoint_flag,
        Withheld_flag,
        Scan_angle,
        User_data,
        Point_source_ID,
        Deleted_flag,
        GPS_time,
        R,
        G,
        B,
        I
      };
    }

    template <typename T, Id::Id id>
    struct Base
    {
      typedef T type;
    };

    typedef Base<double, Id::X> X;
    typedef Base<double, Id::Y> Y;
    typedef Base<double, Id::Z> Z;
    typedef Base<unsigned short, Id::Intensity> Intensity;
    typedef Base<unsigned char, Id::Return_number> Return_number;
    typedef Base<unsigned char, Id::Number_of_returns> Number_of_returns;
    typedef Base<unsigned char, Id::Scan_direction_flag> Scan_direction_flag;
    typedef Base<unsigned char, Id::Edge_of_flight_line> Edge_of_flight_line;
    typedef Base<unsigned char, Id::Classification> Classification;
    typedef Base<unsigned char, Id::Synthetic_flag> Synthetic_flag;
    typedef Base<unsigned char, Id::Keypoint_flag> Keypoint_flag;
    typedef Base<unsigned char, Id::Withheld_flag> Withheld_flag;
    typedef Base<float, Id::Scan_angle> Scan_angle;
    typedef Base<unsigned char, Id::User_data> User_data;
    typedef Base<unsigned short, Id::Point_source_ID> Point_source_ID;
    typedef Base<unsigned int, Id::Deleted_flag> Deleted_flag;
    typedef Base<double, Id::GPS_time> GPS_time;
    typedef Base<unsigned short, Id::R> R;
    typedef Base<unsigned short, Id::G> G;
    typedef Base<unsigned short, Id::B> B;
    typedef Base<unsigned short, Id::I> I;
  }
  /// \endcond


  /**
     \ingroup PkgPointSetProcessing3IOLas

     Generates a %LAS property handler to read 3D points. Points are
     constructed from the input the using 3 %LAS properties
     `LAS_property::X`, `LAS_property::Y` and `LAS_property::Z`.

     \sa `read_las_points_with_properties()`

     \tparam PointMap the property map used to store points.
  */
  template <typename PointMap>
  std::tuple<PointMap,
             typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
             LAS_property::X, LAS_property::Y, LAS_property::Z >
  make_las_point_reader(PointMap point_map)
  {
    return std::make_tuple (point_map, typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3(),
                            LAS_property::X(), LAS_property::Y(), LAS_property::Z());
  }

/// \cond SKIP_IN_MANUAL

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


/// \endcond


/**
   \ingroup PkgPointSetProcessing3IOLas

   Reads user-selected points properties from a .las or .laz stream.
   Potential additional properties are ignored.

   Properties are handled through a variadic list of property
   handlers. A `PropertyHandler` can either be:

   - A `std::pair<PropertyMap, LAS_property::Tag >` if the user wants to
   read a %LAS property as a scalar value `LAS_property::Tag::type` (for
   example, storing an `int` %LAS property into an `int` variable).

   - A `std::tuple<PropertyMap, Constructor,
   LAS_property::Tag...>` if the user wants to use one or several
   %LAS properties to construct a complex object (for example,
   storing 4 `unsigned short` %LAS properties into a %Color object
   that can for example be a `std::array<unsigned short,
   4>`). In that case, the second element of the tuple should be a
   functor that constructs the value type of `PropertyMap` from N
   objects of of type `LAS_property::Tag::type`.

   The %LAS standard defines a fixed set of properties accessible
   through the following tag classes:

   - `LAS_property::X` with type `double`
   - `LAS_property::Y` with type `double`
   - `LAS_property::Z` with type `double`
   - `LAS_property::Intensity` with type `unsigned short`
   - `LAS_property::Return_number` with type `unsigned char`
   - `LAS_property::Number_of_returns` with type `unsigned char`
   - `LAS_property::Scan_direction_flag` with type `unsigned char`
   - `LAS_property::Edge_of_flight_line` with type `unsigned char`
   - `LAS_property::Classification` with type `unsigned char`
   - `LAS_property::Synthetic_flag` with type `unsigned char`
   - `LAS_property::Keypoint_flag` with type `unsigned char`
   - `LAS_property::Withheld_flag` with type `unsigned char`
   - `LAS_property::Scan_angle` with type `double`
   - `LAS_property::User_data` with type `unsigned char`
   - `LAS_property::Point_source_ID` with type `unsigned short`
   - `LAS_property::Deleted_flag` with type `unsigned int`
   - `LAS_property::GPS_time` with type `double`
   - `LAS_property::R` with type `unsigned short`
   - `LAS_property::G` with type `unsigned short`
   - `LAS_property::B` with type `unsigned short`
   - `LAS_property::I` with type `unsigned short`

   \cgalRequiresCPP11

   \sa `make_las_point_reader()`

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.
   \tparam PropertyHandler handlers to recover properties.

   \return `true` on success.
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename ... PropertyHandler>
bool read_las_points_with_properties (std::istream& stream,
                                      OutputIterator output,
                                      PropertyHandler&& ... properties)
{
  typedef OutputIteratorValueType Enriched_point;

  LASreaderLAS lasreader;
  lasreader.open(stream);

  while (lasreader.read_point())
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
bool read_las_points_with_properties (std::istream& stream,
                                      OutputIterator output,
                                      PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<OutputIterator>::type OutputValueType;

  return read_las_points_with_properties<OutputValueType>
    (stream, output, std::forward<PropertyHandler>(properties)...);
}
/// \endcond

/**
   \ingroup PkgPointSetProcessing3IOLas
   Reads points (position only) from a .las or .laz stream.
   Potential additional properties are ignored.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.

   \param stream input stream.
   \param output output iterator over points.

   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return true on success.

   \cgalRequiresCPP11
*/
template < typename OutputIteratorValueType,
           typename OutputIterator,
#ifdef DOXYGEN_RUNNING
           typename NamedParameters
#else
           typename CGAL_BGL_NP_TEMPLATE_PARAMETERS
#endif
>
bool read_las_points(std::istream& stream,
                     OutputIterator output,
#ifdef DOXYGEN_RUNNING
                     const NamedParameters& np)
#else
                     const CGAL_BGL_NP_CLASS& np)
#endif
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::type PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  return read_las_points_with_properties (stream, output,
                                          make_las_point_reader (point_map));
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename OutputIteratorValueType,
          typename OutputIterator>
bool
read_las_points(
  std::istream& stream, ///< input stream.
  OutputIterator output) ///< output iterator over points.
{
  return read_las_points<OutputIteratorValueType>
    (stream, output, CGAL::parameters::all_default());
}

// variant with default output iterator value type
template <typename OutputIterator,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool
read_las_points(
  std::istream& stream, ///< input stream.
  OutputIterator output,
  const CGAL_BGL_NP_CLASS& np)
{
  return read_las_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, np);
}

// variant with default NP and output iterator value type
template <typename OutputIterator>
bool
read_las_points(
  std::istream& stream, ///< input stream.
  OutputIterator output)
{
  return read_las_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, CGAL::parameters::all_default());
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_las_points(), please update your code")
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_las_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map));
}

// deprecated API
template < typename OutputIterator,
           typename PointMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_las_points(), please update your code")
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_las_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


} //namespace CGAL

#endif // CGAL_READ_LAS_POINTS_H
