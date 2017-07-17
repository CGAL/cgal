// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Simon Giraudot

#ifndef CGAL_READ_PLY_POINTS_H
#define CGAL_READ_PLY_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/config.h>
#if defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) || defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES)
#error CGAL PLY reader requires a C++11 compiler
#endif

#include <tuple>

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <sstream>
#include <string>


namespace CGAL {


// PLY types:
// name        type        number of bytes
// ---------------------------------------
// char       character                 1
// uchar      unsigned character        1
// short      short integer             2
// ushort     unsigned short integer    2
// int        integer                   4
// uint       unsigned integer          4
// float      single-precision float    4
// double     double-precision float    8

  /**
     \ingroup PkgPointSetProcessingIOPly
     
     Class used to identify a %PLY property as a type and a name.

     \sa `read_ply_points_with_properties()`
  */
  template <typename T>
  struct PLY_property
  {
    typedef T type;
    const char* name;
    PLY_property (const char* name) : name (name) { }
  };

  /**
     \ingroup PkgPointSetProcessingIOPly
     
     Generates a %PLY property handler to read 3D points. Points are
     constructed from the input using 3 %PLY properties of type
     `double` and named `x`, `y` and `z`.

     \sa `read_ply_points_with_properties()`

     \tparam PointMap the property map used to store points.
  */
  template <typename PointMap>
  std::tuple<PointMap,
             typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
             PLY_property<double>, PLY_property<double>, PLY_property<double> >
   make_ply_point_reader(PointMap point_map)
  {
    return std::make_tuple (point_map, typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3(),
                            PLY_property<double>("x"), PLY_property<double>("y"), PLY_property<double>("z"));
  }

  /**
     \ingroup PkgPointSetProcessingIOPly
     
     Generates a %PLY property handler to read 3D normal
     vectors. Vectors are constructed from the input using 3 PLY
     properties of type `double` and named `nx`, `ny` and `nz`.

     \sa `read_ply_points_with_properties()`

     \tparam VectorMap the property map used to store vectors.
  */
  template <typename VectorMap>
  std::tuple<VectorMap,
             typename Kernel_traits<typename VectorMap::value_type>::Kernel::Construct_vector_3,
             PLY_property<double>, PLY_property<double>, PLY_property<double> >
  make_ply_normal_reader(VectorMap normal_map)
  {
    return std::make_tuple (normal_map, typename Kernel_traits<typename VectorMap::value_type>::Kernel::Construct_vector_3(),
                            PLY_property<double>("nx"), PLY_property<double>("ny"), PLY_property<double>("nz"));
  }

  /// \cond SKIP_IN_MANUAL

namespace internal {

  namespace PLY {
  
  class PLY_read_number
  {
  protected:
    std::string m_name;
    std::size_t m_format;
    
  public:
    PLY_read_number (std::string name, std::size_t format)
      : m_name (name), m_format (format) { }
    virtual ~PLY_read_number() { }

    const std::string& name () const { return m_name; }

    virtual void get (std::istream& stream) const = 0;

    // The two following functions prevent the stream to only extract
    // ONE character (= what the types char imply) by requiring
    // explicitely an integer object when reading the stream
    void read_ascii (std::istream& stream, boost::int8_t& c) const
    {
      short s;
      stream >> s;
      c = static_cast<char>(s);
    }
    void read_ascii (std::istream& stream, boost::uint8_t& c) const
    {
      unsigned short s;
      stream >> s;
      c = static_cast<unsigned char>(s);
    }

    // Default template when Type is not a char type
    template <typename Type>
    void read_ascii (std::istream& stream, Type& t) const
    {
      stream >> t;
    }
    
    template <typename Type>
    Type read (std::istream& stream) const
    {
      if (m_format == 0) // Ascii
        {
          Type t;
          read_ascii (stream, t);
          return t;
        }
      else // Binary (2 = little endian)
        {
          union
          {
            char uChar[sizeof (Type)];
            Type type;
          } buffer;
          
          std::size_t size = sizeof (Type);

          stream.read(buffer.uChar, size);
      
          if (m_format == 2) // Big endian
            {
              for (std::size_t i = 0; i < size / 2; ++ i)
                {
                  unsigned char tmp = buffer.uChar[i];
                  buffer.uChar[i] = buffer.uChar[size - 1 - i];
                  buffer.uChar[size - 1 - i] = tmp;
                }
            }
          return buffer.type;
        }
      return Type();
    }
  };

  template <typename Type>
  class PLY_read_typed_number : public PLY_read_number
  {
    mutable Type m_buffer;
  public:
    PLY_read_typed_number (std::string name, std::size_t format)
      : PLY_read_number (name, format)
    {
    }
    void get (std::istream& stream) const
    {
      m_buffer = (this->read<Type> (stream));
    }
    const Type& buffer() const
    {
      return m_buffer;
    }
  };

  

  class PLY_reader
  {

    std::vector<PLY_read_number*> m_readers;


  public:
    std::size_t m_nb_points;

    PLY_reader () : m_nb_points (0) { }

    const std::vector<PLY_read_number*>& readers() const { return m_readers; }

    template <typename Stream>
    bool init (Stream& stream)
    {
      std::size_t lineNumber = 0; // current line number
      enum Format { ASCII = 0, BINARY_LITTLE_ENDIAN = 1, BINARY_BIG_ENDIAN = 2};
      Format format = ASCII;
    
      std::string line;
      std::istringstream iss;

      // Check the order of the properties of the point set
      bool reading_properties = false;
  
      while (getline (stream,line))
        {
          iss.clear();
          iss.str (line);
          ++ lineNumber;

          // Reads file signature on first line
          if (lineNumber == 1)
            {
              std::string signature;
              if (!(iss >> signature) || (signature != "ply"))
                {
                  // if wrong file format
                  std::cerr << "Error: incorrect file format line " << lineNumber << " of file" << std::endl;
                  return false;
                }
            }

          // Reads format on 2nd line
          else if (lineNumber == 2)
            {
              std::string tag, format_string, version;
              if ( !(iss >> tag >> format_string >> version) )
                {
                  std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                  return false;
                }
              if (format_string == "ascii") format = ASCII;
              else if (format_string == "binary_little_endian") format = BINARY_LITTLE_ENDIAN;
              else if (format_string == "binary_big_endian") format = BINARY_BIG_ENDIAN;
              else
                {
                  std::cerr << "Error: unknown file format \"" << format_string << "\" line " << lineNumber << std::endl;
                  return false;
                }
            }

          // Comments and vertex properties
          else
            {
              std::string keyword;
              if (!(iss >> keyword))
                {
                  std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                  return false;
                }

              if (keyword == "property")
                {
                  if (!reading_properties)
                    continue;

                  std::string type, name;
                  if (!(iss >> type >> name))
                    {
                      std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                      return false;
                    }

                  if (     type == "char"   || type == "int8")
                    m_readers.push_back (new PLY_read_typed_number<boost::int8_t> (name, format));
                  else if (type == "uchar"  || type == "uint8")
                    m_readers.push_back (new PLY_read_typed_number<boost::uint8_t> (name, format));
                  else if (type == "short"  || type == "int16")
                    m_readers.push_back (new PLY_read_typed_number<boost::int16_t> (name, format));
                  else if (type == "ushort" || type == "uint16")
                    m_readers.push_back (new PLY_read_typed_number<boost::uint16_t> (name, format));
                  else if (type == "int"    || type == "int32")
                    m_readers.push_back (new PLY_read_typed_number<boost::int32_t> (name, format));
                  else if (type == "uint"   || type == "uint32")
                    m_readers.push_back (new PLY_read_typed_number<boost::uint32_t> (name, format));
                  else if (type == "float"  || type == "float32")
                    m_readers.push_back (new PLY_read_typed_number<float> (name, format));
                  else if (type == "double" || type == "float64")
                    m_readers.push_back (new PLY_read_typed_number<double> (name, format));
                
                  continue;
                }
              else
                reading_properties = false;
            
              // ignore comments and properties (if not in element
              // vertex - cf below - properties are useless in our case)
              if (keyword == "comment" || keyword == "property")
                continue;

              // When end_header is reached, stop loop and begin reading points
              if (keyword == "end_header")
                break;
            
              if (keyword == "element")
                {
                  std::string type;
                  std::size_t number;
                  if (!(iss >> type >> number))
                    {
                      std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                      return false;
                    }
                
                  if (type == "vertex")
                    {
                      m_nb_points = number;
                      reading_properties = true;
                    }
                  else
                    continue;
                }
            
            }
        }
      return true;
    }

    ~PLY_reader ()
    {
      for (std::size_t i = 0; i < m_readers.size (); ++ i)
        delete m_readers[i];
      m_readers.clear();
    }
  
    template <typename Type>
    bool does_tag_exist (const char* tag)
    {
      return does_tag_exist (tag, Type());
    }

    template <typename Type>
    void assign (Type& t, const char* tag)
    {
      for (std::size_t i = 0; i < m_readers.size (); ++ i)
        if (m_readers[i]->name () == tag)
          {
            PLY_read_typed_number<Type>*
              reader = dynamic_cast<PLY_read_typed_number<Type>*>(m_readers[i]);
            CGAL_assertion (reader != NULL);
            t = reader->buffer();
            return;
          }
    }

    template <typename Type>
    bool does_tag_exist (const char* tag, Type)
    {
      for (std::size_t i = 0; i < m_readers.size (); ++ i)
        if (m_readers[i]->name () == tag)
          return (dynamic_cast<PLY_read_typed_number<Type>*>(m_readers[i]) != NULL);
      return false;
    }
    bool does_tag_exist (const char* tag, double)
    {
      for (std::size_t i = 0; i < m_readers.size (); ++ i)
        if (m_readers[i]->name () == tag)
          return (dynamic_cast<PLY_read_typed_number<double>*>(m_readers[i]) != NULL
                  || dynamic_cast<PLY_read_typed_number<float>*>(m_readers[i]) != NULL);

      return false;
    }
    void assign (double& t, const char* tag)
    {
      for (std::size_t i = 0; i < m_readers.size (); ++ i)
        if (m_readers[i]->name () == tag)
          {
            PLY_read_typed_number<double>*
              reader_double = dynamic_cast<PLY_read_typed_number<double>*>(m_readers[i]);
            if (reader_double == NULL)
              {
                PLY_read_typed_number<float>*
                  reader_float = dynamic_cast<PLY_read_typed_number<float>*>(m_readers[i]);
                CGAL_assertion (reader_float != NULL);
                t = reader_float->buffer();
              }
            else
              t = reader_double->buffer();
          
            return;
          }
    }
  
  };

  template <class Reader, class T>
  void get_value(Reader& r, T& v, PLY_property<T>& wrapper)
  {
    return r.assign(v, wrapper.name);
  }

  
  template <std::size_t N>
  struct Filler
  {
    template <class Reader, class Value_tuple, class PLY_property_tuple>
    static void fill(Reader& r, Value_tuple& values, PLY_property_tuple wrappers)
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
    template <class Reader, class Value_tuple, class PLY_property_tuple>
    static void fill(Reader& r, Value_tuple& values, PLY_property_tuple wrappers)
    {
      get_value(r, std::get<0>(values), std::get<2>(wrappers));
    }
  };
  
  template <typename OutputValueType,
            typename PropertyMap,
            typename Constructor,
            typename ... T>
  void process_properties (PLY_reader& reader, OutputValueType& new_element,
                           std::tuple<PropertyMap, Constructor, PLY_property<T>...>& current)
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
            typename NextPropertyBinder,
            typename ... PropertyMapBinders>
  void process_properties (PLY_reader& reader, OutputValueType& new_element,
                           std::tuple<PropertyMap, Constructor, PLY_property<T>...>& current,
                           NextPropertyBinder& next,
                           PropertyMapBinders&& ... properties)
  {
    typedef typename PropertyMap::value_type PmapValueType;
    std::tuple<T...> values;
    Filler<sizeof...(T)-1>::fill(reader, values, current);
    PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
    put (std::get<0>(current), new_element, new_value);
  
    process_properties (reader, new_element, next, properties...);
  }


  template <typename OutputValueType, typename PropertyMap, typename T>
  void process_properties (PLY_reader& reader, OutputValueType& new_element,
                           std::pair<PropertyMap, PLY_property<T> >& current)
  {
    T new_value = T();
    reader.assign (new_value, current.second.name);
    put (current.first, new_element, new_value);
  }

  template <typename OutputValueType, typename PropertyMap, typename T,
            typename NextPropertyBinder, typename ... PropertyMapBinders>
  void process_properties (PLY_reader& reader, OutputValueType& new_element,
                           std::pair<PropertyMap, PLY_property<T> >& current,
                           NextPropertyBinder& next,
                           PropertyMapBinders&& ... properties)
  {
    T new_value = T();
    reader.assign (new_value, current.second.name);
    put (current.first, new_element, new_value);
    process_properties (reader, new_element, next, properties...);
  }

  } // namespace PLY
  
} // namespace internal

  
  /// \endcond
  

//===================================================================================
/// \ingroup PkgPointSetProcessingIOPly

/// Reads user-selected points properties from a .ply stream (ASCII or
/// binary).
/// Potential additional point properties and faces are ignored.
///
/// Properties are handled through a variadic list of property
/// handlers. A `PropertyHandler` can either be:
///
///  - A `std::pair<PropertyMap, PLY_property<T> >` if the user wants
///  to read a %PLY property as a scalar value T (for example, storing
///  an `int` %PLY property into an `int` variable).
///
///  - A `std::tuple<PropertyMap, Constructor,
///  PLY_property<T>...>` if the user wants to use one or several PLY
///  properties to construct a complex object (for example, storing 3
///  `uchar` %PLY properties into a %Color object that can for example
///  be a `CGAL::cpp11::array<unsigned char, 3>`). In that case, the
///  second element of the tuple should be a functor that constructs
///  the value type of `PropertyMap` from N objects of types `T`.
///
/// @sa `make_ply_point_reader()`
/// @sa `make_ply_normal_reader()`
///
/// @cgalRequiresCPP11
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PropertyHandler handlers to recover properties.
///
/// @return `true` on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename ... PropertyHandler>
bool read_ply_points_with_properties (std::istream& stream,
                                      OutputIterator output,
                                      PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<OutputIterator>::type OutputValueType;
    
  if(!stream)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

  internal::PLY::PLY_reader reader;
  
  if (!(reader.init (stream)))
    return false;
  
  std::size_t points_read = 0;
  
  while (!(stream.eof()) && points_read < reader.m_nb_points)
    {
      for (std::size_t i = 0; i < reader.readers().size (); ++ i)
        reader.readers()[i]->get (stream);

      OutputValueType new_element;

      internal::PLY::process_properties (reader, new_element, properties...);

      *(output ++) = new_element;
      
      ++ points_read;
    }
  // Skip remaining lines

  return (points_read == reader.m_nb_points);
}

/// \cond SKIP_IN_MANUAL
template <typename OutputIterator,
          typename ... PropertyHandler>
bool read_ply_points_with_properties (std::istream& stream,
                                      OutputIterator output,
                                      PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<OutputIterator>::type OutputValueType;

  return read_ply_points_with_properties<OutputValueType>
    (stream, output, properties...);
}
/// \endcond
  
//===================================================================================
/// \ingroup PkgPointSetProcessingIOPly
/// Reads points (positions + normals, if available) from a .ply
/// stream (ASCII or binary).
/// Potential additional point properties and faces are ignored.
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with  value type `CGAL::Point_3`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with value type `CGAL::Vector_3`.
///
/// @return `true` on success.
///
/// @cgalRequiresCPP11

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap,
           typename NormalPMap >
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{

  return read_ply_points_with_properties (stream, output,
                              make_ply_point_reader (point_pmap),
                              make_ply_normal_reader (normal_pmap));
}

/// @cond SKIP_IN_MANUAL
template < typename OutputIterator,
           typename PointPMap,
           typename NormalPMap >
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  // just deduce value_type of OutputIterator
  return read_ply_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap,
                                                       normal_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond


/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename NormalPMap >
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points_and_normals
    <OutputIteratorValueType>(stream,
                              output,
                              make_identity_property_map(OutputIteratorValueType()),
                              normal_pmap);
}

template < typename OutputIterator,
           typename NormalPMap >
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  // just deduce value_type of OutputIterator
  return read_ply_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       normal_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond


//===================================================================================
/// \ingroup PkgPointSetProcessingIOPly
/// Reads points (position only) from a .ply stream (ASCII or binary).
/// Potential additional point properties (including normals) and faces are ignored.
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with  value_type `CGAL::Point_3`.
///        It can be omitted if the value type of `OutputIterator` is convertible to `CGAL::Point_3`.
///
/// @return `true` on success.
///
/// @cgalRequiresCPP11

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap >
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_ply_points_with_properties (stream, output,
                              make_ply_point_reader (point_pmap));
}

/// @cond SKIP_IN_MANUAL
template < typename OutputIterator,
           typename PointPMap >
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  // just deduce value_type of OutputIterator
  return read_ply_points
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond


/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator >
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output) ///< output iterator over points.
{
  return read_ply_points
    <OutputIteratorValueType>(stream,
                              output,
                              make_identity_property_map(OutputIteratorValueType())
                              );
}

template < typename OutputIterator>
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output) ///< output iterator over points.
{
  // just deduce value_type of OutputIterator
  return read_ply_points
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output);
}
//-----------------------------------------------------------------------------------

/// @endcond


} //namespace CGAL

#endif // CGAL_READ_PLY_POINTS_H
