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
// SPDX-License-Identifier: GPL-3.0+
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
#include <CGAL/IO/io.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <sstream>
#include <string>

#define TRY_TO_GENERATE_PROPERTY(STD_TYPE, T_TYPE, TYPE)          \
  if (type == STD_TYPE  || type == T_TYPE)                              \
    m_properties->push_back (new PLY_read_typed_number< TYPE > (name, format))

#define TRY_TO_GENERATE_SIZED_LIST_PROPERTY(STD_SIZE_TYPE, T_SIZE_TYPE, SIZE_TYPE, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE) \
  if ((size_type == STD_SIZE_TYPE  || size_type == T_SIZE_TYPE) &&      \
      (index_type == STD_INDEX_TYPE || index_type == T_INDEX_TYPE))     \
    m_properties->push_back (new PLY_read_typed_list_with_typed_size< SIZE_TYPE , INDEX_TYPE > (name, format))

#define TRY_TO_GENERATE_LIST_PROPERTY(STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE) \
  TRY_TO_GENERATE_SIZED_LIST_PROPERTY("uchar", "uint8", boost::uint8_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("ushort", "uint16", boost::uint16_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("uint", "uint32", boost::uint32_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE)


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
    
    void read_ascii (std::istream& stream, float& t) const
    {
      stream >> iformat(t);
    }

    void read_ascii (std::istream& stream, double& t) const
    {
      stream >> iformat(t);
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

  template <typename Type>
  class PLY_read_typed_list : public PLY_read_number
  {
  protected:
    mutable std::vector<Type> m_buffer;
  public:
    PLY_read_typed_list (std::string name, std::size_t format)
      : PLY_read_number (name, format)
    {
    }
    virtual void get (std::istream& stream) const = 0;

    const std::vector<Type>& buffer() const
    {
      return m_buffer;
    }
  };

  template <typename SizeType, typename IndexType>
  class PLY_read_typed_list_with_typed_size
    : public PLY_read_typed_list<IndexType>
  {

  public:
    PLY_read_typed_list_with_typed_size (std::string name, std::size_t format)
      : PLY_read_typed_list<IndexType> (name, format)
    {
    }
    void get (std::istream& stream) const
    {
      std::size_t size = static_cast<std::size_t>(this->template read<SizeType>(stream));
      this->m_buffer.resize (size);
      for (std::size_t i = 0; i < size; ++ i)
        this->m_buffer[i] = this->template read<IndexType> (stream);
    }
  };


  class PLY_reader
  {
    std::vector<PLY_read_number*> m_point_properties;
    std::vector<PLY_read_number*> m_face_properties;
    std::vector<PLY_read_number*>* m_properties;
    std::string m_comments;

  public:
    std::size_t m_nb_points;
    std::size_t m_nb_faces;

    PLY_reader () : m_properties (&m_point_properties), m_nb_points (0), m_nb_faces(0)  { }

    const std::vector<PLY_read_number*>& readers() const { return *m_properties; }
    void read_faces()
    {
      m_properties = &m_face_properties;
    }

    const std::string& comments() const { return m_comments; }

    template <typename Stream>
    bool init (Stream& stream)
    {
      std::size_t lineNumber = 0; // current line number
      enum Format { ASCII = 0, BINARY_LITTLE_ENDIAN = 1, BINARY_BIG_ENDIAN = 2};
      Format format = ASCII;
    
      std::string line;
      std::istringstream iss;

      // Check the order of the properties of the point set
      bool reading_vertices = false, reading_faces = false;
  
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
                  if (!reading_vertices && !reading_faces)
                    continue;

                  std::string type, name;
                  if (!(iss >> type >> name))
                    {
                      std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                      return false;
                    }

                  
                  if (type == "list") // Special case
                  {
                    std::string size_type = name;
                    std::string index_type;
                    name.clear();
                    if (!(iss >> index_type >> name))
                    {
                      std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                      return false;
                    }
                    
                    TRY_TO_GENERATE_LIST_PROPERTY ("char", "int8", boost::int8_t);
                    else TRY_TO_GENERATE_LIST_PROPERTY ("uchar", "uint8", boost::uint8_t);
                    else TRY_TO_GENERATE_LIST_PROPERTY ("short", "int16", boost::int16_t);
                    else TRY_TO_GENERATE_LIST_PROPERTY ("ushort", "uint16", boost::uint16_t);
                    else TRY_TO_GENERATE_LIST_PROPERTY ("int", "int32", boost::int32_t);
                    else TRY_TO_GENERATE_LIST_PROPERTY ("uint", "uint32", boost::uint32_t);
                    else TRY_TO_GENERATE_LIST_PROPERTY ("float", "float32", float);
                    else TRY_TO_GENERATE_LIST_PROPERTY ("double", "float64", double);
                  }
                  else
                  {
                    TRY_TO_GENERATE_PROPERTY ("char", "int8", boost::int8_t);
                    else TRY_TO_GENERATE_PROPERTY ("uchar", "uint8", boost::uint8_t);
                    else TRY_TO_GENERATE_PROPERTY ("short", "int16", boost::int16_t);
                    else TRY_TO_GENERATE_PROPERTY ("ushort", "uint16", boost::uint16_t);
                    else TRY_TO_GENERATE_PROPERTY ("int", "int32", boost::int32_t);
                    else TRY_TO_GENERATE_PROPERTY ("uint", "uint32", boost::uint32_t);
                    else TRY_TO_GENERATE_PROPERTY ("float", "float32", float);
                    else TRY_TO_GENERATE_PROPERTY ("double", "float64", double);
                  }
                  
                  continue;
                }
              else if (keyword == "comment")
                {
                  std::string str = iss.str();
                  if (str.size() > 8)
                    {
                      std::copy (str.begin() + 8, str.end(), std::back_inserter (m_comments));
                      m_comments += "\n";
                    }
                }
              // When end_header is reached, stop loop and begin reading points
              else if (keyword == "end_header")
                break;
              else if (keyword == "element")
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
                      reading_vertices = true;
                      reading_faces = false;
                    }
                
                  else if (type == "face")
                    {
                      m_nb_faces = number;
                      reading_faces = true;
                      reading_vertices = false;
                      m_properties = &m_face_properties;
                    }
                  else
                    {
                      reading_faces = false;
                      reading_vertices = false;
                      continue;
                    }
                }
            
            }
        }
      m_properties = &m_point_properties;
      return true;
    }

    ~PLY_reader ()
    {
      for (std::size_t i = 0; i < m_point_properties.size (); ++ i)
        delete m_point_properties[i];
      m_point_properties.clear();
    }
  
    template <typename Type>
    bool does_tag_exist (const char* tag)
    {
      return does_tag_exist (tag, Type());
    }

    template <typename Type>
    void assign (Type& t, const char* tag)
    {
      for (std::size_t i = 0; i < m_properties->size (); ++ i)
        if ((*m_properties)[i]->name () == tag)
          {
            PLY_read_typed_number<Type>*
              reader = dynamic_cast<PLY_read_typed_number<Type>*>((*m_properties)[i]);
            CGAL_assertion (reader != NULL);
            t = reader->buffer();
            return;
          }
    }

    template <typename Type>
    void assign (std::vector<Type>& t, const char* tag)
    {
      for (std::size_t i = 0; i < m_properties->size (); ++ i)
        if ((*m_properties)[i]->name () == tag)
          {
            PLY_read_typed_list<Type>*
              reader = dynamic_cast<PLY_read_typed_list<Type>*>((*m_properties)[i]);
            CGAL_assertion (reader != NULL);
            t = reader->buffer();
            return;
          }
    }

    template <typename Type>
    bool does_tag_exist (const char* tag, const std::vector<Type>&)
    {
      for (std::size_t i = 0; i < m_properties->size (); ++ i)
        if ((*m_properties)[i]->name () == tag)
          return (dynamic_cast<PLY_read_typed_list<Type>*>((*m_properties)[i]) != NULL);
      return false;
    }
    
    template <typename Type>
    bool does_tag_exist (const char* tag, Type)
    {
      for (std::size_t i = 0; i < m_properties->size (); ++ i)
        if ((*m_properties)[i]->name () == tag)
          return (dynamic_cast<PLY_read_typed_number<Type>*>((*m_properties)[i]) != NULL);
      return false;
    }
    bool does_tag_exist (const char* tag, double)
    {
      for (std::size_t i = 0; i < m_properties->size (); ++ i)
        if ((*m_properties)[i]->name () == tag)
          return (dynamic_cast<PLY_read_typed_number<double>*>((*m_properties)[i]) != NULL
                  || dynamic_cast<PLY_read_typed_number<float>*>((*m_properties)[i]) != NULL);

      return false;
    }
    void assign (double& t, const char* tag)
    {
      for (std::size_t i = 0; i < m_properties->size (); ++ i)
        if ((*m_properties)[i]->name () == tag)
          {
            PLY_read_typed_number<double>*
              reader_double = dynamic_cast<PLY_read_typed_number<double>*>((*m_properties)[i]);
            if (reader_double == NULL)
              {
                PLY_read_typed_number<float>*
                  reader_float = dynamic_cast<PLY_read_typed_number<float>*>((*m_properties)[i]);
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
                           std::tuple<PropertyMap, Constructor, PLY_property<T>...>&& current)
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
                           std::tuple<PropertyMap, Constructor, PLY_property<T>...>&& current,
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


  template <typename OutputValueType, typename PropertyMap, typename T>
  void process_properties (PLY_reader& reader, OutputValueType& new_element,
                           std::pair<PropertyMap, PLY_property<T> >&& current)
  {
    T new_value = T();
    reader.assign (new_value, current.second.name);
    put (current.first, new_element, new_value);
  }

  template <typename OutputValueType, typename PropertyMap, typename T,
            typename NextPropertyBinder, typename ... PropertyMapBinders>
  void process_properties (PLY_reader& reader, OutputValueType& new_element,
                           std::pair<PropertyMap, PLY_property<T> >&& current,
                           NextPropertyBinder&& next,
                           PropertyMapBinders&& ... properties)
  {
    T new_value = T();
    reader.assign (new_value, current.second.name);
    put (current.first, new_element, new_value);
    process_properties (reader, new_element, std::forward<NextPropertyBinder>(next),
                        std::forward<PropertyMapBinders>(properties)...);
  }

  } // namespace PLY
  
} // namespace internal

  
  /// \endcond
  

/*
  \ingroup PkgPointSetProcessingIOPly

  Reads user-selected points properties from a .ply stream (ASCII or
  binary).
  Potential additional point properties and faces are ignored.

  Properties are handled through a variadic list of property
  handlers. A `PropertyHandler` can either be:

  - A `std::pair<PropertyMap, PLY_property<T> >` if the user wants
  to read a %PLY property as a scalar value T (for example, storing
  an `int` %PLY property into an `int` variable).

  - A `std::tuple<PropertyMap, Constructor,
  PLY_property<T>...>` if the user wants to use one or several PLY
  properties to construct a complex object (for example, storing 3
  `uchar` %PLY properties into a %Color object that can for example
  be a `CGAL::cpp11::array<unsigned char, 3>`). In that case, the
  second element of the tuple should be a functor that constructs
  the value type of `PropertyMap` from N objects of types `T`.

  \sa `make_ply_point_reader()`
  \sa `make_ply_normal_reader()`

  \cgalRequiresCPP11

  \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
  It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
  \tparam OutputIterator iterator over output points.
  \tparam PropertyHandler handlers to recover properties.

  \return `true` on success.
*/
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

      internal::PLY::process_properties (reader, new_element, std::forward<PropertyHandler>(properties)...);

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
    (stream, output, std::forward<PropertyHandler>(properties)...);
}
/// \endcond
  
/**
   \ingroup PkgPointSetProcessingIOPly
   Reads points (positions + normals, if available) from a .ply
   stream (ASCII or binary).
   Potential additional point properties and faces are ignored.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.

   \param stream input stream.
   \param output output iterator over points.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `WritablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
     `geom_traits::Vector_3`. If this parameter is omitted, normals in the input stream are
     ignored.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
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
bool read_ply_points(std::istream& stream,
                     OutputIterator output,
#ifdef DOXYGEN_RUNNING
                     const NamedParameters& np)
#else
                     const CGAL_BGL_NP_CLASS& np)
#endif
{
  using boost::choose_param;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;
  
  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::type NormalMap;

  bool has_normals = !(boost::is_same<NormalMap,
                       typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::NoMap>::value);

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  NormalMap normal_map = choose_param(get_param(np, internal_np::normal_map), NormalMap());

  if (has_normals)
    return read_ply_points_with_properties (stream, output,
                                            make_ply_point_reader (point_map),
                                            make_ply_normal_reader (normal_map));
  // else
  return read_ply_points_with_properties (stream, output,
                                          make_ply_point_reader (point_map));
}


/// \cond SKIP_IN_MANUAL  
// variant with default NP
template <typename OutputIteratorValueType,
          typename OutputIterator>
bool
read_ply_points(
  std::istream& stream, ///< input stream.
  OutputIterator output) ///< output iterator over points.
{
  return read_ply_points<OutputIteratorValueType>
    (stream, output, CGAL::parameters::all_default());
}

// variant with default output iterator value type
template <typename OutputIterator,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool
read_ply_points(
  std::istream& stream, ///< input stream.
  OutputIterator output,
  const CGAL_BGL_NP_CLASS& np)
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, np);
}

// variant with default NP and output iterator value type
template <typename OutputIterator>
bool
read_ply_points(
  std::istream& stream, ///< input stream.
  OutputIterator output)
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, CGAL::parameters::all_default());
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointMap,
           typename NormalMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API
template < typename OutputIterator,
           typename PointMap,
           typename NormalMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename NormalMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API
template < typename OutputIterator,
           typename NormalMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API  
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points(), please update your code")
bool
read_ply_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_ply_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map));
}

// deprecated API
template < typename OutputIterator,
           typename PointMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points(), please update your code")
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


} //namespace CGAL

#undef TRY_TO_GENERATE_POINT_PROPERTY
#undef TRY_TO_GENERATE_SIZED_FACE_PROPERTY
#undef TRY_TO_GENERATE_FACE_PROPERTY

#endif // CGAL_READ_PLY_POINTS_H
