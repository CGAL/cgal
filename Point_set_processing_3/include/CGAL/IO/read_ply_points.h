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

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>

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

namespace internal {

  class Ply_read_number
  {
  protected:
    std::string m_name;
    std::size_t m_format;
    
  public:
    Ply_read_number (std::string name, std::size_t format)
      : m_name (name), m_format (format) { }
    virtual ~Ply_read_number() { }

    const std::string& name () const { return m_name; }

    virtual double operator() (std::istream& stream) const = 0;

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

  class Ply_read_char : public Ply_read_number
  {
  public:
    Ply_read_char (std::string name, std::size_t format) : Ply_read_number (name, format) { }
    double operator() (std::istream& stream) const
    { return static_cast<double> (this->read<boost::int8_t> (stream)); }
  };
  class Ply_read_uchar : public Ply_read_number
  {
  public:
    Ply_read_uchar (std::string name, std::size_t format) : Ply_read_number (name, format) { }
    double operator() (std::istream& stream) const
    { return static_cast<double> (this->read<boost::uint8_t> (stream)); }
  };
  class Ply_read_short : public Ply_read_number
  {
  public:
    Ply_read_short (std::string name, std::size_t format) : Ply_read_number (name, format) { }
    double operator() (std::istream& stream) const
    { return static_cast<double> (this->read<boost::int16_t> (stream)); }
  };
  class Ply_read_ushort : public Ply_read_number
  {
  public:
    Ply_read_ushort (std::string name, std::size_t format) : Ply_read_number (name, format) { }
    double operator() (std::istream& stream) const
    { return static_cast<double> (this->read<boost::uint16_t> (stream)); }
  };
  class Ply_read_int : public Ply_read_number
  {
  public:
    Ply_read_int (std::string name, std::size_t format) : Ply_read_number (name, format) { }
    double operator() (std::istream& stream) const
    { return static_cast<double> (this->read<boost::int32_t> (stream)); }
  };
  class Ply_read_uint : public Ply_read_number
  {
  public:
    Ply_read_uint (std::string name, std::size_t format) : Ply_read_number (name, format) { }
    double operator() (std::istream& stream) const
    { return static_cast<double> (this->read<boost::uint32_t> (stream)); }
  };
  class Ply_read_float : public Ply_read_number
  {
  public:
    Ply_read_float (std::string name, std::size_t format) : Ply_read_number (name, format) { }
    double operator() (std::istream& stream) const
    { return static_cast<double> (this->read<float> (stream)); }
  };
  class Ply_read_double : public Ply_read_number
  {
  public:
    Ply_read_double (std::string name, std::size_t format) : Ply_read_number (name, format) { }
    double operator() (std::istream& stream) const
    { return this->read<float> (stream); }
  };



} //namespace CGAL
  
//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Reads points (positions + normals, if available) from a .ply
/// stream (ASCII or binary).
/// Potential additional point properties and faces are ignored.
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with  value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `OutputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with value type `Vector_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.
///
/// @return true on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap,
           typename NormalPMap,
           typename Kernel >
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_pmap, ///< property map: value_type of OutputIterator -> Vector_3.
                                 const Kernel& /*kernel*/) ///< geometric traits.
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  // typedef typename value_type_traits<OutputIterator>::type Enriched_point;
  typedef OutputIteratorValueType Enriched_point;

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  if(!stream)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // scan points
  std::size_t pointsCount = 0,  // number of items in file
    pointsRead = 0, // current number of points read
    lineNumber = 0; // current line number
  enum Format { ASCII = 0, BINARY_LITTLE_ENDIAN = 1, BINARY_BIG_ENDIAN = 2};
  Format format = ASCII;
    
  std::string line;
  std::istringstream iss;

  // Check the order of the properties of the point set
  bool reading_properties = false;
  std::vector<internal::Ply_read_number*> readers;
  
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
            std::cerr << "Incorrect file format line " << lineNumber << " of file" << std::endl;
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
            std::cerr << "Unknown file format \"" << format_string << "\" line " << lineNumber << std::endl;
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
              readers.push_back (new internal::Ply_read_char (name, format));
            else if (type == "uchar"  || type == "uint8")
              readers.push_back (new internal::Ply_read_uchar (name, format));
            else if (type == "short"  || type == "int16")
              readers.push_back (new internal::Ply_read_short (name, format));
            else if (type == "ushort" || type == "uint16")
              readers.push_back (new internal::Ply_read_ushort (name, format));
            else if (type == "int"    || type == "int32")
              readers.push_back (new internal::Ply_read_int (name, format));
            else if (type == "uint"   || type == "uint32")
              readers.push_back (new internal::Ply_read_uint (name, format));
            else if (type == "float"  || type == "float32")
              readers.push_back (new internal::Ply_read_float (name, format));
            else if (type == "double" || type == "float64")
              readers.push_back (new internal::Ply_read_double (name, format));
                
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
                pointsCount = number;
                reading_properties = true;
              }
            else
              continue;
          }
            
      }
  }
  
  while (!(stream.eof()) && pointsRead < pointsCount)
    {

      FT x = 0., y = 0., z = 0., nx = 0., ny = 0., nz = 0.;
      for (std::size_t i = 0; i < readers.size (); ++ i)
        {
          FT value = (*readers[i])(stream);
          if (readers[i]->name () == "x") x = value;
          else if (readers[i]->name () == "y") y = value;
          else if (readers[i]->name () == "z") z = value;
          else if (readers[i]->name () == "nx") nx = value;
          else if (readers[i]->name () == "ny") ny = value;
          else if (readers[i]->name () == "nz") nz = value;
        }
      Point point(x,y,z);
      Vector normal(nx,ny,nz);
      Enriched_point pwn;
      
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
      put(point_pmap,  &pwn, point);  // point_pmap[&pwn] = point
      put(normal_pmap, &pwn, normal); // normal_pmap[&pwn] = normal
#else
      put(point_pmap,  pwn, point);  // point_pmap[&pwn] = point
      put(normal_pmap, pwn, normal); // normal_pmap[&pwn] = normal
#endif
      *output++ = pwn;
      pointsRead++;
    }
    // Skip remaining lines

  for (std::size_t i = 0; i < readers.size (); ++ i)
    delete readers[i];
  
  return true;
}

/// @cond SKIP_IN_MANUAL
template < typename OutputIterator,
           typename PointPMap,
           typename NormalPMap,
           typename Kernel >
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_pmap, ///< property map: value_type of OutputIterator -> Vector_3.
                                 const Kernel& kernel) ///< geometric traits.
{
  // just deduce value_type of OutputIterator
  return read_ply_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap,
                                                       normal_pmap,
                                                       kernel);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
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
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_ply_points_and_normals
    <OutputIteratorValueType>(stream,
                              output,
                              point_pmap,
                              normal_pmap,
                              Kernel());
}

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
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
                              make_dereference_property_map(output),
#else
                              make_identity_property_map(OutputIteratorValueType()),
#endif
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
/// \ingroup PkgPointSetProcessing
/// Reads points (position only) from a .ply stream (ASCII or binary).
/// Potential additional point properties (including normals) and faces are ignored.
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with  value_type `Point_3<Kernel>`.
///        It can be omitted if the value type of `OutputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from  the value type of `PointPMap`.
///
/// @return `true` on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap,
           typename Kernel >
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                     const Kernel& kernel) ///< geometric traits.
{
  // Calls read_ply_points_and_normals() with a normal property map = boost::dummy_property_map
  return read_ply_points_and_normals
    <OutputIteratorValueType>(stream,
                              output,
                              point_pmap,
                              boost::dummy_property_map(),
                              kernel);
}

/// @cond SKIP_IN_MANUAL
template < typename OutputIterator,
           typename PointPMap,
           typename Kernel >
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                     const Kernel& kernel) ///< geometric traits.
{
  // just deduce value_type of OutputIterator
  return read_ply_points
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap,
                                                       kernel);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap >
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_ply_points
    <OutputIteratorValueType>(stream,
                              output,
                              point_pmap,
                              Kernel());
}

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
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
                              make_dereference_property_map(output)
#else
                              make_identity_property_map(OutputIteratorValueType())
#endif
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
