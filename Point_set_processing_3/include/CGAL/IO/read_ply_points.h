#ifndef CGAL_READ_PLY_POINTS_H
#define CGAL_READ_PLY_POINTS_H

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>

#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <iostream>
#include <sstream>
#include <string>


namespace CGAL {

//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Reads points (positions + normals, if available) from a .ply stream.
/// Faces are ignored.
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

  if(!stream)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // scan points
  std::size_t pointsCount = 0,  // number of items in file
    pointsRead = 0, // current number of points read
    lineNumber = 0; // current line number
  enum Format { ASCII, BINARY_LITTLE_ENDIAN, BINARY_BIG_ENDIAN };
  Format format;
    
  std::string line;
  std::istringstream iss;

  bool in_header = true;
  
  while (getline (stream,line))
  {
    iss.clear();
    iss.str (line);
    ++ lineNumber;

    if (in_header)
      {
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

            // ignore comments and properties (if not in element
            // vertex - cf below - properties are useless in our case)
            if (keyword == "comment" || keyword == "property")
              continue;

            if (keyword == "end_header")
              {
                in_header = false;
                continue;
              }
            
            if (keyword == "element")
              {
                std::string type;
                std::size_t number;
                if (!(iss >> type >> number))
                  {
                    std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                    return false;
                  }
                
                if (type == "face")
                  continue;

                if (type == "vertex")
                  pointsCount = number;
                else
                  {
                    std::cerr << "Unknown element \"" << type << "\" line " << lineNumber << std::endl;
                    return false;
                  }

                while (getline (stream,line))
                  {
                    iss.clear();
                    iss.str (line);
                    ++ lineNumber;
                    std::string property, ftype, name;
                    if ( !(iss >> property >> ftype >> name) || property != "property")
                      {
                        std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                        return false;
                      }
                    
                  }
              }
            
          }
      }


    // Reads 3D points on next lines
    else if (pointsRead < pointsCount)
    {
      break;
//       // Reads position + normal...
//       double x,y,z;
//       double nx,ny,nz;
//       if (iss >> iformat(x) >> iformat(y) >> iformat(z))
//       {
//         Point point(x,y,z);
//         Vector normal = CGAL::NULL_VECTOR;
//         // ... + normal...
//         if (iss >> iformat(nx))
//         {
//           // In case we could read one number, we expect that there are two more
//           if(iss  >> iformat(ny) >> iformat(nz)){
//             normal = Vector(nx,ny,nz);
//           } else {
//             std::cerr << "Error line " << lineNumber << " of file" << std::endl;
//             return false;
//           }
//         }
//         Enriched_point pwn;
// #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
//         put(point_pmap,  &pwn, point);  // point_pmap[&pwn] = point
//         put(normal_pmap, &pwn, normal); // normal_pmap[&pwn] = normal
// #else
//         put(point_pmap,  pwn, point);  // point_pmap[&pwn] = point
//         put(normal_pmap, pwn, normal); // normal_pmap[&pwn] = normal
// #endif
//         *output++ = pwn;
//         pointsRead++;
//       }
      // ...or skip comment line
    }
    // Skip remaining lines
  }


  
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
/// Reads points (position only) from a .ply stream.
/// If the position is followed by the nx ny nz normal, then the normal will be ignored.
/// Faces are ignored.
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
