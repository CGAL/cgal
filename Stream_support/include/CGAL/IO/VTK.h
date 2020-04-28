// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_VTK_H
#define CGAL_IO_VTK_H

#include <CGAL/IO/VTK/VTK_reader.h>
#include <CGAL/IO/VTK/VTK_writer.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#ifdef CGAL_USE_VTK
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkCell.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
namespace CGAL {
namespace IO {
namespace internal {
//append the content of poly_data to a soup.
template <typename PointRange, typename PolygonRange, typename NamedParameters>
bool vtkPointSet_to_polygon_soup(vtkPointSet* poly_data,
                                 PointRange& points,
                                 PolygonRange& polygons,
                                 const NamedParameters&)
{
  typedef typename boost::range_value<PointRange>::type  Point;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  vtkIdType nb_points = poly_data->GetNumberOfPoints();
  vtkIdType nb_cells = poly_data->GetNumberOfCells();
  polygons.reserve(nb_cells);
  std::size_t initial_number_of_pts = points.size();

  // extract points
  points.reserve(initial_number_of_pts+nb_points);
  for(vtkIdType i=0; i<nb_points; ++i)
  {
    double coords[3];
    poly_data->GetPoint(i, coords);
    points.push_back(Point(coords[0], coords[1], coords[2]));
  }

  // extract cells
  for(vtkIdType i=0; i<nb_cells; ++i)
  {
    int cell_type = poly_data->GetCellType(i);
    if(cell_type != 5
       && cell_type != 7
       && cell_type != 9) // only supported cells are triangles, quads and polygons
      continue;

    vtkCell* cell_ptr = poly_data->GetCell(i);

    vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
    if(nb_vertices < 3)
      return false;

    std::vector<std::size_t> ids(nb_vertices);
    for(vtkIdType k=0; k<nb_vertices; ++k)
    {
      vtkIdType id = cell_ptr->GetPointId(k);
      ids[k] = id+initial_number_of_pts;
    }
    polygons.push_back(ids);
  }

  return true;
}
}//end internal
}//end IO


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read
/*!
 * \ingroup IOstreamFunctions
 *
 * reads the content of `is` into `points` and `polygons`, in the VTPformat.
 *
 * \see \ref IOStreamVTK
 */
template <typename PointRange, typename PolygonRange, typename NamedParameters>
bool read_VTP(const char* fname,
              PointRange& points,
              PolygonRange& polygons,
              const NamedParameters& np)
{
  vtkSmartPointer<vtkPointSet> data;
  vtkSmartPointer<IO::internal::ErrorObserverVtk> obs =
      vtkSmartPointer<IO::internal::ErrorObserverVtk>::New();

  data = vtkPolyData::SafeDownCast(IO::internal::read_vtk_file<vtkXMLPolyDataReader>(fname, obs)->GetOutput());

  return  IO::internal::vtkPointSet_to_polygon_soup(data, points, polygons, np);
}

template <typename PointRange, typename PolygonRange, typename NamedParameters>
bool read_VTP(const std::string&  fname,
              PointRange& points,
              PolygonRange& polygons,
              const NamedParameters& np)
{
  return read_VTP(fname.c_str(), points, polygons, np);
}

template <typename PointRange, typename PolygonRange>
bool read_VTP(const char* fname,
              PointRange& points,
              PolygonRange& polygons)
{
  return read_VTP(fname, points, polygons, parameters::all_default());
}

template <typename PointRange, typename PolygonRange>
bool read_VTP(const std::string&  fname,
              PointRange& points,
              PolygonRange& polygons)
{
  return read_VTP(fname.c_str(), points, polygons, parameters::all_default());
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write



namespace IO {
namespace internal {

// writes the polys appended data at the end of the .vtp file
template <typename PolygonRange>
void write_soup_polys(std::ostream& os,
                      const PolygonRange& polygons,
                      const std::vector<std::size_t>& offsets,
                      const std::vector<unsigned char>& cell_type)
{

  std::vector<std::size_t> connectivity_table;
  std::size_t off = 0;
  std::vector<std::size_t> cumul_offsets(offsets);

  for(size_t i = 0; i < polygons.size(); ++i)
  {
    const auto& poly = polygons[i];
    off+=offsets[i];
    cumul_offsets[i]=off;
    for(const std::size_t& i : poly)
      connectivity_table.push_back(i);
  }

  write_vector<std::size_t>(os, connectivity_table);
  write_vector<std::size_t>(os, cumul_offsets);
  write_vector<unsigned char>(os, cell_type);
}

template <typename PointRange>
void write_soup_points_tag(std::ostream& os,
                           const PointRange& points,
                           bool binary,
                           std::size_t& offset)
{
  typedef typename boost::range_value<PointRange>::type     Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel       Gt;
  typedef typename Gt::FT                                   FT;

  std::string format = binary ? "appended" : "ascii";
  std::string type = (sizeof(FT) == 8) ? "Float64" : "Float32";

  os << "    <Points>\n"
     << "      <DataArray type =\"" << type << "\" NumberOfComponents=\"3\" format=\""
     << format;

  if(binary)
  {
    os << "\" offset=\"" << offset << "\"/>\n";
    offset += 3 * points.size() * sizeof(FT) + sizeof(std::size_t);
    // 3 coords per points + length of the encoded data (size_t)
  }
  else
  {
    os << "\">\n";
    for(const Point& p : points)
      os << p.x() << " " << p.y() << " " << p.z() << " ";
    os << "      </DataArray>\n";
  }
  os << "    </Points>\n";
}

template <typename PolygonRange>
void write_soup_polys_tag(std::ostream& os,
                          const PolygonRange& polygons,
                          bool binary,
                          std::size_t& offset,
                          std::vector<std::size_t>& size_map,
                          std::vector<unsigned char>& cell_type)
{
  std::string formatattribute = binary ? " format=\"appended\"" : " format=\"ascii\"";

  std::string typeattribute;
  switch(sizeof(std::size_t))
  {
  case 8: typeattribute = " type=\"UInt64\""; break;
  case 4: typeattribute = " type=\"UInt32\""; break;
  default: CGAL_error_msg("Unknown size of std::size_t");
  }

  // Write connectivity table
  os << "    <Polys>\n"
     << "      <DataArray Name=\"connectivity\""
     << formatattribute << typeattribute;


  //fill a vector of sizes
  size_map.resize(polygons.size());
  cell_type.resize(polygons.size());
  std::size_t total_size=0;
  for(std::size_t i = 0; i < polygons.size(); ++i)
  {
    size_map[i] = polygons[i].size();
    CGAL_assertion(size_map[i]>=3);
    total_size +=size_map[i];
  }

  // if binary output, just write the xml tag
  if(binary)
  {
    os << " offset=\"" << offset << "\"/>\n";
    offset += (total_size + 1) * sizeof(std::size_t);
    // n indices (size_t) per triangle + length of the encoded data (size_t)
  }
  else
  {
    os << ">\n";

    for(const auto& poly : polygons)
    {
      for(const std::size_t& id : poly)
        os << id << " ";
    }
    os << "      </DataArray>\n";
  }

  // Write offsets
  os << "      <DataArray Name=\"offsets\""
     << formatattribute << typeattribute;

  if(binary) {  // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (polygons.size() + 1) * sizeof(std::size_t);
    // 1 offset (size_t) per triangle + length of the encoded data (size_t)
  }
  else
  {
    os << ">\n";
    std::size_t polys_offset = 0;

    for(std::size_t i = 0; i < polygons.size(); ++i)
    {
      polys_offset += size_map[i];
      os << polys_offset << " ";
    }
    os << "      </DataArray>\n";
  }

  // Write cell type (triangle == 5)
  os << "      <DataArray Name=\"types\""
     << formatattribute << " type=\"UInt8\"";

  if(binary)
  {
    os << " offset=\"" << offset << "\"/>\n";
    offset += polygons.size() + sizeof(std::size_t);
    // 1 unsigned char per cell + length of the encoded data (size_t)
  }
  else
  {
    os << ">\n";
    for(std::size_t i = 0; i< polygons.size(); ++i)
    {
      switch(size_map[i]){
      case 3:
        cell_type[i]=5;
        os << "5 ";//triangle
        break;
      case 4:
        cell_type[i]=7;
        os << "7 ";//quad
        break;
      default:
        cell_type[i]=9;
        os << "9 ";//polygon
        break;
      }
    }
    os << "      </DataArray>\n";
  }
  os << "    </Polys>\n";
}

// writes the points appended data at the end of the .vtp file
template <typename PointRange>
void write_soup_polys_points(std::ostream& os,
                             const PointRange& points)
{
  typedef typename boost::range_value<PointRange>::type     Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel       Gt;
  typedef typename Gt::FT                                   FT;

  std::vector<FT> coordinates;

  for(const Point& p : points)
  {
    coordinates.push_back(p.x());
    coordinates.push_back(p.y());
    coordinates.push_back(p.z());
  }

  write_vector<FT>(os, coordinates);
}

}//end internal
}//end IO

/*!
 * \ingroup IOstreamFunctions
 *
 * writes the content of `points` and `polygons` in `out`, in the VTP format.
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{use_binary_mode} a Boolean indicating if the
 *    data should be written in binary (`true`, the default) or in ASCII (`false`).
 *    \cgalParamEnd
 * \cgalNamedParamsEnd
 * \see \ref IOStreamVTK
 */
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_VTP(std::ostream& os,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  os << "<?xml version=\"1.0\"?>\n"
     << "<VTKFile type=\"PolyData\" version=\"0.1\"";

#ifdef CGAL_LITTLE_ENDIAN
  os << " byte_order=\"LittleEndian\"";
#else // CGAL_BIG_ENDIAN
  os << " byte_order=\"BigEndian\"";
#endif

  switch(sizeof(std::size_t))
  {
  case 4: os << " header_type=\"UInt32\""; break;
  case 8: os << " header_type=\"UInt64\""; break;
  default: CGAL_error_msg("Unknown size of std::size_t");
  }

  os << ">\n"
     << "  <PolyData>" << "\n";
  os << "  <Piece NumberOfPoints=\"" << points.size()
     << "\" NumberOfPolys=\"" << polygons.size() << "\">\n";

  std::size_t offset = 0;
  const bool binary = choose_parameter(get_parameter(np, internal_np::use_binary_mode), true);
  std::vector<std::size_t> size_map;
  std::vector<unsigned char> cell_type;
  IO::internal::write_soup_points_tag(os, points, binary, offset);
  IO::internal::write_soup_polys_tag(os, polygons, binary, offset, size_map, cell_type);

  os << "   </Piece>\n"
     << "  </PolyData>\n";
  if(binary)
  {
    os << "<AppendedData encoding=\"raw\">\n_";
    IO::internal::write_soup_polys_points(os, points);
    IO::internal::write_soup_polys(os, polygons,size_map, cell_type);
  }
  os << "</VTKFile>\n";
}

template <typename PointRange, typename PolygonRange>
bool write_VTP(std::ostream& os,
               const PointRange& points,
               const PolygonRange& polygons)
{
  return write_VTP(os, points, polygons, parameters::all_default());
}

/*!
 * \ingroup IOstreamFunctions
 *
 * writes the content of `points` and `polygons` in a file named `fname`, in the VTP format.
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{use_binary_mode} a Boolean indicating if the
 *    data should be written in binary (`true`, the default) or in ASCII (`false`).
 *    \cgalParamEnd
 * \cgalNamedParamsEnd
 * \see \ref IOStreamVTK
 */
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_VTP(const char* fname,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream out(fname);
  return write_VTP(out, points, polygons, np);
}

template <typename PointRange, typename PolygonRange>
bool write_VTP(const char* fname,
               const PointRange& points,
               const PolygonRange& polygons)
{
  return write_VTP(fname, points, polygons, parameters::all_default());
}

template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_VTP(const std::string& fname,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np)
{
  return write_VTP(fname.c_str(), points, polygons, np);
}

template <typename PointRange, typename PolygonRange>
bool write_VTP(const std::string& fname,
               const PointRange& points,
               const PolygonRange& polygons)
{
  return write_VTP(fname, points, polygons, parameters::all_default());
}

}//end CGAL
#endif //CGAL_USE_VTK
#endif // CGAL_IO_VTK_H
