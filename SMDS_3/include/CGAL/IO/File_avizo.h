// Copyright (c) 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_IO_FILE_AVIZO_H
#define CGAL_IO_FILE_AVIZO_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>

#include <iostream>
#include <string>
#include <fstream>
#include <array>
#include <vector>
#include <cstring>

#include <boost/algorithm/string/case_conv.hpp>

#include <CGAL/utility.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

// ---------------------------------
// Output C3t3 to Avizo tetra file format
// ---------------------------------

namespace IO {

namespace internal {

  //a material is composed of an material Id and a material name.
  typedef std::pair<int, std::string> material;

  struct MaterialData
  {
    material innerRegion;
    material outerRegion;
  };

  template <typename MaterialIdType>
  bool get_material_metadata(std::istream& input,
                             std::string& line,
                             material& _material,
                             MaterialIdType material_id)
  {
    std::istringstream iss;
    iss.str(line);

    iss >> _material.second;//name

    while (std::getline(input, line))
    {
      std::string prop; //property
      iss.clear();
      iss.str(line);
      iss >> prop;

      if (prop.compare("Id") == 0)
      {
        int tmp_id;
        iss >> tmp_id;
        _material.first = material_id;

        std::string mat = _material.second;
        boost::algorithm::to_lower(mat);
        if ((0 == mat.compare("exterior"))
          && _material.first != 0)
        {
          std::cerr << "Exterior should have index 0. ";
          std::cerr << "In this file it has index " << _material.first << "." << std::endl;
          std::cerr << "Reader failed, because Meshing will fail to terminate." << std::endl;
          return false;
        }
      }
      else if (prop.compare("}") == 0)
        return true; //end of this material
    }
    return false;
  }

  bool line_starts_with(const std::string& line, const char* cstr)
  {
    const std::size_t fnws = line.find_first_not_of(" \t");
    if (fnws != std::string::npos)
      return (line.compare(fnws, strlen(cstr), cstr) == 0);
    return false;
  }

  template<typename IdType>
  bool treat_surf_materials(std::istream& input,
                            std::vector<material>& materials,
                            IdType& material_id)
  {
    std::string line;
    while (std::getline(input, line))
    {
      if (line_starts_with(line, "}"))
        break;
      else
      {
        material _material;
        if (!get_material_metadata(input, line, _material, material_id++))
          return false;
        materials.push_back(_material);
      }
    }
    return true;
  }

}

/**
 * @ingroup PkgSMDS3ExportFunctions
 * @brief exports a mesh complex to the Avizo (`.am`) file format
 * @tparam C3T3 a class model of `MeshComplex_3InTriangulation_3`
 * @param os the output stream
 * @param c3t3 the mesh complex
 * \see \ref IOStreamAvizo
 */
template <class C3T3>
void
output_to_avizo(std::ostream& os,
                const C3T3& c3t3)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Weighted_point Weighted_point;

  const Tr& tr = c3t3.triangulation();

  CGAL::Unique_hash_map<Vertex_handle, std::size_t> V;

  //-------------------------------------------------------
  // nodes
  //-------------------------------------------------------

  os << std::setprecision(17);
  os << " # Avizo 3D ASCII 2.0\n\n";
  os << "nNodes " << tr.number_of_vertices() << std::endl;
  os << "nTetrahedra " << c3t3.number_of_cells_in_complex() << std::endl;
  os <<  "Parameters {\n"
    "    Materials {\n"
    "        Material3 {\n"
    "            Id 1,\n"
    "            Color 0 0.835294 0.164706\n"
    "        }\n"
    "        Material4 {\n"
    "            Id 2,\n"
    "            Color 0.862745 0.0901961 0.0901961\n"
    "        }\n"
    "        Material5 {\n"
    "            Id 3,\n"
    "            Color 0.94902 0.847059 0.0901961\n"
    "        }\n"
    "        Material6 {\n"
    "            Id 4,\n"
    "            Color 0.8 0.16 0.698646\n"
    "        }\n"
    "        Material7 {\n"
    "            Id 5,\n"
    "            Color 0.494118 0.494118 1\n"
    "        }\n"
    "        Material8 {\n"
    "            Id 6,\n"
    "            Color 0.227451 0.227451 0.968627\n"
    "        }\n"
    "        Material9 {\n"
    "            Id 7,\n"
    "            Color 0.666667 0.666667 0.666667\n"
    "        }\n"
    "    }\n"
    "}\n"
    "Nodes { float[3] Coordinates } @1\n"
    "Tetrahedra { int[4] Nodes } @2\n"
    "TetrahedronData { byte Materials } @3\n"
    "\n"
    "# Data section follows\n"
    "@1\n";

  std::size_t vert_counter = 0;
  for(Finite_vertices_iterator
        vit = tr.finite_vertices_begin(),
        end = tr.finite_vertices_end();
      vit != end; ++vit)
  {
    const Weighted_point& p = tr.point(vit);
    const double x = CGAL::to_double(p.x());
    const double y = CGAL::to_double(p.y());
    const double z = CGAL::to_double(p.z());

    V[vit] = ++vert_counter;

    os << x << " " << y << " " << z << "\n";
  }



  //-------------------------------------------------------
  // Elements
  //-------------------------------------------------------

  os << "\n@2\n";
  for (Cell_iterator
         cit = c3t3.cells_in_complex_begin(),
         end = c3t3.cells_in_complex_end();
       cit != end; ++cit)
  {
    const Cell_handle ch = cit;
    os << V[ch->vertex(0)] ;
    os << " " << V[ch->vertex(1)] ;
    os << " " << V[ch->vertex(2)] ;
    os << " " << V[ch->vertex(3)]  << "\n";
  }

  os << "\n@3\n";
  for (Cell_iterator
         cit = c3t3.cells_in_complex_begin(),
         end = c3t3.cells_in_complex_end();
       cit != end; ++cit)
  {
    os << cit->subdomain_index() << "\n";
  }

} // end output_to_avizo(...)


#ifndef CGAL_NO_DEPRECATED_CODE
using IO::output_to_avizo;
#endif

// ---------------------------------
// Read C3t3 from Avizo file format
// ---------------------------------

namespace internal {
  template<typename Point_3>
  void read_points(std::istream& input,
                   std::vector<Point_3>& points,
                   const int& nb_points,
                   const bool binary)
  {
    points.clear();
    points.reserve(nb_points);

    if (binary)
    {
      std::vector<float> data(3 * nb_points);
      input.read(reinterpret_cast<char*>(&data[0]),
                 3 * nb_points * sizeof(float));
      for (int i = 0; i < nb_points; ++i)
        points.push_back(Point_3(data[3 * i], data[3 * i + 1], data[3 * i + 2]));
    }
    else
    {
      for (int i = 0; i < nb_points; ++i)
      {
        float x, y, z;
        input >> x >> y >> z;
        points.push_back(Point_3(x, y, z));
      }
    }
  }

  void read_tetrahedra(std::istream& input,
                       std::vector<std::array<int, 4> >& tetrahedra,
                       const int& nb_tets,
                       const bool binary)
  {
    tetrahedra.clear();
    tetrahedra.reserve(nb_tets);

    if (binary)
    {
      std::vector<int> data(4 * nb_tets);
      input.read(reinterpret_cast<char*>(&data[0]),
                 4 * nb_tets * sizeof(int));
      for (int i = 0; i < nb_tets; ++i)
      {
        tetrahedra.push_back({ data[4 * i] - 1,
                               data[4 * i + 1] - 1,
                               data[4 * i + 2] - 1,
                               data[4 * i + 3] - 1 });
      }
    }
    else
    {
      for (int n = 0; n < nb_tets; ++n)
      {
        int i, j, k, l;
        input >> i >> j >> k >> l;
        tetrahedra.push_back({ {i - 1,
                                j - 1,
                                k - 1,
                                l - 1} });
      }
    }
  }

  template<typename Material_type>
  void read_materials(std::istream& input,
                      std::vector<Material_type>& materials,
                      const int& nb_tets,
                      const bool binary)
  {
    materials.clear();
    materials.reserve(nb_tets);

    if (binary)
    {
      std::vector<Material_type> data(nb_tets);
      input.read(reinterpret_cast<char*>(&data[0]),
                 nb_tets * sizeof(Material_type));
      for (int i = 0; i < nb_tets; ++i)
        materials.push_back(data[i]);
    }
    else
    {
      for (int i = 0; i < nb_tets; ++i)
      {
        Material_type m;
        input >> m;
        materials.push_back(m);
      }
    }
  }


  struct Data_def
  {
    std::string name;
    std::string type;
    std::string description;
    std::string at_label;
  };

  void go_to_at_label(std::istream& input,
                      std::string& line,
                      const char* at_label)// "@1"
  {
    if (line_starts_with(line, at_label))
      return;

    while (std::getline(input, line))
    {
      if (line_starts_with(line, at_label))
        return;
    }
  }

  bool is_avizo_tetra_format(std::istream& in, const char* binary_or_ascii)
  {
    std::string format(binary_or_ascii);
    boost::algorithm::to_lower(format);

    std::string line;
    while (std::getline(in, line))
    {
      if (line.find("# AmiraMesh") != std::string::npos
        || line.find("# Avizo 3D") != std::string::npos)
      {
        std::cout << "Amira format : " << line << std::endl;
        boost::algorithm::to_lower(line);

        return line.find(format.c_str()) != std::string::npos;
      }
    }
    return false;
  }

  template<typename Tr>
  bool read_tetra_am(std::istream& input, Tr& tr)
  {
    using Material_type = unsigned char;
    using Point_3 = typename Tr::Geom_traits::Point_3;

    int n_nodes, n_tets, n_edges;
    std::vector<Data_def> data;
    std::vector<material> materials; //material = std::pair<int, std::string>

    Material_type material_id = 0;//byte
    bool materials_done = false;

    std::vector<Point_3> points;
    std::vector<std::array<int, 4> > tetrahedra;
    std::vector<Material_type> labels;

    bool binary = true;

    std::string line;
    while (std::getline(input, line))
    {
      if (line.find("# AmiraMesh") != std::string::npos
        || line.find("# Avizo 3D") != std::string::npos)
      {
        std::cout << "Amira format : " << line << std::endl;
        boost::algorithm::to_lower(line);
        if (line.find("ascii") != std::string::npos)
          binary = false;

        continue;
      }
      if (line_starts_with(line, "nNodes"))
      {
        std::string nnodes;
        std::istringstream iss;
        iss.str(line);
        iss >> nnodes >> n_nodes;
      }
      else if (line_starts_with(line, "nTetrahedra"))
      {
        std::string ntetra;
        std::istringstream iss;
        iss.str(line);
        iss >> ntetra >> n_tets;
      }
      else if (line_starts_with(line, "nEdges"))
      {
        std::string nedges;
        std::istringstream iss;
        iss.str(line);
        iss >> nedges >> n_edges;
      }
      else if (!materials_done && line.find("Materials") != std::string::npos)
      {
        if (!IO::internal::treat_surf_materials(input, materials, material_id))
          return false;
        materials_done = true;
      }
      else if (line_starts_with(line, "Nodes"))
      {
        //Nodes { float[3] Coordinates } @1
        Data_def d;
        d.name = "Nodes";
        d.type = "float[3]";
        d.description = "Coordinates";

        std::size_t i = line.find("@");
        d.at_label = line.substr(i);

        data.push_back(d);
      }
      else if (line_starts_with(line, "Tetrahedra"))
      {
        //Tetrahedra{ int[4] Nodes } @3
        Data_def d;
        d.name = "Tetrahedra";
        d.type = "int[4]";
        d.description = "Nodes";

        std::size_t i = line.find("@");
        d.at_label = line.substr(i);

        data.push_back(d);
      }
      else if (line_starts_with(line, "TetrahedronData"))
      {
        //TetrahedronData{ byte Materials } @4
        Data_def d;
        d.name = "TetrahedronData";
        d.type = "byte";
        d.description = "Materials";

        std::size_t i = line.find("@");
        d.at_label = line.substr(i);

        data.push_back(d);
      }
      else if (line.find("Data section follows") != std::string::npos)
      {
        std::getline(input, line);
        for (auto d : data)
        {
          go_to_at_label(input, line, d.at_label.c_str());
          CGAL_assertion(line_starts_with(line, d.at_label.c_str()));

          if (d.name.compare("Nodes") == 0)
            read_points(input, points, n_nodes, binary);
          else if (d.name.compare("Tetrahedra") == 0)
            read_tetrahedra(input, tetrahedra, n_tets, binary);
          else if (d.name.compare("TetrahedronData") == 0)
            read_materials(input, labels, n_tets, binary);
        }
      }
    }

    boost::unordered_map<std::array<int, 3>,
      typename Tr::Cell::Surface_patch_index> empty_facets_map;

    CGAL::SMDS_3::build_triangulation_with_subdomains_range(tr,
      points,
      tetrahedra,
      labels,
      empty_facets_map,
      true,//verbose
      false,//replace subdomain 0
      true);//allow non manifold

    return true;
  }

}// end namespace internal


  template<typename T3>
  bool read_AVIZO_TETRA(std::istream& in, T3& tr)
  {
    if (!in)
    {
      std::cerr << "Cannot open file " << std::endl;
      return false;
    }
    return CGAL::IO::internal::read_tetra_am(in, tr);
  }

} // end namespace IO

} // end namespace CGAL

#endif // CGAL_IO_FILE_AVIZO_H
