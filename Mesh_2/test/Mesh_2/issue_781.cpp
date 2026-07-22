// #define CGAL_MESH_2_VERBOSE 1
#define CGAL_MESH_2_DEBUG_REFINEMENT_POINTS 1
#define CGAL_MESHES_DEBUG_REFINEMENT_POINTS 1
// #define CGAL_MESH_2_DEBUG_BAD_EDGES 1
// #define CGAL_MESH_2_DEBUG_BAD_FACES 1
// #define CGAL_MESH_2_DEBUG_CLUSTERS 1
// #define CGAL_MESH_2_DEBUG_INSERTIONS 1

// main.cxx -- CGAL triangulation test utility
//
// Written by Peter Sadrozinski, started July 2015.
//
// Copyright (C) 2015  Peter Sadrozinski
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
//

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <boost/property_map/function_property_map.hpp>

#include <CGAL/IO/write_VTU.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/boost/graph/graph_traits_Constrained_Delaunay_triangulation_2.h>
#include <CGAL/bisect_failures.h>

// triangulation
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

// mesh refinement
#include <CGAL/Triangulation_simplex_base_with_time_stamp.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>

// IO
#include <CGAL/IO/io.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Triangle_2.h>


using meshTriKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using meshTriPoint = meshTriKernel::Point_2;
using meshTriSegment = meshTriKernel::Segment_2;

using meshTriVertexBase =
    CGAL::Triangulation_simplex_base_with_time_stamp<CGAL::Triangulation_vertex_base_2<meshTriKernel>>;

using Fbb = CGAL::Constrained_triangulation_face_base_2<meshTriKernel>;
using meshTriFaceBase = CGAL::Delaunay_mesh_face_base_2<meshTriKernel, Fbb>;

using meshTriTDS = CGAL::Triangulation_data_structure_2<meshTriVertexBase, meshTriFaceBase>;
using meshTriItag = CGAL::Exact_intersections_tag;
using meshTriCDT = CGAL::Constrained_Delaunay_triangulation_2<meshTriKernel, meshTriTDS, meshTriItag>;

using meshTriEdge = meshTriCDT::Edge;
using meshTriFaceHandle = meshTriCDT::Face_handle;
using meshTriFaceIterator = meshTriCDT::Finite_faces_iterator;
using meshTriangle = CGAL::Triangle_2<meshTriKernel>;

using meshCriteria = CGAL::Delaunay_mesh_size_criteria_2<meshTriCDT>;
using meshRefiner = CGAL::Delaunay_mesher_2<meshTriCDT, meshCriteria>;

struct Data {
  std::vector<meshTriPoint> points;
  std::vector<meshTriSegment> constraints;
};

auto load(std::string filename) -> Data {
  Data data;
  if(filename.empty()) {
    data.points = {
        meshTriPoint{5.4691594172333904, 44.256641611715409},
    };
    data.constraints = {meshTriSegment{meshTriPoint{5.4693788499999929, 44.256578099999999},
                                           meshTriPoint{5.4691178249999917, 44.256653649999997}}};
    return data;
  }
  std::ifstream input_file(filename);
  if(!input_file.is_open()) {
    std::cerr << "Failed to open the " << filename << std::endl;
    std::exit(-1);
  }

  size_t num_points;
  size_t num_constraints;

  input_file >> num_points;
  data.points.reserve(num_points);
  for(unsigned int i = 0; i < num_points; i++) {
    meshTriPoint pt;
    input_file >> pt;
    data.points.push_back(pt);

//    cdt.insert(pt);
  }

  input_file >> num_constraints;
  data.constraints.reserve(num_constraints);
  for(unsigned int i = 0; i < num_constraints; i++) {
    meshTriPoint s, t;
    input_file >> s >> t;
    data.constraints.push_back(meshTriSegment(s, t));

    // cdt.insert_constraint(s, t);
  }
  input_file.close();
  return data;
}

auto save(std::string filename, const Data& data) -> void {
  std::ofstream output_file(filename);
  output_file.precision(17);
  if(!output_file.is_open()) {
    std::cerr << "Failed to open the " << filename << std::endl;
    std::exit(-1);
  }

  output_file << data.points.size() << std::endl;
  for(const auto& p : data.points) {
    output_file << p << std::endl;
  }

  output_file << data.constraints.size() << std::endl;
  for(const auto& s : data.constraints) {
    output_file << s.source() << " " << s.target() << std::endl;
  }
  output_file.close();
}

auto create_cdt(Data data) -> meshTriCDT {
  meshTriCDT cdt;
  for(const auto& p : data.points) {
    cdt.insert(p);
  }
  for(const auto& s : data.constraints) {
    cdt.insert_constraint(s.source(), s.target());
  }
  return cdt;
}

auto write_off(const std::string& filename, const meshTriCDT& cdt) {
  auto vertex_to_point_3 = [](auto vh) {
    const auto& p2 = vh->point();
    return meshTriKernel::Point_3(p2.x(), p2.y(), 0.0);
  };
  auto vpm = boost::make_function_property_map<meshTriCDT::Vertex_handle>(vertex_to_point_3);
  return CGAL::IO::write_OFF(filename, cdt, CGAL::parameters::stream_precision(17).vertex_point_map(vpm));
};

int main(int argc, char* argv[]) {
  std::cerr.precision(17);
  std::string filename = (argc > 1) ? argv[1] : "";
  Data data = load(filename);

  auto get_size_constraints = [](const Data& data) -> std::size_t {
    return data.constraints.size();
  };

  auto get_size_points = [](const Data& data) -> std::size_t {
    return data.points.size();
  };

  auto save = [](const Data& data, const std::string& prefix) {
    std::string out_filename = prefix + ".txt";
    ::save(out_filename, data);
    std::cout << "Saved data with " << data.constraints.size()
              << " constraints to " << out_filename << std::endl;
  };

  auto run = [](Data& data) {
    meshTriCDT cdt = create_cdt(std::move(data));
    CGAL::IO::write_VTU("before_refinement.vtu", cdt, CGAL::IO::ASCII);
    write_off("before_refinement.off", cdt);
    meshRefiner mesher(cdt);
    mesher.set_criteria(meshCriteria(0.125));

    std::cout << "refine mesh" << std::endl;
    mesher.init();
    meshTriPoint p;
    try {
      while(mesher.step_by_step_refine_mesh()) {
        // CGAL_assertion(cdt.is_valid(true) );
        // write_off("current_cdt_2.off", cdt);
      }
    } catch(...) {
      CGAL::IO::write_VTU("current_cdt_2.vtu", cdt, CGAL::IO::ASCII);
      write_off("current_cdt.off", cdt);
      throw;
    }
    std::cout << "complete: " << cdt.number_of_vertices() << " vertices, " << cdt.number_of_faces() << " faces."
              << std::endl;
    CGAL::IO::write_VTU("after_refinement.vtu", cdt, CGAL::IO::ASCII);
    return EXIT_SUCCESS;
  };

  auto simplify_constraints = [](Data& data, int start, int end) -> bool {
    auto begin = data.constraints.begin() + start;
    auto finish = data.constraints.begin() + end;
    data.constraints.erase(begin, finish);
    return true;
  };

  auto simplify_points = [](Data& data, int start, int end) -> bool {
    auto begin = data.points.begin() + start;
    auto finish = data.points.begin() + end;
    data.points.erase(begin, finish);
    return true;
  };

  if(argc > 2 && argv[2] == std::string("points")) {
    return CGAL::bisect_failures(data, get_size_points, simplify_points, run, save);
  }

  if(argc > 2 && argv[2] == std::string("constraints")) {
    return CGAL::bisect_failures(data, get_size_constraints, simplify_constraints, run, save);
  }

  return run(data);
}
