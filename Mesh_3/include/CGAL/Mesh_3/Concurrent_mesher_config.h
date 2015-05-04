// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Clement Jamin
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_3_CONCURRENT_MESHER_CONFIG_H
#define CGAL_MESH_3_CONCURRENT_MESHER_CONFIG_H

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
# include <boost/program_options.hpp>
  namespace po = boost::program_options;
#endif

#include <iostream>
#include <fstream>
#include <CGAL/use.h>

// class Concurrent_mesher_config
/// Singleton containing config
class Concurrent_mesher_config
{
  // Private constructor (singleton)
  Concurrent_mesher_config()
  : locking_grid_num_cells_per_axis(50),
    first_grid_lock_radius(0),
    work_stats_grid_num_cells_per_axis(5),
    num_work_items_per_batch(50),
    refinement_grainsize(10),
    refinement_batch_size(10000),
    min_num_vertices_of_coarse_mesh(100),
    num_vertices_of_coarse_mesh_per_core(3.5f),
    num_pseudo_infinite_vertices_per_core(5.0f),
    m_config_file_loaded(false)
  {}

public:
  static Concurrent_mesher_config &get()
  {
    static Concurrent_mesher_config singleton;
    return singleton;
  }

  static bool load_config_file(const char *filename,
    bool reload_if_already_loaded = false)
  {
    return get().load_file(filename, reload_if_already_loaded);
  }


  //=============== PUBLIC PARAMETERS ==============

  // From config file (or default)
  int     locking_grid_num_cells_per_axis;
  int     first_grid_lock_radius;
  int     work_stats_grid_num_cells_per_axis;
  int     num_work_items_per_batch;
  int     refinement_grainsize;
  int     refinement_batch_size;
  int     min_num_vertices_of_coarse_mesh;
  float   num_vertices_of_coarse_mesh_per_core;
  float   num_pseudo_infinite_vertices_per_core;

  // Others


  //================================================

protected:

  bool load_file(
    const char *filename,
    bool reload_if_already_loaded = false)
  {
    CGAL_USE(reload_if_already_loaded);
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
    if (m_config_file_loaded && reload_if_already_loaded == false)
      return true;

    try
    {
      std::ifstream in(filename);
      if (in.fail())
      {
        std::string err = "could not open file '";
        err += filename;
        err += "'. Using default values.";
        throw std::runtime_error(err);
      }

      // Declare the supported options.
      po::options_description desc("Allowed options");
      desc.add_options()
        ("locking_grid_num_cells_per_axis", po::value<int>(), "")
        ("first_grid_lock_radius", po::value<int>(), "")
        ("work_stats_grid_num_cells_per_axis", po::value<int>(), "")
        ("num_work_items_per_batch", po::value<int>(), "")
        ("refinement_grainsize", po::value<int>(), "")
        ("refinement_batch_size", po::value<int>(), "")
        ("min_num_vertices_of_coarse_mesh", po::value<int>(), "")
        ("num_vertices_of_coarse_mesh_per_core", po::value<float>(), "")
        ("num_pseudo_infinite_vertices_per_core", po::value<float>(), "");


      po::store(po::parse_config_file<char>(in, desc), m_variables_map);
      po::notify(m_variables_map);
    }
    catch (std::exception &e)
    {
      std::cerr << "Concurrency configuration file error: " 
        << e.what() << std::endl;
      return false;
    }

    locking_grid_num_cells_per_axis =
      get_config_file_option_value<int>("locking_grid_num_cells_per_axis");
    first_grid_lock_radius =
      get_config_file_option_value<int>("first_grid_lock_radius");
    work_stats_grid_num_cells_per_axis =
      get_config_file_option_value<int>("work_stats_grid_num_cells_per_axis");
    num_work_items_per_batch =
      get_config_file_option_value<int>("num_work_items_per_batch");
    refinement_grainsize =
      get_config_file_option_value<int>("refinement_grainsize");
    refinement_batch_size =
      get_config_file_option_value<int>("refinement_batch_size");
    min_num_vertices_of_coarse_mesh =
      get_config_file_option_value<int>("min_num_vertices_of_coarse_mesh");
    num_vertices_of_coarse_mesh_per_core =
      get_config_file_option_value<float>("num_vertices_of_coarse_mesh_per_core");
    num_pseudo_infinite_vertices_per_core =
      get_config_file_option_value<float>("num_pseudo_infinite_vertices_per_core");

    m_config_file_loaded = true;

#else // CGAL_USE_BOOST_PROGRAM_OPTIONS not defined
    std::cerr << "Warning: could not load concurrency configuration file '"
      << filename << "'. Default values will be used."
      << std::endl;
#endif // CGAL_USE_BOOST_PROGRAM_OPTIONS
    return true;
  }

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
  template <typename OptionType>
  OptionType get_config_file_option_value(const char *option_name)
  {
    if (m_variables_map.count(option_name))
      return m_variables_map[option_name].as<OptionType>();
    else
      return OptionType();
  }
#endif
  
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
  po::variables_map m_variables_map;
#endif
  bool              m_config_file_loaded;
};

#endif // CGAL_MESH_3_CONCURRENT_MESHER_CONFIG_H
