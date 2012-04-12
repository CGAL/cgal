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

#ifdef CONCURRENT_MESH_3

#ifndef CGAL_MESH_3_CONCURRENT_MESHER_CONFIG_H
#define CGAL_MESH_3_CONCURRENT_MESHER_CONFIG_H

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>

// class Concurrent_mesher_config
/// Singleton containing config
class Concurrent_mesher_config
{
  // Private constructor (singleton)
  Concurrent_mesher_config()
  : m_loaded(false) {}

public:
  static Concurrent_mesher_config &get()
  {
    static Concurrent_mesher_config singleton;
    return singleton;
  }

  static bool load_config_file(const char *filename)
  {
    return get().load_file(filename);
  }

  static void unload_config_file()
  {
    get().unload_file();
  }

  template <typename OptionType>
  static OptionType get_option(const char *option_name)
  {
    return get().get_option_value<OptionType>(option_name);
  }

protected:
  
  bool load_file(const char *filename)
  {
    try
    {
      // Declare the supported options.
      po::options_description desc("Allowed options");
      desc.add_options()
        ("locking_grid_num_cells_per_axis", po::value<int>(), "")
        ("first_grid_lock_radius", po::value<int>(), "")
        ("work_stats_grid_num_cells_per_axis", po::value<int>(), "")
        ("num_work_items_per_batch", po::value<int>(), "")
        ("refinement_grainsize", po::value<int>(), "")
        ("refinement_batch_size", po::value<int>(), "");

      po::store(po::parse_config_file<char>(filename, desc), m_variables_map);
      po::notify(m_variables_map); 
    }
    catch (std::exception &e)
    {
      std::cerr << "Config file error: " << e.what() << std::endl;
      return false;
    }
    m_loaded = true;
    return true;
  }

  void unload_file()
  {
    m_loaded = false;
  }

  template <typename OptionType>
  OptionType get_option_value(const char *option_name)
  {
    if (!m_loaded)
      load_file(CONFIG_FILENAME);
  
    if (m_loaded && m_variables_map.count(option_name))
      return m_variables_map[option_name].as<OptionType>();
    else
      return OptionType();
  }

  bool              m_loaded;
  po::variables_map m_variables_map;
};

#endif // CGAL_MESH_3_CONCURRENT_MESHER_CONFIG_H

#endif // CONCURRENT_MESH_3