// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_DRAW_CT2_H
#define CGAL_DRAW_CT2_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_2/internal/In_domain.h>
#include <CGAL/draw_triangulation_2.h>

namespace CGAL
{

template<class CDT>
struct Graphics_scene_options_constrained_triangulation_2:
    public CGAL::Graphics_scene_options<typename CDT::Triangulation,
                                        typename CDT::Vertex_handle,
                                        typename CDT::Finite_edges_iterator,
                                        typename CDT::Finite_faces_iterator>
{
  using BASET2=typename CDT::Triangulation;

  Graphics_scene_options_constrained_triangulation_2(const CDT& cdt)
  {
    this->colored_edge =
      [&cdt](const BASET2&, typename CDT::Finite_edges_iterator eh) -> bool
      { return cdt.is_constrained(*eh); };

    this->edge_color =
      [&cdt](const BASET2&, typename CDT::Finite_edges_iterator eh) -> CGAL::IO::Color
      { return cdt.is_constrained(*eh)? CGAL::IO::green() : CGAL::IO::black(); };
  };

  template<class InDomainPmap>
  Graphics_scene_options_constrained_triangulation_2(const CDT& cdt, InDomainPmap ipm)
  {
    this->colored_edge =
      [&cdt](const BASET2&, typename CDT::Finite_edges_iterator eh) -> bool
      { return cdt.is_constrained(*eh); };

    this->edge_color =
      [&cdt](const BASET2&, typename CDT::Finite_edges_iterator eh) -> CGAL::IO::Color
      { return cdt.is_constrained(*eh)? CGAL::IO::green() : CGAL::IO::black(); };

    this->colored_face =
      [](const BASET2&, typename CDT::Finite_faces_iterator) -> bool
      { return true; };

    this->face_color =
      [ipm](const BASET2&, typename CDT::Finite_faces_iterator fh) -> CGAL::IO::Color
      { return get(ipm, fh)? CGAL::IO::blue() : CGAL::IO::white(); };
  };
};

// In the Update_cell_from_CDT_3 struct in C3t3_io_plugin.cpp
struct Update_cell_from_CDT_3 {
  typedef Fake_mesh_domain::Surface_patch_index Sp_index;
  template <typename C1, typename C2>
  void operator()(const C1& c1, C2& c2) {
    // Preserve the actual subdomain index based on nesting level
    // instead of setting all to 1
    c2.set_subdomain_index(c1.info().nesting_level + 1);
    
    // Set surface patch indices for constrained facets
    for(int i = 0; i < 4; ++i) {
      if(c1.constrained_facet[i])
        c2.set_surface_patch_index(i, 1); // Use index 1 for constrained facets
      else
        c2.set_surface_patch_index(i, 0); // Use index 0 for non-constrained facets
    }
  }
};

// In the visualization code where colors are defined
// This is likely in the CDT_3 plugin's draw or display method
std::vector<QColor> colors;
colors.resize(2); // For non-constrained (0) and constrained (1) facets
colors[0] = CGAL::IO::black(); 
colors[1] = CGAL::IO::green(); 

std::vector<QColor> domain_colors;
domain_colors.resize(2); // For outside (0) and inside (1) domains
domain_colors[0] = CGAL::IO::white();
domain_colors[1] = CGAL::IO::blue();

// For faces (triangles)
QColor color = colors[cell->surface_patch_index(index)];
f_colors.push_back((float)color.redF());
f_colors.push_back((float)color.greenF());
f_colors.push_back((float)color.blueF());

// For cells (tetrahedra)
QColor cell_color = domain_colors[cell->subdomain_index()];
cell_colors.push_back((float)cell_color.redF());
cell_colors.push_back((float)cell_color.greenF());
cell_colors.push_back((float)cell_color.blueF());

// Specialization of draw function.
#define CGAL_T2_TYPE CGAL::Constrained_triangulation_2<Gt, Tds, Itag>

template <class Gt, class Tds, class Itag, class InDomainPmap>
void add_to_graphics_scene(const CGAL_T2_TYPE& at2, InDomainPmap ipm,
                           CGAL::Graphics_scene& graphics_scene)
{
  Graphics_scene_options_constrained_triangulation_2<CGAL_T2_TYPE> gso(at2, ipm);
  draw_function_for_t2::compute_elements(at2, graphics_scene, gso);
}

template <class Gt, class Tds, class Itag>
void add_to_graphics_scene(const CGAL_T2_TYPE& at2,
                           CGAL::Graphics_scene& graphics_scene)
{
  Graphics_scene_options_constrained_triangulation_2<CGAL_T2_TYPE> gso(at2);
  draw_function_for_t2::compute_elements(at2, graphics_scene, gso);
}

#ifdef CGAL_USE_BASIC_VIEWER

template<class Gt, class Tds, class Itag, class InDomainPmap>
void draw(const CGAL_T2_TYPE& at2, InDomainPmap ipm,
          const char *title="Constrained Triangulation_2 Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(at2, ipm, buffer);
  draw_graphics_scene(buffer, title);
}

template<class Gt, class Tds, class Itag>
void draw(const CGAL_T2_TYPE& at2,
          const char *title="Constrained Triangulation_2 Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(at2, buffer);
  draw_graphics_scene(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_T2_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_CT2_H
