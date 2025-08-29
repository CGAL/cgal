// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//
#ifndef HEXMESHING_CREATE_VERTICES_FOR_TEMPLATES_H
#define HEXMESHING_CREATE_VERTICES_FOR_TEMPLATES_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>

namespace CGAL::internal::Hexmeshing {
  /**
   * @brief Creates new vertices for template refinement by subdividing edges
   * 
   * This function creates new vertices that are needed for the template refinement
   * process. It identifies edges that need to be subdivided based on the marked
   * nodes and creates new vertices at the barycenter of these edges.
   * 
   * The function follows a specific rule for vertex creation:
   * - Two adjacent marked nodes do not produce a new vertex
   * - One marked node adjacent to an unmarked node produces a new vertex
   * 
   * The function performs the following operations:
   * 
   * 1. **Edge Identification**: Iterates through all marked nodes and examines
   *    their incident edges using `lcc.one_dart_per_incident_cell<1, 0>()`
   * 
   * 2. **Edge Filtering**: For each edge incident to a marked node:
   *    - Skips edges that have already been processed (marked with `arrete_done`)
   *    - Skips edges where both endpoints are marked (no vertex creation needed)
   *    - Adds edges to the subdivision list where only one endpoint is marked
   * 
   * 3. **Vertex Creation**: For each edge that needs subdivision:
   *    - Inserts a barycenter in the edge using `lcc.insert_barycenter_in_cell<1>()`
   *    - Calls `thread_number_vertex_in_edge()` to handle any thread-specific
   *      vertex numbering operations
   * 
   * 4. **Cleanup**: Frees the temporary mark used for tracking processed edges
   * 
   * This function is called during the refinement stage to prepare the mesh
   * for template substitution by creating the necessary intermediate vertices.
   * 
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the marked nodes for vertex creation
   */
  template <typename HexData>
  void create_vertices_for_templates(HexData& hdata, RefinementData& rdata)
  {
    // 2 noeuds marqué l'un à coté de l'autre ne produit pas de sommet
    // 1 noeud marqué a coté d'un noeud non marqué produit un sommet

    std::vector<Dart_handle> edges_to_subdivide;
    LCC& lcc = hdata.lcc;

    auto arrete_done = lcc.get_new_mark();

    int vertices_created = 0;
    for (auto dart : rdata.marked_nodes)
    {
      for (auto nit = lcc.one_dart_per_incident_cell<1, 0>(dart).begin(),
                nend = lcc.one_dart_per_incident_cell<1, 0>(dart).end();
          nit != nend;
          nit++)
      {
        if (lcc.is_marked(nit, arrete_done))
          continue;

        // If the node is next to an other marked node, we don't have to create vertices
        if (lcc.is_marked(lcc.beta<1>(nit), hdata.template_mark)){
          lcc.mark_cell<1>(nit, arrete_done);
          continue;
        }

        vertices_created++;
        edges_to_subdivide.push_back(nit);
        lcc.mark_cell<1>(nit, arrete_done);
      }
    }

    for (Dart_handle dart : edges_to_subdivide)
    {
      Dart_handle other_ext = lcc.other_extremity(dart);
      Dart_handle new_node = lcc.insert_barycenter_in_cell<1>(dart);
      thread_number_vertex_in_edge(hdata, new_node, dart, other_ext);
    }

    std::cout << "Vertices created: " << vertices_created << std::endl;

    lcc.free_mark(arrete_done);
  }
}



#endif