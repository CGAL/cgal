// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois

#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>

#include <ostream>
#include <unordered_map>
#include <stack>

namespace CGAL
{
namespace IO
{
/*!
 * @ingroup PkgCDT3IOFunctions
 * @brief outputs a conforming constrained Delaunay triangulation to
 * the MEDIT (`.mesh`) file format.
 *        See \cgalCite{frey:inria-00069921} for a comprehensive description of this
 *        file format.
 * @param os the output stream
 * @param ccdt the conforming constrained Delaunay triangulation to be written
 *
 * \see \ref IOStreamMedit
 */
template <typename Traits, typename Tr>
void write_MEDIT(std::ostream& os,
                 const Conforming_constrained_Delaunay_triangulation_3<Traits, Tr>& ccdt)
{
  using CCDT = Conforming_constrained_Delaunay_triangulation_3<Traits, Tr>;
  using T3 = typename CCDT::Triangulation;
  using Vertex_handle = typename T3::Vertex_handle;
  using Cell_handle = typename T3::Cell_handle;
  using Facet = typename T3::Facet;
  using Point = typename T3::Point;

  const bool all_vertices = true;
  const bool all_cells = true;

  const auto& tr = ccdt.triangulation();

  //-------------------------------------------------------
  // Header
  //-------------------------------------------------------
  os << std::setprecision(17);

  os << "MeshVersionFormatted 1\n"
     << "Dimension 3\n";
  os << "# CGAL::Conforming_constrained_Delaunay_triangulation_3\n";

  //-------------------------------------------------------
  // Vertices

  std::unordered_map<Vertex_handle, std::size_t> V;
  std::size_t inum = 1;
  if(all_vertices || all_cells)
  {
    os << "Vertices\n" << tr.number_of_vertices() << '\n';

    for(auto vit : tr.finite_vertex_handles())
    {
      V[vit] = inum++;
      const Point& p = tr.point(vit);
      os << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << ' ' << CGAL::to_double(p.z()) << ' '
         << "0" << '\n';
    }
  }
  else
  {
    std::ostringstream oss;
    for(Cell_handle c : tr.finite_cell_handles())
    {
      for(int i = 0; i < 4; ++i)
      {
        Vertex_handle vit = c->vertex(i);
        if(V.find(vit) == V.end())
        {
          V[vit] = inum++;
          const Point& p = tr.point(vit);
          oss << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << ' ' << CGAL::to_double(p.z()) << ' '
              << "0" << '\n';
        }
      }
    }
    os << "Vertices\n" << V.size() << "\n";
    os << oss.str();
  }

  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  typename T3::size_type number_of_triangles = ccdt.number_of_constrained_facets();

  os << "Triangles\n" << number_of_triangles << '\n';

  for(Facet f : ccdt.constrained_facets())
  {
   // Get facet vertices in CCW order.
    Vertex_handle vh1 = f.first->vertex((f.second + 1) % 4);
    Vertex_handle vh2 = f.first->vertex((f.second + 2) % 4);
    Vertex_handle vh3 = f.first->vertex((f.second + 3) % 4);

    // Facet orientation also depends on parity.
    if(f.second % 2 != 0)
      std::swap(vh2, vh3);

    os << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3] << ' ';
    os << f.first->ccdt_3_data().face_constraint_index(f.second) + 1 << '\n';
  }

  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  std::unordered_map<Cell_handle, bool> cells_in_domain;
  for(Cell_handle c : tr.all_cell_handles())
    cells_in_domain[c] = true;

  std::stack<Cell_handle> stack;
  stack.push(tr.infinite_cell());
  while(!stack.empty())
  {
    auto ch = stack.top();
    stack.pop();
    cells_in_domain[ch] = false;
    for(int i = 0; i < 4; ++i)
    {
      if(ch->ccdt_3_data().is_facet_constrained(i))
        continue;
      auto n = ch->neighbor(i);
      if(cells_in_domain[n])
        stack.push(n);
    }
  }

  std::vector<std::array<std::size_t, 4>> indexed_tetra;
  indexed_tetra.reserve(tr.number_of_cells());
  for(Cell_handle ch : tr.finite_cell_handles())
  {
    if(cells_in_domain[ch])
    {
      indexed_tetra.push_back({V.at(ch->vertex(0)), V.at(ch->vertex(1)),
                               V.at(ch->vertex(2)), V.at(ch->vertex(3))});
    }
  }

  typename T3::size_type number_of_cells =
      all_cells ? ccdt.triangulation().number_of_finite_cells()
                : indexed_tetra.size();

  os << "Tetrahedra\n" << number_of_cells << '\n';
  for(Cell_handle c : tr.finite_cell_handles())
  {
    for(int i = 0; i < 4; i++)
      os << V[c->vertex(i)] << ' ';

    int subdomain_index = cells_in_domain[c] ? 1 : 0;
    os << subdomain_index << '\n';
  }

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End\n";
}

}// end namespace IO
}// end namespace CGAL
