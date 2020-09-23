// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_TR_INTERNAL_COMPUTE_C3T3_STATISTICS_H
#define CGAL_TR_INTERNAL_COMPUTE_C3T3_STATISTICS_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <limits>
#include <vector>
#include <algorithm>
#include <fstream>

#include <boost/unordered_set.hpp>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
template<typename Triangulation, typename CellSelector>
void compute_statistics(const Triangulation& tr,
                        CellSelector cell_selector,
                        const char* filename = "statistics_c3t3.txt")
{
  typedef Triangulation Tr;
  typedef typename Tr::Geom_traits Gt;
  typedef typename Tr::Cell_handle   Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Gt::Point_3       Point;
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_cells_iterator  Finite_cells_iterator;
  typedef typename Tr::Cell::Subdomain_index  Subdomain_index;

  std::size_t nb_edges = 0;
  double total_edges = 0;
  std::size_t nb_angle = 0;
  double total_angle = 0;

  double min_edges_length = (std::numeric_limits<double>::max)();
  double max_edges_length = 0.;

  double smallest_edge_radius = (std::numeric_limits<double>::max)();
  double smallest_radius_radius = (std::numeric_limits<double>::max)();
  double biggest_v_sma_cube = 0.;
  double max_dihedral_angle = 0.;
  double min_dihedral_angle = 180.;

  for (Finite_facets_iterator fit = tr.finite_facets_begin();
       fit != tr.finite_facets_end(); ++fit)
  {
    const Cell_handle cell = fit->first;
    const int& index = fit->second;
    if (!cell_selector(cell) || !cell_selector(cell->neighbor(index)))
      continue;

    const Point& pa = point(cell->vertex((index + 1) & 3)->point());
    const Point& pb = point(cell->vertex((index + 2) & 3)->point());
    const Point& pc = point(cell->vertex((index + 3) & 3)->point());

    double edges[3];
    edges[0] = (CGAL::sqrt(CGAL::squared_distance(pa, pb)));
    edges[1] = (CGAL::sqrt(CGAL::squared_distance(pa, pc)));
    edges[2] = (CGAL::sqrt(CGAL::squared_distance(pb, pc)));
    for (int i = 0; i < 3; ++i)
    {
      if (edges[i] < min_edges_length){ min_edges_length = edges[i]; }
      if (edges[i] > max_edges_length){ max_edges_length = edges[i]; }
      total_edges += edges[i];
      ++nb_edges;
    }
  }

  double mean_edges_length = total_edges / (double)nb_edges;

  typename Gt::Compute_approximate_dihedral_angle_3 approx_dihedral_angle
    = tr.geom_traits().compute_approximate_dihedral_angle_3_object();

  std::size_t nb_tets = 0;
  boost::unordered_set<Vertex_handle> selected_vertices;
  std::vector<Subdomain_index> sub_ids;
  for (Finite_cells_iterator cit = tr.finite_cells_begin();
       cit != tr.finite_cells_end();
       ++cit)
  {
    const Subdomain_index& si = cit->subdomain_index();
    if (si == Subdomain_index() || !cell_selector(cit))
      continue;

    ++nb_tets;
    if (std::find(sub_ids.begin(), sub_ids.end(), si) == sub_ids.end())
      sub_ids.push_back(cit->subdomain_index());
    for (int i = 0; i < 4; ++i)
      selected_vertices.insert(cit->vertex(i));

    const Point& p0 = point(cit->vertex(0)->point());
    const Point& p1 = point(cit->vertex(1)->point());
    const Point& p2 = point(cit->vertex(2)->point());
    const Point& p3 = point(cit->vertex(3)->point());
    double v = CGAL::abs(tr.tetrahedron(cit).volume());
    if (v == 0.)
    {
      std::cout << "degenerate cell :\n\t";
      std::cout << p0 << "\n\t" << p1 << "\n\t" << p2 << "\n\t" << p3 << std::endl;
    }
    double circumradius = (v == 0.)
                          ? CGAL::sqrt(CGAL::squared_radius(p0, p1, p2))
                          : CGAL::sqrt(CGAL::squared_radius(p0, p1, p2, p3));

    //find shortest edge
    double edges[6];
    edges[0] = CGAL::sqrt(CGAL::squared_distance(p0, p1));
    edges[1] = CGAL::sqrt(CGAL::squared_distance(p0, p2));
    edges[2] = CGAL::sqrt(CGAL::squared_distance(p0, p3));
    edges[3] = CGAL::sqrt(CGAL::squared_distance(p2, p1));
    edges[4] = CGAL::sqrt(CGAL::squared_distance(p2, p3));
    edges[5] = CGAL::sqrt(CGAL::squared_distance(p1, p3));

    double min_edge = edges[0];
    for (int i = 1; i < 6; ++i)
    {
      if (edges[i] < min_edge)
        min_edge = edges[i];
    }

    double sumar = CGAL::sqrt(CGAL::squared_area(p0, p1, p2))
                   + CGAL::sqrt(CGAL::squared_area(p1, p2, p3))
                   + CGAL::sqrt(CGAL::squared_area(p2, p3, p0))
                   + CGAL::sqrt(CGAL::squared_area(p3, p1, p0));
    double inradius = 3. * v / sumar;
    double smallest_edge_radius_ = min_edge / circumradius*CGAL::sqrt(6.) / 4.;//*sqrt(6)/4 so that the perfect tet ratio is 1
    double smallest_radius_radius_ = inradius / circumradius * 3.; //*3 so that the perfect tet ratio is 1 instead of 1/3
    double biggest_v_sma_cube_ = v / std::pow(min_edge, 3) * 6. * CGAL::sqrt(2.);//*6*sqrt(2) so that the perfect tet ratio is 1 instead

    if (smallest_edge_radius_ < smallest_edge_radius)
      smallest_edge_radius = smallest_edge_radius_;

    if (smallest_radius_radius_ < smallest_radius_radius)
      smallest_radius_radius = smallest_radius_radius_;

    if (biggest_v_sma_cube_ > biggest_v_sma_cube)
      biggest_v_sma_cube = biggest_v_sma_cube_;

    double a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p1, p2, p3)));
    if (a < min_dihedral_angle) { min_dihedral_angle = a; }
    if (a > max_dihedral_angle) { max_dihedral_angle = a; }
    total_angle += a;
    ++nb_angle;
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p2, p1, p3)));
    if (a < min_dihedral_angle) { min_dihedral_angle = a; }
    if (a > max_dihedral_angle) { max_dihedral_angle = a; }
    total_angle += a;
    ++nb_angle;
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p3, p1, p2)));
    if (a < min_dihedral_angle) { min_dihedral_angle = a; }
    if (a > max_dihedral_angle) { max_dihedral_angle = a; }
    total_angle += a;
    ++nb_angle;
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p1, p2, p0, p3)));
    if (a < min_dihedral_angle) { min_dihedral_angle = a; }
    if (a > max_dihedral_angle) { max_dihedral_angle = a; }
    total_angle += a;
    ++nb_angle;
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p1, p3, p0, p2)));
    if (a < min_dihedral_angle) { min_dihedral_angle = a; }
    if (a > max_dihedral_angle) { max_dihedral_angle = a; }
    total_angle += a;
    ++nb_angle;
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p2, p3, p0, p1)));
    if (a < min_dihedral_angle) { min_dihedral_angle = a; }
    if (a > max_dihedral_angle) { max_dihedral_angle = a; }
    total_angle += a;
    ++nb_angle;
  }

  std::size_t nb_subdomains = sub_ids.size();
  //std::size_t nb_vertices = d->c3t3.number_of_vertices_in_complex();

  std::ofstream ofs(filename);
  if (!ofs)
    return;

  ofs << "Nb subdomains               : " << nb_subdomains << std::endl;
  ofs << "Total number of vertices    : " << tr.number_of_vertices() << std::endl;
  ofs << "Number of selected cells    : " << nb_tets << std::endl;
  ofs << "Number of selected vertices : " << selected_vertices.size() << std::endl;
  ofs << std::endl;
  ofs << "Min dihedral angle : " << min_dihedral_angle << std::endl;
  ofs << "Max dihedral angle : " << max_dihedral_angle << std::endl;
  ofs << std::endl;
  ofs << "Shortest edge       : " << min_edges_length << std::endl;
  ofs << "Longest edge        : " << max_edges_length << std::endl;
  ofs << "Average edge length : " << mean_edges_length << std::endl;

  ofs.close();
}

}//end namespace internal
}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif // CGAL_TR_INTERNAL_COMPUTE_C3T3_STATISTICS_H
