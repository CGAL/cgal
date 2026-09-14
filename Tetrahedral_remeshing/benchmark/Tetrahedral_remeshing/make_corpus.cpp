// Build a volumetric corpus from surface meshes with make_mesh_3.
//
// Purpose: estimate ParMmg's speedup over sequential Mmg over a BROAD set of
// meshes. Both the existing cdt and mesh3 datasets are mostly unusable for that
// -- Mmg's MMG5_minQualCheck rejects a mesh outright if any element has zero
// quality, which killed 5/9 cdt and 6/9 mesh3 meshes. Mmg's criterion is
// VOLUME-based, so a min-dihedral screen does not predict it: mesh3_bear has a
// min dihedral of 0.21 degrees and is still rejected.
//
// The fix is exudation. make_mesh_3's optimisers (perturb, then exude) exist
// precisely to remove slivers, and without them make_mesh_3 output retains a
// few near-degenerate cells. We run both, so the corpus should load in Mmg.
//
// Sizing is relative to each model's bounding-box diagonal, so meshes of very
// different scales get comparable element counts and the corpus spans a range
// of shapes rather than a range of units.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/Bbox_3.h>

#include <fstream>
#include <iostream>
#include <string>

using K          = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Polyhedron_3<K>;
using Domain     = CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>;
using Tr         = CGAL::Mesh_triangulation_3<Domain>::type;
using C3t3       = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
using Criteria   = CGAL::Mesh_criteria_3<Tr>;

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <input.off> <output.mesh> [size_divisor=25]\n";
    return 2;
  }
  const std::string in = argv[1], out = argv[2];
  const double div = (argc > 3) ? std::stod(argv[3]) : 25.0;

  Polyhedron poly;
  {
    std::ifstream is(in);
    if (!is || !(is >> poly) || poly.empty())
    {
      std::cerr << "SKIP " << in << " : cannot read polyhedron\n";
      return 1;
    }
  }
  if (!poly.is_closed())
  {
    std::cerr << "SKIP " << in << " : surface is not closed\n";
    return 1;
  }

  CGAL::Bbox_3 bb;
  for (auto v = poly.points_begin(); v != poly.points_end(); ++v)
    bb += v->bbox();
  const double dx = bb.xmax() - bb.xmin(), dy = bb.ymax() - bb.ymin(), dz = bb.zmax() - bb.zmin();
  const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double h = diag / div;

  try
  {
    Domain domain(poly);
    Criteria criteria(CGAL::parameters::facet_angle(25)
                                       .facet_size(h)
                                       .facet_distance(h / 10.)
                                       .cell_radius_edge_ratio(3)
                                       .cell_size(h));

    // perturb + exude are what remove the slivers that make Mmg refuse the mesh
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::perturb().exude());

    if (c3t3.number_of_cells_in_complex() == 0)
    {
      std::cerr << "SKIP " << in << " : empty complex\n";
      return 1;
    }
    std::ofstream os(out);
    CGAL::IO::write_MEDIT(os, c3t3.triangulation());
    std::cout << "OK " << out << "  cells=" << c3t3.number_of_cells_in_complex()
              << "  h=" << h << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cerr << "SKIP " << in << " : exception: " << e.what() << "\n";
    return 1;
  }
  catch (...)
  {
    std::cerr << "SKIP " << in << " : unknown exception\n";
    return 1;
  }
  return 0;
}
