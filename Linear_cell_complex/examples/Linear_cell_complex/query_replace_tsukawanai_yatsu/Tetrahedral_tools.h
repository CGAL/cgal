// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#include <CGAL/Triangulation_3_to_lcc.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include "lcc_to_face_graph.h"
#include "lcc_to_tetgen_io.h"
#include "lcc_triangulate_faces.h"
#include "MeshComplex_3InTriangulation_3_to_lcc.h"
#include <tetgen.h>
////////////////////////////////////////////////////////////////////////////////
/* Code trying to use CGAL mesher, problem to fix parameters.
template<typename LCC>
void tetrahedralize(LCC &lcc)
{
  typedef typename LCC::Traits K;
  typedef typename CGAL::Surface_mesh<typename K::Point_3> Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;
  typedef CGAL::Sequential_tag Concurrency_tag; // Or CGAL::Parallel_tag
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default,
      Concurrency_tag>::type Tr;
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3T3;

  Polyhedron p;
  lcc_to_face_graph(lcc, p);
  CGAL::Polygon_mesh_processing::triangulate_faces(p);

  double cursize, sumsize=0., minsize, maxsize;
  std::size_t nbf=0;
  for(typename boost::graph_traits<Polyhedron>::face_descriptor f: faces(p))
  {
    auto he1=CGAL::halfedge(f, p);
    auto he2=CGAL::next(he1, p);
    auto he3=CGAL::next(he2, p);
    cursize=std::sqrt(CGAL::squared_area(p.point(CGAL::source(he1, p)),
                                         p.point(CGAL::source(he2, p)),
                                         p.point(CGAL::source(he3, p))));
    sumsize+=cursize;
    if (nbf==0 || cursize<minsize) { minsize=cursize; }
    if (nbf==0 || cursize>maxsize) { maxsize=cursize; }
    ++nbf;
  }

  double myvol=CGAL::Polygon_mesh_processing::volume(p);

  std::cout<<"minsize="<<minsize<<"  maxsize="<<maxsize<<"  meansize="
          <<sumsize/nbf<<"   myvol="<<myvol<<std::endl;
  std::size_t nbcells=10000; // Wanted number of cells

  Mesh_domain domain(p);
  using namespace CGAL::parameters;
  double fs=10*sumsize/nbf, fd=fs/100., es=fs;
  std::cout<<"make_mesh_3: facet_size="<<fs<<", edge_size="<<es
          <<", facet_distance="<<fd<<std::endl;

  double s=pow((myvol*10 / * /nbcells * /)*6*sqrt(2), 1.0/3.0);
  fs=s*s*sqrt(3)/4.;
  fd=fs/10;
  es=fs;

  std::cout<<"make_mesh_3: facet_size="<<fs<<", edge_size="<<es
          <<", facet_distance="<<fd<<std::endl;

  Mesh_criteria criteria(cell_size=myvol, ///nbcells,
                         facet_size=fs,
                         edge_size=es,
                         facet_distance=fd);//facet_size=fs, edge_size=es, facet_distance=fd);
                         // cell_radius_edge_ratio=3, facet_angle=30);
        //facet_angle=25, facet_size=0.15, facet_distance=0.008,
         //                cell_radius_edge_ratio=3); // TODO compute parameters from the mesh

  C3T3 c3t3 = CGAL::make_mesh_3<C3T3>(domain, criteria);//, manifold());
  / * edge_size = 8,
                                      facet_angle = 25, facet_size = 0.1, facet_distance = 0.2,
                                      cell_radius_edge_ratio = 3, cell_size = 0.1);* /
  // CGAL::refine_mesh_3(c3t3, domain, criteria);

  lcc.clear();
  import_from_complex_3_in_triangulation_3(lcc, c3t3);
}*/
////////////////////////////////////////////////////////////////////////////////
inline void tetgen_options(tetgenbehavior& b)
{
  // Required option
  b.plc=1;           // '-p', 0. Tetrahedralizes a piecewise linear complex (PLC).
  b.zeroindex=1;     // '-z', 0. Numbers all output items starting from zero.
  b.neighout=1;      // '-n', 0. Outputs tetrahedra neighbors to .neigh file.

  // Required to be able to tetrahedralize only some volumes
  b.nobisect=1;      // '-Y', 0. Preserves the input surface mesh (does not modify it).

  // Optimisation
  b.quality=1;       // '-q', 0. Refines mesh (to improve mesh quality).

  // To debug
  // b.diagnose=1;      // '-d', 0. Detects self-intersections of facets of the PLC.

  /*
  b.psc;             // '-s', 0.
  b.refine;          // '-r', 0. Reconstructs a previously generated mesh.
  b.coarsen;         // '-R', 0. Mesh coarsening (to reduce the mesh elements).
  b.weighted;        // '-w', 0. Generates weighted Delaunay (regular) triangulation.
  b.brio_hilbert;    // '-b', 1.
  b.incrflip;        // '-l', 0.
  b.flipinsert;      // '-L', 0.
  b.metric;          // '-m', 0. Applies a mesh sizing function.
  b.varvolume;       // '-a', 0. Applies a maximum tetrahedron volume constraint.
  b.fixedvolume;     // '-a', 0.
  b.regionattrib;    // '-A', 0. Assigns attributes to tetrahedra in different regions.
  b.conforming;      // '-D', 0.
  b.insertaddpoints; // '-i', 0. Inserts a list of additional points.
  b.convex;          // '-c', 0. Retains the convex hull of the PLC.
  b.nomergefacet;    // '-M', 0. No merge of coplanar facets or very close vertices.
  b.nomergevertex;   // '-M', 0.
  b.noexact;         // '-X', 0. Suppresses use of exact arithmetic.
  b.nostaticfilter;  // '-X', 0.
  b.facesout;        // '-f', 0. Outputs all faces to .face file.
  b.edgesout;        // '-e', 0. Outputs all edges to .edge file.
  b.voroout;         // '-v', 0. Outputs Voronoi diagram to files.
  b.meditview;       // '-g', 0. Outputs mesh to .mesh file for viewing by Medit.
  b.vtkview;         // '-k', 0. Outputs mesh to .vtk file for viewing by Paraview.
  b.nobound;         // '-B', 0. Suppresses output of boundary information.
  b.nonodewritten;   // '-N', 0. Suppresses output of .node file.
  b.noelewritten;    // '-E', 0. Suppresses output of .ele file.
  b.nofacewritten;   // '-F', 0. Suppresses output of .face and .edge file.
  b.noiterationnum;  // '-I', 0. Suppresses mesh iteration numbers.
  b.nojettison;      // '-J', 0. No jettison of unused vertices from output .node file.
  b.reversetetori;   // '-R', 0.
  b.docheck;         // '-C', 0. Checks the consistency of the final mesh.
  b.quiet;           // '-Q', 0. Quiet: No terminal output except errors.
  b.verbose;         // '-V', 0. Verbose: Detailed information, more terminal output.

  b.vertexperblock;     // '-x', 4092.
  b.tetrahedraperblock; // '-x', 8188.
  b.shellfaceperblock;  // '-x', 2044.
  b.nobisect_param;     // '-Y', 2.
  b.addsteiner_algo;    // '-Y/', 1.
  b.coarsen_param;      // '-R', 0. Mesh coarsening (to reduce the mesh elements).
  b.weighted_param;     // '-w', 0. Generates weighted Delaunay (regular) triangulation.
  b.fliplinklevel;      // -1.
  b.flipstarsize;       // -1.
  b.fliplinklevelinc;   //  1.
  b.reflevel;           // '-D', 3.
  b.optlevel;           // '-O', 2. Specifies the level of mesh optimization.
  b.optscheme;          // '-O', 7.
  b.delmaxfliplevel;    // 1.
  b.order;              // '-o', 1.
  b.steinerleft;        // '-S', 0. Specifies maximum number of added points.
  b.no_sort;            // 0.
  b.hilbert_order;      // '-b///', 52.
  b.hilbert_limit;      // '-b//'  8.
  b.brio_threshold;     // '-b' 64.
  b.brio_ratio;         // '-b/' 0.125.
  b.facet_ang_tol;      // '-p', 179.9.
  b.maxvolume;          // '-a', -1.0. Applies a maximum tetrahedron volume constraint.
  b.minratio;           // '-q', 0.0.
  b.mindihedral;        // '-q', 5.0.
  b.optmaxdihedral;     // 165.0.
  b.optminsmtdihed;     // 179.0.
  b.optminslidihed;     // 179.0.
  b.epsilon;            // '-T', 1.0e-8. Sets a tolerance for coplanar test (default 10−8).
  b.minedgelength;      // 0.0.
  b.coarsen_percent;    // -R1/#, 1.0.
*/
}
////////////////////////////////////////////////////////////////////////////////
/// Tetrahedralize all marked volumes of the lcc.
///   (a volume is considered marked if one its dart is marked)
/// At the end, no more dart is marked.
template <typename LCC>
void tetrahedralize_with_tetgen(LCC& lcc, typename LCC::size_type amark)
{
  // 1) Triangulate all marked faces
  triangulate_marked_faces(lcc, amark);

  auto newdart=lcc.get_new_mark();
  lcc.negate_mark(newdart); // Old darts are marked
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if(lcc.is_marked(it, amark))
    {
      // 2) Fill tetgen surface
      tetgenio mysurface;
      tetgenbehavior b;
      lcc_to_tetgen_one_volume_(lcc, it, mysurface);

      // 3) Tetrahedralize
      tetgenio result;
      tetgen_options(b);
      tetrahedralize(&b, &mysurface, &result); // char *switches, tetgenio *in, tetgenio *out, tetgenio *addin = NULL, tetgenio *bgmin = NULL

      // 4) Transfert tetra into lcc
      lcc.template remove_cell<3>(it);
      tetgen_to_lcc(result, lcc);
    }
  }

  lcc.negate_mark(newdart); // Now only new darts are marked
  lcc.sew3_same_facets(newdart); // 3-sew faces of the new tetrahedra with old faces

  lcc.free_mark(newdart);
}
////////////////////////////////////////////////////////////////////////////////
/// Tetrahedralize all volumes of the lcc.
template<typename LCC>
void tetrahedralize_with_tetgen(LCC &lcc)
{
  // 1) Triangulate all faces
  triangulate_all_faces(lcc);

  // 2) Fill tetgen surface
  tetgenio mysurface;
  tetgenbehavior b;
  lcc_to_tetgen(lcc, mysurface);

  // 3) Tetrahedralize
  tetgenio result;
  tetgen_options(b);
  tetrahedralize(&b, &mysurface, &result); // char *switches, tetgenio *in, tetgenio *out, tetgenio *addin = NULL, tetgenio *bgmin = NULL

  // 4) Transfert tetra into lcc
  lcc.clear();
  tetgen_to_lcc(result, lcc);
}
////////////////////////////////////////////////////////////////////////////////
