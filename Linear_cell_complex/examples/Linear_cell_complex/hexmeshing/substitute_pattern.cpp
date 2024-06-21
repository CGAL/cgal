#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>

#include "../query_replace/cmap_query_replace.h"
#include "../query_replace/init_to_preserve_for_query_replace.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT       FT;
typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3> LCC;
typedef typename LCC::Dart_handle Dart_handle;
typedef typename LCC::Vertex_attribute_handle Vertex_handle;
typedef typename LCC::size_type   size_type;

const size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

void render(LCC &lcc, size_type mark)
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(lcc, buffer);

  // Hacking a bit here, we are going to draw colored points for marked vertices just for debugging purposes
  auto mono = const_cast<CGAL::Buffer_for_vao &>(buffer.get_buffer_for_mono_points());
  mono.clear();

  for (auto it = lcc.one_dart_per_cell<0>().begin(), end = lcc.one_dart_per_cell<0>().end(); it != end; it++)
  {

    if (lcc.is_marked(it, mark))
    {
      buffer.add_point(
          lcc.point(it),
          CGAL::Color(0, 255, 0));
      continue;
    }

    buffer.add_point(lcc.point(it));
  }
  draw_graphics_scene(buffer);
  }


void refine_hexes_with_1_2_4_templates(LCC &lcc, Pattern_substituer<LCC> &ps, size_type toprocess_mark)
{
    /**
   * Il se pourrait qu'on ait pas besoins d'identifier les patterns avec leur marking lors du remplacage surfacique/volumique
   * car on a déjà donné la structure à l'étape précedente. Il faudrait donner quel remplacement utiliser de manière directe
   * pour éviter d'essayer tous les patterns 
   */

  // Replace all volumes that has a matching surfacic pattern
  std::vector<Dart_handle> marked_cells;


  for (auto dit = lcc.one_dart_per_cell<2>().begin(),
            dend = lcc.one_dart_per_cell<2>().end();
       dit != dend;
       dit++)
  {
    Dart_handle dart = dit;
    marked_cells.push_back(dart);
  }


  int nbsub = 0;
  for (auto& dart : marked_cells)
  {
    if (ps.query_replace_one_face(lcc, dart, toprocess_mark) != SIZE_T_MAX) nbsub++;
  }
  
  std::cout << nbsub << " Faces substitution was made" << std::endl;

  // Refine the resulting volumes
  marked_cells.clear();

  for (auto dit = lcc.one_dart_per_cell<3>().begin(),
            dend = lcc.one_dart_per_cell<3>().end();
       dit != dend;
       dit++)
  {
    Dart_handle dart = dit;
    marked_cells.push_back(dart);
  }

  nbsub = 0;
  for (auto& dart : marked_cells)
  {
    if (ps.query_replace_one_volume(lcc, dart) != SIZE_T_MAX) nbsub++;
  }

  std::cout << nbsub << " volumic substitution was made" << std::endl;
}

void create_vertices_for_templates(LCC &lcc, std::vector<Dart_handle>& marked_cells, size_type toprocess_mark, size_type arrete_done)
{

  // 2 noeuds marqué l'un à coté de l'autre ne produit pas de sommet
  // 1 noeud marqué a coté d'un noeud non marqué produit un sommet

  // TODO A changer pour une itération sur les arretes ou pas selon comment l'algo va s'executer   

  std::vector<Dart_handle> edges_to_subdivide;

  for (auto dart : marked_cells)
  {
    for (auto nit = lcc.one_dart_per_incident_cell<1, 0>(dart).begin(),
              nend = lcc.one_dart_per_incident_cell<1, 0>(dart).end();
         nit != nend;
         nit++)
    {
      if (lcc.is_marked(nit, arrete_done))
        continue;

      // If the node is next to an other marked node, we don't have to create vertices 
      if (lcc.is_marked(lcc.beta<1>(nit), toprocess_mark)){
        lcc.mark_cell<1>(nit, arrete_done);
        continue;
      }

      edges_to_subdivide.push_back(nit);
      lcc.mark_cell<1>(nit, arrete_done);
    }

    dart = lcc.beta<1>(dart);
  }

  for (Dart_handle dart : edges_to_subdivide)
  {
    lcc.insert_barycenter_in_cell<1>(dart);
  }
}

void mark_1template_face(LCC& lcc, size_type mark){
  auto dart = lcc.one_dart_per_cell<3>().begin();
  lcc.mark_cell<0>(lcc.beta(dart, 0, 2, 0, 0), mark); 
}

void mark_2template_face(LCC& lcc, size_type mark){
  auto dart = lcc.one_dart_per_cell<3>().begin();
  lcc.mark_cell<0>(dart, mark); 
  lcc.mark_cell<0>(lcc.beta(dart, 1), mark); 
}


int main()
{
  LCC lcc;
  Pattern_substituer<LCC> ps;
  auto d1=
    lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                        Point(5,5,0), Point(0,5,0),
                        Point(0,5,5), Point(0,0,5),
                        Point(5,0,5), Point(5,5,5));

  auto d2=
    lcc.make_hexahedron(Point(5,0,0), Point(10,0,0),
                        Point(10,5,0), Point(5,5,0),
                        Point(5,5,5), Point(5,0,5),
                        Point(10,0,5), Point(10,5,5));

  auto d3=
    lcc.make_hexahedron(Point(0,0,-5), Point(5,0,-5),
                        Point(5,5,-5), Point(0,5,-5),
                        Point(0,5,0), Point(0,0,0),
                        Point(5,0,0), Point(5,5,0));

  lcc.sew<3>(lcc.beta(d1, 1, 1, 2), lcc.beta(d2, 2));
  lcc.sew<3>(lcc.beta(d1, 0, 2), lcc.beta(d3, 1, 2));

  ps.load_vpatterns(CGAL::data_file_path("hexmeshing/vpattern"));

  // BUG Manque des opérateurs de déplacements pour Pattern, le reserve est un fix temporaire
  // Pour pouvoir charger les patterns correctement sans réallocation
  ps.m_fpatterns.reserve(10);
  ps.load_additional_fpattern(CGAL::data_file_path("hexmeshing/fpattern/pattern1-face.moka"), mark_1template_face);
  ps.load_additional_fpattern(CGAL::data_file_path("hexmeshing/fpattern/pattern2-face.moka"), mark_2template_face);
  ps.fpattern(0).display_characteristics(std::cout);
  // ps.fpattern(1).display_characteristics(std::cout);

  auto toprocess_mark = lcc.get_new_mark();
  auto arrete_done = lcc.get_new_mark();
  // auto corner_mark = lcc.get_new_mark();
  // lcc.negate_mark(corner_mark);

  lcc.mark_cell<0>(d1, toprocess_mark);
  auto md1 = lcc.beta(d1, 1);
  lcc.mark_cell<0>(md1, toprocess_mark);
  auto md2 = lcc.beta(d2, 1);
  lcc.mark_cell<0>(md2, toprocess_mark);
  auto md3 = lcc.beta(d2, 1, 1);
  lcc.mark_cell<0>(md3, toprocess_mark);
  auto md4 = d3;
  lcc.mark_cell<0>(md4, toprocess_mark);
  
  // lcc.mark_cell<0>(md3, toprocess_mark);

  std::vector<Dart_handle> marked_cells; 
  
  // We assume we already filled out marked_cells earlier

  marked_cells.push_back(d1);
  marked_cells.push_back(md1);
  marked_cells.push_back(md2);
  marked_cells.push_back(md3);
  marked_cells.push_back(md4);

  create_vertices_for_templates(lcc, marked_cells, toprocess_mark, arrete_done);
  refine_hexes_with_1_2_4_templates(lcc, ps, toprocess_mark);

  render(lcc, toprocess_mark);
  return EXIT_SUCCESS;
}
