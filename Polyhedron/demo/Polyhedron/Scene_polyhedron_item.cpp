#include "Scene_polyhedron_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/AABB_intersections.h>
#include "Kernel_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_writer_wavefront.h>
#include <CGAL/IO/generic_copy_OFF.h>
#include <CGAL/IO/OBJ_reader.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/property_map.h>
#include <CGAL/statistics_helpers.h>

#include <list>
#include <queue>
#include <iostream>
#include <limits>

#include <QVariant>
#include <QDebug>
#include <QDialog>
#include <QApplication>

#include <boost/foreach.hpp>
#include "triangulate_primitive.h"

namespace PMP = CGAL::Polygon_mesh_processing;
typedef Polyhedron::Traits Traits;
typedef Polyhedron::Facet Facet;
typedef Polyhedron::Traits	                           Kernel;
typedef Kernel::Point_3	                                   Point;
typedef Kernel::Vector_3	                           Vector;
typedef Polyhedron::Halfedge_around_facet_circulator       HF_circulator;
typedef boost::graph_traits<Polyhedron>::face_descriptor   face_descriptor;
//Used to triangulate the AABB_Tree
class Primitive
{
public:
  // types
  typedef Polyhedron::Facet_iterator Id; // Id type
  typedef Kernel::Point_3 Point; // point type
  typedef Kernel::Triangle_3 Datum; // datum type

private:
  // member data
  Id m_it; // iterator
  Datum m_datum; // 3D triangle

  // constructor
public:
  Primitive() {}
  Primitive(Datum triangle, Id it)
    : m_it(it), m_datum(triangle)
  {
  }
public:
  Id& id() { return m_it; }
  const Id& id() const { return m_it; }
  Datum& datum() { return m_datum; }
  const Datum& datum() const { return m_datum; }

  /// Returns a point on the primitive
  Point reference_point() const { return m_datum.vertex(0); }
};

typedef CGAL::AABB_traits<Kernel, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Input_facets_AABB_tree;

struct Scene_polyhedron_item_priv{
  typedef Polyhedron::Facet_iterator Facet_iterator;
  typedef std::vector<QColor> Color_vector;
  Scene_polyhedron_item_priv(Scene_polyhedron_item* item)
    : item(item), poly(new Polyhedron)
  {
    init_default_values();
  }
  Scene_polyhedron_item_priv(const Polyhedron& poly_, Scene_polyhedron_item* item)
    : item(item), poly(new Polyhedron(poly_))
  {
    init_default_values();
  }

  Scene_polyhedron_item_priv(Polyhedron* const poly_, Scene_polyhedron_item* item)
    : item(item), poly(poly_)
  {
    init_default_values();
  }

  void init_default_values() {
    show_only_feature_edges_m = false;
    show_feature_edges_m = false;
    facet_picking_m = false;
    erase_next_picked_facet_m = false;
    plugin_has_set_color_vector_m = false;
    nb_facets = 0;
    nb_lines = 0;
    nb_f_lines = 0;
    is_multicolor = false;
    targeted_id = NULL;
    all_ids_displayed = false;
    invalidate_stats();
  }


  void compute_normals_and_vertices(const bool colors_only = false) const;
  bool isFacetConvex(Facet_iterator, const Polyhedron::Traits::Vector_3&)const;
  template<typename VertexNormalPmap>
  void triangulate_convex_facet(Facet_iterator,
                                const Polyhedron::Traits::Vector_3&,
                                const VertexNormalPmap&,
                                const bool)const;
  template<typename VertexNormalPmap>
  void triangulate_facet(Scene_polyhedron_item::Facet_iterator,
                         const Traits::Vector_3& normal,
                         const VertexNormalPmap&,
                         const bool colors_only) const;
  void init();
  void invalidate_stats();
  void destroy()
  {
    delete poly;
    delete targeted_id;
  }
  void* get_aabb_tree();
  QList<Kernel::Triangle_3> triangulate_primitive(Polyhedron::Facet_iterator fit,
                                                  Traits::Vector_3 normal);
  Color_vector colors_;
  bool show_only_feature_edges_m;
  bool show_feature_edges_m;
  bool facet_picking_m;
  bool erase_next_picked_facet_m;
  //the following variable is used to indicate if the color vector must not be automatically updated.
  // If this boolean is true, the item's color will be independant of the color's wheel, if it is false
  // changing the color in the wheel will change the color of the item, even if it is multicolor.
  bool plugin_has_set_color_vector_m;
  bool is_multicolor;

  Scene_polyhedron_item* item;
  Polyhedron *poly;
  double volume, area;
  QVector<QColor> colors;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_feature_lines;
  mutable std::vector<float> positions_facets;
  mutable std::vector<float> normals_flat;
  mutable std::vector<float> normals_gouraud;
  mutable std::vector<float> color_lines;
  mutable std::vector<float> color_facets;
  mutable std::size_t nb_facets;
  mutable std::size_t nb_lines;
  mutable std::size_t nb_f_lines;
  mutable QOpenGLShaderProgram *program;
  unsigned int number_of_null_length_edges;
  unsigned int number_of_degenerated_faces;
  bool self_intersect;
  int m_min_patch_id; // the min value of the patch ids initialized in init()
  mutable bool all_ids_displayed;
  mutable QList<double> text_ids;
  mutable TextItem* targeted_id;
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer = 0) const;
  enum VAOs {
    Facets=0,
    Edges,
    Feature_edges,
    Gouraud_Facets,
    NbOfVaos
  };
  enum VBOs {
    Facets_vertices = 0,
    Facets_normals_flat,
    Facets_color,
    Edges_vertices,
    Feature_edges_vertices,
    Edges_color,
    Facets_normals_gouraud,
    NbOfVbos
  };
  // Initialization
};



const char* aabb_property_name = "Scene_polyhedron_item aabb tree";


QList<Kernel::Triangle_3> Scene_polyhedron_item_priv::triangulate_primitive(Polyhedron::Facet_iterator fit,
                                                Traits::Vector_3 normal)
{
  typedef FacetTriangulator<Polyhedron, Polyhedron::Traits, boost::graph_traits<Polyhedron>::vertex_descriptor> FT;
  //The output list
  QList<Kernel::Triangle_3> res;
  //check if normal contains NaN values
  if (normal.x() != normal.x() || normal.y() != normal.y() || normal.z() != normal.z())
  {
    qDebug()<<"Warning in triangulation of the selection item: normal contains NaN values and is not valid.";
    return QList<Kernel::Triangle_3>();
  }
  double diagonal;
  if(item->diagonalBbox() != std::numeric_limits<double>::infinity())
    diagonal = item->diagonalBbox();
  else
    diagonal = 0.0;
  FT triangulation(fit,normal,poly,diagonal);
  //iterates on the internal faces to add the vertices to the positions
  //and the normals to the appropriate vectors
  for( FT::CDT::Finite_faces_iterator
      ffit = triangulation.cdt->finite_faces_begin(),
      end = triangulation.cdt->finite_faces_end();
      ffit != end; ++ffit)
  {
    if(ffit->info().is_external)
      continue;


    res << Kernel::Triangle_3(ffit->vertex(0)->point(),
                              ffit->vertex(1)->point(),
                              ffit->vertex(2)->point());

  }
  return res;
}



void* Scene_polyhedron_item_priv::get_aabb_tree()
{
  QVariant aabb_tree_property = item->property(aabb_property_name);
  if(aabb_tree_property.isValid()) {
    void* ptr = aabb_tree_property.value<void*>();
    return static_cast<Input_facets_AABB_tree*>(ptr);
  }
  else {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Polyhedron* poly = item->polyhedron();
    if(poly) {

      Input_facets_AABB_tree* tree =
          new Input_facets_AABB_tree();
      typedef Polyhedron::Traits	    Kernel;
      int index =0;
      Q_FOREACH( Polyhedron::Facet_iterator f, faces(*poly))
      {
        if(!f->is_triangle())
        {
          Traits::Vector_3 normal = f->plane().orthogonal_vector(); //initialized in compute_normals_and_vertices
          index +=3;
          Q_FOREACH(Kernel::Triangle_3 triangle, triangulate_primitive(f,normal))
          {
            Primitive primitive(triangle, f);
            tree->insert(primitive);
          }
        }
        else
        {
          Kernel::Triangle_3 triangle(
                f->halfedge()->vertex()->point(),
                f->halfedge()->next()->vertex()->point(),
                f->halfedge()->next()->next()->vertex()->point()
                );
          Primitive primitive(triangle, f);
          tree->insert(primitive);
        }
      }
      item->setProperty(aabb_property_name,
                        QVariant::fromValue<void*>(tree));
      QApplication::restoreOverrideCursor();
      return tree;
    }
    else return 0;
  }
}

void delete_aabb_tree(Scene_polyhedron_item* item)
{
    QVariant aabb_tree_property = item->property(aabb_property_name);
    if(aabb_tree_property.isValid()) {
        void* ptr = aabb_tree_property.value<void*>();
        Input_facets_AABB_tree* tree = static_cast<Input_facets_AABB_tree*>(ptr);
        if(tree) {
            delete tree;
            tree = 0;
        }
        item->setProperty(aabb_property_name, QVariant());
    }
}

template<typename TypeWithXYZ, typename ContainerWithPushBack>
void push_back_xyz(const TypeWithXYZ& t,
                   ContainerWithPushBack& vector)
{
  vector.push_back(t.x());
  vector.push_back(t.y());
  vector.push_back(t.z());
}


template<typename TypeWithRGB, typename ContainerWithPushBack>
void push_back_rgb(const TypeWithRGB& t,
                   ContainerWithPushBack& vector)
{
  vector.push_back(t.redF());
  vector.push_back(t.greenF());
  vector.push_back(t.blueF());
}

bool Scene_polyhedron_item_priv::isFacetConvex(Facet_iterator f, const Polyhedron::Traits::Vector_3& N) const
{
  typedef Polyhedron::Traits::Vector_3 Vector;
  typedef Polyhedron::Traits::Orientation Orientation;
  Orientation orientation;
  Vector normal = N;
  Facet::Halfedge_around_facet_circulator
          he = f->facet_begin(),
          he_end(he);
  bool normal_is_ok;
  do{
    normal_is_ok = true;

    //Initializes the facet orientation

    Polyhedron::Traits::Point_3 S,T;
    S = source(he, *poly)->point();
    T = target(he, *poly)->point();
    Vector V1 = Vector((T-S).x(), (T-S).y(), (T-S).z());
    S = source(he->next(), *poly)->point();
    T = target(he->next(), *poly)->point();
    Vector V2 = Vector((T-S).x(), (T-S).y(), (T-S).z());

    if(normal == Vector(0,0,0))
      normal_is_ok = false;
    if(normal_is_ok)
    {
      orientation = Polyhedron::Traits::Orientation_3()(V1, V2, normal);
      if( orientation == CGAL::COPLANAR )
        normal_is_ok = false;
    }
    //Checks if the normal is good : if the normal is null
    // or if it is coplanar to the facet, we need another one.
    if(!normal_is_ok)
    {
      normal = CGAL::cross_product(V1,V2);
    }

  }while( ++he != he_end && !normal_is_ok);
  //if no good normal can be found, stop here.
  if (!normal_is_ok)
    return false;

  //computes convexness

  //re-initializes he_end;
  he = f->facet_begin(),
  he_end = he;
  do
  {
    Polyhedron::Traits::Point_3 S,T;
    S = source(he, *poly)->point();
    T = target(he, *poly)->point();
    Vector V1 = Vector((T-S).x(), (T-S).y(), (T-S).z());
    S = source(he->next(), *poly)->point();
    T = target(he->next(), *poly)->point();
    Vector V2 = Vector((T-S).x(), (T-S).y(), (T-S).z());
    Orientation res = Polyhedron::Traits::Orientation_3()(V1, V2, normal) ;

    if(res!= orientation && res != CGAL::ZERO)
      return false;
  }while( ++he != he_end);
  return true;
}

template<typename VertexNormalPmap>
void Scene_polyhedron_item_priv::triangulate_convex_facet(Facet_iterator f,
                                                     const Traits::Vector_3& normal,
                                                     const VertexNormalPmap& vnmap,
                                                     const bool colors_only)const
{
  Polyhedron::Traits::Point_3 p0,p1,p2;
  Facet::Halfedge_around_facet_circulator
      he = f->facet_begin(),
      he_end(he);
  const int this_patch_id = f->patch_id();
  while( next(he, *poly) != prev(he_end, *poly))
  {
    ++he;
    if (is_multicolor)
    {
      for (int i = 0; i<3; ++i)
      {
        color_facets.push_back(colors_[this_patch_id-m_min_patch_id].redF());
        color_facets.push_back(colors_[this_patch_id-m_min_patch_id].greenF());
        color_facets.push_back(colors_[this_patch_id-m_min_patch_id].blueF());

        color_facets.push_back(colors_[this_patch_id-m_min_patch_id].redF());
        color_facets.push_back(colors_[this_patch_id-m_min_patch_id].greenF());
        color_facets.push_back(colors_[this_patch_id-m_min_patch_id].blueF());
      }
    }
    if (colors_only)
      continue;


    p0 = he_end->vertex()->point();
    p1 = he->vertex()->point();
    p2 = next(he, *poly)->vertex()->point();
    push_back_xyz(p0, positions_facets);
    positions_facets.push_back(1.0);
    push_back_xyz(p1, positions_facets);
    positions_facets.push_back(1.0);
    push_back_xyz(p2, positions_facets);
    positions_facets.push_back(1.0);

    push_back_xyz(normal, normals_flat);
    push_back_xyz(normal, normals_flat);
    push_back_xyz(normal, normals_flat);

    Traits::Vector_3 ng = get(vnmap, he_end->vertex());
    push_back_xyz(ng, normals_gouraud);

    ng = get(vnmap, he->vertex());
    push_back_xyz(ng, normals_gouraud);

    ng = get(vnmap, next(he, *poly)->vertex());
    push_back_xyz(ng, normals_gouraud);
  }


}
//Make sure all the facets are triangles
template<typename VertexNormalPmap>
void
Scene_polyhedron_item_priv::triangulate_facet(Scene_polyhedron_item::Facet_iterator fit,
                                         const Traits::Vector_3& normal,
                                         const VertexNormalPmap& vnmap,
                                         const bool colors_only) const
{
  typedef FacetTriangulator<Polyhedron, Polyhedron::Traits, boost::graph_traits<Polyhedron>::vertex_descriptor> FT;
  double diagonal;
  if(item->diagonalBbox() != std::numeric_limits<double>::infinity())
    diagonal = item->diagonalBbox();
  else
    diagonal = 0.0;
  FT triangulation(fit,normal,poly,diagonal);

  if(triangulation.cdt->dimension() != 2 )
  {
    qDebug()<<"Warning : cdt not right. Facet not displayed";
    return;
  }

  //iterates on the internal faces to add the vertices to the positions
  //and the normals to the appropriate vectors
  const int this_patch_id = fit->patch_id();

  for(FT::CDT::Finite_faces_iterator
      ffit = triangulation.cdt->finite_faces_begin(),
      end = triangulation.cdt->finite_faces_end();
      ffit != end; ++ffit)
  {
    if(ffit->info().is_external)
      continue;

    if (is_multicolor)
    {
     for (int i = 0; i<3; ++i)
     {
      push_back_rgb(colors_[this_patch_id-m_min_patch_id], color_facets);
     }
    }
    if (colors_only)
     continue;

    push_back_xyz(ffit->vertex(0)->point(), positions_facets);
    positions_facets.push_back(1.0);

    push_back_xyz(ffit->vertex(1)->point(), positions_facets);
    positions_facets.push_back(1.0);

    push_back_xyz(ffit->vertex(2)->point(), positions_facets);
    positions_facets.push_back(1.0);

    push_back_xyz(normal, normals_flat);
    push_back_xyz(normal, normals_flat);
    push_back_xyz(normal, normals_flat);

    Traits::Vector_3 ng = get(vnmap, triangulation.v2v[ffit->vertex(0)]);
    push_back_xyz(ng, normals_gouraud);

    ng = get(vnmap, triangulation.v2v[ffit->vertex(1)]);
    push_back_xyz(ng, normals_gouraud);

    ng = get(vnmap, triangulation.v2v[ffit->vertex(2)]);
    push_back_xyz(ng, normals_gouraud);
  }
}


#include <QObject>
#include <QMenu>
#include <QAction>


void
Scene_polyhedron_item_priv::initialize_buffers(CGAL::Three::Viewer_interface* viewer) const
{
    //vao containing the data for the facets
    {
        program = item->getShaderProgram(Scene_polyhedron_item::PROGRAM_WITH_LIGHT, viewer);
        program->bind();
        //flat
        item->vaos[Facets]->bind();
        item->buffers[Facets_vertices].bind();
        item->buffers[Facets_vertices].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[Facets_vertices].release();



        item->buffers[Facets_normals_flat].bind();
        item->buffers[Facets_normals_flat].allocate(normals_flat.data(),
                            static_cast<int>(normals_flat.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        item->buffers[Facets_normals_flat].release();

        if(is_multicolor)
        {
            item->buffers[Facets_color].bind();
            item->buffers[Facets_color].allocate(color_facets.data(),
                                static_cast<int>(color_facets.size()*sizeof(float)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            item->buffers[Facets_color].release();
        }
        else
        {
          program->disableAttributeArray("colors");
        }
        item->vaos[Facets]->release();
        //gouraud
        item->vaos[Gouraud_Facets]->bind();
        item->buffers[Facets_vertices].bind();
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[Facets_vertices].release();

        item->buffers[Facets_normals_gouraud].bind();
        item->buffers[Facets_normals_gouraud].allocate(normals_gouraud.data(),
                            static_cast<int>(normals_gouraud.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        item->buffers[Facets_normals_gouraud].release();
        if(is_multicolor)
        {
            item->buffers[Facets_color].bind();
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            item->buffers[Facets_color].release();
        }
        else
        {
            program->disableAttributeArray("colors");
        }
        item->vaos[Gouraud_Facets]->release();

        program->release();

    }
    //vao containing the data for the lines
    {
        program = item->getShaderProgram(Scene_polyhedron_item::PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        item->vaos[Edges]->bind();

        item->buffers[Edges_vertices].bind();
        item->buffers[Edges_vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[Edges_vertices].release();

       if(is_multicolor)
       {
         item->buffers[Edges_color].bind();
         item->buffers[Edges_color].allocate(color_lines.data(),
                            static_cast<int>(color_lines.size()*sizeof(float)));
         program->enableAttributeArray("colors");
         program->setAttributeBuffer("colors",GL_FLOAT,0,3);
         item->buffers[Edges_color].release();
       }
       else
       {
           program->disableAttributeArray("colors");
       }
        program->release();

        item->vaos[Edges]->release();

    }
  //vao containing the data for the feature_edges
  {
      program = item->getShaderProgram(Scene_polyhedron_item::PROGRAM_NO_SELECTION, viewer);
      program->bind();
      item->vaos[Feature_edges]->bind();

      item->buffers[Feature_edges_vertices].bind();
      item->buffers[Feature_edges_vertices].allocate(positions_feature_lines.data(),
                          static_cast<int>(positions_feature_lines.size()*sizeof(float)));
      program->enableAttributeArray("vertex");
      program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
      item->buffers[Feature_edges_vertices].release();
      program->disableAttributeArray("colors");
      program->release();

      item->vaos[Feature_edges]->release();

  }
    nb_f_lines = positions_feature_lines.size();
    positions_feature_lines.resize(0);
    std::vector<float>(positions_feature_lines).swap(positions_feature_lines);
    nb_lines = positions_lines.size();
    positions_lines.resize(0);
    std::vector<float>(positions_lines).swap(positions_lines);
    nb_facets = positions_facets.size();
    positions_facets.resize(0);
    std::vector<float>(positions_facets).swap(positions_facets);


    color_lines.resize(0);
    std::vector<float>(color_lines).swap(color_lines);
    color_facets.resize(0);
    std::vector<float>(color_facets).swap(color_facets);
    normals_flat.resize(0);
    std::vector<float>(normals_flat).swap(normals_flat);
    normals_gouraud.resize(0);
    std::vector<float>(normals_gouraud).swap(normals_gouraud);

    CGAL::set_halfedgeds_items_id(*poly);
    if (viewer->hasText())
        item->printPrimitiveIds(viewer);
    item->are_buffers_filled = true;

}

void
Scene_polyhedron_item_priv::compute_normals_and_vertices(const bool colors_only) const
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    positions_facets.resize(0);
    positions_lines.resize(0);
    positions_feature_lines.resize(0);
    normals_flat.resize(0);
    normals_gouraud.resize(0);
    color_lines.resize(0);
    color_facets.resize(0);

    //Facets
    typedef Polyhedron::Traits	    Kernel;
    typedef Kernel::Point_3	    Point;
    typedef Kernel::Vector_3	    Vector;
    typedef Polyhedron::Facet_iterator Facet_iterator;
    typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;
    typedef boost::graph_traits<Polyhedron>::face_descriptor   face_descriptor;
    typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

    CGAL::Unique_hash_map<face_descriptor, Vector> face_normals_map;
    boost::associative_property_map<CGAL::Unique_hash_map<face_descriptor, Vector> >
      nf_pmap(face_normals_map);
    CGAL::Unique_hash_map<vertex_descriptor, Vector> vertex_normals_map;
    boost::associative_property_map< CGAL::Unique_hash_map<vertex_descriptor, Vector> >
      nv_pmap(vertex_normals_map);

    PMP::compute_normals(*poly, nv_pmap, nf_pmap);

    Facet_iterator f = poly->facets_begin();
    for(f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {
      if (f == boost::graph_traits<Polyhedron>::null_face())
        continue;
      Vector nf = get(nf_pmap, f);
      f->plane() = Kernel::Plane_3(f->halfedge()->vertex()->point(), nf);
      if(is_triangle(f->halfedge(),*poly))
      {
          const int this_patch_id = f->patch_id();
          HF_circulator he = f->facet_begin();
          HF_circulator end = he;
          CGAL_For_all(he,end)
          {
            if (item->isItemMulticolor())
            {
              color_facets.push_back(colors_[this_patch_id-m_min_patch_id].redF());
              color_facets.push_back(colors_[this_patch_id-m_min_patch_id].greenF());
              color_facets.push_back(colors_[this_patch_id-m_min_patch_id].blueF());
            }
            if (colors_only)
              continue;

            // If Flat shading:1 normal per polygon added once per vertex
            push_back_xyz(nf, normals_flat);

            //// If Gouraud shading: 1 normal per vertex
            Vector nv = get(nv_pmap, he->vertex());
            push_back_xyz(nv, normals_gouraud);

            //position
            const Point& p = he->vertex()->point();
            push_back_xyz(p, positions_facets);
            positions_facets.push_back(1.0);
         }
      }
      else if (is_quad(f->halfedge(), *poly))
      {
        if (item->isItemMulticolor())
        {
          const int this_patch_id = f->patch_id();
          for (unsigned int i = 0; i < 6; ++i)
          { //6 "halfedges" for the quad, because it is 2 triangles
            color_facets.push_back(colors_[this_patch_id-m_min_patch_id].redF());
            color_facets.push_back(colors_[this_patch_id-m_min_patch_id].greenF());
            color_facets.push_back(colors_[this_patch_id-m_min_patch_id].blueF());
          }
        }
        if (colors_only)
          continue;

        //1st half-quad
        Point p0 = f->halfedge()->vertex()->point();
        Point p1 = f->halfedge()->next()->vertex()->point();
        Point p2 = f->halfedge()->next()->next()->vertex()->point();

        push_back_xyz(p0, positions_facets);
        positions_facets.push_back(1.0);

        push_back_xyz(p1, positions_facets);
        positions_facets.push_back(1.0);

        push_back_xyz(p2, positions_facets);
        positions_facets.push_back(1.0);

        push_back_xyz(nf, normals_flat);
        push_back_xyz(nf, normals_flat);
        push_back_xyz(nf, normals_flat);

        Vector nv = get(nv_pmap, f->halfedge()->vertex());
        push_back_xyz(nv, normals_gouraud);

        nv = get(nv_pmap, f->halfedge()->next()->vertex());
        push_back_xyz(nv, normals_gouraud);

        nv = get(nv_pmap, f->halfedge()->next()->next()->vertex());
        push_back_xyz(nv, normals_gouraud);

        //2nd half-quad
        p0 = f->halfedge()->next()->next()->vertex()->point();
        p1 = f->halfedge()->prev()->vertex()->point();
        p2 = f->halfedge()->vertex()->point();

        push_back_xyz(p0, positions_facets);
        positions_facets.push_back(1.0);

        push_back_xyz(p1, positions_facets);
        positions_facets.push_back(1.0);

        push_back_xyz(p2, positions_facets);
        positions_facets.push_back(1.0);

        push_back_xyz(nf, normals_flat);
        push_back_xyz(nf, normals_flat);
        push_back_xyz(nf, normals_flat);

        nv = get(nv_pmap, f->halfedge()->next()->next()->vertex());
        push_back_xyz(nv, normals_gouraud);

        nv = get(nv_pmap, f->halfedge()->prev()->vertex());
        push_back_xyz(nv, normals_gouraud);

        nv = get(nv_pmap, f->halfedge()->vertex());
        push_back_xyz(nv, normals_gouraud);
      }
      else
      {
        if(isFacetConvex(f, nf))
        {
          triangulate_convex_facet(f, nf, nv_pmap, colors_only);
        }
        else
        {
          this->triangulate_facet(f, nf, nv_pmap, colors_only);
        }
      }

    }
    //Lines
    typedef Kernel::Point_3		Point;
    typedef Polyhedron::Edge_iterator	Edge_iterator;
    Edge_iterator he;
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        if ( he->is_feature_edge())
        {
          if (colors_only)
            continue;

          push_back_xyz(a, positions_feature_lines);
          positions_feature_lines.push_back(1.0);

          push_back_xyz(b, positions_feature_lines);
          positions_feature_lines.push_back(1.0);
        }
        else
        {
          if (item->isItemMulticolor())
          {
            color_lines.push_back(item->color().lighter(50).redF());
            color_lines.push_back(item->color().lighter(50).greenF());
            color_lines.push_back(item->color().lighter(50).blueF());

            color_lines.push_back(item->color().lighter(50).redF());
            color_lines.push_back(item->color().lighter(50).greenF());
            color_lines.push_back(item->color().lighter(50).blueF());
          }
          if (colors_only)
            continue;

          push_back_xyz(a, positions_lines);
          positions_lines.push_back(1.0);

          push_back_xyz(b, positions_lines);
          positions_lines.push_back(1.0);
        }
    }
    QApplication::restoreOverrideCursor();
}

Scene_polyhedron_item::Scene_polyhedron_item()
    : Scene_item(Scene_polyhedron_item_priv::NbOfVbos,Scene_polyhedron_item_priv::NbOfVaos),
      d(new Scene_polyhedron_item_priv(this))
{
    cur_shading=FlatPlusEdges;
    is_selected = true;
    textItems = new TextListItem(this);

}

Scene_polyhedron_item::Scene_polyhedron_item(Polyhedron* const p)
    : Scene_item(Scene_polyhedron_item_priv::NbOfVbos,Scene_polyhedron_item_priv::NbOfVaos),
      d(new Scene_polyhedron_item_priv(p,this))
{
    cur_shading=FlatPlusEdges;
    is_selected = true;
    textItems = new TextListItem(this);
}

Scene_polyhedron_item::Scene_polyhedron_item(const Polyhedron& p)
    : Scene_item(Scene_polyhedron_item_priv::NbOfVbos,Scene_polyhedron_item_priv::NbOfVaos),
      d(new Scene_polyhedron_item_priv(p,this))
{
    cur_shading=FlatPlusEdges;
    is_selected=true;
    textItems = new TextListItem(this);
}

Scene_polyhedron_item::~Scene_polyhedron_item()
{
    delete_aabb_tree(this);
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    if(viewer)
    {
      CGAL::Three::Viewer_interface* v = qobject_cast<CGAL::Three::Viewer_interface*>(viewer);

      //Clears the targeted Id
      v->textRenderer->removeText(d->targeted_id);
      //Remove textitems
      v->textRenderer->removeTextList(textItems);
      delete textItems;
    }

    d->destroy();
    delete d;
}

#include "Color_map.h"

void
Scene_polyhedron_item_priv::
init()
{
  typedef Polyhedron::Facet_iterator Facet_iterator;

  if ( !plugin_has_set_color_vector_m )
  {
    // Fill indices map and get max subdomain value
    int max = 0;
    int min = (std::numeric_limits<int>::max)();
    for(Facet_iterator fit = poly->facets_begin(), end = poly->facets_end() ;
        fit != end; ++fit)
    {
      max = (std::max)(max, fit->patch_id());
      min = (std::min)(min, fit->patch_id());
    }
    
    colors_.clear();
    compute_color_map(item->color(), (std::max)(0, max + 1 - min),
                      std::back_inserter(colors_));
    m_min_patch_id=min;
  }
  else
    m_min_patch_id=0;
  invalidate_stats();
}

void
Scene_polyhedron_item_priv::
invalidate_stats()
{
  number_of_degenerated_faces = (unsigned int)(-1);
  number_of_null_length_edges = (unsigned int)(-1);
  volume = -std::numeric_limits<double>::infinity();
  area = -std::numeric_limits<double>::infinity();
  self_intersect = false;

}

Scene_polyhedron_item*
Scene_polyhedron_item::clone() const {
    return new Scene_polyhedron_item(*(d->poly));}

// Load polyhedron from .OFF file
bool
Scene_polyhedron_item::load(std::istream& in)
{


    in >> *(d->poly);

    if ( in && !isEmpty() )
    {
        invalidateOpenGLBuffers();
        return true;
    }
    return false;
}
// Load polyhedron from .obj file
bool
Scene_polyhedron_item::load_obj(std::istream& in)
{
  typedef Polyhedron::Vertex::Point Point;
  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > faces;
  bool failed = !CGAL::read_OBJ(in,points,faces);

  if(CGAL::Polygon_mesh_processing::orient_polygon_soup(points,faces)){
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points,faces,*(d->poly));
  }else{
    std::cerr << "not orientable"<< std::endl;
    return false;
  }
    if ( (! failed) && !isEmpty() )
    {
        invalidateOpenGLBuffers();
        return true;
    }
    return false;
}

// Write polyhedron to .OFF file
bool
Scene_polyhedron_item::save(std::ostream& out) const
{
  out.precision(17);
    out << *(d->poly);
    return (bool) out;
}

bool
Scene_polyhedron_item::save_obj(std::ostream& out) const
{
  CGAL::File_writer_wavefront  writer;
  CGAL::generic_print_polyhedron(out, *(d->poly), writer);
  return out.good();
}


QString
Scene_polyhedron_item::toolTip() const
{
    if(!d->poly)
        return QString();

  QString str =
         QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                     "Number of facets: %4")
            .arg(this->name())
            .arg(d->poly->size_of_vertices())
            .arg(d->poly->size_of_halfedges()/2)
            .arg(d->poly->size_of_facets())
            .arg(this->renderingModeName())
            .arg(this->color().name());
  str += QString("<br />Number of isolated vertices : %1<br />").arg(getNbIsolatedvertices());
  return str;
}

QMenu* Scene_polyhedron_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_polyhedron_item.";

  QMenu* menu = Scene_item::contextMenu();

  // Use dynamic properties:
  // http://doc.qt.io/qt-5/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if(!menuChanged) {

    QAction* actionShowOnlyFeatureEdges =
        menu->addAction(tr("Show Only &Feature Edges"));
    actionShowOnlyFeatureEdges->setCheckable(true);
    actionShowOnlyFeatureEdges->setChecked(d->show_only_feature_edges_m);
    actionShowOnlyFeatureEdges->setObjectName("actionShowOnlyFeatureEdges");
    connect(actionShowOnlyFeatureEdges, SIGNAL(toggled(bool)),
            this, SLOT(show_only_feature_edges(bool)));

    QAction* actionShowFeatureEdges =
        menu->addAction(tr("Show Feature Edges"));
    actionShowFeatureEdges->setCheckable(true);
    actionShowFeatureEdges->setChecked(d->show_feature_edges_m);
    actionShowFeatureEdges->setObjectName("actionShowFeatureEdges");
    connect(actionShowFeatureEdges, SIGNAL(toggled(bool)),
            this, SLOT(show_feature_edges(bool)));

    QAction* actionPickFacets =
        menu->addAction(tr("Facets Picking"));
    actionPickFacets->setCheckable(true);
    actionPickFacets->setObjectName("actionPickFacets");
    connect(actionPickFacets, SIGNAL(toggled(bool)),
            this, SLOT(enable_facets_picking(bool)));

    QAction* actionEraseNextFacet =
        menu->addAction(tr("Erase Next Picked Facet"));
    actionEraseNextFacet->setCheckable(true);
    actionEraseNextFacet->setObjectName("actionEraseNextFacet");
    connect(actionEraseNextFacet, SIGNAL(toggled(bool)),
            this, SLOT(set_erase_next_picked_facet(bool)));
    menu->setProperty(prop_name, true);
  }

  QAction* action = menu->findChild<QAction*>("actionShowOnlyFeatureEdges");
  if(action) action->setChecked(d->show_only_feature_edges_m);
  action = menu->findChild<QAction*>("actionShowFeatureEdges");
  if(action) action->setChecked(d->show_feature_edges_m);
  action = menu->findChild<QAction*>("actionPickFacets");
  if(action) action->setChecked(d->facet_picking_m);
  action = menu->findChild<QAction*>("actionEraseNextFacet");
  if(action) action->setChecked(d->erase_next_picked_facet_m);
  return menu;
}

void Scene_polyhedron_item::show_only_feature_edges(bool b)
{
    d->show_only_feature_edges_m = b;
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
}

void Scene_polyhedron_item::show_feature_edges(bool b)
{
  d->show_feature_edges_m = b;
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

void Scene_polyhedron_item::enable_facets_picking(bool b)
{
    d->facet_picking_m = b;
}

void Scene_polyhedron_item::set_erase_next_picked_facet(bool b)
{
    if(b) { d->facet_picking_m = true; } // automatically activate facet_picking
    d->erase_next_picked_facet_m = b;
}

void Scene_polyhedron_item::draw(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        d->compute_normals_and_vertices();
        d->initialize_buffers(viewer);
        compute_bbox();
    }

    if(renderingMode() == Flat || renderingMode() == FlatPlusEdges)
        vaos[Scene_polyhedron_item_priv::Facets]->bind();
    else
    {
        vaos[Scene_polyhedron_item_priv::Gouraud_Facets]->bind();
    }
    attribBuffers(viewer, PROGRAM_WITH_LIGHT);
    d->program = getShaderProgram(PROGRAM_WITH_LIGHT);
    d->program->bind();
    if(!d->is_multicolor)
    {
            d->program->setAttributeValue("colors", this->color());
    }
    if(is_selected)
            d->program->setUniformValue("is_selected", true);
    else
            d->program->setUniformValue("is_selected", false);
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->nb_facets/4));
    d->program->release();
    if(renderingMode() == Flat || renderingMode() == FlatPlusEdges)
        vaos[Scene_polyhedron_item_priv::Facets]->release();
    else
        vaos[Scene_polyhedron_item_priv::Gouraud_Facets]->release();
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
    if (!are_buffers_filled)
    {
        d->compute_normals_and_vertices();
        d->initialize_buffers(viewer);
        compute_bbox();
    }

    if(!d->show_only_feature_edges_m)
    {
        vaos[Scene_polyhedron_item_priv::Edges]->bind();

        attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
        d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        d->program->bind();
        //draw the edges
        if(!d->is_multicolor)
        {
            d->program->setAttributeValue("colors", this->color().lighter(50));
            if(is_selected)
                d->program->setUniformValue("is_selected", true);
            else
                d->program->setUniformValue("is_selected", false);
        }
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->nb_lines/4));
        d->program->release();
        vaos[Scene_polyhedron_item_priv::Edges]->release();
    }

    //draw the feature edges
    vaos[Scene_polyhedron_item_priv::Feature_edges]->bind();
    attribBuffers(viewer, PROGRAM_NO_SELECTION);
    d->program = getShaderProgram(PROGRAM_NO_SELECTION);
    d->program->bind();
    if(d->show_feature_edges_m || d->show_only_feature_edges_m)
        d->program->setAttributeValue("colors", Qt::red);
    else
    {
        if(!is_selected)
            d->program->setAttributeValue("colors", this->color().lighter(50));
        else
            d->program->setAttributeValue("colors",QColor(0,0,0));
    }
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->nb_f_lines/4));
    d->program->release();
    vaos[Scene_polyhedron_item_priv::Feature_edges]->release();
    }

void
Scene_polyhedron_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        d->compute_normals_and_vertices();
        d->initialize_buffers(viewer);
        compute_bbox();
    }

    vaos[Scene_polyhedron_item_priv::Edges]->bind();
    attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
    d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    d->program->bind();
    //draw the points
    d->program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(d->nb_lines/4));
    // Clean-up
    d->program->release();
    vaos[Scene_polyhedron_item_priv::Edges]->release();
}

Polyhedron*
Scene_polyhedron_item::polyhedron()       { return d->poly; }
const Polyhedron*
Scene_polyhedron_item::polyhedron() const { return d->poly; }

bool
Scene_polyhedron_item::isEmpty() const {
    return (d->poly == 0) || d->poly->empty();
}

void Scene_polyhedron_item::compute_bbox() const {
    const Kernel::Point_3& p = *(d->poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedron::Point_iterator it = d->poly->points_begin();
        it != d->poly->points_end();
        ++it) {
        bbox = bbox + it->bbox();
    }
    _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}


void
Scene_polyhedron_item::
invalidateOpenGLBuffers()
{
  Q_EMIT item_is_about_to_be_changed();
    delete_aabb_tree(this);
    d->init();
    Base::invalidateOpenGLBuffers();
    are_buffers_filled = false;

    d->invalidate_stats();
}

void
Scene_polyhedron_item::selection_changed(bool p_is_selected)
{
    if(p_is_selected != is_selected)
    {
        is_selected = p_is_selected;
    }

}

void
Scene_polyhedron_item::setColor(QColor c)
{
  // reset patch ids
  if (d->colors_.size()>2 || d->plugin_has_set_color_vector_m)
  {
    d->colors_.clear();
    d->is_multicolor = false;
  }
  Scene_item::setColor(c);
  if(d->is_multicolor)
    invalidateOpenGLBuffers();
}

void
Scene_polyhedron_item::select(double orig_x,
                              double orig_y,
                              double orig_z,
                              double dir_x,
                              double dir_y,
                              double dir_z)
{

    void* vertex_to_emit = 0;
    if(d->facet_picking_m)
    {
        typedef Input_facets_AABB_tree Tree;
        typedef Tree::Object_and_primitive_id Object_and_primitive_id;

        Tree* aabb_tree = static_cast<Tree*>(d->get_aabb_tree());
        if(aabb_tree)
        {
            const Kernel::Point_3 ray_origin(orig_x, orig_y, orig_z);
            const Kernel::Vector_3 ray_dir(dir_x, dir_y, dir_z);
            const Kernel::Ray_3 ray(ray_origin, ray_dir);
            typedef std::list<Object_and_primitive_id> Intersections;
            Intersections intersections;
            aabb_tree->all_intersections(ray, std::back_inserter(intersections));
            Intersections::iterator closest = intersections.begin();
            if(closest != intersections.end())
            {
                const Kernel::Point_3* closest_point =
                        CGAL::object_cast<Kernel::Point_3>(&closest->first);
                for(Intersections::iterator
                    it = boost::next(intersections.begin()),
                    end = intersections.end();
                    it != end; ++it)
                {
                    if(! closest_point) {
                        closest = it;
                    }
                    else {
                        const Kernel::Point_3* it_point =
                                CGAL::object_cast<Kernel::Point_3>(&it->first);
                        if(it_point &&
                                (ray_dir * (*it_point - *closest_point)) < 0)
                        {
                            closest = it;
                            closest_point = it_point;
                        }
                    }
                }
                if(closest_point) {
                    Polyhedron::Facet_handle selected_fh = closest->second;

                    // The computation of the nearest vertex may be costly.  Only
                    // do it if some objects are connected to the signal
                    // 'selected_vertex'.
                    if(QObject::receivers(SIGNAL(selected_vertex(void*))) > 0)
                    {
                        Polyhedron::Halfedge_around_facet_circulator
                                he_it = selected_fh->facet_begin(),
                                around_end = he_it;

                        Polyhedron::Vertex_handle v = he_it->vertex(), nearest_v = v;

                        Kernel::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                                    v->point());
                        while(++he_it != around_end) {
                            v = he_it->vertex();
                            Kernel::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                                            v->point());
                            if(new_sq_dist < sq_dist) {
                                sq_dist = new_sq_dist;
                                nearest_v = v;
                            }
                        }
                        //bottleneck
                        vertex_to_emit = (void*)(&*nearest_v);
                    }

                    if(QObject::receivers(SIGNAL(selected_edge(void*))) > 0
                            || QObject::receivers(SIGNAL(selected_halfedge(void*))) > 0)
                    {
                        Polyhedron::Halfedge_around_facet_circulator
                                he_it = selected_fh->facet_begin(),
                                around_end = he_it;

                        Polyhedron::Halfedge_handle nearest_h = he_it;
                        Kernel::FT sq_dist =
                                CGAL::squared_distance(*closest_point,
                                                       Kernel::Segment_3(he_it->vertex()->point(),
                                                                         he_it->opposite()->
                                                                         vertex()->
                                                                         point()));

                        while(++he_it != around_end)
                        {
                            Kernel::FT new_sq_dist =
                                    CGAL::squared_distance(*closest_point,
                                                           Kernel::Segment_3(he_it->vertex()->point(),
                                                                             he_it->opposite()->
                                                                             vertex()->
                                                                             point()));
                            if(new_sq_dist < sq_dist) {
                                sq_dist = new_sq_dist;
                                nearest_h = he_it;
                            }
                        }

                        Q_EMIT selected_halfedge((void*)(&*nearest_h));
                        Q_EMIT selected_edge((void*)(std::min)(&*nearest_h, &*nearest_h->opposite()));
                    }
                    Q_EMIT selected_vertex(vertex_to_emit);
                    Q_EMIT selected_facet((void*)(&*selected_fh));
                    if(d->erase_next_picked_facet_m) {
                        polyhedron()->erase_facet(selected_fh->halfedge());
                        polyhedron()->normalize_border();
                        //set_erase_next_picked_facet(false);
                        invalidateOpenGLBuffers();

                        Q_EMIT itemChanged();
                    }
                }
            }
        }
    }
    Base::select(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z);
    Q_EMIT selection_done();
}

void Scene_polyhedron_item::update_vertex_indices()
{
    std::size_t id=0;
    for (Polyhedron::Vertex_iterator vit = polyhedron()->vertices_begin(),
         vit_end = polyhedron()->vertices_end(); vit != vit_end; ++vit)
    {
        vit->id()=id++;
    }
}
void Scene_polyhedron_item::update_facet_indices()
{
    std::size_t id=0;
    for (Polyhedron::Facet_iterator  fit = polyhedron()->facets_begin(),
         fit_end = polyhedron()->facets_end(); fit != fit_end; ++fit)
    {
        fit->id()=id++;
    }
}
void Scene_polyhedron_item::update_halfedge_indices()
{
    std::size_t id=0;
    for (Polyhedron::Halfedge_iterator hit = polyhedron()->halfedges_begin(),
         hit_end = polyhedron()->halfedges_end(); hit != hit_end; ++hit)
    {
        hit->id()=id++;
    }
}
void Scene_polyhedron_item::invalidate_aabb_tree()
{
  delete_aabb_tree(this);
}
QString Scene_polyhedron_item::computeStats(int type)
{
  double minl, maxl, meanl, midl;
  switch (type)
  {
  case MIN_LENGTH:
  case MAX_LENGTH:
  case MID_LENGTH:
  case MEAN_LENGTH:
  case NB_NULL_LENGTH:
    d->poly->normalize_border();
    edges_length(d->poly, minl, maxl, meanl, midl, d->number_of_null_length_edges);
  }

  double mini, maxi, ave;
  switch (type)
  {
  case MIN_ANGLE:
  case MAX_ANGLE:
  case MEAN_ANGLE:
    angles(d->poly, mini, maxi, ave);
  }

  switch(type)
  {
  case NB_VERTICES:
    return QString::number(d->poly->size_of_vertices());

  case NB_FACETS:
    return QString::number(d->poly->size_of_facets());
  
  case NB_CONNECTED_COMPOS:
  {
    typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
    int i = 0;
    BOOST_FOREACH(face_descriptor f, faces(*(d->poly))){
      f->id() = i++;
    }
    boost::vector_property_map<int,
      boost::property_map<Polyhedron, boost::face_index_t>::type>
      fccmap(get(boost::face_index, *(d->poly)));
    return QString::number(PMP::connected_components(*(d->poly), fccmap));
  }
  case NB_BORDER_EDGES:
    d->poly->normalize_border();
    return QString::number(d->poly->size_of_border_halfedges());

  case NB_EDGES:
    return QString::number(d->poly->size_of_halfedges() / 2);

  case NB_DEGENERATED_FACES:
  {
    if (d->poly->is_pure_triangle())
    {
      if (d->number_of_degenerated_faces == (unsigned int)(-1))
        d->number_of_degenerated_faces = nb_degenerate_faces(d->poly, get(CGAL::vertex_point, *(d->poly)));
      return QString::number(d->number_of_degenerated_faces);
    }
    else
      return QString("n/a");
  }
  case AREA:
  {
    if (d->poly->is_pure_triangle())
    {
      if(d->area == -std::numeric_limits<double>::infinity())
        d->area = CGAL::Polygon_mesh_processing::area(*(d->poly));
      return QString::number(d->area);
    }
    else
      return QString("n/a");
  }
  case VOLUME:
  {
    if (d->poly->is_pure_triangle() && d->poly->is_closed())
    {
      if (d->volume == -std::numeric_limits<double>::infinity())
        d->volume = CGAL::Polygon_mesh_processing::volume(*(d->poly));
      return QString::number(d->volume);
    }
    else
      return QString("n/a");
  }
  case SELFINTER:
  {
    //todo : add a test about cache validity
    if (d->poly->is_pure_triangle())
      d->self_intersect = CGAL::Polygon_mesh_processing::does_self_intersect(*(d->poly));
    if (d->self_intersect)
      return QString("Yes");
    else if (d->poly->is_pure_triangle())
      return QString("No");
    else
      return QString("n/a");
  }
  case MIN_LENGTH:
    return QString::number(minl);
  case MAX_LENGTH:
    return QString::number(maxl);
  case MID_LENGTH:
    return QString::number(midl);
  case MEAN_LENGTH:
    return QString::number(meanl);
  case NB_NULL_LENGTH:
    return QString::number(d->number_of_null_length_edges);

  case MIN_ANGLE:
    return QString::number(mini);
  case MAX_ANGLE:
    return QString::number(maxi);
  case MEAN_ANGLE:
    return QString::number(ave);

  case HOLES:
    return QString::number(nb_holes(d->poly));
  }
  return QString();
}

CGAL::Three::Scene_item::Header_data Scene_polyhedron_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  //categories
  data.categories.append(std::pair<QString,int>(QString("Properties"),9));
  data.categories.append(std::pair<QString,int>(QString("Edges"),6));
  data.categories.append(std::pair<QString,int>(QString("Angles"),3));


  //titles
  data.titles.append(QString("#Vertices"));
  data.titles.append(QString("#Facets"));
  data.titles.append(QString("#Connected Components"));
  data.titles.append(QString("#Border Edges"));
  data.titles.append(QString("#Degenerated Faces"));
  data.titles.append(QString("Connected Components of the Boundary"));
  data.titles.append(QString("Area"));
  data.titles.append(QString("Volume"));
  data.titles.append(QString("Self-Intersecting"));
  data.titles.append(QString("#Edges"));
  data.titles.append(QString("Minimum Length"));
  data.titles.append(QString("Maximum Length"));
  data.titles.append(QString("Median Length"));
  data.titles.append(QString("Mean Length"));
  data.titles.append(QString("#Null Length"));
  data.titles.append(QString("Minimum"));
  data.titles.append(QString("Maximum"));
  data.titles.append(QString("Average"));
  return data;
}


void Scene_polyhedron_item::printPrimitiveId(QPoint point, CGAL::Three::Viewer_interface *viewer)
{
  TextRenderer *renderer = viewer->textRenderer;
  renderer->getLocalTextItems().removeAll(d->targeted_id);
  renderer->removeTextList(textItems);
  textItems->clear();
  QFont font;
  font.setBold(true);

  typedef Input_facets_AABB_tree Tree;
  typedef Tree::Intersection_and_primitive_id<Kernel::Ray_3>::Type Intersection_and_primitive_id;

  Tree* aabb_tree = static_cast<Input_facets_AABB_tree*>(d->get_aabb_tree());
  if(aabb_tree) {
    //find clicked facet
    bool found = false;
    const Kernel::Point_3 ray_origin(viewer->camera()->position().x, viewer->camera()->position().y, viewer->camera()->position().z);
    qglviewer::Vec point_under = viewer->camera()->pointUnderPixel(point,found);
    qglviewer::Vec dir = point_under - viewer->camera()->position();
    const Kernel::Vector_3 ray_dir(dir.x, dir.y, dir.z);
    const Kernel::Ray_3 ray(ray_origin, ray_dir);
    typedef std::list<Intersection_and_primitive_id> Intersections;
    Intersections intersections;
    aabb_tree->all_intersections(ray, std::back_inserter(intersections));

    if(!intersections.empty()) {
      Intersections::iterator closest = intersections.begin();
      const Kernel::Point_3* closest_point =
          boost::get<Kernel::Point_3>(&closest->first);
      for(Intersections::iterator
          it = boost::next(intersections.begin()),
          end = intersections.end();
          it != end; ++it)
      {
        if(! closest_point) {
          closest = it;
        }
        else {
          const Kernel::Point_3* it_point =
              boost::get<Kernel::Point_3>(&it->first);
          if(it_point &&
             (ray_dir * (*it_point - *closest_point)) < 0)
          {
            closest = it;
            closest_point = it_point;
          }
        }
      }
      if(closest_point) {
        Polyhedron::Facet_handle selected_fh = closest->second;
        //Test spots around facet to find the closest to point

        double min_dist = (std::numeric_limits<double>::max)();
        TextItem text_item;
        Kernel::Point_3 pt_under(point_under.x, point_under.y, point_under.z);

        // test the vertices of the closest face
        Q_FOREACH(Polyhedron::Vertex_handle vh, vertices_around_face(selected_fh->halfedge(), *d->poly))
        {
          Kernel::Point_3 test=vh->point();
          double dist = CGAL::squared_distance(test, pt_under);
          if( dist < min_dist){
            min_dist = dist;
            text_item = TextItem(test.x(), test.y(), test.z(), QString("%1").arg(vh->id()), true, font, Qt::red);
          }
        }
        // test the midpoint of edges of the closest face
        Q_FOREACH(boost::graph_traits<Polyhedron>::halfedge_descriptor e, halfedges_around_face(selected_fh->halfedge(), *d->poly))
        {
          Kernel::Point_3 test=CGAL::midpoint(source(e, *d->poly)->point(),target(e, *d->poly)->point());
          double dist = CGAL::squared_distance(test, pt_under);
          if(dist < min_dist){
            min_dist = dist;
            text_item = TextItem(test.x(), test.y(), test.z(), QString("%1").arg(e->id()/2), true, font, Qt::green);
          }
        }

        // test the centroid of the closest face
        double x(0), y(0), z(0);
        int total(0);
        Q_FOREACH(Polyhedron::Vertex_handle vh, vertices_around_face(selected_fh->halfedge(), *d->poly))
        {
          x+=vh->point().x();
          y+=vh->point().y();
          z+=vh->point().z();
          ++total;
        }

        Kernel::Point_3 test(x/total, y/total, z/total);
        double dist = CGAL::squared_distance(test, pt_under);
        if(dist < min_dist){
          min_dist = dist;
          text_item = TextItem(test.x(), test.y(), test.z(), QString("%1").arg(selected_fh->id()), true, font, Qt::blue);
        }

        TextItem* former_targeted_id=d->targeted_id;
        if (d->targeted_id == NULL || d->targeted_id->position() != text_item.position() )
        {
          d->targeted_id = new TextItem(text_item);
          textItems->append(d->targeted_id);
          renderer->addTextList(textItems);
        }
        else
          d->targeted_id=NULL;
        if(former_targeted_id != NULL) renderer->removeText(former_targeted_id);
      }
    }
  }
}

void Scene_polyhedron_item::printPrimitiveIds(CGAL::Three::Viewer_interface *viewer) const 
{
  TextRenderer *renderer = viewer->textRenderer;


  if(!d->all_ids_displayed)
  {
    QFont font;
    font.setBold(true);

    //fills textItems
    Q_FOREACH(Polyhedron::Vertex_const_handle vh, vertices(*d->poly))
    {
      const Point& p = vh->point();
      textItems->append(new TextItem((float)p.x(), (float)p.y(), (float)p.z(), QString("%1").arg(vh->id()), true, font, Qt::red));

    }

    Q_FOREACH(boost::graph_traits<Polyhedron>::edge_descriptor e, edges(*d->poly))
    {
      const Point& p1 = source(e, *d->poly)->point();
      const Point& p2 = target(e, *d->poly)->point();
      textItems->append(new TextItem((float)(p1.x() + p2.x()) / 2, (float)(p1.y() + p2.y()) / 2, (float)(p1.z() + p2.z()) / 2, QString("%1").arg(e.halfedge()->id() / 2), true, font, Qt::green));
    }

    Q_FOREACH(Polyhedron::Facet_handle fh, faces(*d->poly))
    {
      double x(0), y(0), z(0);
      int total(0);
      Q_FOREACH(Polyhedron::Vertex_handle vh, vertices_around_face(fh->halfedge(), *d->poly))
      {
        x += vh->point().x();
        y += vh->point().y();
        z += vh->point().z();
        ++total;
      }

      textItems->append(new TextItem((float)x / total, (float)y / total, (float)z / total, QString("%1").arg(fh->id()), true, font, Qt::blue));
    }
    //add the QList to the render's pool
    renderer->addTextList(textItems);
    if(textItems->size() > static_cast<std::size_t>(renderer->getMax_textItems()))
      d->all_ids_displayed = !d->all_ids_displayed;
  }
  if(d->all_ids_displayed)
  {
    //clears TextItems
    textItems->clear();
    renderer->removeTextList(textItems);
    if(d->targeted_id)
    {
      textItems->append(d->targeted_id);
      renderer->addTextList(textItems);
    }
  }
  d->all_ids_displayed = !d->all_ids_displayed;
}

bool Scene_polyhedron_item::testDisplayId(double x, double y, double z, CGAL::Three::Viewer_interface* viewer)
{
  Kernel::Point_3 src(x,y,z);
  Kernel::Point_3 dest(viewer->camera()->position().x, viewer->camera()->position().y,viewer->camera()->position().z);
  Kernel::Vector_3 v(src,dest);
  v = 0.01*v;
  Kernel::Point_3 point = src;
  point = point + v;
  Kernel::Segment_3 query(point, dest);
  return !static_cast<Input_facets_AABB_tree*>(d->get_aabb_tree())->do_intersect(query);

}

std::vector<QColor>& Scene_polyhedron_item::color_vector() {return d->colors_;}
void Scene_polyhedron_item::set_color_vector_read_only(bool on_off) {d->plugin_has_set_color_vector_m=on_off;}
bool Scene_polyhedron_item::is_color_vector_read_only() { return d->plugin_has_set_color_vector_m;}
int Scene_polyhedron_item::getNumberOfNullLengthEdges(){return d->number_of_null_length_edges;}
int Scene_polyhedron_item::getNumberOfDegeneratedFaces(){return d->number_of_degenerated_faces;}
bool Scene_polyhedron_item::triangulated(){return d->poly->is_pure_triangle();}
bool Scene_polyhedron_item::self_intersected(){return !(d->self_intersect);}
void Scene_polyhedron_item::setItemIsMulticolor(bool b){ d->is_multicolor = b;}
bool Scene_polyhedron_item::isItemMulticolor(){ return d->is_multicolor;}

