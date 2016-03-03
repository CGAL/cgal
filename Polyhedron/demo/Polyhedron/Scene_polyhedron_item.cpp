#include "Scene_polyhedron_item.h"
#include <CGAL/AABB_intersections.h>
#include "Kernel_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
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

#include <boost/foreach.hpp>
#include <boost/container/flat_map.hpp>

namespace PMP = CGAL::Polygon_mesh_processing;


typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Input_facets_AABB_tree;
const char* aabb_property_name = "Scene_polyhedron_item aabb tree";

Input_facets_AABB_tree* get_aabb_tree(Scene_polyhedron_item* item)
{
    QVariant aabb_tree_property = item->property(aabb_property_name);
    if(aabb_tree_property.isValid()) {
        void* ptr = aabb_tree_property.value<void*>();
        return static_cast<Input_facets_AABB_tree*>(ptr);
    }
    else {
        Polyhedron* poly = item->polyhedron();
        if(poly) {
            Input_facets_AABB_tree* tree =
                    new Input_facets_AABB_tree(faces(*poly).first,
                                               faces(*poly).second,
                                               *poly);
            item->setProperty(aabb_property_name,
                              QVariant::fromValue<void*>(tree));
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

typedef Polyhedron::Traits Traits;
typedef Polyhedron::Facet Facet;
typedef CGAL::Triangulation_2_filtered_projection_traits_3<Traits>   P_traits;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
struct Face_info {
    Polyhedron::Halfedge_handle e[3];
    bool is_external;
};
typedef CGAL::Triangulation_vertex_base_with_info_2<Halfedge_handle,
P_traits>        Vb;
typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
P_traits>          Fb1;
typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>   Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                  TDS;
typedef CGAL::No_intersection_tag                                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits,
TDS,
Itag>             CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase>              CDT;

//Make sure all the facets are triangles

template<typename FaceNormalPmap, typename VertexNormalPmap>
void
Scene_polyhedron_item::triangulate_facet(Facet_iterator fit,
                                         const FaceNormalPmap& fnmap,
                                         const VertexNormalPmap& vnmap,
                                         const bool colors_only) const
{
    //Computes the normal of the facet
    Traits::Vector_3 normal = get(fnmap, fit);

    //check if normal contains NaN values
    if (normal.x() != normal.x() || normal.y() != normal.y() || normal.z() != normal.z())
    {
        qDebug()<<"Warning : normal is not valid. Facet not displayed";
        return;
    }
    P_traits cdt_traits(normal);
    CDT cdt(cdt_traits);

    Facet::Halfedge_around_facet_circulator
            he_circ = fit->facet_begin(),
            he_circ_end(he_circ);

    // Iterates on the vector of facet handles
    typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
    boost::container::flat_map<CDT::Vertex_handle, vertex_descriptor> v2v;
    CDT::Vertex_handle previous, first;
    do {
        CDT::Vertex_handle vh = cdt.insert(he_circ->vertex()->point());
        v2v.insert(std::make_pair(vh, he_circ->vertex()));
        if(first == 0) {
            first = vh;
        }
        vh->info() = he_circ;
        if(previous != 0 && previous != vh) {
            cdt.insert_constraint(previous, vh);
        }
        previous = vh;
    } while( ++he_circ != he_circ_end );
    cdt.insert_constraint(previous, first);
    // sets mark is_external
    for(CDT::All_faces_iterator
        fit2 = cdt.all_faces_begin(),
        end = cdt.all_faces_end();
        fit2 != end; ++fit2)
    {
        fit2->info().is_external = false;
    }
    //check if the facet is external or internal
    std::queue<CDT::Face_handle> face_queue;
    face_queue.push(cdt.infinite_vertex()->face());
    while(! face_queue.empty() ) {
        CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(fh->info().is_external) continue;
        fh->info().is_external = true;
        for(int i = 0; i <3; ++i) {
            if(!cdt.is_constrained(std::make_pair(fh, i)))
            {
                face_queue.push(fh->neighbor(i));
            }
        }
    }
    //iterates on the internal faces to add the vertices to the positions
    //and the normals to the appropriate vectors
    const int this_patch_id = fit->patch_id();
    for(CDT::Finite_faces_iterator
        ffit = cdt.finite_faces_begin(),
        end = cdt.finite_faces_end();
        ffit != end; ++ffit)
    {
        if(ffit->info().is_external)
            continue;

        if (!is_monochrome)
        {
          for (int i = 0; i<3; ++i)
          {
            color_facets.push_back(colors_[this_patch_id].redF());
            color_facets.push_back(colors_[this_patch_id].greenF());
            color_facets.push_back(colors_[this_patch_id].blueF());

            color_facets.push_back(colors_[this_patch_id].redF());
            color_facets.push_back(colors_[this_patch_id].greenF());
            color_facets.push_back(colors_[this_patch_id].blueF());
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

        Traits::Vector_3 ng = get(vnmap, v2v[ffit->vertex(0)]);
        push_back_xyz(ng, normals_gouraud);

        ng = get(vnmap, v2v[ffit->vertex(1)]);
        push_back_xyz(ng, normals_gouraud);

        ng = get(vnmap, v2v[ffit->vertex(2)]);
        push_back_xyz(ng, normals_gouraud);
    }
}


#include <QObject>
#include <QMenu>
#include <QAction>


void
Scene_polyhedron_item::initialize_buffers(CGAL::Three::Viewer_interface* viewer) const
{
    //vao containing the data for the facets
    {
        program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
        program->bind();
        //flat
        vaos[Facets]->bind();
        buffers[Facets_vertices].bind();
        buffers[Facets_vertices].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[Facets_vertices].release();



        buffers[Facets_normals_flat].bind();
        buffers[Facets_normals_flat].allocate(normals_flat.data(),
                            static_cast<int>(normals_flat.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[Facets_normals_flat].release();

        if(!is_monochrome)
        {
            buffers[Facets_color].bind();
            buffers[Facets_color].allocate(color_facets.data(),
                                static_cast<int>(color_facets.size()*sizeof(float)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            buffers[Facets_color].release();
        }
        vaos[Facets]->release();
        //gouraud
        vaos[Gouraud_Facets]->bind();
        buffers[Facets_vertices].bind();
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[Facets_vertices].release();

        buffers[Facets_normals_gouraud].bind();
        buffers[Facets_normals_gouraud].allocate(normals_gouraud.data(),
                            static_cast<int>(normals_gouraud.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[Facets_normals_gouraud].release();
        if(!is_monochrome)
        {
            buffers[Facets_color].bind();
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            buffers[Facets_color].release();
        }
        else
        {
            program->disableAttributeArray("colors");
        }
        vaos[Gouraud_Facets]->release();

        program->release();

    }
    //vao containing the data for the lines
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[Edges]->bind();

        buffers[Edges_vertices].bind();
        buffers[Edges_vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[Edges_vertices].release();

        buffers[Edges_color].bind();
        buffers[Edges_color].allocate(color_lines.data(),
                            static_cast<int>(color_lines.size()*sizeof(float)));
       if(!is_monochrome)
       {
           program->enableAttributeArray("colors");
           program->setAttributeBuffer("colors",GL_FLOAT,0,3);
           buffers[Edges_color].release();
       }
       else
       {
           program->disableAttributeArray("colors");
       }
        program->release();

        vaos[Edges]->release();

    }
  //vao containing the data for the feature_edges
  {
      program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
      program->bind();
      vaos[Feature_edges]->bind();

      buffers[Feature_edges_vertices].bind();
      buffers[Feature_edges_vertices].allocate(positions_feature_lines.data(),
                          static_cast<int>(positions_feature_lines.size()*sizeof(float)));
      program->enableAttributeArray("vertex");
      program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
      buffers[Feature_edges_vertices].release();
      program->disableAttributeArray("colors");
      program->release();

      vaos[Feature_edges]->release();

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
    are_buffers_filled = true;
}

void
Scene_polyhedron_item::compute_normals_and_vertices(const bool colors_only) const
{
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

    boost::container::flat_map<face_descriptor, Vector> face_normals_map;
    boost::associative_property_map< boost::container::flat_map<face_descriptor, Vector> >
      nf_pmap(face_normals_map);
    boost::container::flat_map<vertex_descriptor, Vector> vertex_normals_map;
    boost::associative_property_map< boost::container::flat_map<vertex_descriptor, Vector> >
      nv_pmap(vertex_normals_map);

    PMP::compute_normals(*poly, nv_pmap, nf_pmap);

    Facet_iterator f = poly->facets_begin();
    for(f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {
      if (f == boost::graph_traits<Polyhedron>::null_face())
        continue;

      if(is_triangle(f->halfedge(),*poly))
      {
          const int this_patch_id = f->patch_id();
          Vector n = get(nf_pmap, f);
          HF_circulator he = f->facet_begin();
          HF_circulator end = he;
          CGAL_For_all(he,end)
          {
            if (!is_monochrome)
            {
              color_facets.push_back(colors_[this_patch_id].redF());
              color_facets.push_back(colors_[this_patch_id].greenF());
              color_facets.push_back(colors_[this_patch_id].blueF());
            }
            if (colors_only)
              continue;

            // If Flat shading:1 normal per polygon added once per vertex
            push_back_xyz(n, normals_flat);

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
        if (!is_monochrome)
        {
          const int this_patch_id = f->patch_id();
          for (unsigned int i = 0; i < 6; ++i)
          { //6 "halfedges" for the quad, because it is 2 triangles
            color_facets.push_back(colors_[this_patch_id].redF());
            color_facets.push_back(colors_[this_patch_id].greenF());
            color_facets.push_back(colors_[this_patch_id].blueF());
          }
        }
        if (colors_only)
          continue;

        Vector nf = get(nf_pmap, f);

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
        triangulate_facet(f, nf_pmap, nv_pmap, colors_only);
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
          if (!is_monochrome)
          {
            color_lines.push_back(this->color().lighter(50).redF());
            color_lines.push_back(this->color().lighter(50).greenF());
            color_lines.push_back(this->color().lighter(50).blueF());

            color_lines.push_back(this->color().lighter(50).redF());
            color_lines.push_back(this->color().lighter(50).greenF());
            color_lines.push_back(this->color().lighter(50).blueF());
          }
          if (colors_only)
            continue;

          push_back_xyz(a, positions_lines);
          positions_lines.push_back(1.0);

          push_back_xyz(b, positions_lines);
          positions_lines.push_back(1.0);
        }
    }
}

Scene_polyhedron_item::Scene_polyhedron_item()
    : Scene_item(NbOfVbos,NbOfVaos),
      poly(new Polyhedron),
      show_only_feature_edges_m(false),
      show_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    cur_shading=FlatPlusEdges;
    is_selected = true;
    nb_facets = 0;
    nb_lines = 0;
    nb_f_lines = 0;
    init();
}

Scene_polyhedron_item::Scene_polyhedron_item(Polyhedron* const p)
    : Scene_item(NbOfVbos,NbOfVaos),
      poly(p),
      show_only_feature_edges_m(false),
      show_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    cur_shading=FlatPlusEdges;
    is_selected = true;
    nb_facets = 0;
    nb_lines = 0;
    nb_f_lines = 0;
    init();
    invalidateOpenGLBuffers();
}

Scene_polyhedron_item::Scene_polyhedron_item(const Polyhedron& p)
    : Scene_item(NbOfVbos,NbOfVaos),
      poly(new Polyhedron(p)),
      show_only_feature_edges_m(false),
      show_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    //setItemIsMulticolor(true);
    cur_shading=FlatPlusEdges;
    is_selected=true;
    init();
    nb_facets = 0;
    nb_lines = 0;
    nb_f_lines = 0;
    invalidateOpenGLBuffers();
}

Scene_polyhedron_item::~Scene_polyhedron_item()
{
    delete_aabb_tree(this);
    delete poly;
}

#include "Color_map.h"

void
Scene_polyhedron_item::
init()
{
  typedef Polyhedron::Facet_iterator Facet_iterator;

  if ( !plugin_has_set_color_vector_m )
  {
    // Fill indices map and get max subdomain value
    int max = 0;
    for(Facet_iterator fit = poly->facets_begin(), end = poly->facets_end() ;
        fit != end; ++fit)
    {
      max = (std::max)(max, fit->patch_id());
    }

    colors_.resize(0);
    compute_color_map(this->color(), max + 1,
                      std::back_inserter(colors_));
  }
  invalidate_stats();
}

void
Scene_polyhedron_item::
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
    return new Scene_polyhedron_item(*poly);}

// Load polyhedron from .OFF file
bool
Scene_polyhedron_item::load(std::istream& in)
{


    in >> *poly;

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
  
  Point p;
  std::string line;
  std::size_t i, j, k;
  bool failed = false;
  while(getline(in, line)){
    if(line[0] == 'v'){
      std::istringstream iss(line.substr(1));
      iss >> p;
      if(! iss) failed = true;
      points.push_back(p);
    } else if(line[0] == 'f'){
      std::istringstream iss(line.substr(1));
      iss >> i >> j >> k;
      if(! iss) failed = true;
      std::vector<std::size_t> face;
      face.push_back(i-1);
      face.push_back(j-1);
      face.push_back(k-1);
      faces.push_back(face);
    }
  }
  if(CGAL::Polygon_mesh_processing::orient_polygon_soup(points,faces)){
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points,faces,*poly);
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
    out << *poly;
    return (bool) out;
}

QString
Scene_polyhedron_item::toolTip() const
{
    if(!poly)
        return QString();

  QString str =
         QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                     "Number of facets: %4")
            .arg(this->name())
            .arg(poly->size_of_vertices())
            .arg(poly->size_of_halfedges()/2)
            .arg(poly->size_of_facets())
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
    actionShowOnlyFeatureEdges->setChecked(show_only_feature_edges_m);
    actionShowOnlyFeatureEdges->setObjectName("actionShowOnlyFeatureEdges");
    connect(actionShowOnlyFeatureEdges, SIGNAL(toggled(bool)),
            this, SLOT(show_only_feature_edges(bool)));

    QAction* actionShowFeatureEdges =
        menu->addAction(tr("Show Feature Edges"));
    actionShowFeatureEdges->setCheckable(true);
    actionShowFeatureEdges->setChecked(show_feature_edges_m);
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
  if(action) action->setChecked(show_only_feature_edges_m);
  action = menu->findChild<QAction*>("actionShowFeatureEdges");
  if(action) action->setChecked(show_feature_edges_m);
  action = menu->findChild<QAction*>("actionPickFacets");
  if(action) action->setChecked(facet_picking_m);
  action = menu->findChild<QAction*>("actionEraseNextFacet");
  if(action) action->setChecked(erase_next_picked_facet_m);
  return menu;
}

void Scene_polyhedron_item::show_only_feature_edges(bool b)
{
    show_only_feature_edges_m = b;
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
}

void Scene_polyhedron_item::show_feature_edges(bool b)
{
  show_feature_edges_m = b;
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

void Scene_polyhedron_item::enable_facets_picking(bool b)
{
    facet_picking_m = b;
}

void Scene_polyhedron_item::set_erase_next_picked_facet(bool b)
{
    if(b) { facet_picking_m = true; } // automatically activate facet_picking
    erase_next_picked_facet_m = b;
}

void Scene_polyhedron_item::draw(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        compute_normals_and_vertices();
        initialize_buffers(viewer);
        compute_bbox();
    }

    if(renderingMode() == Flat || renderingMode() == FlatPlusEdges)
        vaos[Facets]->bind();
    else
    {
        vaos[Gouraud_Facets]->bind();
    }
    attrib_buffers(viewer, PROGRAM_WITH_LIGHT);
    program = getShaderProgram(PROGRAM_WITH_LIGHT);
    program->bind();
    if(is_monochrome)
    {
            program->setAttributeValue("colors", this->color());
    }
    if(is_selected)
            program->setUniformValue("is_selected", true);
    else
            program->setUniformValue("is_selected", false);
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_facets/4));
    program->release();
    if(renderingMode() == Flat || renderingMode() == FlatPlusEdges)
        vaos[Facets]->release();
    else
        vaos[Gouraud_Facets]->release();
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::draw_edges(CGAL::Three::Viewer_interface* viewer) const
{
    if (!are_buffers_filled)
    {
        compute_normals_and_vertices();
        initialize_buffers(viewer);
        compute_bbox();
    }

    if(!show_only_feature_edges_m)
    {
        vaos[Edges]->bind();

        attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        program->bind();
        //draw the edges
        if(is_monochrome)
        {
            program->setAttributeValue("colors", this->color().lighter(50));
            if(is_selected)
                program->setUniformValue("is_selected", true);
            else
                program->setUniformValue("is_selected", false);
        }
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_lines/4));
        program->release();
        vaos[Edges]->release();
    }

    //draw the feature edges
    vaos[Feature_edges]->bind();
    attrib_buffers(viewer, PROGRAM_NO_SELECTION);
    program = getShaderProgram(PROGRAM_NO_SELECTION);
    program->bind();
    if(show_feature_edges_m || show_only_feature_edges_m)
        program->setAttributeValue("colors", Qt::red);
    else
    {
        if(!is_selected)
            program->setAttributeValue("colors", this->color().lighter(50));
        else
            program->setAttributeValue("colors",QColor(0,0,0));
    }
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_f_lines/4));
    program->release();
    vaos[Feature_edges]->release();
    }

void
Scene_polyhedron_item::draw_points(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        compute_normals_and_vertices();
        initialize_buffers(viewer);
        compute_bbox();
    }

    vaos[Edges]->bind();
    attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    program->bind();
    //draw the points
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(nb_lines/4));
    // Clean-up
    program->release();
    vaos[Edges]->release();
}

Polyhedron*
Scene_polyhedron_item::polyhedron()       { return poly; }
const Polyhedron*
Scene_polyhedron_item::polyhedron() const { return poly; }

bool
Scene_polyhedron_item::isEmpty() const {
    return (poly == 0) || poly->empty();
}

void Scene_polyhedron_item::compute_bbox() const {
    const Kernel::Point_3& p = *(poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedron::Point_iterator it = poly->points_begin();
        it != poly->points_end();
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
    init();
    Base::invalidateOpenGLBuffers();
    are_buffers_filled = false;

    invalidate_stats();
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
  if (colors_.size()>2 || plugin_has_set_color_vector_m)
  {
    BOOST_FOREACH(Polyhedron::Facet_handle fh, faces(*poly))
      fh->set_patch_id(1);
    colors_[1]=c;
  }
  Scene_item::setColor(c);
}

void
Scene_polyhedron_item::select(double orig_x,
                              double orig_y,
                              double orig_z,
                              double dir_x,
                              double dir_y,
                              double dir_z)
{
    if(facet_picking_m) {
        typedef Input_facets_AABB_tree Tree;
        typedef Tree::Object_and_primitive_id Object_and_primitive_id;

        Tree* aabb_tree = get_aabb_tree(this);
        if(aabb_tree) {
            const Kernel::Point_3 ray_origin(orig_x, orig_y, orig_z);
            const Kernel::Vector_3 ray_dir(dir_x, dir_y, dir_z);
            const Kernel::Ray_3 ray(ray_origin, ray_dir);
            typedef std::list<Object_and_primitive_id> Intersections;
            Intersections intersections;
            aabb_tree->all_intersections(ray, std::back_inserter(intersections));
            Intersections::iterator closest = intersections.begin();
            if(closest != intersections.end()) {
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
            Q_EMIT selected_vertex((void*)(&*nearest_v));
                    }

                    if(QObject::receivers(SIGNAL(selected_edge(void*))) > 0
                            || QObject::receivers(SIGNAL(selected_halfedge(void*))) > 0)
                    {
                        Polyhedron::Halfedge_around_facet_circulator
                                he_it = selected_fh->facet_begin(),
                                around_end = he_it;

                        Polyhedron::Halfedge_handle nearest_h = he_it;
                        Kernel::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                                    Kernel::Segment_3(he_it->vertex()->point(), he_it->opposite()->vertex()->point()));

                        while(++he_it != around_end) {
                            Kernel::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                                            Kernel::Segment_3(he_it->vertex()->point(), he_it->opposite()->vertex()->point()));
                            if(new_sq_dist < sq_dist) {
                                sq_dist = new_sq_dist;
                                nearest_h = he_it;
                            }
                        }

            Q_EMIT selected_halfedge((void*)(&*nearest_h));
            Q_EMIT selected_edge((void*)(std::min)(&*nearest_h, &*nearest_h->opposite()));
                    }

          Q_EMIT selected_facet((void*)(&*selected_fh));
                    if(erase_next_picked_facet_m) {
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
QString Scene_polyhedron_item::compute_stats(int type)
{
  double minl, maxl, meanl, midl;
  switch (type)
  {
  case MIN_LENGTH:
  case MAX_LENGTH:
  case MID_LENGTH:
  case MEAN_LENGTH:
  case NB_NULL_LENGTH:
    poly->normalize_border();
    edges_length(poly, minl, maxl, meanl, midl, number_of_null_length_edges);
  }

  double mini, maxi, ave;
  switch (type)
  {
  case MIN_ANGLE:
  case MAX_ANGLE:
  case MEAN_ANGLE:
    angles(poly, mini, maxi, ave);
  }

  switch(type)
  {
  case NB_VERTICES:
    return QString::number(poly->size_of_vertices());

  case NB_FACETS:
    return QString::number(poly->size_of_facets());
  
  case NB_CONNECTED_COMPOS:
  {
    typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
    int i = 0;
    BOOST_FOREACH(face_descriptor f, faces(*poly)){
      f->id() = i++;
    }
    boost::vector_property_map<int,
      boost::property_map<Polyhedron, boost::face_index_t>::type>
      fccmap(get(boost::face_index, *poly));
    return QString::number(PMP::connected_components(*poly, fccmap));
  }
  case NB_BORDER_EDGES:
    poly->normalize_border();
    return QString::number(poly->size_of_border_halfedges());

  case NB_EDGES:
    return QString::number(poly->size_of_halfedges() / 2);

  case NB_DEGENERATED_FACES:
  {
    if (poly->is_pure_triangle())
    {
      if (number_of_degenerated_faces == (unsigned int)(-1))
        number_of_degenerated_faces = nb_degenerate_faces(poly, get(CGAL::vertex_point, *poly));
      return QString::number(number_of_degenerated_faces);
    }
    else
      return QString("n/a");
  }
  case AREA:
  {
    if (poly->is_pure_triangle())
    {
      if(area == -std::numeric_limits<double>::infinity())
        area = CGAL::Polygon_mesh_processing::area(*poly);
      return QString::number(area);
    }
    else
      return QString("n/a");
  }
  case VOLUME:
  {
    if (poly->is_pure_triangle() && poly->is_closed())
    {
      if (volume == -std::numeric_limits<double>::infinity())
        volume = CGAL::Polygon_mesh_processing::volume(*poly);
      return QString::number(volume);
    }
    else
      return QString("n/a");
  }
  case SELFINTER:
  {
    //todo : add a test about cache validity
    if (poly->is_pure_triangle())
      self_intersect = CGAL::Polygon_mesh_processing::does_self_intersect(*poly);
    if (self_intersect)
      return QString("Yes");
    else if (poly->is_pure_triangle())
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
    return QString::number(number_of_null_length_edges);

  case MIN_ANGLE:
    return QString::number(mini);
  case MAX_ANGLE:
    return QString::number(maxi);
  case MEAN_ANGLE:
    return QString::number(ave);

  case HOLES:
    return QString::number(nb_holes(poly));
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
