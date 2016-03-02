#include "Scene_polyhedron_selection_item.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>





#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <boost/container/flat_map.hpp>


void Scene_polyhedron_selection_item::initialize_buffers(CGAL::Three::Viewer_interface *viewer)const
{
    //vao containing the data for the unselected facets
    {
        program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
        program->bind();

        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[0].release();



        buffers[1].bind();
        buffers[1].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(float)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[1].release();

        vaos[0]->release();
        program->release();

    }
    //vao containing the data for the unselected lines
    {
        program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
        program->bind();
        vaos[1]->bind();

        buffers[2].bind();
        buffers[2].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[2].release();

        program->release();

        vaos[1]->release();

    }
    //vao containing the data for the points
    {
        program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
        program->bind();
        vaos[2]->bind();

        buffers[3].bind();
        buffers[3].allocate(positions_points.data(),
                            static_cast<int>(positions_points.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[3].release();

        buffers[6].release();
        program->release();

        vaos[2]->release();
    }
    nb_facets = positions_facets.size();
    positions_facets.resize(0);
    std::vector<float>(positions_facets).swap(positions_facets);

    normals.resize(0);
    std::vector<float>(normals).swap(normals);

    nb_lines = positions_lines.size();
    positions_lines.resize(0);
    std::vector<float>(positions_lines).swap(positions_lines);

    nb_points = positions_points.size();
    positions_points.resize(0);
    std::vector<float>(positions_points).swap(positions_points);

    are_buffers_filled = true;
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
typedef Polyhedron::Traits	    Kernel;
typedef Kernel::Point_3	            Point;
typedef Kernel::Vector_3	    Vector;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;
typedef boost::graph_traits<Polyhedron>::face_descriptor   face_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;


template<typename FaceNormalPmap>
void
Scene_polyhedron_selection_item::triangulate_facet(Facet_handle fit,
                                         const FaceNormalPmap& fnmap) const
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
    for(CDT::Finite_faces_iterator
        ffit = cdt.finite_faces_begin(),
        end = cdt.finite_faces_end();
        ffit != end; ++ffit)
    {
        if(ffit->info().is_external)
            continue;

        push_back_xyz(ffit->vertex(0)->point(), positions_facets);
        push_back_xyz(ffit->vertex(1)->point(), positions_facets);
        push_back_xyz(ffit->vertex(2)->point(), positions_facets);

        push_back_xyz(normal, normals);
        push_back_xyz(normal, normals);
        push_back_xyz(normal, normals);
    }
}


void Scene_polyhedron_selection_item::compute_elements()const
{
    positions_facets.clear();
    positions_lines.clear();
    positions_points.clear();
    normals.clear();
    //The facet
    boost::container::flat_map<face_descriptor, Vector> face_normals_map;
    boost::associative_property_map< boost::container::flat_map<face_descriptor, Vector> >
      nf_pmap(face_normals_map);
    boost::container::flat_map<vertex_descriptor, Vector> vertex_normals_map;
    boost::associative_property_map< boost::container::flat_map<vertex_descriptor, Vector> >
      nv_pmap(vertex_normals_map);
if(!poly)
  return;
    PMP::compute_normals(*poly, nv_pmap, nf_pmap);
    for(Selection_set_facet::iterator
        it = selected_facets.begin(),
        end = selected_facets.end();
        it != end; it++)
    {
      Facet_handle f = (*it);
      if (f == boost::graph_traits<Polyhedron>::null_face())
        continue;

      if(is_triangle(f->halfedge(),*poly))
      {
        const Kernel::Vector_3 n =
            CGAL::Polygon_mesh_processing::compute_face_normal(f, *this->poly_item->polyhedron());

        normals.push_back(n.x());
        normals.push_back(n.y());
        normals.push_back(n.z());

        normals.push_back(n.x());
        normals.push_back(n.y());
        normals.push_back(n.z());

        normals.push_back(n.x());
        normals.push_back(n.y());
        normals.push_back(n.z());


        Polyhedron::Halfedge_around_facet_circulator
            he = f->facet_begin(),
            cend = he;

        CGAL_For_all(he,cend)
        {
          const Kernel::Point_3& p = he->vertex()->point();
          positions_facets.push_back(p.x());
          positions_facets.push_back(p.y());
          positions_facets.push_back(p.z());
        }
      }
      else if (is_quad(f->halfedge(), *poly))
      {
        Vector nf = get(nf_pmap, f);

        //1st half-quad
        Point p0 = f->halfedge()->vertex()->point();
        Point p1 = f->halfedge()->next()->vertex()->point();
        Point p2 = f->halfedge()->next()->next()->vertex()->point();

        push_back_xyz(p0, positions_facets);
        push_back_xyz(p1, positions_facets);
        push_back_xyz(p2, positions_facets);

        push_back_xyz(nf, normals);
        push_back_xyz(nf, normals);
        push_back_xyz(nf, normals);

        //2nd half-quad
        p0 = f->halfedge()->next()->next()->vertex()->point();
        p1 = f->halfedge()->prev()->vertex()->point();
        p2 = f->halfedge()->vertex()->point();

        push_back_xyz(p0, positions_facets);
        push_back_xyz(p1, positions_facets);
        push_back_xyz(p2, positions_facets);

        push_back_xyz(nf, normals);
        push_back_xyz(nf, normals);
        push_back_xyz(nf, normals);
      }
      else
      {
        triangulate_facet(f, nf_pmap);
      }
    }

    //The Lines
    {

        for(Selection_set_edge::iterator it = selected_edges.begin(); it != selected_edges.end(); ++it) {
            const Kernel::Point_3& a = (it->halfedge())->vertex()->point();
            const Kernel::Point_3& b = (it->halfedge())->opposite()->vertex()->point();
            positions_lines.push_back(a.x());
            positions_lines.push_back(a.y());
            positions_lines.push_back(a.z());

            positions_lines.push_back(b.x());
            positions_lines.push_back(b.y());
            positions_lines.push_back(b.z());
        }

    }

    //The points
    {
        for(Selection_set_vertex::iterator
            it = selected_vertices.begin(),
            end = selected_vertices.end();
            it != end; ++it)
        {
            const Kernel::Point_3& p = (*it)->point();
            positions_points.push_back(p.x());
            positions_points.push_back(p.y());
            positions_points.push_back(p.z());
        }
    }
}

void Scene_polyhedron_selection_item::draw(CGAL::Three::Viewer_interface* viewer) const
{

    if(!are_buffers_filled)
    {
        compute_elements();
        initialize_buffers(viewer);
    }

    draw_points(viewer);
    GLfloat offset_factor;
    GLfloat offset_units;
    viewer->glGetFloatv( GL_POLYGON_OFFSET_FACTOR, &offset_factor);
    viewer->glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
    glPolygonOffset(-1.f, 1.f);

    vaos[0]->bind();
    program = getShaderProgram(PROGRAM_WITH_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITH_LIGHT);
    program->bind();
    program->setAttributeValue("colors",this->color());
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_facets/3));
    program->release();
    vaos[0]->release();
    glPolygonOffset(offset_factor, offset_units);
    draw_edges(viewer);


}

void Scene_polyhedron_selection_item::draw_edges(CGAL::Three::Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
    {
        compute_elements();
        initialize_buffers(viewer);
    }

    viewer->glLineWidth(3.f);
    vaos[1]->bind();
    program = getShaderProgram(PROGRAM_NO_SELECTION);
    attrib_buffers(viewer,PROGRAM_NO_SELECTION);
    program->bind();

    program->setAttributeValue("colors",edge_color);
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_lines/3));
    program->release();
    vaos[1]->release();
    viewer->glLineWidth(1.f);
}

void Scene_polyhedron_selection_item::draw_points(CGAL::Three::Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
    {
        compute_elements();
        initialize_buffers(viewer);
    }
    viewer->glPointSize(5.f);
    vaos[2]->bind();
    program = getShaderProgram(PROGRAM_NO_SELECTION);
    attrib_buffers(viewer,PROGRAM_NO_SELECTION);
    program->bind();
    program->setAttributeValue("colors",vertex_color);
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(nb_points/3));
    program->release();
    vaos[2]->release();
    viewer->glPointSize(1.f);

}

void Scene_polyhedron_selection_item::inverse_selection()
{
  switch(k_ring_selector.active_handle_type)
  {
  case Active_handle::VERTEX:
  {
    Selection_set_vertex temp_select = selected_vertices;
    select_all();
    Q_FOREACH(Vertex_handle vh, temp_select)
    {
      selected_vertices.erase(vh);
    }
    break;
  }
  case Active_handle::EDGE:
  {
    Selection_set_edge temp_select = selected_edges;
    select_all();
    Q_FOREACH(edge_descriptor ed , temp_select)
      selected_edges.erase(ed);
    break;
  }
  default:
  {
    Selection_set_facet temp_select = selected_facets;
    select_all();
    Q_FOREACH(Facet_handle fh, temp_select)
      selected_facets.erase(fh);
    break;
  }
  }
  invalidateOpenGLBuffers();
  QGLViewer* v = *QGLViewer::QGLViewerPool().begin();
  v->update();
}
