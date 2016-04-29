#include "Scene_polyhedron_selection_item.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>


void Scene_polyhedron_selection_item::initialize_buffers(CGAL::Three::Viewer_interface *viewer)const
{
  //vao containing the data for the facets
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
  //vao containing the data for the  lines
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

void Scene_polyhedron_selection_item::initialize_temp_buffers(CGAL::Three::Viewer_interface *viewer)const
{
  //vao containing the data for the temp facets
  {
    program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
    program->bind();

    vaos[3]->bind();
    buffers[4].bind();
    buffers[4].allocate(positions_temp_facets.data(),
                        static_cast<int>(positions_temp_facets.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[4].release();



    buffers[5].bind();
    buffers[5].allocate(temp_normals.data(),
                        static_cast<int>(temp_normals.size()*sizeof(float)));
    program->enableAttributeArray("normals");
    program->setAttributeBuffer("normals",GL_FLOAT,0,3);
    buffers[5].release();

    vaos[3]->release();
    program->release();

  }
  //vao containing the data for the temp lines
  {
    program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
    program->bind();
    vaos[4]->bind();

    buffers[6].bind();
    buffers[6].allocate(positions_temp_lines.data(),
                        static_cast<int>(positions_temp_lines.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[6].release();

    program->release();

    vaos[4]->release();

  }
  //vao containing the data for the temp points
  {
    program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
    program->bind();
    vaos[5]->bind();

    buffers[7].bind();
    buffers[7].allocate(positions_temp_points.data(),
                        static_cast<int>(positions_temp_points.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[7].release();

    program->release();

    vaos[5]->release();
  }
  nb_temp_facets = positions_temp_facets.size();
  positions_temp_facets.resize(0);
  std::vector<float>(positions_temp_facets).swap(positions_temp_facets);

  temp_normals.resize(0);
  std::vector<float>(temp_normals).swap(temp_normals);

  nb_temp_lines = positions_temp_lines.size();
  positions_temp_lines.resize(0);
  std::vector<float>(positions_temp_lines).swap(positions_temp_lines);

  nb_temp_points = positions_temp_points.size();
  positions_temp_points.resize(0);
  std::vector<float>(positions_temp_points).swap(positions_temp_points);

  are_temp_buffers_filled = true;
}
void Scene_polyhedron_selection_item::initialize_HL_buffers(CGAL::Three::Viewer_interface *viewer)const
{
  //vao containing the data for the temp facets
  {
    program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
    program->bind();

    vaos[6]->bind();
    buffers[8].bind();
    buffers[8].allocate(positions_HL_facets.data(),
                        static_cast<int>(positions_HL_facets.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[8].release();



    buffers[9].bind();
    buffers[9].allocate(HL_normals.data(),
                        static_cast<int>(HL_normals.size()*sizeof(float)));
    program->enableAttributeArray("normals");
    program->setAttributeBuffer("normals",GL_FLOAT,0,3);
    buffers[9].release();

    vaos[6]->release();
    program->release();

  }
  //vao containing the data for the temp lines
  {
    program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
    program->bind();
    vaos[7]->bind();

    buffers[10].bind();
    buffers[10].allocate(positions_HL_lines.data(),
                        static_cast<int>(positions_HL_lines.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[10].release();

    program->release();

    vaos[7]->release();

  }
  //vao containing the data for the temp points
  {
    program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
    program->bind();
    vaos[8]->bind();

    buffers[11].bind();
    buffers[11].allocate(positions_HL_points.data(),
                        static_cast<int>(positions_HL_points.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[11].release();

    program->release();

    vaos[8]->release();
  }
  are_HL_buffers_filled = true;
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
typedef CGAL::Triangulation_2_projection_traits_3<Traits>   P_traits;
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
typedef CGAL::Exact_predicates_tag                                    Itag;
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
Scene_polyhedron_selection_item::triangulate_facet(Facet_handle fit,const FaceNormalPmap& fnmap,
                                                   std::vector<float> &p_facets,std::vector<float> &p_normals ) const
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

        push_back_xyz(ffit->vertex(0)->point(), p_facets);
        push_back_xyz(ffit->vertex(1)->point(), p_facets);
        push_back_xyz(ffit->vertex(2)->point(), p_facets);

        push_back_xyz(normal, p_normals);
        push_back_xyz(normal, p_normals);
        push_back_xyz(normal, p_normals);
    }
}


void Scene_polyhedron_selection_item::compute_any_elements(std::vector<float>& p_facets, std::vector<float>& p_lines, std::vector<float>& p_points, std::vector<float>& p_normals,
                                                           const Selection_set_vertex& p_sel_vertices, const Selection_set_facet& p_sel_facets, const Selection_set_edge& p_sel_edges)const
{
    p_facets.clear();
    p_lines.clear();
    p_points.clear();
    p_normals.clear();
    //The facet

if(!poly)
  return;
    for(Selection_set_facet::iterator
        it = p_sel_facets.begin(),
        end = p_sel_facets.end();
        it != end; it++)
    {
      Facet_handle f = (*it);
      if (f == boost::graph_traits<Polyhedron>::null_face())
        continue;

      if(is_triangle(f->halfedge(),*poly))
      {
        const Kernel::Vector_3 n =
            CGAL::Polygon_mesh_processing::compute_face_normal(f, *this->poly_item->polyhedron());
        p_normals.push_back(n.x());
        p_normals.push_back(n.y());
        p_normals.push_back(n.z());

        p_normals.push_back(n.x());
        p_normals.push_back(n.y());
        p_normals.push_back(n.z());

        p_normals.push_back(n.x());
        p_normals.push_back(n.y());
        p_normals.push_back(n.z());


        Polyhedron::Halfedge_around_facet_circulator
            he = f->facet_begin(),
            cend = he;

        CGAL_For_all(he,cend)
        {
          const Kernel::Point_3& p = he->vertex()->point();
          p_facets.push_back(p.x());
          p_facets.push_back(p.y());
          p_facets.push_back(p.z());
        }
      }
      else if (is_quad(f->halfedge(), *poly))
      {
        Vector nf = get(nf_pmap, f);
        //1st half-quad
        Point p0 = f->halfedge()->vertex()->point();
        Point p1 = f->halfedge()->next()->vertex()->point();
        Point p2 = f->halfedge()->next()->next()->vertex()->point();

        push_back_xyz(p0, p_facets);
        push_back_xyz(p1, p_facets);
        push_back_xyz(p2, p_facets);

        push_back_xyz(nf, p_normals);
        push_back_xyz(nf, p_normals);
        push_back_xyz(nf, p_normals);

        //2nd half-quad
        p0 = f->halfedge()->next()->next()->vertex()->point();
        p1 = f->halfedge()->prev()->vertex()->point();
        p2 = f->halfedge()->vertex()->point();

        push_back_xyz(p0, p_facets);
        push_back_xyz(p1, p_facets);
        push_back_xyz(p2, p_facets);

        push_back_xyz(nf, p_normals);
        push_back_xyz(nf, p_normals);
        push_back_xyz(nf, p_normals);
      }
      else
      {
        triangulate_facet(f, nf_pmap, p_facets, p_normals);
      }
    }

    //The Lines
    {

        for(Selection_set_edge::iterator it = p_sel_edges.begin(); it != p_sel_edges.end(); ++it) {
            const Kernel::Point_3& a = (it->halfedge())->vertex()->point();
            const Kernel::Point_3& b = (it->halfedge())->opposite()->vertex()->point();
            p_lines.push_back(a.x());
            p_lines.push_back(a.y());
            p_lines.push_back(a.z());

            p_lines.push_back(b.x());
            p_lines.push_back(b.y());
            p_lines.push_back(b.z());
        }

    }
    //The points
    {
        for(Selection_set_vertex::iterator
            it = p_sel_vertices.begin(),
            end = p_sel_vertices.end();
            it != end; ++it)
        {
            const Kernel::Point_3& p = (*it)->point();
            p_points.push_back(p.x());
            p_points.push_back(p.y());
            p_points.push_back(p.z());
        }
    }
}
void Scene_polyhedron_selection_item::compute_elements()const
{
  compute_any_elements(positions_facets, positions_lines, positions_points, normals,
                       selected_vertices, selected_facets, selected_edges);
}
void Scene_polyhedron_selection_item::compute_temp_elements()const
{
  compute_any_elements(positions_temp_facets, positions_temp_lines, positions_temp_points, temp_normals,
                       temp_selected_vertices, temp_selected_facets, temp_selected_edges);
}

void Scene_polyhedron_selection_item::compute_HL_elements()const
{
  compute_any_elements(positions_HL_facets, positions_HL_lines, positions_HL_points, HL_normals,
                       HL_selected_vertices, HL_selected_facets, HL_selected_edges);
}

void Scene_polyhedron_selection_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  GLfloat offset_factor;
  GLfloat offset_units;

  if(!are_HL_buffers_filled)
  {
    compute_HL_elements();
    initialize_HL_buffers(viewer);
  }

  viewer->glGetFloatv(GL_POLYGON_OFFSET_FACTOR, &offset_factor);
  viewer->glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
  glPolygonOffset(-0.1f, 0.2f);
  vaos[6]->bind();
  program = getShaderProgram(PROGRAM_WITH_LIGHT);
  attrib_buffers(viewer,PROGRAM_WITH_LIGHT);
  program->bind();
  program->setAttributeValue("colors",QColor(255,153,51));
  viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(positions_HL_facets.size())/3);
  program->release();
  vaos[6]->release();


  if(!are_temp_buffers_filled)
  {
      compute_temp_elements();
      initialize_temp_buffers(viewer);
  }
  vaos[3]->bind();
  program = getShaderProgram(PROGRAM_WITH_LIGHT);
  attrib_buffers(viewer,PROGRAM_WITH_LIGHT);
  program->bind();
  program->setAttributeValue("colors",QColor(0,255,0));
  viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_temp_facets/3));
  program->release();
  vaos[3]->release();

    if(!are_buffers_filled)
    {
        compute_elements();
        initialize_buffers(viewer);
    }

    draw_points(viewer);
    viewer->glGetFloatv( GL_POLYGON_OFFSET_FACTOR, &offset_factor);
    viewer->glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);

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

  viewer->glLineWidth(3.f);

  if(!are_HL_buffers_filled)
  {
    compute_HL_elements();
    initialize_HL_buffers(viewer);
  }

  vaos[7]->bind();
  program = getShaderProgram(PROGRAM_NO_SELECTION);
  attrib_buffers(viewer,PROGRAM_NO_SELECTION);
  program->bind();

  program->setAttributeValue("colors",QColor(255,153,51));
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_HL_lines.size()/3));
  program->release();
  vaos[7]->release();


  if(!are_temp_buffers_filled)
  {
    compute_temp_elements();
    initialize_temp_buffers(viewer);
  }

  vaos[4]->bind();
  program = getShaderProgram(PROGRAM_NO_SELECTION);
  attrib_buffers(viewer,PROGRAM_NO_SELECTION);
  program->bind();

  program->setAttributeValue("colors",QColor(0,200,0));
  viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_temp_lines/3));
  program->release();
  vaos[4]->release();
  viewer->glLineWidth(3.0f);
  if(!are_buffers_filled)
  {
    compute_elements();
    initialize_buffers(viewer);
  }

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

  viewer->glPointSize(5.5f);

  if(!are_HL_buffers_filled)
  {
    compute_HL_elements();
    initialize_HL_buffers(viewer);
  }
  vaos[8]->bind();
  program = getShaderProgram(PROGRAM_NO_SELECTION);
  attrib_buffers(viewer,PROGRAM_NO_SELECTION);
  program->bind();
  program->setAttributeValue("colors",QColor(255,153,51));
  viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(positions_HL_points.size()/3));
  program->release();
  vaos[8]->release();


  if(!are_temp_buffers_filled)
  {
    compute_temp_elements();
    initialize_temp_buffers(viewer);
  }
  vaos[5]->bind();
  program = getShaderProgram(PROGRAM_NO_SELECTION);
  attrib_buffers(viewer,PROGRAM_NO_SELECTION);
  program->bind();
  program->setAttributeValue("colors",QColor(0,50,0));
  viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(nb_temp_points/3));
  program->release();
  vaos[5]->release();

  viewer->glPointSize(5.5f);
  if(!are_buffers_filled)
  {
    compute_elements();
    initialize_buffers(viewer);
  }
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

void Scene_polyhedron_selection_item::set_operation_mode(int mode)
{
  k_ring_selector.setEditMode(true);
  Q_EMIT updateInstructions(QString("SHIFT + left click to apply operation."));
  switch(mode)
  {
  case -1:
    //restore original selection_type
    set_active_handle_type(original_sel_mode);
    clearHL();
    k_ring_selector.setEditMode(false);
    break;
    //Join vertex
  case 0:
    Q_EMIT updateInstructions("Select the edge with extremities you want to join. (1/2)");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Split vertex
  case 1:
    Q_EMIT updateInstructions("Select the vertex you want to split. (1/3)");
    //set the selection type to Vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
    //Split edge
  case 2:
    Q_EMIT updateInstructions("Select the edge you want to split.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Join face
  case 3:
    Q_EMIT updateInstructions("Select the edge separating the faces you want to join.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Split face
  case 4:
    Q_EMIT updateInstructions("Select the facet you want to split (degree >= 4). (1/3)");
    //set the selection type to Facet
    set_active_handle_type(static_cast<Active_handle::Type>(1));
    break;
    //Collapse edge
  case 5:
    Q_EMIT updateInstructions("Select the edge you want to collapse.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Flip edge
  case 6:
    Q_EMIT updateInstructions("Select the edge you want to flip.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Add center vertex
  case 7:
    Q_EMIT updateInstructions("Select a facet.");
    //set the selection type to Facet
    set_active_handle_type(static_cast<Active_handle::Type>(1));
    break;
    //Remove center vertex
  case 8:
    Q_EMIT updateInstructions("Select the vertex you want to remove.");
    //set the selection type to vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
    //Add vertex and face to border
  case 9:
    Q_EMIT updateInstructions("Select a border edge. (1/2)");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Add face to border
  case 10:
    Q_EMIT updateInstructions("Select a border edge. (1/2)");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
  default:
    break;
  }
  operation_mode = mode;
}
template<typename HandleRange>
bool Scene_polyhedron_selection_item::treat_classic_selection(const HandleRange& selection)
{
  typedef typename HandleRange::value_type HandleType;
  Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);
  bool any_change = false;
  if(is_insert) {
    BOOST_FOREACH(HandleType h, selection)
        any_change |= tr.container().insert(h).second;
  }
  else{
    BOOST_FOREACH(HandleType h, selection)
        any_change |= (tr.container().erase(h)!=0);
  }
  if(any_change) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
  return any_change;
}

bool Scene_polyhedron_selection_item::treat_selection(const std::set<Polyhedron::Vertex_handle>& selection)
{
  if(!is_treated)
  {
    Vertex_handle vh = *selection.begin();
    Selection_traits<Vertex_handle, Scene_polyhedron_selection_item> tr(this);
    switch(operation_mode)
    {
    //classic selection
    case -1:
    {
      return treat_classic_selection(selection);
      break;
    }
      //Join vertex
    case 0:
    {
      bool belong = false;
      Halfedge_handle target = halfedge(to_join_ed, *polyhedron());
      if(halfedge(to_join_ed, *polyhedron())->vertex() == vh)
        belong = true;
      if(halfedge(to_join_ed, *polyhedron())->opposite()->vertex() == vh)
      {
        belong = true;
        target = halfedge(to_join_ed, *polyhedron())->opposite();
      }
      if(!belong)
      {
        tempInstructions("Vertices not joined : the vertex must belong to the selected edge.",
                         "Select the vertex that will remain. (2/2)");
      }
      else
      {

        polyhedron()->join_vertex(target);

        temp_selected_edges.clear();
        //set to select edge
        set_active_handle_type(static_cast<Active_handle::Type>(2));
        tempInstructions("Vertices joined.",
                         "Select the edge with extremities you want to join. (1/2)");
        invalidateOpenGLBuffers();
        polyhedron_item()->invalidateOpenGLBuffers();
      }
      break;
    }
      //Split vertex
    case 1:
    {
      //save VH
      to_split_vh = vh;
      temp_selected_vertices.insert(to_split_vh);
      //set to select facet
      set_active_handle_type(static_cast<Active_handle::Type>(1));
      invalidateOpenGLBuffers();
      Q_EMIT updateInstructions("Select first facet. (2/3)");
      break;
    }
      //Split face
    case 4:
    {
      static Vertex_handle s;
      static Polyhedron::Halfedge_handle h1,h2;
      static bool found_h1(false), found_h2(false);
      if(!first_selected)
      {
          //Is the vertex on the face ?
          Polyhedron::Halfedge_around_facet_circulator hafc = to_split_fh->facet_begin();
          Polyhedron::Halfedge_around_facet_circulator end = hafc;
          CGAL_For_all(hafc, end)
          {
            if(hafc->vertex()==vh)
            {
              h1 = hafc;
              s = vh;
              found_h1 = true;
                break;
            }
          }
          if(!found_h1)
          {
            tempInstructions("Vertex not selected : The vertex is not on the face.",
                             "Select the first vertex. (2/3)");
          }
          else
          {
            first_selected = true;
            temp_selected_vertices.insert(s);
            invalidateOpenGLBuffers();
            Q_EMIT updateInstructions("Select the second vertex (3/3)");
          }
      }
      else
      {
        bool is_same(false), are_next(false);
        Polyhedron::Halfedge_around_facet_circulator hafc = to_split_fh->facet_begin();
        Polyhedron::Halfedge_around_facet_circulator end = hafc;
        for(int i=0; i<1; i++) //seems useless but allow the use of break.
        {
          //Is the vertex on the face ?
          CGAL_For_all(hafc, end)

              if(hafc->vertex()==vh)
          {
            h2 = hafc;
            found_h2 = true;
            break;
          }
          if(!found_h2)
          {
            break;
          }
          //Are they different ?
          if(h1 == h2)
          {
            is_same = true;
            break;
          }
          is_same = false;
          //Are they directly following each other?
          if(next(h1, *polyhedron()) == h2 ||
             next(h2, *polyhedron()) == h1)
          {
            are_next = true;
            break;
          }
          are_next = false;
        }
        if(!found_h2)
          tempInstructions("Vertex not selected : The vertex is not on the face.",
                           "Select the second vertex (3/3).");
        else if(is_same)
          tempInstructions("Vertex not selected : The vertices must be different.",
                           "Select the second vertex (3/3).");
        else if(are_next)
          tempInstructions("Vertex not selected : The vertices must not directly follow each other.",
                           "Select the second vertex (3/3).");
        else
        {
          CGAL::Euler::split_face(h1,h2, *polyhedron());
          first_selected = false;
          temp_selected_vertices.clear();
          temp_selected_facets.clear();
          invalidateOpenGLBuffers();
          //reset selection type to Facet
          set_active_handle_type(static_cast<Active_handle::Type>(1));
          tempInstructions("Face split.",
                           "Select a facet (1/3).");
          polyhedron_item()->invalidateOpenGLBuffers();
        }
      }
      break;
    }
      //Remove center vertex
    case 8:

        bool has_hole = false;
        Polyhedron::Halfedge_around_vertex_circulator hc = vh->vertex_begin();
        Polyhedron::Halfedge_around_vertex_circulator end(hc);
        CGAL_For_all(hc, end)
        {
          if(hc->is_border())
          {
            has_hole = true;
            break;
          }
        }
        if(!has_hole)
        {
          CGAL::Euler::remove_center_vertex(vh->halfedge(),*polyhedron());
          polyhedron_item()->invalidateOpenGLBuffers();
        }
        else
        {
          tempInstructions("Vertex not selected : There must be no hole incident to the selection.",
                           "Select the vertex you want to remove.");
        }
      break;

    }
  }
  is_treated = true;
  return false;
}

//returns true if halfedge's facet's degree >= degree

int facet_degree(Halfedge_handle h)
{
  if(h->is_border())
  {
    Halfedge_handle it = h;
    int deg =0;
    do
    {
      deg ++;
      it=it->next();
    }
    while(it != h);
    return deg;
  }
  else
    return h->facet()->facet_degree();
}

bool Scene_polyhedron_selection_item:: treat_selection(const std::set<edge_descriptor>& selection)
{
  edge_descriptor ed =  *selection.begin();
  if(!is_treated)
  {
    Selection_traits<edge_descriptor, Scene_polyhedron_selection_item> tr(this);
    switch(operation_mode)
    {
    //classic selection
    case -1:
    {
      return treat_classic_selection(selection);
      break;
    }
      //Join vertex
    case 0:
        if(facet_degree(halfedge(ed, *polyhedron())) < 4
           ||
           facet_degree(halfedge(ed, *polyhedron())->opposite())< 4)
        {
          tempInstructions("Edge not selected: the incident facets must have a degree of at least 4.",
                           "Select the edge with extremities you want to join.(1/2)");
        }
        else
        {
          to_join_ed = ed;
          temp_selected_edges.insert(to_join_ed);
          invalidateOpenGLBuffers();
          //set to select vertex
          set_active_handle_type(static_cast<Active_handle::Type>(0));
          Q_EMIT updateInstructions("Select the vertex that will remain.");
        }
      break;
      //Split edge
    case 2:
    {
        Polyhedron::Point_3 a(halfedge(ed, *polyhedron())->vertex()->point()),b(halfedge(ed, *polyhedron())->opposite()->vertex()->point());
        Polyhedron::Halfedge_handle hhandle = polyhedron()->split_edge(halfedge(ed, *polyhedron()));
        Polyhedron::Point_3 p((b.x()+a.x())/2.0, (b.y()+a.y())/2.0,(b.z()+a.z())/2.0);

        hhandle->vertex()->point() = p;
        selected_vertices.insert(hhandle->vertex());
        invalidateOpenGLBuffers();
        poly_item->invalidateOpenGLBuffers();
      tempInstructions("Edge splitted.",
                       "Select the edge you want to split.");
      break;
    }
      //Join face
    case 3:
        if(out_degree(source(halfedge(ed,*polyhedron()),*polyhedron()),*polyhedron())<3 ||
           out_degree(target(halfedge(ed,*polyhedron()),*polyhedron()),*polyhedron())<3)
          tempInstructions("Faces not joined : the two ends of the edge must have a degree of at least 3.",
                           "Select the edge separating the faces you want to join.");
        else
        {
          polyhedron()->join_facet(halfedge(ed, *polyhedron()));
          poly_item->invalidateOpenGLBuffers();
        }
      break;
      //Collapse edge
    case 5:
        if(!is_triangle_mesh(*polyhedron()))
        {
          tempInstructions("Edge not collapsed : the graph must be triangulated.",
                           "Select the edge you want to collapse.");
        }
        else if(!CGAL::Euler::does_satisfy_link_condition(ed, *polyhedron()))
        {
          tempInstructions("Edge not collapsed : link condition not satidfied.",
                           "Select the edge you want to collapse.");
        }
        else
        {
          CGAL::Euler::collapse_edge(ed, *polyhedron());
          polyhedron_item()->invalidateOpenGLBuffers();
          tempInstructions("Edge collapsed.",
                           "Select the edge you want to collapse.");
        }
      break;
      //Flip edge
    case 6:

        //check preconditions
        if(facet_degree(halfedge(ed, *polyhedron())) == 3 && facet_degree(halfedge(ed, *polyhedron())->opposite()) == 3)
        {
          CGAL::Euler::flip_edge(halfedge(ed, *polyhedron()), *polyhedron());
          polyhedron_item()->invalidateOpenGLBuffers();
        }
        else
        {
          tempInstructions("Edge not selected : incident facets must be triangles.",
                           "Select the edge you want to flip.");
        }

      break;
      //Add vertex and face to border
    case 9:
    {
      static Halfedge_handle t;
      if(!first_selected)
      {
          bool found = false;
          Halfedge_handle hc = halfedge(ed, *polyhedron());
          if(hc->is_border())
          {
            t = hc;
            found = true;
          }
          else if(hc->opposite()->is_border())
          {
            t = hc->opposite();
            found = true;
          }
          if(found)
          {
            first_selected = true;
            temp_selected_edges.insert(edge(t, *polyhedron()));
            temp_selected_vertices.insert(t->vertex());
            invalidateOpenGLBuffers();
            Q_EMIT updateInstructions("Select second edge. (2/2)");
          }
          else
          {
            tempInstructions("Edge not selected : no border found.",
                             "Select a border edge. (1/2)");
          }
      }
      else
      {
        bool found(false), is_equal(true), is_border(false);
          Halfedge_handle hc = halfedge(ed, *polyhedron());
          //seems strange but allows to use break efficiently
          for(int i=0; i< 2; i++)
          {
            //if the selected halfedge is not a border, stop and signal it.
            if(hc->is_border())
              is_border = true;
            else if(hc->opposite()->is_border())
            {
              hc = halfedge(ed, *polyhedron())->opposite();
              is_border = true;
            }
            if(!is_border)
              break;
            //if the halfedges are the same, stop and signal it.
            if(hc == t)
            {
              is_equal = true;
              break;
            }
            is_equal = false;
            //if the halfedges are not on the same border, stop and signal it.
            boost::graph_traits<Polyhedron>::halfedge_descriptor iterator = next(t, *polyhedron());
            while(iterator != t)
            {
              if(iterator == hc)
              {
                found = true;
                Halfedge_handle res = CGAL::Euler::add_vertex_and_face_to_border(t,hc, *polyhedron());
                //res seems to be the opposite of what it is said to be in the doc.
                if(res->vertex() == hc->vertex())
                  res = res->opposite();
                //create and add a point point

                Polyhedron::Point_3 a = t->vertex()->point();
                Polyhedron::Point_3 b = t->opposite()->vertex()->point();
                Polyhedron::Point_3 c = t->opposite()->next()->vertex()->point();
                double x = b.x()+a.x()-c.x() ;
                double y = b.y()+a.y()-c.y() ;
                double z = b.z()+a.z()-c.z() ;
                res->vertex()->point() = Polyhedron::Point_3(x,y,z);
                break;
              }
              iterator = next(iterator, *polyhedron());
            }
          }
        if(is_equal)
        {
          tempInstructions("Edge not selected : halfedges must be different.",
                           "Select the second edge.");
        }
        else if(!is_border || !found)
        {
          tempInstructions("Edge not selected : no shared border found.",
                           "Select the second edge.");
        }
        else
        {
          first_selected = false;


          temp_selected_edges.clear();
          temp_selected_vertices.clear();
          invalidateOpenGLBuffers();
          polyhedron_item()->invalidateOpenGLBuffers();
          tempInstructions("Face and vertex added.",
                           "Select a border edge. (1/2)");
        }
      }
      break;
    }
      //Add face to border
    case 10:
    {
      static Halfedge_handle t;
      if(!first_selected)
      {
          bool found = false;
          Halfedge_handle hc = halfedge(ed, *polyhedron());
          if(hc->is_border())
          {
            t = hc;
            found = true;
          }
          else if(hc->opposite()->is_border())
          {
            t = hc->opposite();
            found = true;
          }
          if(found)
          {
            first_selected = true;
            temp_selected_edges.insert(edge(t, *polyhedron()));
            temp_selected_vertices.insert(t->vertex());
            invalidateOpenGLBuffers();
            Q_EMIT updateInstructions("Select second edge. (2/2)");
            set_active_handle_type(static_cast<Active_handle::Type>(2));
          }
          else
          {
            tempInstructions("Edge not selected : no border found.",
                             "Select a border edge. (1/2)");
          }
      }
      else
      {
        bool found(false), is_equal(true), is_next(true), is_border(false);
          Halfedge_handle hc = halfedge(ed, *polyhedron());
          //seems strange but allows to use break efficiently
          for(int i= 0; i<1; i++)
          {
            //if the selected halfedge is not a border, stop and signal it.
            if(hc->is_border())
              is_border = true;
            else if(hc->opposite()->is_border())
            {
              hc = hc->opposite();
              is_border = true;
            }
            if(!is_border)
              break;
            //if the halfedges are the same, stop and signal it.
            if(hc == t)
            {
              is_equal = true;
              break;
            }
            is_equal = false;
            //if the halfedges are adjacent, stop and signal it.
            if(next(t, *polyhedron()) == hc || next(hc, *polyhedron()) == t)
            {
              is_next = true;
              break;
            }
            is_next = false;
            //if the halfedges are not on the same border, stop and signal it.
            boost::graph_traits<Polyhedron>::halfedge_descriptor iterator = next(t, *polyhedron());
            while(iterator != t)
            {
              if(iterator == hc)
              {
                found = true;
                CGAL::Euler::add_face_to_border(t,hc, *polyhedron());
                break;
              }
              iterator = next(iterator, *polyhedron());
            }
          }
        if(is_equal)
        {
          tempInstructions("Edge not selected : halfedges must be different.",
                           "Select the second edge. (2/2)");
        }

        else if(is_next)
        {
          tempInstructions("Edge not selected : halfedges must not be adjacent.",
                           "Select the second edge. (2/2)");
        }
        else if(!is_border || !found)
        {
          tempInstructions("Edge not selected : no shared border found.",
                           "Select the second edge. (2/2)");
        }
        else
        {
          first_selected = false;
          temp_selected_vertices.clear();
          temp_selected_edges.clear();
          invalidateOpenGLBuffers();
          polyhedron_item()->invalidateOpenGLBuffers();
          tempInstructions("Face added.",
                           "Select a border edge. (1/2)");
        }
      }
      break;
    }
    }
  }
  is_treated = true;
  return false;
}

bool Scene_polyhedron_selection_item::treat_selection(const std::vector<Polyhedron::Facet_handle>& selection)
{
  return treat_classic_selection(selection);
}

bool Scene_polyhedron_selection_item::treat_selection(const std::set<Polyhedron::Facet_handle>& selection)
{
  if(!is_treated)
  {
    Facet_handle fh = *selection.begin();
    Selection_traits<Facet_handle, Scene_polyhedron_selection_item> tr(this);
    switch(operation_mode)
    {
    //classic selection
    case -1:
    {
      return treat_classic_selection(selection);
      break;
    }
    //Split vertex
    case 1:
    {
      static Polyhedron::Halfedge_handle h1;
      //stores first fh and emit change label
      if(!first_selected)
      {
          bool found = false;
          //test preco
          Polyhedron::Halfedge_around_facet_circulator hafc = fh->facet_begin();
          Polyhedron::Halfedge_around_facet_circulator end = hafc;
          CGAL_For_all(hafc, end)
          {
            if(hafc->vertex()==to_split_vh)
            {
              h1 = hafc;
              found = true;
              break;
            }
          }
          if(found)
          {
            first_selected = true;
            temp_selected_facets.insert(fh);
            invalidateOpenGLBuffers();
            Q_EMIT updateInstructions("Select the second facet. (3/3)");
          }
          else
            tempInstructions("Facet not selected : no valid halfedge",
                             "Select first facet. (2/3)");
      }
      //call the function with point and facets.
      else
      {
          //get the right halfedges
          Polyhedron::Halfedge_handle h2;
          bool found = false;
          Polyhedron::Halfedge_around_facet_circulator hafc = fh->facet_begin();
          Polyhedron::Halfedge_around_facet_circulator end = hafc;
          CGAL_For_all(hafc, end)
          {
            if(hafc->vertex()==to_split_vh)
            {
              h2 = hafc;
              found = true;
              break;
            }
          }

          if(found &&(h1 != h2))
          {
            Polyhedron::Halfedge_handle hhandle = CGAL::Euler::split_vertex(h1,h2,*polyhedron());

            temp_selected_facets.clear();
            Polyhedron::Point_3 p1t = h1->vertex()->point();
            Polyhedron::Point_3 p1s = h1->opposite()->vertex()->point();
            double x =  p1t.x() + 0.01 * (p1s.x() - p1t.x());
            double y =  p1t.y() + 0.01 * (p1s.y() - p1t.y());
            double z =  p1t.z() + 0.01 * (p1s.z() - p1t.z());
            hhandle->opposite()->vertex()->point() = Polyhedron::Point_3(x,y,z);;
            first_selected = false;
            temp_selected_vertices.clear();
            invalidateOpenGLBuffers();
            //reset selection mode
            set_active_handle_type(static_cast<Active_handle::Type>(0));
            poly_item->invalidateOpenGLBuffers();
            tempInstructions("Vertex splitted.", "Select the vertex you want splitted. (1/3)");
          }
          else if(h1 == h2)
          {
             tempInstructions("Facet not selected : same as the first.", "Select the second facet. (3/3)");
          }
          else
          {
            tempInstructions("Facet not selected : no valid halfedge.", "Select the second facet. (3/3)");
          }
      }
      break;
    }
      //Split face
    case 4:
      if(is_triangle(fh->halfedge(), *poly))
      {
        tempInstructions("Facet not selected : Facet must not be a triangle.",
                         "Select the facet you want to split (degree >= 4). (1/3)");
      }
      else
      {
        to_split_fh = fh;
        temp_selected_facets.insert(to_split_fh);
        invalidateOpenGLBuffers();
        //set to select vertex
        set_active_handle_type(static_cast<Active_handle::Type>(0));
        Q_EMIT updateInstructions("Select first vertex. (2/3)");
      }
      break;
      //Add center vertex
    case 7:
        if(fh->halfedge()->is_border())
        {
          tempInstructions("Facet not selected : Facet must not be null.",
                           "Select a Facet.");
        }
        else
        {
          Polyhedron::Halfedge_around_facet_circulator hafc = fh->facet_begin();
          Polyhedron::Halfedge_around_facet_circulator end = hafc;

          double x(0), y(0), z(0);
          int total(0);
          CGAL_For_all(hafc, end)
          {
            x+=hafc->vertex()->point().x(); y+=hafc->vertex()->point().y(); z+=hafc->vertex()->point().z();
            total++;
          }
          Polyhedron::Halfedge_handle hhandle = CGAL::Euler::add_center_vertex(fh->facet_begin(), *polyhedron());
          if(total !=0)
            hhandle->vertex()->point() = Polyhedron::Point_3(x/(double)total, y/(double)total, z/(double)total);
          poly_item->invalidateOpenGLBuffers();

        }
      break;
    }
  }
  is_treated = true;
  return false;
}

void Scene_polyhedron_selection_item::tempInstructions(QString s1, QString s2)
{
  m_temp_instructs = s2;
  Q_EMIT updateInstructions(QString("<font color='red'>%1</font>").arg(s1));
  QTimer timer;
  timer.singleShot(5500, this, SLOT(emitTempInstruct()));
}
void Scene_polyhedron_selection_item::emitTempInstruct()
{
  Q_EMIT updateInstructions(QString("<font color='black'>%1</font>").arg(m_temp_instructs));
}

void Scene_polyhedron_selection_item::on_Ctrlz_pressed()
{
  first_selected = false;
  temp_selected_vertices.clear();
  temp_selected_edges.clear();
  temp_selected_facets.clear();
  are_temp_buffers_filled = false;
  set_operation_mode(operation_mode);
  Q_EMIT itemChanged();
}

void Scene_polyhedron_selection_item::compute_normal_maps()
{

  face_normals_map.clear();
  vertex_normals_map.clear();
  nf_pmap = boost::associative_property_map< boost::container::flat_map<boost::graph_traits<Polyhedron>::face_descriptor, Kernel::Vector_3> >(face_normals_map);
  nv_pmap = boost::associative_property_map< boost::container::flat_map<boost::graph_traits<Polyhedron>::vertex_descriptor, Kernel::Vector_3> >(vertex_normals_map);
  PMP::compute_normals(*poly, nv_pmap, nf_pmap);
}
