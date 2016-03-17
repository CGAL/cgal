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

void Scene_polyhedron_selection_item::set_operation_mode(int mode)
{
  k_ring_selector.setEditMode(true);
  Q_EMIT updateInstructions(QString("SHIFT + left click to apply operation."));
  switch(mode)
  {
  case -1:
    //restore original selection_type
    set_active_handle_type(original_sel_mode);
    k_ring_selector.setEditMode(false);
    break;
    //Join vertex
  case 0:
    Q_EMIT updateInstructions("Select the edge with extremities you want to join.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Split vertex
  case 1:
    //set the selection type to Vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    Q_EMIT updateInstructions("Select the vertex you want to split.");
    break;
    //Split edge
  case 2:
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
    //set the selection type to Vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
    //Add edge
  case 5:
    //set the selection type to Vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
    //Collapse edge
  case 6:
    //set the selection type to Edge
    Q_EMIT updateInstructions("Select the edge you want to collapse.");
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Flip edge
  case 7:
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    Q_EMIT updateInstructions("Select the edge you want to flip.");
    break;
    //Add center vertex
  case 8:
    //set the selection type to Facet
    set_active_handle_type(static_cast<Active_handle::Type>(1));
    break;
    //Remove center vertex
  case 9:
    Q_EMIT updateInstructions("Select the vertex you want to remove.");
    //set the selection type to vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
    //Add vertex and face to border
  case 10:
    //set the selection type to vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
    //Add face to border
  case 11:
    Q_EMIT updateInstructions("Select a vertex on a border.");
    //set the selection type to vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
  default:
    break;
  }
  operation_mode = mode;
}

bool Scene_polyhedron_selection_item::treat_selection(const std::set<Polyhedron::Vertex_handle>& selection)
{
  if(!is_treated)
  {
    static bool first_selected = false;
    static Vertex_handle s;
    static Polyhedron::Halfedge_handle t;
    Selection_traits<Vertex_handle, Scene_polyhedron_selection_item> tr(this);
    switch(operation_mode)
    {
    //classic selection
    case -1:
    {

      bool any_change = false;
      if(is_insert) {
        BOOST_FOREACH(Vertex_handle vh, selection)
            any_change |= tr.container().insert(vh).second;
      }
      else{
        BOOST_FOREACH(Vertex_handle vh, selection)
            any_change |= (tr.container().erase(vh)!=0);
      }
      if(any_change) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
      return any_change;
      break;
    }
      //Split vertex
    case 1:
    {
      //save VH
      BOOST_FOREACH(Vertex_handle vh, selection)
          to_split_vh = vh;
      selected_vertices.insert(to_split_vh);
      //set to select facet
      set_active_handle_type(static_cast<Active_handle::Type>(1));
      Q_EMIT updateInstructions("Select first facet.");
      break;
    }
      //Split face
    case 4:
      if(!first_selected)
      {
        BOOST_FOREACH(Vertex_handle vh, selection)
        {
          s = vh;
        }
        first_selected = true;
        selected_vertices.insert(s);
      }
      else
      {
        BOOST_FOREACH(Vertex_handle vh, selection)
        {
          CGAL::Euler::split_face(s->halfedge(),vh->halfedge(), *poly);
        }
        first_selected = false;
        selected_vertices.erase(s);
      }
      break;
      //Add edge
    case 5:
      if(!first_selected)
      {
        BOOST_FOREACH(Vertex_handle vh, selection)
        {
          s = vh;
          first_selected = true;
        }
        selected_vertices.insert(s);
      }
      else
      {
        BOOST_FOREACH(Vertex_handle vh, selection)
        {
          CGAL::Euler::add_edge(s,vh,*poly);
        }
        first_selected = false;
        selected_vertices.erase(s);
      }
      break;
      //Remove center vertex
    case 9:
      BOOST_FOREACH(Vertex_handle vh, selection)
      {
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
          CGAL::Euler::remove_center_vertex(vh->halfedge(),*poly);
        else
        {
          tempInstructions("Vertex not selected : There must be no hole incident to the selection.",
                           "Select the vertex you want to remove.");
        }
      }
      //Avoids a segfault in Scene_polyhedron_item::select()
      Q_EMIT skipEmits(true);
      break;
      //Add vertex and face to border
    case 10:
      if(!first_selected)
      {
        BOOST_FOREACH(Vertex_handle vh, selection)
        {
          Polyhedron::Halfedge_around_vertex_circulator hc = vh->vertex_begin();
          Polyhedron::Halfedge_around_vertex_circulator end(hc);
          CGAL_For_all(hc, end)
          {
            if(hc->is_border())
            {
              t = hc;
              break;
            }
          }
          first_selected = true;
          selected_vertices.insert(t->vertex());
        }
      }
      else
      {
        BOOST_FOREACH(Vertex_handle vh, selection)
        {

            Q_FOREACH( boost::graph_traits<Polyhedron>::halfedge_descriptor hc,
                         CGAL::halfedges_around_source(vh->halfedge(), *polyhedron()))
            {
            if(hc->is_border())
            {
              CGAL::Euler::add_vertex_and_face_to_border(t,hc, *poly);
              break;
            }
          }
        }
        first_selected = false;
        selected_vertices.erase(t->vertex());
      }
      break;
      //Add face to border
    case 11:
      if(!first_selected)
      {
        BOOST_FOREACH(Vertex_handle vh, selection)
        {
          bool found = false;
          Q_FOREACH( boost::graph_traits<Polyhedron>::halfedge_descriptor hc,
                     CGAL::halfedges_around_target(vh, *polyhedron()))
          {
            if(hc->is_border())
            {
              t = hc;
              found = true;
              break;
            }
          }
          if(found)
          {
            first_selected = true;
            selected_vertices.insert(t->vertex());
            selected_edges.insert(edge(t,*polyhedron()));
            Q_EMIT updateInstructions("Select second facet.");
          }
          else
          {
            tempInstructions("Vertex not selected : no border found.",
                             "Select a vertex on a border.");
          }
        }
      }
      else
      {

        bool found(false), is_equal(true), is_next(true), is_border(false);
        BOOST_FOREACH(Vertex_handle vh, selection)
        {
          Q_FOREACH( boost::graph_traits<Polyhedron>::halfedge_descriptor hc,
                     CGAL::halfedges_around_target(vh, *polyhedron()))
          {
            //if the selected halfedge is not a border, stop and signal it.
            if(!hc->is_border())
            {
              is_border = false;
              continue;
            }
            is_border = true;
            //if the halfedges are the same, stop and signal it.
            if(hc == t)
            {
              is_equal = true;
              continue;
            }
            is_equal = false;
            //if the halfedges are adjacent, stop and signal it.
            if(next(t, *polyhedron()) == hc || next(hc, *polyhedron()) == t)
            {
              is_next = true;
              continue;
            }
            is_next = false;
            //if the halfedges are not on the same border, stop and signal it.
            boost::graph_traits<Polyhedron>::halfedge_descriptor iterator = next(t, *polyhedron());
            while(iterator != t)
            {
              if(iterator == hc)
              {
                found = true;
                CGAL::Euler::add_face_to_border(t,hc, *poly);
                break;
              }
              iterator = next(iterator, *polyhedron());
            }
            //if an  halfedges pass all the testes, stop and keep it.
            if(found)
              break;
          }
        }
        if(is_equal)
        {
          tempInstructions("Vertex not selected : halfedges must be different.",
                           "Select the second vertex.");
        }

        else if(is_next)
        {
          tempInstructions("Vertex not selected : halfedges must not be adjacent.",
                           "Select the second vertex.");
        }
        else if(!is_border or !found)
        {
          tempInstructions("Vertex not selected : no shared border found.",
                           "Select the second vertex.");
        }
        else
        {
          first_selected = false;
          selected_vertices.erase(t->vertex());
          selected_edges.erase(edge(t,*polyhedron()));
          tempInstructions("Face added.",
                           "Select a vertex on a border.");
        }
      }
      break;
    }
    polyhedron_item()->invalidateOpenGLBuffers();
    invalidateOpenGLBuffers();
  }
  is_treated = true;
}

bool Scene_polyhedron_selection_item:: treat_selection(const std::set<edge_descriptor>& selection)
{
  if(!is_treated)
  {
    Selection_traits<edge_descriptor, Scene_polyhedron_selection_item> tr(this);
    switch(operation_mode)
    {
    //classic selection
    case -1:
    {

      bool any_change = false;
      if(is_insert) {
        BOOST_FOREACH(edge_descriptor ed, selection)
            any_change |= tr.container().insert(ed).second;
      }
      else{
        BOOST_FOREACH(edge_descriptor ed, selection)
            any_change |= (tr.container().erase(ed)!=0);
      }
      if(any_change) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
      return any_change;
      break;
    }
      //Join vertex
    case 0:
      BOOST_FOREACH(edge_descriptor ed, selection)
      {
        if(halfedge(ed, *polyhedron())->facet()->facet_degree()<4
          ||
          opposite(halfedge(ed, *polyhedron()), *polyhedron())->facet()->facet_degree()<4)
        {
          tempInstructions("Vertices not joined : the incident facets must have a degree of at least 4.",
                           "Select the edge with extremities you want to join.");
        }
        else
          polyhedron()->join_vertex(halfedge(ed, *polyhedron()));
      }
      break;
      //Split edge
    case 2:
      BOOST_FOREACH(edge_descriptor ed, selection)
      {
        Polyhedron::Point_3 a(halfedge(ed, *poly)->vertex()->point()),b(halfedge(ed, *poly)->opposite()->vertex()->point());
        Polyhedron::Halfedge_handle hhandle = polyhedron()->split_edge(halfedge(ed, *poly));
        Polyhedron::Point_3 p((b.x()+a.x())/2.0, (b.y()+a.y())/2.0,(b.z()+a.z())/2.0);

        hhandle->vertex()->point() = p;

      }
      break;
      //Join face
    case 3:
      BOOST_FOREACH(edge_descriptor ed, selection)
      {
        if(out_degree(source(halfedge(ed,*polyhedron()),*polyhedron()),*polyhedron())<3 ||
           out_degree(target(halfedge(ed,*polyhedron()),*polyhedron()),*polyhedron())<3)
          tempInstructions("Faces not joined : the two ends of the edge must have a degree of at least 3.",
                           "Select the edge separating the faces you want to join.");
        else
          polyhedron()->join_facet(halfedge(ed, *polyhedron()));
      }
      break;
      //Collapse edge
    case 6:
      BOOST_FOREACH(edge_descriptor ed, selection)
      {
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
          CGAL::Euler::collapse_edge(ed, *poly);
          tempInstructions("Edge collapsed.",
                                   "Select the edge you want to collapse.");
        }
      }
      break;
      //Flip edge
    case 7:
      BOOST_FOREACH(edge_descriptor ed, selection)
      {
      //check preconditions
        if(halfedge(ed, *polyhedron())->facet()->facet_degree() == 3 && halfedge(ed, *polyhedron())->opposite()->facet()->facet_degree() == 3)
        {
          CGAL::Euler::flip_edge(halfedge(ed, *polyhedron()), *polyhedron());
        }
        else
        {
          tempInstructions("Edge not selected : incident facets must be triangles.",
                                   "Select the edge you want to flip.");
        }
      }
      break;
    }
    polyhedron_item()->invalidateOpenGLBuffers();
    invalidateOpenGLBuffers();
  }
  is_treated = true;
}

bool Scene_polyhedron_selection_item::treat_selection(const std::set<Polyhedron::Facet_handle>& selection)
{
  if(!is_treated)
  {
    Selection_traits<Facet_handle, Scene_polyhedron_selection_item> tr(this);
    switch(operation_mode)
    {
    //classic selection
    case -1:
    {

      bool any_change = false;
      if(is_insert) {
        BOOST_FOREACH(Facet_handle fh, selection)
            any_change |= tr.container().insert(fh).second;
      }
      else{
        BOOST_FOREACH(Facet_handle fh, selection)
            any_change |= (tr.container().erase(fh)!=0);
      }
      if(any_change) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
      return any_change;
      break;
    }
    case 1:
    {
      static bool first_selected = false;
      static Polyhedron::Halfedge_handle h1;
      //stores first fh and emit change label
      if(!first_selected)
      {
        BOOST_FOREACH(Facet_handle fh, selection)
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
            selected_facets.insert(fh);
            Q_EMIT updateInstructions("Select second facet.");
          }
          else
            tempInstructions("Facet not selected : no valid halfedge",
                             "Select first facet");
        }
      }
      //call the function with point and facets.
      else
      {
        BOOST_FOREACH(Facet_handle fh, selection)
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
            Polyhedron::Halfedge_handle hhandle = polyhedron()->split_vertex(h1,h2);

            selected_facets.erase(h1->facet());
            hhandle->vertex()->point() = to_split_vh->point();
            first_selected = false;
            selected_vertices.erase(to_split_vh);
            //reset selection mode
            set_active_handle_type(static_cast<Active_handle::Type>(0));
            tempInstructions("Point splitted.", "Select the vertex you want splitted");
          }
          else if(h1 == h2)
          {
             tempInstructions("Facet not selected : same as the first.", "Select the second facet.");
          }
          else
          {
            tempInstructions("Facet not selected : no valid halfedge.", "Select the second facet.");
          }
        }
      }
      break;
    }
      //Add center vertex
    case 8:
      BOOST_FOREACH(Facet_handle fh, selection)
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
        Polyhedron::Halfedge_handle hhandle = CGAL::Euler::add_center_vertex(fh->facet_begin(), *poly);
        if(total !=0)
          hhandle->vertex()->point() = Polyhedron::Point_3(x/(double)total, y/(double)total, z/(double)total);


      }
      break;
    }
    polyhedron_item()->invalidateOpenGLBuffers();
    invalidateOpenGLBuffers();
  }
  is_treated = true;
}

void Scene_polyhedron_selection_item::tempInstructions(QString s1, QString s2)
{
  m_temp_instructs = s2;
  Q_EMIT updateInstructions(s1);
  QTimer timer;
  timer.singleShot(3000, this, SLOT(emitTempInstruct()));
}
void Scene_polyhedron_selection_item::emitTempInstruct()
{
  Q_EMIT updateInstructions(m_temp_instructs);
}

