//#define CGAL_PMP_REMESHING_VERBOSE

#include "opengl_tools.h"
#include "Scene_edit_polyhedron_item.h"
#include "Scene_spheres_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <boost/foreach.hpp>
#include <algorithm>
#include <QTime>

#include <QApplication>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

struct Scene_edit_polyhedron_item_priv
{
  Scene_edit_polyhedron_item_priv(Scene_polyhedron_item* poly_item,  Ui::DeformMesh* ui_widget, QMainWindow* mw, Scene_edit_polyhedron_item* parent)
    : ui_widget(ui_widget),
      poly_item(poly_item),
      is_rot_free(true),
      own_poly_item(true),
      k_ring_selector(poly_item, mw, Scene_polyhedron_item_k_ring_selection::Active_handle::VERTEX, true)
  {
    item = parent;
    nb_ROI = 0;
    nb_control = 0;
    nb_axis = 0;
    nb_bbox = 0;
    spheres = NULL;
    spheres_ctrl = NULL;
    need_change = false;

  }
  ~Scene_edit_polyhedron_item_priv()
  {
    delete deform_mesh;
    if (own_poly_item)
      delete poly_item;
  }
  void remesh();
  void expand_or_reduce(int);
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_normals_and_vertices(void);
  void compute_bbox(const CGAL::Three::Scene_interface::Bbox&);
  void reset_drawing_data();
  enum Buffer
  {
      Facet_vertices =0,
      Facet_normals,
      Roi_vertices,
      Control_vertices,
      Control_color,
      Bbox_vertices,
      Axis_vertices,
      Axis_colors,
      Frame_vertices,
      NumberOfBuffers
  };
  enum Vao
  {
      Facets=0,
      Roi_points,
      Edges,
      BBox,
      Control_points,
      Axis,
      Frame_plane,
      NumberOfVaos
  };

  Ui::DeformMesh* ui_widget;
  Scene_polyhedron_item* poly_item;
  bool need_change;
  // For drawing
  mutable std::vector<GLdouble> positions;
  mutable std::vector<unsigned int> tris;
  mutable std::vector<unsigned int> edges;
  mutable std::vector<GLdouble> color_lines;
  mutable std::vector<GLdouble> color_bbox;
  mutable std::vector<GLdouble> ROI_points;
  mutable std::vector<GLdouble> control_points;
  mutable std::vector<GLdouble> ROI_color;
  mutable std::vector<GLdouble> control_color;
  mutable std::vector<GLdouble> normals;
  mutable std::vector<GLdouble> pos_bbox;
  mutable std::vector<GLdouble> pos_axis;
  mutable std::vector<GLdouble> pos_frame_plane;
  mutable QOpenGLShaderProgram *program;
  mutable QOpenGLShaderProgram bbox_program;
  mutable std::size_t nb_ROI;
  mutable std::size_t nb_control;
  mutable std::size_t nb_axis;
  mutable std::size_t nb_bbox;
  mutable Scene_spheres_item* spheres;
  mutable Scene_spheres_item* spheres_ctrl;
  mutable QOpenGLBuffer *in_bu;

  Deform_mesh* deform_mesh;
  typedef std::list<Control_vertices_data> Ctrl_vertices_group_data_list;
  Ctrl_vertices_group_data_list::iterator active_group;
  Ctrl_vertices_group_data_list ctrl_vertex_frame_map; // keep list of group of control vertices with assoc data

  double length_of_axis; // for drawing axis at a group of control vertices

  // by interleaving 'viewer's events (check constructor), keep followings:
  Mouse_keyboard_state_deformation state;

  //For constraint rotation
  qglviewer::LocalConstraint rot_constraint;
  bool is_rot_free;

  bool own_poly_item; //indicates if the poly_item should be deleted by the destructor
  Scene_polyhedron_item_k_ring_selection k_ring_selector;
  Scene_edit_polyhedron_item *item;

}; //end Scene_edit_polyhedron_item_priv

Scene_edit_polyhedron_item::Scene_edit_polyhedron_item
(Scene_polyhedron_item* poly_item,
 Ui::DeformMesh* ui_widget,
 QMainWindow* mw)
  : Scene_group_item("unnamed",Scene_edit_polyhedron_item_priv::NumberOfBuffers,Scene_edit_polyhedron_item_priv::NumberOfVaos)

{
  mw->installEventFilter(this);
  d = new Scene_edit_polyhedron_item_priv(poly_item, ui_widget, mw, this);
  // bind vertex picking
  connect(&d->k_ring_selector, SIGNAL(selected(const std::set<Polyhedron::Vertex_handle>&)), this,
          SLOT(selected(const std::set<Polyhedron::Vertex_handle>&)));

  d->poly_item->set_color_vector_read_only(true); // to prevent recomputation of color vector in invalidateOpenGLBuffers()
  d->poly_item->update_vertex_indices();

  d->deform_mesh = new Deform_mesh(*(poly_item->polyhedron()),
                                Deform_mesh::Vertex_index_map(),
                                Deform_mesh::Hedge_index_map(),
                                Array_based_vertex_point_map(&d->positions));

  d->length_of_axis = (CGAL::sqrt((bbox().xmax()-bbox().xmin())*(bbox().xmax()-bbox().xmin()) + (bbox().ymax()-bbox().ymin())*(bbox().ymax()-bbox().ymin()) + (bbox().zmax()-bbox().zmin())*(bbox().zmax()-bbox().zmin()))) / 15.0;

  // interleave events of viewer (there is only one viewer) 
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->installEventFilter(this);
    
  // create an empty group of control vertices for starting
  create_ctrl_vertices_group();

  // Required for drawing functionality
  d->reset_drawing_data();

    //Generates an integer which will be used as ID for each buffer

    const char vertex_shader_source_bbox[] =
    {
        "#version 120 \n"
        "attribute highp vec3 vertex; \n"
        "attribute highp vec3 colors; \n"

        "uniform highp mat4 mvp_matrix; \n"
        "uniform highp mat4 rotations; \n"
        "uniform highp vec3 translation; \n"
        "uniform highp vec3 translation_2; \n"
        "varying highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = colors; \n"
        "   gl_Position = mvp_matrix * (rotations *(vec4(translation_2,0.0)+vec4(vertex,1.0) )+ vec4(translation,0.0)) ; \n"
        "} \n"
    };
    const char fragment_shader_source[]=
    {
        "#version 120 \n"
        "varying vec3 fColors; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        " gl_FragColor = vec4(fColors, 1.0); \n"
        "} \n"
    };
    d->bbox_program.addShaderFromSourceCode(QOpenGLShader::Vertex,vertex_shader_source_bbox);
    d->bbox_program.addShaderFromSourceCode(QOpenGLShader::Fragment,fragment_shader_source);
    d->bbox_program.link();

    d->ui_widget->remeshing_iterations_spinbox->setValue(1);

    d->ui_widget->remeshing_edge_length_spinbox->setValue(d->length_of_axis);
    d->ui_widget->remeshing_edge_length_spinbox->setDisabled(true);
    d->ui_widget->remeshingEdgeLengthInput_checkBox->setChecked(false);
    connect(ui_widget->remeshingEdgeLengthInput_checkBox, SIGNAL(toggled(bool)),
            ui_widget->remeshing_edge_length_spinbox, SLOT(setEnabled(bool)));

    invalidateOpenGLBuffers();
}

Scene_edit_polyhedron_item::~Scene_edit_polyhedron_item()
{
  while(is_there_any_ctrl_vertices_group())
  {
    delete_ctrl_vertices_group(false);
  }

  delete d;
}
/////////////////////////////
/// For the Shader gestion///
void Scene_edit_polyhedron_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer =0) const
{
    //vao for the facets
    {
        std::vector<GLdouble> vertices;
        std::vector<GLdouble> *vertices_ptr;
        const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
        if(offset.norm() !=0)
        {
          vertices.resize(positions.size());
          for(std::size_t i=0; i<positions.size(); ++i)
          {
            (vertices)[i] = positions[i]+offset[i%3];
          }
          vertices_ptr = &vertices;
        }
        else
        {
          vertices_ptr = &positions;
        }
        program = item->getShaderProgram(Scene_edit_polyhedron_item::PROGRAM_WITH_LIGHT, viewer);
        program->bind();

        item->vaos[Facets]->bind();
        item->buffers[Facet_vertices].bind();
        item->buffers[Facet_vertices].allocate(vertices_ptr->data(),
                            static_cast<int>(vertices_ptr->size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Facet_vertices].release();

        item->buffers[Facet_normals].bind();
        item->buffers[Facet_normals].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(double)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_DOUBLE,0,3);
        item->buffers[Facet_normals].release();
        item->vaos[Facets]->release();
        program->release();
    }
    //vao for the ROI points
    {   program = item->getShaderProgram(Scene_edit_polyhedron_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();
        item->vaos[Roi_points]->bind();
        item->buffers[Roi_vertices].bind();
        item->buffers[Roi_vertices].allocate(ROI_points.data(),
                            static_cast<int>(ROI_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Roi_vertices].release();
        item->vaos[Roi_points]->release();
        program->release();

        nb_ROI = ROI_points.size();
        ROI_points.clear();
        ROI_points.swap(ROI_points);
    }
   //vao for the edges
    {
        program = item->getShaderProgram(Scene_edit_polyhedron_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();
        item->vaos[Edges]->bind();
        item->buffers[Facet_vertices].bind();
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Facet_vertices].release();
        item->vaos[Edges]->release();
        program->release();
    }
    //vao for the BBOX
    {
        bbox_program.bind();
        item->vaos[BBox]->bind();
        item->buffers[Bbox_vertices].bind();
        item->buffers[Bbox_vertices].allocate(pos_bbox.data(),
                             static_cast<int>(pos_bbox.size()*sizeof(double)));
        bbox_program.enableAttributeArray("vertex");
        bbox_program.setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Bbox_vertices].release();

        item->vaos[BBox]->release();
        nb_bbox = pos_bbox.size();
        pos_bbox.resize(0);
        std::vector<double>(pos_bbox).swap(pos_bbox);
        color_bbox.resize(0);
        std::vector<double>(color_bbox).swap(color_bbox);
        bbox_program.release();
    }
    //vao for the control points
    {
        program = item->getShaderProgram(Scene_edit_polyhedron_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();
        item->vaos[Control_points]->bind();
        item->buffers[Control_vertices].bind();
        item->buffers[Control_vertices].allocate(control_points.data(),
                             static_cast<int>(control_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Control_vertices].release();
        item->buffers[Control_color].bind();
        item->buffers[Control_color].allocate(control_color.data(),
                                                 static_cast<int>(control_color.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        item->buffers[Control_color].release();


        item->vaos[Control_points]->release();
        program->release();
        nb_control = control_points.size();
        control_points.clear();
        control_points.swap(control_points);
    }
    //vao for the axis
    {
        program = item->getShaderProgram(Scene_edit_polyhedron_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();
        item->vaos[Axis]->bind();
        item->buffers[Axis_vertices].bind();
        item->buffers[Axis_vertices].allocate(pos_axis.data(),
                             static_cast<int>(pos_axis.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Axis_vertices].release();
        item->buffers[Axis_colors].bind();
        item->buffers[Axis_colors].allocate(color_lines.data(),
                             static_cast<int>(color_lines.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        item->buffers[Axis_colors].release();
        item->vaos[Axis]->release();
        program->release();
        nb_axis = pos_axis.size();
        pos_axis.resize(0);
        std::vector<double>(pos_axis).swap(pos_axis);
        color_lines.resize(0);
        std::vector<double>(color_lines).swap(color_lines);
    }
    //vao for the frame plane
    {
        program = item->getShaderProgram(Scene_edit_polyhedron_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();
        item->vaos[Frame_plane]->bind();
        item->buffers[Frame_vertices].bind();
        item->buffers[Frame_vertices].allocate(pos_frame_plane.data(),
                             static_cast<int>(pos_frame_plane.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Frame_vertices].release();
        program->disableAttributeArray("colors");
        item->vaos[Frame_plane]->release();
        program->release();
    }
    item->are_buffers_filled = true;
}

void Scene_edit_polyhedron_item_priv::reset_drawing_data()
{
  positions.clear();
  positions.resize(num_vertices(*item->polyhedron()) * 3);

  normals.clear();
  normals.resize(positions.size());

  std::size_t counter = 0;
  BOOST_FOREACH(vertex_descriptor vb, vertices(*item->polyhedron()))
  {
    positions[counter * 3] = vb->point().x();
    positions[counter * 3 + 1] = vb->point().y();
    positions[counter * 3 + 2] = vb->point().z();

    const Polyhedron::Traits::Vector_3& n =
      CGAL::Polygon_mesh_processing::compute_vertex_normal(vb, deform_mesh->halfedge_graph());
    normals[counter * 3] = n.x();
    normals[counter * 3 + 1] = n.y();
    normals[counter * 3 + 2] = n.z();

    ++counter;
  }

  tris.clear();
  tris.resize(item->polyhedron()->size_of_facets() * 3);
  counter = 0;
  BOOST_FOREACH(face_descriptor fb, faces(*item->polyhedron()))
  {
    tris[counter * 3] = static_cast<unsigned int>(fb->halfedge()->vertex()->id());
    tris[counter * 3 + 1] = static_cast<unsigned int>(fb->halfedge()->next()->vertex()->id());
    tris[counter * 3 + 2] = static_cast<unsigned int>(fb->halfedge()->prev()->vertex()->id());
    ++counter;
  }

  edges.clear();
  edges.resize(item->polyhedron()->size_of_halfedges());
  counter = 0;
  for (Polyhedron::Edge_iterator eb = item->polyhedron()->edges_begin();
       eb != item->polyhedron()->edges_end(); ++eb, ++counter)
  {
    edges[counter * 2] = static_cast<unsigned int>(eb->vertex()->id());
    edges[counter * 2 + 1] = static_cast<unsigned int>(eb->opposite()->vertex()->id());
  }
}

void Scene_edit_polyhedron_item_priv::compute_normals_and_vertices(void)
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    ROI_points.resize(0);
    control_points.resize(0);
    control_color.resize(0);
    pos_frame_plane.resize(0);
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();

    BOOST_FOREACH(vertex_descriptor vd, deform_mesh->roi_vertices())
    {
        if(!deform_mesh->is_control_vertex(vd))
        {
            ROI_points.push_back(vd->point().x()+offset.x);
            ROI_points.push_back(vd->point().y()+offset.y);
            ROI_points.push_back(vd->point().z()+offset.z);

            if(spheres)
            {
              CGAL::Color c(0,255,0);
              Kernel::Point_3 point(vd->point().x()+offset.x, vd->point().y()+offset.y, vd->point().z()+offset.z);
              spheres->add_sphere(Kernel::Sphere_3(point, length_of_axis/15.0*length_of_axis/15.0), c);
            }
        }

    }
    ROI_color.assign(ROI_points.size(),0);
    for(std::size_t i=0; i<ROI_color.size()/3; i++)
      ROI_color[3*i+1]=1.0;

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
    {
        if(hgb_data->frame == viewer->manipulatedFrame())
        {
            if(!ui_widget->ActivatePivotingCheckBox->isChecked())
            {
                // draw bbox
                compute_bbox(hgb_data->bbox);
            }
        }

        const double r=hgb_data == active_group?1:0;
        const double b=hgb_data == active_group?0:1;

        for(std::vector<vertex_descriptor>::const_iterator hb = hgb_data->ctrl_vertices_group.begin(); hb != hgb_data->ctrl_vertices_group.end(); ++hb)
        {
            control_points.push_back((*hb)->point().x()+offset.x);
            control_points.push_back((*hb)->point().y()+offset.y);
            control_points.push_back((*hb)->point().z()+offset.z);
            control_color.push_back(r);
            control_color.push_back(0);
            control_color.push_back(b);

            if(spheres_ctrl)
            {
              CGAL::Color c(255*r,0,255*b);

              Kernel::Point_3 center((*hb)->point().x()+offset.x,
                                     (*hb)->point().y()+offset.y,
                                     (*hb)->point().z()+offset.z);
              spheres_ctrl->add_sphere(Kernel::Sphere_3(center, length_of_axis/15.0*length_of_axis/15.0), c);
            }
        }
    }

    //The box color
    color_bbox.resize(pos_bbox.size());
    for(int i =0; i< (int)pos_bbox.size(); i++)
        color_bbox[i]=0.0;

    for(int i =0; i< (int)pos_bbox.size(); i+=3)
        color_bbox[i]=1.0;

    //The axis

    pos_axis.resize(18);
    for(int i =0; i< 18; i++)
        pos_axis[i]=0.0;
    pos_axis[3] = length_of_axis; pos_axis[10] = length_of_axis; pos_axis[17] = length_of_axis;
    color_lines.resize(18);
    for(int i =0; i< 18; i++)
        color_lines[i]=0.0;

    color_lines[2] = 1.0; color_lines[5] = 1.0;
    color_lines[6] = 1.0; color_lines[9] = 1.0;
    color_lines[13] = 1.0; color_lines[16] = 1.0;

    if(ui_widget->ActivateFixedPlaneCheckBox->isChecked())
      item->draw_frame_plane(viewer);
    QApplication::restoreOverrideCursor();
}

/////////////////////////////////////////////////////////
/////////// Most relevant functions lie here ///////////
void Scene_edit_polyhedron_item::deform()
{
  if(!is_there_any_ctrl_vertices()) { return; }

  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it)
  { it->set_target_positions(); }
  d->deform_mesh->deform();

  d->poly_item->invalidate_aabb_tree(); // invalidate the AABB-tree of the poly_item
  Q_EMIT itemChanged();
}

struct ROI_faces_pmap
{
  typedef face_descriptor                    key_type;
  typedef std::size_t                        value_type;
  typedef std::size_t&                       reference;
  typedef boost::read_write_property_map_tag category;

  friend value_type get(const ROI_faces_pmap&, const key_type& f)
  { 
    /*magic number 12345*/
    if (f->patch_id() == 12345) return 1;
    else                        return 0;
  }
  friend void put(ROI_faces_pmap&, const key_type& f, const value_type b)
  {
    if(b != 0)  f->set_patch_id(12345);
    else        f->set_patch_id(0);
  }
};

struct ROI_border_pmap
{
  std::set<edge_descriptor>* m_set_ptr;

  typedef edge_descriptor                    key_type;
  typedef bool                               value_type;
  typedef bool                               reference;
  typedef boost::read_write_property_map_tag category;

  ROI_border_pmap() : m_set_ptr(NULL) {}
  ROI_border_pmap(std::set<edge_descriptor>* set_)
    : m_set_ptr(set_)
  {}
  friend bool get(const ROI_border_pmap& map, const key_type& k)
  {
    CGAL_assertion(map.m_set_ptr != NULL);
    return map.m_set_ptr->count(k);
  }
  friend void put(ROI_border_pmap& map, const key_type& k, const value_type b)
  {
    CGAL_assertion(map.m_set_ptr != NULL);
    if (b)              map.m_set_ptr->insert(k);
    else if(get(map,k)) map.m_set_ptr->erase(k);
  }
};

struct halfedge2edge
{
  halfedge2edge(const Polyhedron& m, std::set<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.insert(edge(h, m_mesh));
  }
  const Polyhedron& m_mesh;
  std::set<edge_descriptor>& m_edges;
};

void Scene_edit_polyhedron_item::remesh()
{
  d->remesh();
}

struct Is_constrained_map
{
  boost::unordered_set<vertex_descriptor, CGAL::Handle_hash_function>* m_set_ptr;

  typedef vertex_descriptor                  key_type;
  typedef bool                               value_type;
  typedef bool                               reference;
  typedef boost::read_write_property_map_tag category;

  Is_constrained_map()
    : m_set_ptr(NULL)
  {}
  Is_constrained_map( boost::unordered_set<vertex_descriptor, CGAL::Handle_hash_function>* set_)
    : m_set_ptr(set_)
  {}
  friend bool get(const Is_constrained_map& map, const key_type& k)
  {
    CGAL_assertion(map.m_set_ptr != NULL);
    return map.m_set_ptr->count(k);
  }
  friend void put(Is_constrained_map& map, const key_type& k, const value_type b)
  {
    CGAL_assertion(map.m_set_ptr != NULL);
    if (b)  map.m_set_ptr->insert(k);
    else    map.m_set_ptr->erase(k);
  }
};

void Scene_edit_polyhedron_item_priv::remesh()
{
  if(deform_mesh->roi_vertices().empty())
    return;
  boost::unordered_set<vertex_descriptor, CGAL::Handle_hash_function> constrained_set;
  std::vector<std::vector<vertex_descriptor> > control_groups;
  const Polyhedron& g = deform_mesh->halfedge_graph();
  Array_based_vertex_point_map vpmap(&positions);

  std::set<face_descriptor> roi_facets;
  std::set<vertex_descriptor> roi_vertices(
    deform_mesh->roi_vertices().begin(),deform_mesh->roi_vertices().end());
  QApplication::setOverrideCursor(Qt::WaitCursor);
  for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
  {
    std::vector<vertex_descriptor> group;
    BOOST_FOREACH(vertex_descriptor vd, hgb_data->ctrl_vertices_group)
        group.push_back(vd);
    control_groups.push_back(group);
  }

  ROI_faces_pmap roi_faces_pmap;
  BOOST_FOREACH(vertex_descriptor v, deform_mesh->roi_vertices())
  {
    if(deform_mesh->is_control_vertex(v))
      constrained_set.insert(v);

    BOOST_FOREACH(face_descriptor fv, CGAL::faces_around_target(halfedge(v, g), g))
    {
      if(fv == boost::graph_traits<Polyhedron>::null_face())
        continue;
      bool add_face=true;
      BOOST_FOREACH(vertex_descriptor vfd, CGAL::vertices_around_face(halfedge(fv,g),g))
        if (roi_vertices.count(vfd)==0)
          add_face=false;
      if(add_face)
      {
        roi_facets.insert(fv);
        put(roi_faces_pmap, fv, 1/*true*/);
      }
    }
  }

  if (roi_facets.empty())
  {
    std::cout << "Remeshing canceled (there is no facet with "
              << "its 3 vertices in the ROI)." << std::endl;
    QApplication::restoreOverrideCursor();
    return;
  }
  // set face_index map needed for border_halfedges and isotropic_remeshing
  boost::property_map<Polyhedron, CGAL::face_index_t>::type fim
    = get(CGAL::face_index, *item->polyhedron());
  unsigned int id = 0;

  // estimate the target_length using the perimeter of the region to remesh
  bool automatic_target_length = !ui_widget->remeshingEdgeLengthInput_checkBox->isChecked();
  double estimated_target_length = 0.;
  
    BOOST_FOREACH(face_descriptor f, faces(*item->polyhedron()))
      put(fim, f, id++);

    std::set<edge_descriptor> roi_border;
    CGAL::Polygon_mesh_processing::border_halfedges(roi_facets, g,
      boost::make_function_output_iterator(halfedge2edge(g, roi_border)));

  if (automatic_target_length)
  {
    double sum_len = 0.;
    BOOST_FOREACH(edge_descriptor e, roi_border)
    {
      halfedge_descriptor h = halfedge(e, g);
      sum_len += CGAL::sqrt(CGAL::squared_distance(
                    get(vpmap, source(h, g)), get(vpmap, target(h, g))));
    }
    if (sum_len==0) automatic_target_length = false;
    else
      estimated_target_length = sum_len / (0. + roi_border.size());
  }

  double target_length = automatic_target_length
    ? estimated_target_length
    : ui_widget->remeshing_edge_length_spinbox->value();

  unsigned int nb_iter = ui_widget->remeshing_iterations_spinbox->value();

  std::cout << "Remeshing (target edge length = " << target_length <<")...";

  ROI_border_pmap border_pmap(&roi_border);
  CGAL::Polygon_mesh_processing::isotropic_remeshing(
      roi_facets
    , target_length
    , *item->polyhedron()
    , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
    .protect_constraints(false)
    .vertex_point_map(vpmap)
    .edge_is_constrained_map(border_pmap)
    .face_patch_map(roi_faces_pmap)
    .vertex_is_constrained_map(Is_constrained_map(&constrained_set))
    );
  std::cout << "done." << std::endl;
  poly_item->update_vertex_indices();
  poly_item->update_facet_indices();
  poly_item->update_halfedge_indices();
  //reset ROI from its outside border roi_border
  item->clear_roi();
  do{
    item->delete_ctrl_vertices_group(false);
  }
  while(!ctrl_vertex_frame_map.empty());

  poly_item->update_vertex_indices();
  poly_item->update_halfedge_indices();
  delete deform_mesh;
  deform_mesh = new Deform_mesh(*(poly_item->polyhedron()),
                                Deform_mesh::Vertex_index_map(),
                                Deform_mesh::Hedge_index_map(),
                                vpmap);

  for(std::size_t i=0; i<control_groups.size() ; i++)
  {
    item->create_ctrl_vertices_group();
    BOOST_FOREACH(vertex_descriptor vd, control_groups[i]){
      item->insert_control_vertex(vd);
    }
  }

  BOOST_FOREACH(face_descriptor f, faces(g))
  {
    if (get(roi_faces_pmap, f) == 0/*false*/)
      continue;
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(f, g), g))
    {
      vertex_descriptor v = target(h, g);
      item->insert_roi_vertex(v);
    }
    put(roi_faces_pmap, f, 0/*false*/); //reset ids
  }

  reset_drawing_data();
  compute_normals_and_vertices();

  poly_item->invalidate_aabb_tree(); // invalidate the AABB tree
  QApplication::restoreOverrideCursor();
  Q_EMIT item->itemChanged();
}

void Scene_edit_polyhedron_item::updateDeform()
{
  if(d->need_change)
  {
    // just handle deformation - paint like selection is handled in eventFilter()
    invalidateOpenGLBuffers();
    if(!d->ui_widget->ActivatePivotingCheckBox->isChecked()) {
      deform();
    }
    else {
      Q_EMIT itemChanged(); // for redraw while Pivoting (since we close signals of manipulatedFrames while pivoting,
      // for now redraw with timer)
    }
    d->need_change = 0;
  }
}

void Scene_edit_polyhedron_item::change()
{
  d->need_change = true;
  QTimer::singleShot(0, this, SLOT(updateDeform()));
}
struct Is_selected_property_map{
  std::vector<bool>* is_selected_ptr;
  Is_selected_property_map()
    : is_selected_ptr(NULL) {}
  Is_selected_property_map(std::vector<bool>& is_selected)
    : is_selected_ptr( &is_selected) {}

  std::size_t id(vertex_descriptor v){ return v->id(); }

  friend bool get(Is_selected_property_map map, vertex_descriptor v)
  {
    CGAL_assertion(map.is_selected_ptr!=NULL);
    return (*map.is_selected_ptr)[map.id(v)];
  }

  friend void put(Is_selected_property_map map, vertex_descriptor v, bool b)
  {
    CGAL_assertion(map.is_selected_ptr!=NULL);
    (*map.is_selected_ptr)[map.id(v)]=b;
  }
};

void Scene_edit_polyhedron_item_priv::expand_or_reduce(int steps)
{

  std::vector<bool> mark(poly_item->polyhedron()->size_of_vertices(),false);



  bool ctrl_active = ui_widget->CtrlVertRadioButton->isChecked();
  if (ctrl_active && active_group == ctrl_vertex_frame_map.end() ) return;
  std::size_t original_size;
  if(ctrl_active)
    original_size = active_group->ctrl_vertices_group.size();
  else
    original_size = deform_mesh->roi_vertices().size();
  BOOST_FOREACH(vertex_descriptor v,deform_mesh->roi_vertices())
  {
    if(ctrl_active)
    {
      if(deform_mesh->is_control_vertex(v))
        mark[v->id()]=true;
    }
    else
    {
      if(!deform_mesh->is_control_vertex(v))
        mark[v->id()]=true;
    }
  }
  std::vector<bool> mask = mark;
  if(steps > 0)
  {
    if(ctrl_active)
      expand_vertex_selection(active_group->ctrl_vertices_group, *poly_item->polyhedron(), steps, Is_selected_property_map(mark),
                            CGAL::Emptyset_iterator());
    else
      expand_vertex_selection(deform_mesh->roi_vertices(), *poly_item->polyhedron(), steps, Is_selected_property_map(mark),
                            CGAL::Emptyset_iterator());
  }
  else
  {
    if(ctrl_active)
      reduce_vertex_selection(active_group->ctrl_vertices_group, *poly_item->polyhedron(), -steps, Is_selected_property_map(mark),
                              CGAL::Emptyset_iterator());
    else
      reduce_vertex_selection(deform_mesh->roi_vertices(), *poly_item->polyhedron(), -steps, Is_selected_property_map(mark),
                              CGAL::Emptyset_iterator());
  }

  for(Polyhedron::Vertex_iterator it = poly_item->polyhedron()->vertices_begin() ; it != poly_item->polyhedron()->vertices_end(); ++it)
  {
    if(ctrl_active)
    {
      if(mark[it->id()] && !mask[it->id()])
        item->insert_control_vertex(it);
      else if(!mark[it->id()] && mask[it->id()])
        item->erase_control_vertex(it);
    }
    else
    {
      if(mark[it->id()] && !mask[it->id()])
        item->insert_roi_vertex(it);
      else if(!mark[it->id()] && mask[it->id()])
        item->erase_roi_vertex(it);
    }
  }
  if(active_group->ctrl_vertices_group.empty() && ctrl_vertex_frame_map.size()>1)
    item->delete_ctrl_vertices_group(false);

  if(
     (!ctrl_active && deform_mesh->roi_vertices().size() != original_size)
     || (ctrl_active && active_group->ctrl_vertices_group.size() != original_size)
     )
  { item->invalidateOpenGLBuffers(); Q_EMIT item->itemChanged(); }
}

bool Scene_edit_polyhedron_item::eventFilter(QObject* /*target*/, QEvent *event)
{
  // This filter is both filtering events from 'viewer' and 'main window'
  Mouse_keyboard_state_deformation old_state = d->state;
  ////////////////// TAKE EVENTS /////////////////////
  // key events
  if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease) 
  {
    QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
    Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

    d->state.ctrl_pressing = modifiers.testFlag(Qt::ControlModifier);
    d->state.shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
  }
  // mouse events
  if(event->type() == QEvent::Wheel
     &&d->state.shift_pressing)
  {
    QWheelEvent *w_event = static_cast<QWheelEvent*>(event);
    int steps = w_event->delta() / 120;
    d->expand_or_reduce(steps);
  }
  if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease)
  {
    QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
    if(mouse_event->button() == Qt::LeftButton) {
      d->state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
    }
    if(mouse_event->button() == Qt::RightButton) {
      d->state.right_button_pressing = event->type() == QEvent::MouseButtonPress;
    }    
  }
  ////////////////// //////////////// /////////////////////

  if(!d->poly_item->visible()) { return false; } // if not visible just update event state but don't do any action

  // check state changes between old and current state
  bool ctrl_pressed_now = d->state.ctrl_pressing && !old_state.ctrl_pressing;
  bool ctrl_released_now = !d->state.ctrl_pressing && old_state.ctrl_pressing;
  if(ctrl_pressed_now || ctrl_released_now || event->type() == QEvent::HoverMove)
  {// activate a handle manipulated frame
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    const QPoint& p = viewer->mapFromGlobal(QCursor::pos());
    bool need_repaint = activate_closest_manipulated_frame(p.x(), p.y());

    if (!d->ui_widget->ActivatePivotingCheckBox->isChecked() &&
        ctrl_released_now && d->ui_widget->RemeshingCheckBox->isChecked())
    {
      remesh();
    }
    if(need_repaint)
     invalidateOpenGLBuffers();

    need_repaint |= d->state.left_button_pressing || d->state.right_button_pressing;
    if(need_repaint) { Q_EMIT itemChanged(); }
  }

  return false;
}

#include "opengl_tools.h"
void Scene_edit_polyhedron_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
      d->compute_normals_and_vertices();
      d->initializeBuffers(viewer);
    }
    vaos[Scene_edit_polyhedron_item_priv::Edges]->bind();
    d->program = getShaderProgram(PROGRAM_NO_SELECTION);
    attribBuffers(viewer,PROGRAM_NO_SELECTION);
    d->program->bind();
    d->program->setAttributeValue("colors", QColor(0,0,0));
    viewer->glDrawElements(GL_LINES, (GLsizei) d->edges.size(), GL_UNSIGNED_INT, d->edges.data());
    d->program->release();
    vaos[Scene_edit_polyhedron_item_priv::Edges]->release();


    vaos[Scene_edit_polyhedron_item_priv::Frame_plane]->bind();
    d->program = getShaderProgram(PROGRAM_NO_SELECTION);
    attribBuffers(viewer,PROGRAM_NO_SELECTION);
    d->program->bind();
    d->program->setAttributeValue("colors", QColor(0,0,0));
    viewer->glDrawArrays(GL_LINE_LOOP, 0, (GLsizei)d->pos_frame_plane.size()/3);
    d->program->release();
    vaos[Scene_edit_polyhedron_item_priv::Frame_plane]->release();


  if(rendering_mode == Wireframe) {
        draw_ROI_and_control_vertices(viewer);
  }
}
void Scene_edit_polyhedron_item::draw(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
      d->compute_normals_and_vertices();
      d->initializeBuffers(viewer);
    }
    vaos[Scene_edit_polyhedron_item_priv::Facets]->bind();
    d->program = getShaderProgram(PROGRAM_WITH_LIGHT);
    attribBuffers(viewer,PROGRAM_WITH_LIGHT);
    d->program->bind();
    QColor color = this->color();
    d->program->setAttributeValue("colors", color);
    viewer->glDrawElements(GL_TRIANGLES, (GLsizei) d->tris.size(), GL_UNSIGNED_INT, d->tris.data());
    d->program->release();
    vaos[Scene_edit_polyhedron_item_priv::Facets]->release();
    drawEdges(viewer);
    draw_ROI_and_control_vertices(viewer);

}

void Scene_edit_polyhedron_item::draw_frame_plane(QGLViewer* ) const
{
    d->pos_frame_plane.resize(15);
    for(Scene_edit_polyhedron_item_priv::Ctrl_vertices_group_data_list::const_iterator hgb_data = d->ctrl_vertex_frame_map.begin(); hgb_data != d->ctrl_vertex_frame_map.end(); ++hgb_data)
    {
          const double diag = scene_diag();
          qglviewer::Vec base1(1,0,0);
          qglviewer::Vec base2(0,1,0);

          qglviewer::Quaternion orientation=hgb_data->frame->orientation();
          base1=orientation.rotate(base1);
          base2=orientation.rotate(base2);

          qglviewer::Vec center = hgb_data->calculate_initial_center();
          qglviewer::Vec p1 = center - diag*base1 - diag*base2;
          qglviewer::Vec p2 = center + diag*base1 - diag*base2;
          qglviewer::Vec p3 = center + diag*base1 + diag*base2;
          qglviewer::Vec p4 = center - diag*base1 + diag*base2;

          d->pos_frame_plane[0] = p1.x ; d->pos_frame_plane[1] = p1.y; d->pos_frame_plane[2] =p1.z ;
          d->pos_frame_plane[3] = p2.x ; d->pos_frame_plane[4] = p2.y; d->pos_frame_plane[5] =p2.z ;
          d->pos_frame_plane[6] = p3.x ; d->pos_frame_plane[7] = p3.y; d->pos_frame_plane[8] =p3.z ;
          d->pos_frame_plane[9] = p4.x ; d->pos_frame_plane[10]= p4.y; d->pos_frame_plane[11] =p4.z ;
          d->pos_frame_plane[12] = p1.x ; d->pos_frame_plane[13]= p1.y; d->pos_frame_plane[14] =p1.z ;
    }
}


void Scene_edit_polyhedron_item::draw_ROI_and_control_vertices(CGAL::Three::Viewer_interface* viewer) const {

  CGAL::GL::Point_size point_size; point_size.set_point_size(5);

  //Draw the points
  if(d->ui_widget->ShowROICheckBox->isChecked()) {

        if(!d->ui_widget->ShowAsSphereCheckBox->isChecked() || !viewer->extension_is_found) {

            vaos[Scene_edit_polyhedron_item_priv::Roi_points]->bind();
            d->program = getShaderProgram(PROGRAM_NO_SELECTION);
            attribBuffers(viewer,PROGRAM_NO_SELECTION);
            d->program->bind();
            d->program->setAttributeValue("colors", QColor(0,255,0));
            viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(d->nb_ROI/3));
            d->program->release();
            vaos[Scene_edit_polyhedron_item_priv::Roi_points]->release();
        }
        else
        {
          d->spheres->setVisible(true);
          Scene_group_item::draw(viewer);
        }
  }
  else
  {
    if(d->ui_widget->ShowAsSphereCheckBox->isChecked() && viewer->extension_is_found) {
      d->spheres->setVisible(false);
      Scene_group_item::draw(viewer);
    }
  }

    if(!d->ui_widget->ShowAsSphereCheckBox->isChecked() || !viewer->extension_is_found) {
        vaos[Scene_edit_polyhedron_item_priv::Control_points]->bind();
        d->program = getShaderProgram(PROGRAM_NO_SELECTION);
        attribBuffers(viewer,PROGRAM_NO_SELECTION);
        d->program->bind();
        //d->program->setAttributeValue("colors", QColor(255,0,0));
        viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(d->nb_control/3));
        d->program->release();
        vaos[Scene_edit_polyhedron_item_priv::Control_points]->release();
    }
    // Draw the axis
    QGLViewer* viewerB = *QGLViewer::QGLViewerPool().begin();
    for(Scene_edit_polyhedron_item_priv::Ctrl_vertices_group_data_list::const_iterator hgb_data = d->ctrl_vertex_frame_map.begin(); hgb_data != d->ctrl_vertex_frame_map.end(); ++hgb_data)
    {
         if(hgb_data->frame == viewerB->manipulatedFrame())
         {
              GLfloat f_matrix[16];
              for(int i =0; i<16; i++)
                  f_matrix[i] = hgb_data->frame->matrix()[i];
              QMatrix4x4 f_mat;
                  for(int i=0; i<16; i++)
                      f_mat.data()[i] = (float)f_matrix[i];
              vaos[Scene_edit_polyhedron_item_priv::Axis]->bind();
              d->program = getShaderProgram(PROGRAM_NO_SELECTION);
              attribBuffers(viewer, PROGRAM_NO_SELECTION);
              d->program->bind();
              d->program->setUniformValue("f_matrix", f_mat);
              viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->nb_axis/3));
              d->program->release();
              vaos[Scene_edit_polyhedron_item_priv::Axis]->release();

              //QGLViewer::drawAxis(length_of_axis);
              // draw bbox
              if(!d->ui_widget->ActivatePivotingCheckBox->isChecked())
              {
                   GLfloat f_matrix[16];
                   GLfloat trans[3];
                   GLfloat trans2[3];

                   trans[0] = hgb_data->frame->position().x;
                   trans[1] = hgb_data->frame->position().y;
                   trans[2] = hgb_data->frame->position().z;

                   trans2[0] = -hgb_data->frame_initial_center.x;
                   trans2[1] = -hgb_data->frame_initial_center.y;
                   trans2[2] = -hgb_data->frame_initial_center.z;

                   for(int i =0; i<16; i++)
                       f_matrix[i] = hgb_data->frame->orientation().matrix()[i];
                   QMatrix4x4 f_mat;
                   QMatrix4x4 mvp_mat;

                   QVector3D vec(trans[0], trans[1], trans[2]);
                   QVector3D vec2(trans2[0], trans2[1], trans2[2]);
                   for(int i=0; i<16; i++)
                       f_mat.data()[i] = (float)f_matrix[i];
                   GLdouble temp_mat[16];
                   viewer->camera()->getModelViewProjectionMatrix(temp_mat);
                   for(int i=0; i<16; i++)
                       mvp_mat.data()[i] = (float)temp_mat[i];
                   vaos[Scene_edit_polyhedron_item_priv::BBox]->bind();
                   d->bbox_program.bind();
                   d->bbox_program.setUniformValue("rotations", f_mat);
                   d->bbox_program.setUniformValue("translation", vec);
                   d->bbox_program.setUniformValue("translation_2", vec2);
                   d->bbox_program.setUniformValue("mvp_matrix", mvp_mat);
                   d->program->setAttributeValue("colors", QColor(255,0,0));
                   viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->nb_bbox/3));
                   d->bbox_program.release();
                   vaos[Scene_edit_polyhedron_item_priv::BBox]->release();
              }
         }
    }

}


void Scene_edit_polyhedron_item_priv::compute_bbox(const CGAL::Three::Scene_interface::Bbox& bb){
    pos_bbox.resize(24*3);

    pos_bbox[0]=bb.xmin(); pos_bbox[1]=bb.ymin(); pos_bbox[2]=bb.zmin();
    pos_bbox[3]=bb.xmax(); pos_bbox[4]=bb.ymin(); pos_bbox[5]=bb.zmin();
    pos_bbox[6]=bb.xmin(); pos_bbox[7]=bb.ymin(); pos_bbox[8]=bb.zmin();
    pos_bbox[9]=bb.xmin(); pos_bbox[10]=bb.ymax(); pos_bbox[11]=bb.zmin();
    
    pos_bbox[12]=bb.xmin(); pos_bbox[13]=bb.ymin(); pos_bbox[14]=bb.zmin();
    pos_bbox[15]=bb.xmin(); pos_bbox[16]=bb.ymin(); pos_bbox[17]=bb.zmax();
    pos_bbox[18]= bb.xmax(); pos_bbox[19]=bb.ymin(); pos_bbox[20]=bb.zmin();
    pos_bbox[21]= bb.xmax(); pos_bbox[22]=bb.ymax(); pos_bbox[23]=bb.zmin();
    
    pos_bbox[24]= bb.xmax(); pos_bbox[25]=bb.ymin(); pos_bbox[26]=bb.zmin();
    pos_bbox[27]= bb.xmax(); pos_bbox[28]=bb.ymin(); pos_bbox[29]=bb.zmax();
    pos_bbox[30]=bb.xmin(); pos_bbox[31]=bb.ymax(); pos_bbox[32]=bb.zmin();
    pos_bbox[33]=bb.xmax(); pos_bbox[34]=bb.ymax(); pos_bbox[35]=bb.zmin();
    
    pos_bbox[36]=bb.xmin(); pos_bbox[37]=bb.ymax(); pos_bbox[38]=bb.zmin();
    pos_bbox[39]=bb.xmin(); pos_bbox[40]=bb.ymax(); pos_bbox[41]=bb.zmax();
    pos_bbox[42]=bb.xmin(); pos_bbox[43]=bb.ymin(); pos_bbox[44]=bb.zmax();
    pos_bbox[45]=bb.xmax(); pos_bbox[46]=bb.ymin(); pos_bbox[47]=bb.zmax();

    pos_bbox[48]=bb.xmin(); pos_bbox[49]=bb.ymin(); pos_bbox[50]=bb.zmax();
    pos_bbox[51]=bb.xmin(); pos_bbox[52]=bb.ymax(); pos_bbox[53]=bb.zmax();
    pos_bbox[54]=bb.xmax(); pos_bbox[55]=bb.ymax(); pos_bbox[56]=bb.zmax();
    pos_bbox[57]=bb.xmin(); pos_bbox[58]=bb.ymax(); pos_bbox[59]=bb.zmax();

    pos_bbox[60]=bb.xmax(); pos_bbox[61]=bb.ymax(); pos_bbox[62]=bb.zmax();
    pos_bbox[63]=bb.xmax(); pos_bbox[64]=bb.ymin(); pos_bbox[65]=bb.zmax();
    pos_bbox[66]=bb.xmax(); pos_bbox[67]=bb.ymax(); pos_bbox[68]=bb.zmax();
    pos_bbox[69]=bb.xmax(); pos_bbox[70]=bb.ymax(); pos_bbox[71]=bb.zmin();
    
}

void Scene_edit_polyhedron_item::invalidateOpenGLBuffers()
{
    if(d->spheres)
      d->spheres->clear_spheres();
    if(d->spheres_ctrl)
      d->spheres_ctrl->clear_spheres();
    update_normals();
    compute_bbox();
    are_buffers_filled = false;
    if(d->spheres)
      d->spheres->invalidateOpenGLBuffers();
    if(d->spheres_ctrl)
      d->spheres_ctrl->invalidateOpenGLBuffers();
}

Scene_polyhedron_item* Scene_edit_polyhedron_item::to_polyhedron_item() {
  Scene_polyhedron_item* poly_item_tmp = d->poly_item;
  d->poly_item->set_color_vector_read_only(false);
  d->own_poly_item=false;
  poly_item_tmp->invalidateOpenGLBuffers();
  return poly_item_tmp;
}

Polyhedron* Scene_edit_polyhedron_item::polyhedron()       
{ return d->poly_item->polyhedron(); }
const Polyhedron* Scene_edit_polyhedron_item::polyhedron() const 
{ return d->poly_item->polyhedron(); }
QString Scene_edit_polyhedron_item::toolTip() const
{
  if(!d->poly_item->polyhedron())
    return QString();

  return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
    .arg(this->name())
    .arg(d->poly_item->polyhedron()->size_of_vertices())
    .arg(d->poly_item->polyhedron()->size_of_halfedges()/2)
    .arg(d->poly_item->polyhedron()->size_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}
bool Scene_edit_polyhedron_item::isEmpty() const {
  return d->poly_item->isEmpty();
}
void Scene_edit_polyhedron_item::compute_bbox() const {
  _bbox = d->poly_item->bbox();
}

void Scene_edit_polyhedron_item::setVisible(bool b) {
  d->poly_item->setVisible(b);
  Scene_item::setVisible(b);
  if(!b) {
    (*QGLViewer::QGLViewerPool().begin())->setManipulatedFrame(NULL);
  }
}
void Scene_edit_polyhedron_item::setColor(QColor c) {
  d->poly_item->setColor(c);
  Scene_item::setColor(c);
}
void Scene_edit_polyhedron_item::setName(QString n) {
  Scene_item::setName(n);
  n.replace(" (edit)", "");
  d->poly_item->setName(n);
}
void Scene_edit_polyhedron_item::setRenderingMode(RenderingMode m) {
  d->poly_item->setRenderingMode(m);
  Scene_item::setRenderingMode(m);
}
Scene_edit_polyhedron_item* Scene_edit_polyhedron_item::clone() const {
  return 0;
}
void Scene_edit_polyhedron_item::select(
          double orig_x,
          double orig_y,
          double orig_z,
          double dir_x,
          double dir_y,
          double dir_z)
{
  Scene_item::select(orig_x,
                     orig_y,
                     orig_z,
                     dir_x,
                     dir_y,
                     dir_z);
  d->poly_item->select(orig_x,
                       orig_y,
                       orig_z,
                       dir_x,
                       dir_y,
                       dir_z);
}

bool Scene_edit_polyhedron_item::keyPressEvent(QKeyEvent* e)
{
  //setting/unsetting rotation constraints
  if (e->key()==Qt::Key_R && !d->state.ctrl_pressing)
  {
    d->is_rot_free = !d->is_rot_free;
    d->rot_constraint.setRotationConstraintType( d->is_rot_free?
        qglviewer::AxisPlaneConstraint::FREE:
        qglviewer::AxisPlaneConstraint::AXIS);
    return true;
  }

  return false;
}


//#include "Scene_edit_polyhedron_item.moc"
void Scene_edit_polyhedron_item::update_frame_plane()
{
  for(Ctrl_vertices_group_data_list::iterator hgb_data = d->ctrl_vertex_frame_map.begin(); hgb_data != d->ctrl_vertex_frame_map.end(); ++hgb_data)
  {
      hgb_data->refresh();
  }
}


struct Reset_spheres_ctrl
{
  Scene_edit_polyhedron_item_priv* d;
  Reset_spheres_ctrl(Scene_edit_polyhedron_item_priv* d) : d(d) {}
  void operator()() const { d->spheres_ctrl = NULL; }
};

void Scene_edit_polyhedron_item::ShowAsSphere(bool b)
{
  if(b)
  {
    if(!d->spheres)
    {
      d->spheres = new Scene_spheres_item(this, false);
      d->spheres->setName("ROI spheres");
      d->spheres->setRenderingMode(Gouraud);
      connect(d->spheres, SIGNAL(destroyed()), this, SLOT(reset_spheres()));
      scene->setSelectedItem(scene->item_id(this));
      scene->addItem(d->spheres);
      scene->changeGroup(d->spheres, this);
      lockChild(d->spheres);
      invalidateOpenGLBuffers();
    }
    if(!d->spheres_ctrl)
    {
      d->spheres_ctrl = new Scene_spheres_item(this, false);
      d->spheres_ctrl->setName("Control spheres");
      d->spheres_ctrl->setRenderingMode(Gouraud);
      connect(d->spheres_ctrl, &QObject::destroyed, this, Reset_spheres_ctrl(d) );
      scene->setSelectedItem(scene->item_id(this));
      scene->addItem(d->spheres_ctrl);
      scene->changeGroup(d->spheres_ctrl, this);
      lockChild(d->spheres_ctrl);
      invalidateOpenGLBuffers();
    }
  }
  else if(!b )
  {
    if(d->spheres!=0)
    {
      unlockChild(d->spheres);
      removeChild(d->spheres);
      scene->erase(scene->item_id(d->spheres));
    }
    if(d->spheres_ctrl!=0)
    {
      unlockChild(d->spheres_ctrl);
      removeChild(d->spheres_ctrl);
      scene->erase(scene->item_id(d->spheres_ctrl));
    }
  }
}
bool Scene_edit_polyhedron_item::is_there_any_ctrl_vertices_group(Ctrl_vertices_group_data_list::iterator& hgb, Ctrl_vertices_group_data_list::iterator& hge)
{
  hgb = d->ctrl_vertex_frame_map.begin(); hge = d->ctrl_vertex_frame_map.end();
  return hgb != hge;
}
int Scene_edit_polyhedron_item::get_k_ring()       { return d->k_ring_selector.k_ring; }
void Scene_edit_polyhedron_item::set_k_ring(int v) { d->k_ring_selector.k_ring = v; }

void Scene_edit_polyhedron_item::reset_spheres()
{
  d->spheres = NULL;
}

void Scene_edit_polyhedron_item::selected(const std::set<Polyhedron::Vertex_handle>& m)
{
  bool any_changes = false;
  for(std::set<vertex_descriptor>::const_iterator it = m.begin(); it != m.end(); ++it)
  {
    vertex_descriptor vh = *it;
    bool changed = false;
    if(d->ui_widget->ROIRadioButton->isChecked()) {
      if(d->ui_widget->InsertRadioButton->isChecked()) { changed = insert_roi_vertex(vh);}
      else          { changed = erase_roi_vertex(vh);  }
    }
    else {
      if(d->ui_widget->InsertRadioButton->isChecked()) { changed = insert_control_vertex(vh); }
      else          { changed = erase_control_vertex(vh);  }
    }
    any_changes |= changed;
  }
  if(any_changes) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
}

bool Scene_edit_polyhedron_item::insert_control_vertex(vertex_descriptor v)
{
  if(!is_there_any_ctrl_vertices_group()) {
    std::cerr<<"There is no group of control vertices, create one!\n";
    return false;
  } // no group of control vertices to insert

  bool inserted = d->deform_mesh->insert_control_vertex(v);
  if(inserted) {
    d->active_group->ctrl_vertices_group.push_back(v);
    d->active_group->refresh();
  }
  return inserted;
}

bool Scene_edit_polyhedron_item::insert_roi_vertex(vertex_descriptor v)
{
  return d->deform_mesh->insert_roi_vertex(v);
}

bool Scene_edit_polyhedron_item::erase_control_vertex(vertex_descriptor v)
{
  if(d->deform_mesh->erase_control_vertex(v)) // API should be safe enough to do that (without checking empty group of control vertices etc.)
  {
    refresh_all_group_centers(); // since we don't know which group of control vertices v is erased from, refresh all
    return true;
  }

  print_message("Selected vertex is not a control vertex!");
  return false;
}

bool Scene_edit_polyhedron_item::erase_roi_vertex(vertex_descriptor v)
{
  erase_control_vertex(v); // erase control vertex
  return d->deform_mesh->erase_roi_vertex(v);
}

void Scene_edit_polyhedron_item::clear_roi()
{
  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it)
  {
    delete it->frame;
  }
  d->ctrl_vertex_frame_map.clear();
  d->deform_mesh->clear_roi_vertices();

  create_ctrl_vertices_group(); // create one new group of control vertices
}

void Scene_edit_polyhedron_item::create_ctrl_vertices_group()
{
  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it) {
    if(it->ctrl_vertices_group.empty()) {
      d->active_group = it;
      return;
    }
  }

  // No empty group of control vertices
  const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
  qglviewer::ManipulatedFrame* new_frame = new qglviewer::ManipulatedFrame();
  new_frame->setPosition(offset);
  new_frame->setRotationSensitivity(2.0f);
  connect(new_frame, SIGNAL(manipulated()), this, SLOT(change()));

  Control_vertices_data hgd(d->deform_mesh, new_frame);
  d->ctrl_vertex_frame_map.push_back(hgd);
  hgd.refresh();

  d->active_group = --d->ctrl_vertex_frame_map.end();

  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();

  print_message("A new empty group of control vertices is created.");
}

void Scene_edit_polyhedron_item::delete_ctrl_vertices_group(bool create_new)
{
  if(!is_there_any_ctrl_vertices_group()) {
    print_message("There is no group of control vertices to be deleted!");
    return;
  } // no group of control vertices

  // delete group representative
  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it)
  {
    if(it == d->active_group)
    {
      delete it->frame;
      for(std::vector<vertex_descriptor>::iterator v_it = it->ctrl_vertices_group.begin(); v_it != it->ctrl_vertices_group.end(); ++v_it) {
        d->deform_mesh->erase_control_vertex(*v_it);
      }
      d->ctrl_vertex_frame_map.erase(it);
      break;
    }
  }

  // assign another ctrl_vertices_group to active_group
  Ctrl_vertices_group_data_list::iterator hgb, hge;
  if( is_there_any_ctrl_vertices_group(hgb, hge) )
  {
    d->active_group = hgb;
  } // no group of control vertices
  else if(create_new)
  {
    create_ctrl_vertices_group();
  }
}

void Scene_edit_polyhedron_item::prev_ctrl_vertices_group()
{
  Ctrl_vertices_group_data_list::iterator hgb, hge;
  if( !is_there_any_ctrl_vertices_group(hgb, hge) ) {
    print_message("There is no group of control vertices to iterate on!");
    return;
  }
  // shift
  if(hgb == d->active_group) { d->active_group = --hge; }
  else                    {--d->active_group; }
}

void Scene_edit_polyhedron_item::next_ctrl_vertices_group()
{
  Ctrl_vertices_group_data_list::iterator hgb, hge;
  if( !is_there_any_ctrl_vertices_group(hgb, hge) ) {
    print_message("There is no group of control vertices to iterate on!");
    return;
  }
  // shift
  if(--hge == d->active_group) { d->active_group = hgb; }
  else                      {++d->active_group; }
}

void Scene_edit_polyhedron_item::pivoting_end()
{
  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it)
  {
    //update constraint rotation vector, set only for the last group
    it->rot_direction = it->frame->rotation().rotate( qglviewer::Vec(0.,0.,1.) );
    //translate center of the frame
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    qglviewer::Vec vec= it->frame->position();
    it->refresh();
    it->frame_initial_center = vec-offset;
    it->frame->setPosition(vec);
  }
  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it)
  {
    it->frame->blockSignals(false);
  }
}

void Scene_edit_polyhedron_item::pivoting_begin()
{
  d->is_rot_free=true;
  d->rot_constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
  d->rot_constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);

  // just block signals to prevent deformation
  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it)
  {
    it->frame->blockSignals(true);
  }
}

void Scene_edit_polyhedron_item::save_roi(const char* file_name) const
{
  std::ofstream out(file_name);
  // save roi
  out << d->deform_mesh->roi_vertices().size() << std::endl;
  BOOST_FOREACH(vertex_descriptor vd, d->deform_mesh->roi_vertices())
  {
    out << vd->id() << " ";
  }
  out << std::endl;
  // save control vertices

  out << d->ctrl_vertex_frame_map.size() << std::endl; // control vertices count
  for(Ctrl_vertices_group_data_list::const_iterator hgb = d->ctrl_vertex_frame_map.begin(); hgb != d->ctrl_vertex_frame_map.end(); ++hgb) {

    out << hgb->ctrl_vertices_group.size() << std::endl;
    for(std::vector<vertex_descriptor>::const_iterator hb = hgb->ctrl_vertices_group.begin(); hb != hgb->ctrl_vertices_group.end(); ++hb)
    {
      out << (*hb)->id() << " ";
    }
    out << std::endl;
  }
}

void Scene_edit_polyhedron_item::read_roi(const char* file_name)
{
  clear_roi();
  delete_ctrl_vertices_group(false);

  // put vertices to vector
  std::vector<vertex_descriptor> all_vertices;
  all_vertices.reserve(num_vertices(d->deform_mesh->halfedge_graph()));
  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(d->deform_mesh->halfedge_graph()); vb != ve; ++vb) {
    all_vertices.push_back(*vb);
  }
  // read roi
  std::ifstream in(file_name);
  int roi_size;
  in >> roi_size;
  while(roi_size-- > 0)
  {
    std::size_t v_id;
    in >> v_id;
    insert_roi_vertex(all_vertices[v_id]);
  }
  // read control vertices
  int ctrl_vertices_group_size;
  in >> ctrl_vertices_group_size;
  while(ctrl_vertices_group_size-- > 0)
  {
    create_ctrl_vertices_group();
    int ctrl_size;
    in >> ctrl_size;
    while(ctrl_size-- > 0)
    {
      std::size_t v_id;
      in >> v_id;
      insert_control_vertex(all_vertices[v_id]);
    }
  }
}

void Scene_edit_polyhedron_item::overwrite_deform_object()
{
  d->deform_mesh->overwrite_initial_geometry();

  refresh_all_group_centers();
}

void Scene_edit_polyhedron_item::reset_deform_object()
{
  d->deform_mesh->reset();
  refresh_all_group_centers();
}


boost::optional<std::size_t> Scene_edit_polyhedron_item::get_minimum_isolated_component() {
  Travel_isolated_components::Minimum_visitor visitor;
  Travel_isolated_components().travel<Vertex_handle>
    (vertices(*polyhedron()).first, vertices(*polyhedron()).second,
     polyhedron()->size_of_vertices(), Is_selected(d->deform_mesh), visitor);
  return visitor.minimum;
}



boost::optional<std::size_t> Scene_edit_polyhedron_item::select_isolated_components(std::size_t threshold) {
  typedef boost::function_output_iterator<Select_roi_output> Output_iterator;
  Output_iterator out(d->deform_mesh);

  Travel_isolated_components::Selection_visitor<Output_iterator> visitor(threshold, out);
  Travel_isolated_components().travel<Vertex_handle>
    (vertices(*polyhedron()).first, vertices(*polyhedron()).second,
    polyhedron()->size_of_vertices(), Is_selected(d->deform_mesh), visitor);

  if(visitor.any_inserted) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
  return visitor.minimum_visitor.minimum;
}

bool Scene_edit_polyhedron_item::is_there_any_ctrl_vertices_group()
{
  Ctrl_vertices_group_data_list::iterator hgb, hge;
  return is_there_any_ctrl_vertices_group(hgb, hge);
}

bool Scene_edit_polyhedron_item::is_there_any_ctrl_vertices()
{
  Ctrl_vertices_group_data_list::iterator hgb, hge;
  if(!is_there_any_ctrl_vertices_group(hgb, hge)) { return false; } // there isn't any group of control vertices

  for(; hgb != hge; ++hgb) // check inside groups of control vertices
  {
    if(!hgb->ctrl_vertices_group.empty()) { return true; }
  }
  return false;
}

void Scene_edit_polyhedron_item::refresh_all_group_centers()
{
  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it)
  { it->refresh(); }
}

bool Scene_edit_polyhedron_item::activate_closest_manipulated_frame(int x, int y)
{
  if(d->state.ctrl_pressing && (d->state.left_button_pressing || d->state.right_button_pressing) )
  { // user is deforming currently don't change the state
    return false;
  }
  if(d->ctrl_vertex_frame_map.empty()) { return false; }

  d->rot_constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
  d->rot_constraint.setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);

  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  qglviewer::Camera* camera = viewer->camera();

  if(!d->state.ctrl_pressing)
  {
    if(viewer->manipulatedFrame() == NULL)
    { return false;}
    viewer->setManipulatedFrame(NULL);
    return true;
  }

  // now find closest frame and make it active manipulated frame
  Ctrl_vertices_group_data_list::iterator min_it = d->ctrl_vertex_frame_map.begin();
  const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(min_it->frame->position());
  float min_dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);

  for(Ctrl_vertices_group_data_list::iterator it = d->ctrl_vertex_frame_map.begin(); it != d->ctrl_vertex_frame_map.end(); ++it)
  {
    const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(it->frame->position());
    float dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);
    if(dist < min_dist) {
      min_dist = dist;
      min_it = it;
    }
  }

  //set rotation constraint for the manipulated frame
  if (!d->is_rot_free){
    d->rot_constraint.setRotationConstraintDirection(min_it->rot_direction);
    d->rot_constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
    min_it->frame->setConstraint(&d->rot_constraint);
  }
  else
  {
    if( d->ui_widget->ActivateFixedPlaneCheckBox->isChecked())
    {
      // the constraint is local to the frame
      d->rot_constraint.setTranslationConstraint(qglviewer::AxisPlaneConstraint::PLANE,qglviewer::Vec(0,0,1));
       if(!d->ui_widget->ActivatePivotingCheckBox->isChecked()){
           d->rot_constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
       }
       min_it->frame->setConstraint(&d->rot_constraint);
    }
  }

  if(viewer->manipulatedFrame() == min_it->frame)
  { return false; }
  viewer->setManipulatedFrame(min_it->frame);
  d->active_group = min_it;
  return true;
}

void Scene_edit_polyhedron_item::update_normals() {
  BOOST_FOREACH(vertex_descriptor vd, d->deform_mesh->roi_vertices())
  {
    std::size_t id = vd->id();
    const Polyhedron::Traits::Vector_3& n =
      CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, d->deform_mesh->halfedge_graph());
    d->normals[id*3] = n.x();
    d->normals[id*3+1] = n.y();
    d->normals[id*3+2] = n.z();

  }
}

