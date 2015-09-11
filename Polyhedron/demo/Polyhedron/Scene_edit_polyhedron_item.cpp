#include "opengl_tools.h"
#include "create_sphere.h"
#include "Scene_edit_polyhedron_item.h"
#include <boost/foreach.hpp>
#include <algorithm>
#include <QTime>
#include <CGAL/gl_render.h>
Scene_edit_polyhedron_item::Scene_edit_polyhedron_item
(Scene_polyhedron_item* poly_item,
 Ui::DeformMesh* ui_widget,
 QMainWindow* mw)
  : Scene_item(19,8),
    ui_widget(ui_widget),
    poly_item(poly_item),
    deform_mesh(*(poly_item->polyhedron()), Deform_mesh::Vertex_index_map(), Deform_mesh::Hedge_index_map(), Array_based_vertex_point_map(&positions)),
    is_rot_free(true),
    own_poly_item(true),
    k_ring_selector(poly_item, mw, Scene_polyhedron_item_k_ring_selection::Active_handle::VERTEX, true)
{
  nb_ROI = 0;
  nb_sphere = 0;
  nb_control = 0;
  nb_axis = 0;
  nb_bbox = 0;
  mw->installEventFilter(this);
  // bind vertex picking
  connect(&k_ring_selector, SIGNAL(selected(const std::set<Polyhedron::Vertex_handle>&)), this,
          SLOT(selected(const std::set<Polyhedron::Vertex_handle>&)));

  poly_item->set_color_vector_read_only(true); // to prevent recomputation of color vector in invalidate_buffers()
  poly_item->update_vertex_indices();

  length_of_axis = bbox().diagonal_length() / 15.0;

  // interleave events of viewer (there is only one viewer) 
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->installEventFilter(this);
    
  // create an empty group of control vertices for starting
  create_ctrl_vertices_group();
   
  // start QObject's timer for continuous effects 
  // (deforming mesh while mouse not moving)
  startTimer(0);

  // Required for drawing functionality
  positions.resize(num_vertices(*polyhedron())*3);
  normals.resize(positions.size());
  Polyhedron::Vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(vb=polyhedron()->vertices_begin(), ve = polyhedron()->vertices_end();vb != ve; ++vb, ++counter) {
    positions[counter*3] = vb->point().x();
    positions[counter*3+1] = vb->point().y();
    positions[counter*3+2] = vb->point().z();

    const Polyhedron::Traits::Vector_3& n = 
      CGAL::Polygon_mesh_processing::compute_vertex_normal(vb, deform_mesh.halfedge_graph());

    normals[counter*3] = n.x();
    normals[counter*3+1] = n.y();
    normals[counter*3+2] = n.z();
  }
  tris.resize(polyhedron()->size_of_facets()*3);
  counter = 0;
  for(Polyhedron::Facet_handle fb = polyhedron()->facets_begin(); fb != polyhedron()->facets_end(); ++fb, ++counter) {
    tris[counter*3] =  static_cast<unsigned int>(fb->halfedge()->vertex()->id());
    tris[counter*3+1] = static_cast<unsigned int>(fb->halfedge()->next()->vertex()->id());
    tris[counter*3+2] = static_cast<unsigned int>(fb->halfedge()->prev()->vertex()->id());
  }

  edges.resize(polyhedron()->size_of_halfedges());
  counter = 0;
  for(Polyhedron::Edge_iterator eb = polyhedron()->edges_begin(); eb != polyhedron()->edges_end(); ++eb, ++counter) {
    edges[counter*2] = static_cast<unsigned int>(eb->vertex()->id());
    edges[counter*2+1] = static_cast<unsigned int>(eb->opposite()->vertex()->id());
  }
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
    bbox_program.addShaderFromSourceCode(QOpenGLShader::Vertex,vertex_shader_source_bbox);
    bbox_program.addShaderFromSourceCode(QOpenGLShader::Fragment,fragment_shader_source);
    bbox_program.link();

    //the spheres :
    create_Sphere(length_of_axis/15.0);
    invalidate_buffers();
}

Scene_edit_polyhedron_item::~Scene_edit_polyhedron_item()
{
    setVisible(false);
  while(is_there_any_ctrl_vertices_group())
  {
    delete_ctrl_vertices_group(false);
  }
  if (own_poly_item) delete poly_item;

}
/////////////////////////////
/// For the Shader gestion///
void Scene_edit_polyhedron_item::initialize_buffers(Viewer_interface *viewer =0) const
{
    //vao for the facets
    {
        program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
        program->bind();

        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions.data(),
                            static_cast<int>(positions.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[0].release();

        buffers[1].bind();
        buffers[1].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(double)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_DOUBLE,0,3);
        buffers[1].release();
        vaos[0]->release();
        program->release();
    }
    //vao for the ROI points
    {   program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[1]->bind();
        buffers[2].bind();
        buffers[2].allocate(ROI_points.data(),
                            static_cast<int>(ROI_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[2].release();

        buffers[3].bind();
        buffers[3].allocate(ROI_color.data(),
                            static_cast<int>(ROI_color.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[3].release();
        vaos[1]->release();

        program->release();
    }

   //vao for the edges
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[2]->bind();
        buffers[4].bind();
        buffers[4].allocate(positions.data(),
                            static_cast<int>(positions.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[4].release();

        vaos[2]->release();
        program->release();
    }
    //vao for the ROI spheres
    {
        program = getShaderProgram(PROGRAM_INSTANCED, viewer);
        program->bind();
        vaos[3]->bind();
        buffers[5].bind();
        buffers[5].allocate(pos_sphere.data(),
                            static_cast<int>(pos_sphere.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[5].release();

        buffers[6].bind();
        buffers[6].allocate(normals_sphere.data(),
                            static_cast<int>(normals_sphere.size()*sizeof(double)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_DOUBLE,0,3);
        buffers[6].release();

        buffers[7].bind();
        buffers[7].allocate(ROI_color.data(),
                            static_cast<int>(ROI_color.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[7].release();

        buffers[8].bind();
        buffers[8].allocate(ROI_points.data(),
                            static_cast<int>(ROI_points.size()*sizeof(double)));
        program->enableAttributeArray("center");
        program->setAttributeBuffer("center",GL_DOUBLE,0,3);
        buffers[8].release();

        if(viewer->extension_is_found)
        {
            viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
            viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);
        }
        vaos[3]->release();
        ROI_color.resize(0);
        std::vector<double>(ROI_color).swap(ROI_color);
        nb_ROI = ROI_points.size();
        ROI_points.resize(0);
        std::vector<double>(ROI_points).swap(ROI_points);
    }
    //vao for the BBOX
    {
        bbox_program.bind();
        vaos[4]->bind();
        buffers[9].bind();
        buffers[9].allocate(pos_bbox.data(),
                             static_cast<int>(pos_bbox.size()*sizeof(double)));
        bbox_program.enableAttributeArray("vertex");
        bbox_program.setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[9].release();

        buffers[10].bind();
        buffers[10].allocate(color_bbox.data(),
                             static_cast<int>(color_bbox.size()*sizeof(double)));
        bbox_program.enableAttributeArray("colors");
        bbox_program.setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[10].release();
        vaos[4]->release();
        nb_bbox = pos_bbox.size();
        pos_bbox.resize(0);
        std::vector<double>(pos_bbox).swap(pos_bbox);
        color_bbox.resize(0);
        std::vector<double>(color_bbox).swap(color_bbox);
        bbox_program.release();
    }
    //vao for the control points
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[5]->bind();
        buffers[11].bind();
        buffers[11].allocate(control_points.data(),
                             static_cast<int>(control_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[11].release();

        buffers[12].bind();
        buffers[12].allocate(control_color.data(),
                             static_cast<int>(control_color.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[12].release();
        vaos[5]->release();
        program->release();
    }
    //vao for the control spheres
    {
        program = getShaderProgram(PROGRAM_INSTANCED, viewer);
        program->bind();
        vaos[6]->bind();
        buffers[13].bind();
        buffers[13].allocate(pos_sphere.data(),
                             static_cast<int>(pos_sphere.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[13].release();

        buffers[14].bind();
        buffers[14].allocate(normals_sphere.data(),
                             static_cast<int>(normals_sphere.size()*sizeof(double)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_DOUBLE,0,3);
        buffers[14].release();

        buffers[15].bind();
        buffers[15].allocate(control_color.data(),
                             static_cast<int>(control_color.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[15].release();

        buffers[16].bind();
        buffers[16].allocate(control_points.data(),
                             static_cast<int>(control_points.size()*sizeof(double)));
        program->enableAttributeArray("center");
        program->setAttributeBuffer("center",GL_DOUBLE,0,3);
        buffers[16].release();

        if(viewer->extension_is_found)
        {
            viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
            viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);
        }
        vaos[6]->release();
        nb_sphere = pos_sphere.size();
        pos_sphere.resize(0);
        std::vector<double>(pos_sphere).swap(pos_sphere);
        normals_sphere.resize(0);
        std::vector<double>(normals_sphere).swap(normals_sphere);
        control_color.resize(0);
        std::vector<double>(control_color).swap(control_color);
        nb_control = control_points.size();
        control_points.resize(0);
        std::vector<double>(control_points).swap(control_points);
    }
    //vao for the axis
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[7]->bind();
        buffers[17].bind();
        buffers[17].allocate(pos_axis.data(),
                             static_cast<int>(pos_axis.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[17].release();

        buffers[18].bind();
        buffers[18].allocate(color_lines.data(),
                             static_cast<int>(color_lines.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[18].release();
        vaos[7]->release();
        program->release();
        nb_axis = pos_axis.size();
        pos_axis.resize(0);
        std::vector<double>(pos_axis).swap(pos_axis);
        color_lines.resize(0);
        std::vector<double>(color_lines).swap(color_lines);
    }
    are_buffers_filled = true;
}

void Scene_edit_polyhedron_item::compute_normals_and_vertices(void)
{
    ROI_points.resize(0);
    control_points.resize(0);
    control_color.resize(0);
    BOOST_FOREACH(vertex_descriptor vd, deform_mesh.roi_vertices())
    {
        if(!deform_mesh.is_control_vertex(vd))
        {//gl_draw_point( vd->point() );
            ROI_points.push_back(vd->point().x());
            ROI_points.push_back(vd->point().y());
            ROI_points.push_back(vd->point().z());
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
            control_points.push_back((*hb)->point().x());
            control_points.push_back((*hb)->point().y());
            control_points.push_back((*hb)->point().z());
            control_color.push_back(r);
            control_color.push_back(0);
            control_color.push_back(b);
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

}

/////////////////////////////////////////////////////////
/////////// Most relevant functions lie here ///////////
void Scene_edit_polyhedron_item::deform()
{
  if(!is_there_any_ctrl_vertices()) { return; }

  for(Ctrl_vertices_group_data_list::iterator it = ctrl_vertex_frame_map.begin(); it != ctrl_vertex_frame_map.end(); ++it)
  { it->set_target_positions(); }
  deform_mesh.deform();

  poly_item->invalidate_aabb_tree(); // invalidate the AABB-tree of the poly_item
  Q_EMIT itemChanged();
}

void Scene_edit_polyhedron_item::timerEvent(QTimerEvent* /*event*/)
{ // just handle deformation - paint like selection is handled in eventFilter()
  if(state.ctrl_pressing && (state.left_button_pressing || state.right_button_pressing)) {
    if(!ui_widget->ActivatePivotingCheckBox->isChecked()) {
      invalidate_buffers();
      deform();
    }
    else {
      Q_EMIT itemChanged(); // for redraw while Pivoting (since we close signals of manipulatedFrames while pivoting, 
                          // for now redraw with timer)
    }
  }
}
bool Scene_edit_polyhedron_item::eventFilter(QObject* /*target*/, QEvent *event)
{
  // This filter is both filtering events from 'viewer' and 'main window'
  Mouse_keyboard_state_deformation old_state = state;
  ////////////////// TAKE EVENTS /////////////////////
  // key events
  if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease) 
  {
    QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
    Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

    state.ctrl_pressing = modifiers.testFlag(Qt::ControlModifier);
    state.shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
  }
  // mouse events
  if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease)
	{
    QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
    if(mouse_event->button() == Qt::LeftButton) {
      state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
    }
    if(mouse_event->button() == Qt::RightButton) {
      state.right_button_pressing = event->type() == QEvent::MouseButtonPress;
    }    
  }
  ////////////////// //////////////// /////////////////////

  if(!poly_item->visible()) { return false; } // if not visible just update event state but don't do any action

  // check state changes between old and current state
  bool ctrl_pressed_now = state.ctrl_pressing && !old_state.ctrl_pressing;
  bool ctrl_released_now = !state.ctrl_pressing && old_state.ctrl_pressing;
  if(ctrl_pressed_now || ctrl_released_now || event->type() == QEvent::HoverMove) 
  {// activate a handle manipulated frame
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    const QPoint& p = viewer->mapFromGlobal(QCursor::pos());
    bool need_repaint = activate_closest_manipulated_frame(p.x(), p.y());

    if(need_repaint) { Q_EMIT itemChanged(); }
  }

  return false;
}

#include "opengl_tools.h"
void Scene_edit_polyhedron_item::draw_edges(Viewer_interface* viewer) const {
    if(!are_buffers_filled)
        initialize_buffers(viewer);
    vaos[2]->bind();
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    program->setAttributeValue("colors", QColor(0,0,0));
    viewer->glDrawElements(GL_LINES, (GLsizei) edges.size(), GL_UNSIGNED_INT, edges.data());
    program->release();
    vaos[2]->release();

  if(rendering_mode == Wireframe) {
        draw_ROI_and_control_vertices(viewer);
  }
}
void Scene_edit_polyhedron_item::draw(Viewer_interface* viewer) const {
    if(!are_buffers_filled)
        initialize_buffers(viewer);
    vaos[0]->bind();
    program = getShaderProgram(PROGRAM_WITH_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITH_LIGHT);
    program->bind();
    QColor color = this->color();
    program->setAttributeValue("colors", color);
    viewer->glDrawElements(GL_TRIANGLES, (GLsizei) tris.size(), GL_UNSIGNED_INT, tris.data());
    program->release();
    vaos[0]->release();
    draw_edges(viewer);
    draw_ROI_and_control_vertices(viewer);


}

void Scene_edit_polyhedron_item::draw_ROI_and_control_vertices(Viewer_interface* viewer) const {

  CGAL::GL::Color color;
  CGAL::GL::Point_size point_size; point_size.set_point_size(5);

  color.set_rgb_color(0, 1.f, 0);
  if(ui_widget->ShowROICheckBox->isChecked()) {

        if(!ui_widget->ShowAsSphereCheckBox->isChecked() || !viewer->extension_is_found) {

            vaos[1]->bind();
            program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
            attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
            program->bind();
            viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(nb_ROI/3));
            program->release();
            vaos[1]->release();
        }
        else{
            vaos[3]->bind();
            program = getShaderProgram(PROGRAM_INSTANCED);
            attrib_buffers(viewer,PROGRAM_INSTANCED);
            program->bind();
            viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                        static_cast<GLsizei>(nb_sphere/3),
                                        static_cast<GLsizei>(nb_ROI/3));
            program->release();
            vaos[3]->release();
    }
  }

    if(!ui_widget->ShowAsSphereCheckBox->isChecked() || !viewer->extension_is_found) {
        vaos[5]->bind();
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
        program->bind();
        viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(nb_control/3));
        program->release();
        vaos[5]->release();
    }
    else{
        vaos[6]->bind();
        program = getShaderProgram(PROGRAM_INSTANCED);
        attrib_buffers(viewer,PROGRAM_INSTANCED);
        program->bind();
        viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                    static_cast<GLsizei>(nb_sphere/3),
                                    static_cast<GLsizei>(nb_control/3));
        program->release();
        vaos[6]->release();
    }

    QGLViewer* viewerB = *QGLViewer::QGLViewerPool().begin();
  for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
  {
        if(hgb_data->frame == viewerB->manipulatedFrame())
    {      
            GLfloat f_matrix[16];
            for(int i =0; i<16; i++)
                f_matrix[i] = hgb_data->frame->matrix()[i];
            QMatrix4x4 f_mat;
                for(int i=0; i<16; i++)
                    f_mat.data()[i] = (float)f_matrix[i];
            vaos[7]->bind();
            program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
            attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
            program->bind();
            program->setUniformValue("f_matrix", f_mat);
            viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_axis/3));
            program->release();
            vaos[7]->release();

            //QGLViewer::drawAxis(length_of_axis);
      // draw bbox
      if(!ui_widget->ActivatePivotingCheckBox->isChecked())
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
                vaos[4]->bind();
                bbox_program.bind();
                bbox_program.setUniformValue("rotations", f_mat);
                bbox_program.setUniformValue("translation", vec);
                bbox_program.setUniformValue("translation_2", vec2);
                bbox_program.setUniformValue("mvp_matrix", mvp_mat);
                viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_bbox/3));
                bbox_program.release();
                vaos[4]->release();
    }
    }
  }

  } 


void Scene_edit_polyhedron_item::compute_bbox(const Scene_interface::Bbox& bb){
    pos_bbox.resize(24*3);

    pos_bbox[0]=bb.xmin; pos_bbox[1]=bb.ymin; pos_bbox[2]=bb.zmin;
    pos_bbox[3]=bb.xmax; pos_bbox[4]=bb.ymin; pos_bbox[5]=bb.zmin;
    pos_bbox[6]=bb.xmin; pos_bbox[7]=bb.ymin; pos_bbox[8]=bb.zmin;
    pos_bbox[9]=bb.xmin; pos_bbox[10]=bb.ymax; pos_bbox[11]=bb.zmin;
    
    pos_bbox[12]=bb.xmin; pos_bbox[13]=bb.ymin; pos_bbox[14]=bb.zmin;
    pos_bbox[15]=bb.xmin; pos_bbox[16]=bb.ymin; pos_bbox[17]=bb.zmax;
    pos_bbox[18]= bb.xmax; pos_bbox[19]=bb.ymin; pos_bbox[20]=bb.zmin;
    pos_bbox[21]= bb.xmax; pos_bbox[22]=bb.ymax; pos_bbox[23]=bb.zmin;
    
    pos_bbox[24]= bb.xmax; pos_bbox[25]=bb.ymin; pos_bbox[26]=bb.zmin;
    pos_bbox[27]= bb.xmax; pos_bbox[28]=bb.ymin; pos_bbox[29]=bb.zmax;
    pos_bbox[30]=bb.xmin; pos_bbox[31]=bb.ymax; pos_bbox[32]=bb.zmin;
    pos_bbox[33]=bb.xmax; pos_bbox[34]=bb.ymax; pos_bbox[35]=bb.zmin;
    
    pos_bbox[36]=bb.xmin; pos_bbox[37]=bb.ymax; pos_bbox[38]=bb.zmin;
    pos_bbox[39]=bb.xmin; pos_bbox[40]=bb.ymax; pos_bbox[41]=bb.zmax;
    pos_bbox[42]=bb.xmin; pos_bbox[43]=bb.ymin; pos_bbox[44]=bb.zmax;
    pos_bbox[45]=bb.xmax; pos_bbox[46]=bb.ymin; pos_bbox[47]=bb.zmax;

    pos_bbox[48]=bb.xmin; pos_bbox[49]=bb.ymin; pos_bbox[50]=bb.zmax;
    pos_bbox[51]=bb.xmin; pos_bbox[52]=bb.ymax; pos_bbox[53]=bb.zmax;
    pos_bbox[54]=bb.xmax; pos_bbox[55]=bb.ymax; pos_bbox[56]=bb.zmax;
    pos_bbox[57]=bb.xmin; pos_bbox[58]=bb.ymax; pos_bbox[59]=bb.zmax;

    pos_bbox[60]=bb.xmax; pos_bbox[61]=bb.ymax; pos_bbox[62]=bb.zmax;
    pos_bbox[63]=bb.xmax; pos_bbox[64]=bb.ymin; pos_bbox[65]=bb.zmax;
    pos_bbox[66]=bb.xmax; pos_bbox[67]=bb.ymax; pos_bbox[68]=bb.zmax;
    pos_bbox[69]=bb.xmax; pos_bbox[70]=bb.ymax; pos_bbox[71]=bb.zmin;
    
}

void Scene_edit_polyhedron_item::invalidate_buffers()
{
    compute_normals_and_vertices();
    update_normals();
    are_buffers_filled = false;
}

Scene_polyhedron_item* Scene_edit_polyhedron_item::to_polyhedron_item() {
  Scene_polyhedron_item* poly_item_tmp = poly_item;
  poly_item->set_color_vector_read_only(false);
  own_poly_item=false;
  return poly_item_tmp;
}

Polyhedron* Scene_edit_polyhedron_item::polyhedron()       
{ return poly_item->polyhedron(); }
const Polyhedron* Scene_edit_polyhedron_item::polyhedron() const 
{ return poly_item->polyhedron(); }
QString Scene_edit_polyhedron_item::toolTip() const
{
  if(!poly_item->polyhedron())
    return QString();

  return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
    .arg(this->name())
    .arg(poly_item->polyhedron()->size_of_vertices())
    .arg(poly_item->polyhedron()->size_of_halfedges()/2)
    .arg(poly_item->polyhedron()->size_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}
bool Scene_edit_polyhedron_item::isEmpty() const {
  return poly_item->isEmpty();
}
Scene_edit_polyhedron_item::Bbox Scene_edit_polyhedron_item::bbox() const {
  return poly_item->bbox();
}

void Scene_edit_polyhedron_item::setVisible(bool b) {
  poly_item->setVisible(b);
  Scene_item::setVisible(b);
  if(!b) {
    (*QGLViewer::QGLViewerPool().begin())->setManipulatedFrame(NULL);
  }
}
void Scene_edit_polyhedron_item::setColor(QColor c) {
  poly_item->setColor(c);
  Scene_item::setColor(c);
}
void Scene_edit_polyhedron_item::setName(QString n) {
  Scene_item::setName(n);
  n.replace(" (edit)", "");
  poly_item->setName(n);
}
void Scene_edit_polyhedron_item::setRenderingMode(RenderingMode m) {
  poly_item->setRenderingMode(m);
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
  poly_item->select(orig_x,
                       orig_y,
                       orig_z,
                       dir_x,
                       dir_y,
                       dir_z);
}

bool Scene_edit_polyhedron_item::keyPressEvent(QKeyEvent* e)
{
  //setting/unsetting rotation constraints
  if (e->key()==Qt::Key_R && !state.ctrl_pressing)
  {
    is_rot_free = !is_rot_free;
    rot_constraint.setRotationConstraintType( is_rot_free?
        qglviewer::AxisPlaneConstraint::FREE:
        qglviewer::AxisPlaneConstraint::AXIS);
    return true;
  }

  return false;
}

void Scene_edit_polyhedron_item::create_Sphere(double R)
{
  create_flat_sphere(R, pos_sphere, normals_sphere);
}

//#include "Scene_edit_polyhedron_item.moc"
