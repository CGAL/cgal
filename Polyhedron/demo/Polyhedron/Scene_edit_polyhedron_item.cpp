#include "opengl_tools.h"
#include "Scene_edit_polyhedron_item.h"
#include <boost/foreach.hpp>
#include <algorithm>


#include <CGAL/gl_render.h>
struct light_info
{
    //position
    GLfloat position[4];

    //ambient
    GLfloat ambient[4];

    //diffuse
    GLfloat diffuse[4];

    //specular
    GLfloat specular[4];
};
Scene_edit_polyhedron_item::Scene_edit_polyhedron_item
(Scene_polyhedron_item* poly_item,
 Ui::DeformMesh* ui_widget,
 QMainWindow* mw)
    : ui_widget(ui_widget),
      poly_item(poly_item),
      deform_mesh(*(poly_item->polyhedron()), Deform_mesh::Vertex_index_map(), Deform_mesh::Hedge_index_map(), Array_based_vertex_point_map(&positions)),
      is_rot_free(true),
      own_poly_item(true),
      ROI_points(0),
      control_points(0),
      control_color(0),
      ROI_color(0),
      k_ring_selector(poly_item, mw, Scene_polyhedron_item_k_ring_selection::Active_handle::VERTEX, true),
      quadric(gluNewQuadric())
{
    mw->installEventFilter(this);
    gluQuadricNormals(quadric, GLU_SMOOTH);
    // bind vertex picking
    connect(&k_ring_selector, SIGNAL(selected(const std::set<Polyhedron::Vertex_handle>&)), this,
            SLOT(selected(const std::set<Polyhedron::Vertex_handle>&)));

    poly_item->set_color_vector_read_only(true); // to prevent recomputation of color vector in changed()
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
                compute_vertex_normal<Polyhedron::Vertex, Polyhedron::Traits>(*vb);

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
    glGenVertexArrays(3, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(20, buffer);
    compile_shaders();
    changed();
}

Scene_edit_polyhedron_item::~Scene_edit_polyhedron_item()
{
    glDeleteBuffers(4, buffer);
    glDeleteVertexArrays(1, vao);
    glDeleteProgram(rendering_program_facets);
    glDeleteProgram(rendering_program_lines);
    while(is_there_any_ctrl_vertices_group())
    {
        delete_ctrl_vertices_group(false);
    }
    gluDeleteQuadric(quadric);
    if (own_poly_item) delete poly_item;

}
/////////////////////////////
/// For the Shader gestion///
void Scene_edit_polyhedron_item::initialize_buffers()
{
    glBindVertexArray(vao[0]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions.size())*sizeof(double),
                 positions.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0, //number of the buffer
                          3, //number of floats to be taken
                          GL_DOUBLE, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[1]);
    glBufferData(GL_ARRAY_BUFFER,
                 (normals.size())*sizeof(double),
                 normals.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[2]);
    glBufferData(GL_ARRAY_BUFFER,
                 (ROI_points.size())*sizeof(double),
                 ROI_points.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(2,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[4]);
    glBufferData(GL_ARRAY_BUFFER,
                 (ROI_color.size())*sizeof(double),
                 ROI_color.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(3,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[8]);
    glBufferData(GL_ARRAY_BUFFER,
                 (color_edges.size())*sizeof(double),
                 color_edges.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(4,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(4);



    // Clean-up
    glBindVertexArray(0);

    glBindVertexArray(vao[1]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[6]);
    glBufferData(GL_ARRAY_BUFFER,
                 (pos_bbox.size())*sizeof(double),
                 pos_bbox.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0, //number of the buffer
                          3, //number of floats to be taken
                          GL_DOUBLE, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, buffer[3]);
    glBufferData(GL_ARRAY_BUFFER,
                 (control_points.size())*sizeof(double),
                 control_points.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(2,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[5]);
    glBufferData(GL_ARRAY_BUFFER,
                 (control_color.size())*sizeof(double),
                 control_color.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(3,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[9]);
    glBufferData(GL_ARRAY_BUFFER,
                 (color_bbox.size())*sizeof(double),
                 color_bbox.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(4,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(4);
    glBindVertexArray(0);



    glBindVertexArray(vao[2]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[7]);
    glBufferData(GL_ARRAY_BUFFER,
                 (pos_axis.size())*sizeof(double),
                 pos_axis.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[10]);
    glBufferData(GL_ARRAY_BUFFER,
                 (color_lines.size())*sizeof(double),
                 color_lines.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(4,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(4);
    glBindVertexArray(0);




}

void Scene_edit_polyhedron_item::compile_shaders(void)
{
    //fill the vertex shader
    static const GLchar* vertex_shader_source[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec3 positions; \n"
        "layout (location = 1) in vec3 vNormals; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 mv_matrix; \n"
        "uniform vec3 u_color;  \n"
        "uniform int is_two_side; \n"
        "uniform vec3 light_pos;  \n"
        "uniform vec3 light_diff; \n"
        "uniform vec3 light_spec; \n"
        "uniform vec3 light_amb;  \n"
        "float spec_power = 128.0; \n"
        "vec4 positions_facets = vec4(positions, 1.0); \n"
        "out highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   vec4 P = mv_matrix * positions_facets; \n"
        "   vec3 N = mat3(mv_matrix)* vNormals; \n"
        "   vec3 L = light_pos - P.xyz; \n"
        "   vec3 V = -P.xyz; \n"

        "   N = normalize(N); \n"
        "   L = normalize(L); \n"
        "   V = normalize(V); \n"

        "   vec3 R = reflect(-L, N); \n"
        "   vec3 diffuse; \n"
        "   if(is_two_side == 1) \n"
        "       diffuse = abs(dot(N,L)) * light_diff * u_color; \n"
        "   else \n"
        "       diffuse = max(dot(N,L), 0.0) * light_diff * u_color; \n"
        "   vec3 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

        "   fColors = light_amb * u_color + diffuse + specular ; \n"

        "   gl_Position =  mvp_matrix *positions_facets; \n"
        "} \n"
    };
    //fill the fragment shader
    static const GLchar* fragment_shader_source[]=
    {
        "#version 300 es \n"
        " \n"
        "precision mediump float; \n"
        "in vec3 fColors; \n"

        "out vec3 color; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        " color = fColors; \n"
        "} \n"
    };

    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source, NULL);
    glCompileShader(vertex_shader);
    GLuint fragment_shader =	glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    //creates the program, attaches and links the shaders
    GLuint program= glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    //Clean-up
    glDeleteShader(vertex_shader);

    rendering_program_facets = program;

    //For the edges
    static const GLchar* vertex_shader_source_lines[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec3 positions_lines; \n"
        "layout (location = 4) in vec3 color_lines; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 rotations; \n"
        "uniform vec3 translation; \n"
        "uniform vec3 translation_2; \n"
        "out highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color_lines; \n"
        "   gl_Position = mvp_matrix * (rotations *(vec4(translation_2,0.0)+vec4(positions_lines,1.0) )+ vec4(translation,0.0)) ; \n"
        "} \n"
    };

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_lines, NULL);
    glCompileShader(vertex_shader);

    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    //Clean-up
    glDeleteShader(vertex_shader);
    rendering_program_lines = program;

    //For the points
    static const GLchar* vertex_shader_source_points[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 2) in vec3 positions; \n"
        "layout (location = 3) in vec3 color; \n"

        "uniform mat4 mvp_matrix; \n"

        "out highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color; \n"
        "   gl_Position = mvp_matrix * vec4(positions,1.0); \n"
        "} \n"
    };

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_points, NULL);
    glCompileShader(vertex_shader);

    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    //Clean-up
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    rendering_program_points = program;

}

void Scene_edit_polyhedron_item::compute_normals_and_vertices(void)
{
    ROI_points.clear();
    control_points.clear();
    BOOST_FOREACH(vertex_descriptor vd, deform_mesh.roi_vertices())
    {
        if(!deform_mesh.is_control_vertex(vd))
        {//gl_draw_point( vd->point() );
            ROI_points.push_back(vd->point().x());
            ROI_points.push_back(vd->point().y());
            ROI_points.push_back(vd->point().z());
            ROI_color.push_back(0.0);
            ROI_color.push_back(1.0);
            ROI_color.push_back(0);
        }
    }
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
    {
        if(hgb_data->frame == viewer->manipulatedFrame())
        {
            // draw axis

            if(ui_widget->ActivatePivotingCheckBox->isChecked())
            {
                // draw bbox
                compute_bbox(hgb_data->bbox);
            }
        }
        // draw control vertices
        if(hgb_data == active_group)
        {
            //set color to red
            control_color.push_back(1.0);
            control_color.push_back(0.0);
            control_color.push_back(0.0);
        }
        else
        {
            //set color to blue
            control_color.push_back(0.0);
            control_color.push_back(0.0);
            control_color.push_back(1.0);
        }
        for(std::vector<vertex_descriptor>::const_iterator hb = hgb_data->ctrl_vertices_group.begin(); hb != hgb_data->ctrl_vertices_group.end(); ++hb)
        {
            control_points.push_back((*hb)->point().x());
            control_points.push_back((*hb)->point().y());
            control_points.push_back((*hb)->point().z());

        }
    }

    //The edges color
    color_edges.resize(edges.size());
    for(int i =0; i< edges.size(); i++)
        color_edges[i]=0.0;

    //The box color
    color_bbox.resize(pos_bbox.size());
    for(int i =0; i< pos_bbox.size(); i++)
        color_bbox[i]=0.0;

    for(int i =0; i< pos_bbox.size(); i+=3)
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

    location[0] = glGetUniformLocation(rendering_program_facets, "mvp_matrix");
    location[1] = glGetUniformLocation(rendering_program_facets, "mv_matrix");
    location[2] = glGetUniformLocation(rendering_program_facets, "light_pos");
    location[3] = glGetUniformLocation(rendering_program_facets, "light_diff");
    location[4] = glGetUniformLocation(rendering_program_facets, "light_spec");
    location[5] = glGetUniformLocation(rendering_program_facets, "light_amb");
    location[6] = glGetUniformLocation(rendering_program_facets, "is_two_side");
    location[8] = glGetUniformLocation(rendering_program_facets, "u_color");

    location[7] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
    location[11] = glGetUniformLocation(rendering_program_lines, "rotations");
    location[13] = glGetUniformLocation(rendering_program_lines, "translation");
    location[14] = glGetUniformLocation(rendering_program_lines, "translation_2");

    location[10] = glGetUniformLocation(rendering_program_points, "mvp_matrix");



}

void Scene_edit_polyhedron_item::uniform_attrib(Viewer_interface* viewer, int mode) const
{

    light_info light;
    GLint is_both_sides = 0;
    GLfloat mvp_mat[16];
    GLfloat mv_mat[16];

    //fills the MVP and MV matrices.

    GLdouble d_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i){
        mvp_mat[i] = GLfloat(d_mat[i]);
    }

    viewer->camera()->getModelViewMatrix(d_mat);
    for (int i=0; i<16; ++i)
        mv_mat[i] = GLfloat(d_mat[i]);

    glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE, &is_both_sides);


    //Gets lighting info :

    //position
    glGetLightfv(GL_LIGHT0, GL_POSITION, light.position);

    //ambient
    glGetLightfv(GL_LIGHT0, GL_AMBIENT, light.ambient);


    //specular
    glGetLightfv(GL_LIGHT0, GL_SPECULAR, light.specular);

    //diffuse
    glGetLightfv(GL_LIGHT0, GL_DIFFUSE, light.diffuse);


    if(mode ==0)
    {
        GLfloat color[3];
        color[0] = this->color().redF();
        color[1] = this->color().greenF();
        color[2] = this->color().blueF();

        glUseProgram(rendering_program_facets);

        glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
        glUniformMatrix4fv(location[1], 1, GL_FALSE, mv_mat);

        glUniform3fv(location[2], 1, light.position);
        glUniform3fv(location[3], 1, light.diffuse);
        glUniform3fv(location[4], 1, light.specular);
        glUniform3fv(location[5], 1, light.ambient);
        glUniform1i(location[6], is_both_sides);
        glUniform3fv(location[8], 1, color);

    }
    else if(mode ==1)
    {
        glUseProgram(rendering_program_lines);
        glUniformMatrix4fv(location[7], 1, GL_FALSE, mvp_mat);
    }
    else if(mode ==2)
    {
        glUseProgram(rendering_program_points);
        glUniformMatrix4fv(location[10], 1, GL_FALSE, mvp_mat);
    }

}


/////////////////////////////////////////////////////////
/////////// Most relevant functions lie here ///////////
void Scene_edit_polyhedron_item::deform()
{
    if(!is_there_any_ctrl_vertices()) { return; }

    for(Ctrl_vertices_group_data_list::iterator it = ctrl_vertex_frame_map.begin(); it != ctrl_vertex_frame_map.end(); ++it)
    { it->set_target_positions(); }
    deform_mesh.deform();

    poly_item->changed(); // now we need to call poly_item changed to delete AABB tree
    emit itemChanged();
}

void Scene_edit_polyhedron_item::timerEvent(QTimerEvent* /*event*/)
{ // just handle deformation - paint like selection is handled in eventFilter()
    if(state.ctrl_pressing && (state.left_button_pressing || state.right_button_pressing)) {
        if(!ui_widget->ActivatePivotingCheckBox->isChecked()) {
            deform();
        }
        else {
            emit itemChanged(); // for redraw while Pivoting (since we close signals of manipulatedFrames while pivoting,
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

        if(need_repaint) { emit itemChanged(); }
    }

    return false;
}

#include "opengl_tools.h"
void Scene_edit_polyhedron_item::draw_edges(Viewer_interface* viewer) const {
    GLfloat vec[3];
    for(int i=0; i< 3; i++)
        vec[i]=0.0;
    GLfloat f_matrix[16];
    for(int i=0; i<16; i++)
        f_matrix[i]=0.0;
    f_matrix[0]=1.0; f_matrix[5]=1.0; f_matrix[10]=1.0; f_matrix[15]=1.0;

    glBindVertexArray(vao[0]);
    glUseProgram(rendering_program_lines);
    glUniform3fv(location[13],1,vec);
    glUniform3fv(location[14],1,vec);
    glUniformMatrix4fv(location[11], 1, GL_FALSE, f_matrix);
    uniform_attrib(viewer,1);
    glDrawElements(GL_LINES, (GLsizei) edges.size(), GL_UNSIGNED_INT, edges.data());
    glUseProgram(0);
    glBindVertexArray(0);

    if(rendering_mode == Wireframe) {
        draw_ROI_and_control_vertices(viewer);
    }
}
void Scene_edit_polyhedron_item::draw(Viewer_interface* viewer) const {

    glBindVertexArray(vao[0]);
    glUseProgram(rendering_program_facets);
    uniform_attrib(viewer,0);
    glDrawElements(GL_TRIANGLES, (GLsizei) tris.size(), GL_UNSIGNED_INT, tris.data());
    glUseProgram(0);
    glBindVertexArray(0);
    draw_edges(viewer);
    draw_ROI_and_control_vertices(viewer);


}

void Scene_edit_polyhedron_item::draw_ROI_and_control_vertices(Viewer_interface* viewer) const {

    GLboolean enable_back_lighting = glIsEnabled(GL_LIGHTING);
    (GL_LIGHTING);

    CGAL::GL::Color color;
    CGAL::GL::Point_size point_size; point_size.set_point_size(5);

    color.set_rgb_color(0, 1.f, 0);
    if(!ui_widget->ShowROICheckBox->isChecked()) {


        glBindVertexArray(vao[0]);
        glUseProgram(rendering_program_points);
        uniform_attrib(viewer,2);
        glDrawArrays(GL_POINTS, 0, ROI_points.size()/3);
        glUseProgram(0);

        glBindVertexArray(0);
    }


    glBindVertexArray(vao[1]);
    glUseProgram(rendering_program_points);
    uniform_attrib(viewer,2);
    glDrawArrays(GL_POINTS, 0, control_points.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);

    QGLViewer* viewerB = *QGLViewer::QGLViewerPool().begin();
    for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
    {
        if(hgb_data->frame == viewerB->manipulatedFrame())
        {
            //Draw the axis
            GLfloat vec[3];
            for(int i=0; i< 3; i++)
                vec[i]=0.0;
            GLfloat f_matrix[16];
            for(int i =0; i<16; i++)
                f_matrix[i] = hgb_data->frame->matrix()[i];

            glBindVertexArray(vao[2]);
            glUseProgram(rendering_program_lines);
            glUniform3fv(location[13], 1, vec);
            glUniform3fv(location[14], 1, vec);
            glUniformMatrix4fv(location[11], 1, GL_FALSE, f_matrix);
            uniform_attrib(viewer,1);
            glDrawArrays(GL_LINES, 0, pos_axis.size()/3);
            glUseProgram(0);
            glBindVertexArray(0);

            //QGLViewer::drawAxis(length_of_axis);
            // draw bbox
            if(!ui_widget->ActivatePivotingCheckBox->isChecked())
            {
                GLfloat colors[3];
                GLfloat f_matrix[16];
                GLfloat trans[3];
                GLfloat trans2[3];
                colors[0]=1.0;
                colors[1]=0.0;
                colors[2]=0.0;

                trans[0] = hgb_data->frame->position().x;
                trans[1] = hgb_data->frame->position().y;
                trans[2] = hgb_data->frame->position().z;

                trans2[0] = -hgb_data->frame_initial_center.x;
                trans2[1] = -hgb_data->frame_initial_center.y;
                trans2[2] = -hgb_data->frame_initial_center.z;

                for(int i =0; i<16; i++)
                    f_matrix[i] = hgb_data->frame->orientation().matrix()[i];

                glBindVertexArray(vao[1]);
                glUseProgram(rendering_program_lines);
                glUniform3fv(location[13], 1, trans);
                glUniform3fv(location[14], 1, trans2);
                glUniformMatrix4fv(location[11], 1, GL_FALSE, f_matrix);
                uniform_attrib(viewer,1);
                glDrawArrays(GL_LINES, 0, pos_bbox.size()/3);
                glUseProgram(0);
                glBindVertexArray(0);
            }
        }
    }

    /* // draw ROI

     if(ui_widget->ShowROICheckBox->isChecked()) {
        BOOST_FOREACH(vertex_descriptor vd, deform_mesh.roi_vertices())
        {
            if(!deform_mesh.is_control_vertex(vd))
                gl_draw_point( vd->point() );
        }
    }
    // draw control vertices related things
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();

    for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
    {
        if(hgb_data->frame == viewer->manipulatedFrame())
        {
            // draw axis
            ::glPushMatrix();
            ::glMultMatrixd(hgb_data->frame->matrix());
            QGLViewer::drawAxis(length_of_axis);
            ::glPopMatrix();
            // draw bbox
            if(!ui_widget->ActivatePivotingCheckBox->isChecked())
            {
                color.set_rgb_color(1.0f, 0, 0);
                ::glPushMatrix();
                ::glTranslated(hgb_data->frame->position().x, hgb_data->frame->position().y, hgb_data->frame->position().z);
                ::glMultMatrixd(hgb_data->frame->orientation().matrix());
                ::glTranslated(-hgb_data->frame_initial_center.x, -hgb_data->frame_initial_center.y, -hgb_data->frame_initial_center.z);
                draw_bbox(hgb_data->bbox);
                ::glPopMatrix();
            }
        }
        // draw control vertices
        if(hgb_data == active_group) { color.set_rgb_color(1.0f, 0, 0); }
        else                    { color.set_rgb_color(0, 0, 1.0f); }
        for(std::vector<vertex_descriptor>::const_iterator hb = hgb_data->ctrl_vertices_group.begin(); hb != hgb_data->ctrl_vertices_group.end(); ++hb)
        {  gl_draw_point( (*hb)->point() );
        }
    }

    if(enable_back_lighting) { glEnable(GL_LIGHTING); }*/
}
void Scene_edit_polyhedron_item::gl_draw_point(const Point& p) const
{
    if(!ui_widget->ShowAsSphereCheckBox->isChecked()) {
        ::glBegin(GL_POINTS);
        ::glVertex3d(p.x(), p.y(), p.z());
        ::glEnd();
    }
    else {
        GLint shading;
        ::glGetIntegerv(GL_SHADE_MODEL, &shading);
        ::glShadeModel(GL_SMOOTH);

        ::glPushMatrix();
        ::glTranslated(p.x(), p.y(), p.z());
        ::gluSphere(quadric, length_of_axis/15, 8, 8);
        ::glPopMatrix();

        ::glShadeModel(shading);
    }
}
//////////////////////////////////////////////////////////

/////////////// from trivial_plugin //////////////////////
void Scene_edit_polyhedron_item::draw_bbox(const Scene_interface::Bbox &bb)const{

    ::glBegin(GL_LINES);
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                 bb.xmax, bb.ymin, bb.zmin);
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                 bb.xmin, bb.ymax, bb.zmin);
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                 bb.xmin, bb.ymin, bb.zmax);

    gl_draw_edge(bb.xmax, bb.ymin, bb.zmin,
                 bb.xmax, bb.ymax, bb.zmin);
    gl_draw_edge(bb.xmax, bb.ymin, bb.zmin,
                 bb.xmax, bb.ymin, bb.zmax);

    gl_draw_edge(bb.xmin, bb.ymax, bb.zmin,
                 bb.xmax, bb.ymax, bb.zmin);
    gl_draw_edge(bb.xmin, bb.ymax, bb.zmin,
                 bb.xmin, bb.ymax, bb.zmax);

    gl_draw_edge(bb.xmin, bb.ymin, bb.zmax,
                 bb.xmax, bb.ymin, bb.zmax);
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmax,
                 bb.xmin, bb.ymax, bb.zmax);

    gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                 bb.xmin, bb.ymax, bb.zmax);
    gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                 bb.xmax, bb.ymin, bb.zmax);
    gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                 bb.xmax, bb.ymax, bb.zmin);


    ::glEnd();

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

void Scene_edit_polyhedron_item::gl_draw_edge(double px, double py, double pz,
                                              double qx, double qy, double qz) const
{
    ::glVertex3d(px,py,pz);
    ::glVertex3d(qx,qy,qz);
}
/////////////////////////////////////////////////////////////

void Scene_edit_polyhedron_item::changed()
{
    compute_normals_and_vertices();
    initialize_buffers();
    update_normals();
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

#include "Scene_edit_polyhedron_item.moc"
