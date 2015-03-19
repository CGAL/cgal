#include "Scene_polyhedron_transform_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

Scene_polyhedron_transform_item::Scene_polyhedron_transform_item(const qglviewer::Vec& pos,const Scene_polyhedron_item* poly_item_,const Scene_interface*):
    poly_item(poly_item_),
    manipulable(false),
    frame(new ManipulatedFrame()),
    positions_lines(0),
    poly(poly_item->polyhedron()),
    center_(pos) {
    frame->setPosition(pos);

    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(1, buffer);
    compile_shaders();

}

Scene_polyhedron_transform_item::~Scene_polyhedron_transform_item()
{

    glDeleteBuffers(1, buffer);
    glDeleteVertexArrays(1, vao);
    glDeleteProgram(rendering_program);
    delete frame;
    emit killed();
}

void Scene_polyhedron_transform_item::initialize_buffers()
{

    glBindVertexArray(vao[0]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_lines.size())*sizeof(float),
                 positions_lines.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0, //number of the buffer
                          3, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(0);
}
void Scene_polyhedron_transform_item::compile_shaders()
{
    //fill the vertex shader
    static const GLchar* vertex_shader_source_lines[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec3 positions; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 f_matrix; \n"
        "uniform vec3 color_lines; \n"
        "vec4 positions_lines = vec4(positions, 1.0); \n"
        "out highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color_lines; \n"
        "   gl_Position = mvp_matrix * f_matrix * positions_lines; \n"
        "} \n"
    };

    //fill the fragment shader
    static const GLchar* fragment_shader_source[]=
    {
        "#version 300 es \n"
        " \n"
        "in highp vec3 fColors; \n"

        "out highp vec3 color; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        " color = fColors; \n"
        "} \n"
    };

    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_lines, NULL);
    glCompileShader(vertex_shader);
    GLuint fragment_shader =	glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    GLuint program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    rendering_program = program;

}
void Scene_polyhedron_transform_item::uniform_attrib(Viewer_interface* viewer) const
{
    GLfloat mvp_mat[16];
    GLfloat f_mat[16];
    GLdouble d_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    for (int i=0; i<16; ++i){
        mvp_mat[i] = GLfloat(d_mat[i]);
        f_mat[i] = frame->matrix()[i];
    }

    GLfloat colors[3];
    colors[0] = this->color().redF();
    colors[1] = this->color().greenF();
    colors[2] = this->color().blueF();
    glUseProgram(rendering_program);
    glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
    glUniform3fv(location[1], 1, colors);
    glUniformMatrix4fv(location[2], 1, GL_FALSE, f_mat);
}
void Scene_polyhedron_transform_item::compute_elements()
{
     positions_lines.clear();
    typedef Kernel::Point_3		Point;
    typedef Polyhedron::Edge_const_iterator	Edge_iterator;

    Edge_iterator he;
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        positions_lines.push_back(a.x()-center_.x);
        positions_lines.push_back(a.y()-center_.y);
        positions_lines.push_back(a.z()-center_.z);

        positions_lines.push_back(b.x()-center_.x);
        positions_lines.push_back(b.y()-center_.y);
        positions_lines.push_back(b.z()-center_.z);

    }

    location[0] = glGetUniformLocation(rendering_program, "mvp_matrix");
    location[1] = glGetUniformLocation(rendering_program, "color_lines");
    location[2] = glGetUniformLocation(rendering_program, "f_matrix");
}


void Scene_polyhedron_transform_item::draw() const{
    glPushMatrix();
    glMultMatrixd(frame->matrix());
    direct_draw_edges();
    //Scene_item_with_display_list::draw();
    glPopMatrix();
}
void Scene_polyhedron_transform_item::draw_edges(Viewer_interface* viewer) const
{
    glBindVertexArray(vao[0]);
    glUseProgram(rendering_program);
    uniform_attrib(viewer);
    glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);

}
void Scene_polyhedron_transform_item::direct_draw_edges() const {
    typedef Kernel::Point_3		Point;
    typedef Polyhedron::Edge_const_iterator	Edge_iterator;

    ::glDisable(GL_LIGHTING);
    ::glBegin(GL_LINES);
    Edge_iterator he;
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        ::glVertex3d(a.x()-center_.x,a.y()-center_.y,a.z()-center_.z);
        ::glVertex3d(b.x()-center_.x,b.y()-center_.y,b.z()-center_.z);
    }
    ::glEnd();
    ::glEnable(GL_LIGHTING);
}
QString Scene_polyhedron_transform_item::toolTip() const {
    return QObject::tr("<p>Affine transformation of <b>%1</b></p>"
                       "<p>Keep <b>Ctrl</b> pressed and use the arcball to define an affine transformation.<br />"
                       "Press <b>S</b> to apply the affine transformation to a copy of <b>%1</b>.</p>")
            .arg(getBase()->name());
}
bool Scene_polyhedron_transform_item::keyPressEvent(QKeyEvent* e){
    if (e->key()==Qt::Key_S){
        emit stop();
        return true;
    }
    return false;
}

Scene_polyhedron_transform_item::Bbox
Scene_polyhedron_transform_item::bbox() const {
    const Kernel::Point_3& p = *(poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedron::Point_const_iterator it = poly->points_begin();
        it != poly->points_end();
        ++it) {
        bbox = bbox + it->bbox();
    }
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}


void Scene_polyhedron_transform_item::changed()
{
    compute_elements();
    initialize_buffers();
}
#include "Scene_polyhedron_transform_item.moc"

