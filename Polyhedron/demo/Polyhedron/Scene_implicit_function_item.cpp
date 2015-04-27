#include "Scene_implicit_function_item.h"
#include <QColor>
#include <map>
#include <CGAL/gl.h>
#include <CGAL/Simple_cartesian.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

#include "Color_ramp.h"
#include <Viewer_interface.h>

#include <CGAL/double.h>

inline
bool is_nan(double d)
{
    return !CGAL::Is_valid<double>()( d );
}

void Scene_implicit_function_item::initialize_buffers()
{
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_tex_quad.size())*sizeof(float),
                 positions_tex_quad.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0, //number of the buffer
                          4, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[1]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_cube.size())*sizeof(float),
                 positions_cube.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(1, //number of the buffer
                          4, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(1);



    glBindBuffer(GL_ARRAY_BUFFER, buffer[2]);
    glBufferData(GL_ARRAY_BUFFER,
                 (texture_map.size())*sizeof(float),
                 texture_map.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(2,
                          2,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[3]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_grid.size())*sizeof(float),
                 positions_grid.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(3,
                          4,
                          GL_FLOAT,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(3);


    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGB,
                 texture->getWidth(),
                 texture->getHeight(),
                 0,
                 GL_RGB,
                 GL_UNSIGNED_BYTE,
                 texture->getData());
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Clean-up
    glBindVertexArray(0);
}

void Scene_implicit_function_item::compile_shaders(void)
{
    //fill the vertex shader
    static const GLchar* vertex_shader_source[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec4 positions_tex_quad; \n"
        "layout (location = 2) in vec2 v_texCoord; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 m_matrix; \n"

        "out highp vec2 f_texCoord; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   f_texCoord = v_texCoord; \n"
        "   gl_Position =  mvp_matrix *m_matrix * positions_tex_quad; \n"
        "} \n"
    };
    //fill the fragment shader
    static const GLchar* fragment_shader_source[]=
    {
        "#version 300 es \n"
        " \n"
        "in highp vec2 f_texCoord; \n"
        "uniform sampler2D s_texture; \n"
        "out highp vec4 color; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        " color = texture2D(s_texture, f_texCoord); \n"
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
    GLint result;
    glGetShaderiv(vertex_shader,GL_COMPILE_STATUS,&result);
    if(result == GL_TRUE){
        std::cout<<"Vertex compilation OK"<<std::endl;
    } else {
        int maxLength;
        int length;
        glGetShaderiv(vertex_shader,GL_INFO_LOG_LENGTH,&maxLength);
        char* log = new char[maxLength];
        glGetShaderInfoLog(vertex_shader,maxLength,&length,log);
        std::cout<<"link error : Length = "<<length<<", log ="<<log<<std::endl;
    }
    glGetShaderiv(fragment_shader,GL_COMPILE_STATUS,&result);
    if(result == GL_TRUE){
        std::cout<<"Fragment compilation OK"<<std::endl;
    } else {
        int maxLength;
        int length;
        glGetShaderiv(fragment_shader,GL_INFO_LOG_LENGTH,&maxLength);
        char* log = new char[maxLength];
        glGetShaderInfoLog(fragment_shader,maxLength,&length,log);
        std::cout<<"link error : Length = "<<length<<", log ="<<log<<std::endl;
    }
    //Clean-up
    glDeleteShader(vertex_shader);

    rendering_program_tex_quad = program;

    //For the cube
    static const GLchar* vertex_shader_source_cube[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 1) in vec4 positions_lines; \n"
        "uniform mat4 mvp_matrix; \n"
        "out highp vec3 fColors; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   fColors = vec3(0.6, 0.6, 0.6); \n"
        "   gl_Position = mvp_matrix * positions_lines; \n"
        "} \n"
    };
    static const GLchar* fragment_shader_source_lines[]=
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
    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_cube, NULL);
    glCompileShader(vertex_shader);


    glShaderSource(fragment_shader, 1, fragment_shader_source_lines, NULL);
    glCompileShader(fragment_shader);

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    rendering_program_cube = program;


    //Clean-up
    glDeleteShader(vertex_shader);

    //For the grid
    static const GLchar* vertex_shader_source_grid[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 3) in vec4 positions_lines; \n"
        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 m_matrix; \n"
        "out highp vec3 fColors; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   fColors = vec3(0.6, 0.6, 0.6); \n"
        "   gl_Position = mvp_matrix * m_matrix* positions_lines; \n"
        "} \n"
    };

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_grid, NULL);
    glCompileShader(vertex_shader);


    glShaderSource(fragment_shader, 1, fragment_shader_source_lines, NULL);
    glCompileShader(fragment_shader);

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    rendering_program_grid = program;
    //Clean-up
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);

}

void Scene_implicit_function_item::uniform_attrib(Viewer_interface* viewer, int mode) const
{

    GLfloat mvp_mat[16];
    GLfloat m_mat[16];

    //fills the MVP matrix.
    GLdouble d_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i){
        mvp_mat[i] = GLfloat(d_mat[i]);
    }

    for (int i=0; i<16; ++i){
        m_mat[i] = GLfloat(frame_->matrix()[i]);
    }
    if(mode ==0)
    {
        glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
        glUniformMatrix4fv(location[2], 1, GL_FALSE, m_mat);
        glUniform1i(sampler_location, 0);

    }
    else if(mode ==1)
    {

        glUniformMatrix4fv(location[1], 1, GL_FALSE, mvp_mat);

    }
    else if (mode ==2)
    {
        glUniformMatrix4fv(location[3], 1, GL_FALSE, mvp_mat);
        glUniformMatrix4fv(location[4], 1, GL_FALSE, m_mat);
    }

}
void Scene_implicit_function_item::compute_vertices_and_texmap(void)
{
    positions_tex_quad.clear();
    positions_cube.clear();
    positions_grid.clear();
    texture_map.clear();

    const Bbox& b = bbox();
    float x,y,z;
    z = 0;
    x = (b.xmax-b.xmin)/10.0;
    y = (b.ymax-b.ymin)/10.0;
    // The Quad
    {


        //A
        positions_tex_quad.push_back(b.xmin);
        positions_tex_quad.push_back(b.ymin);
        positions_tex_quad.push_back(z);
        positions_tex_quad.push_back(1.0);

        //B
        positions_tex_quad.push_back(b.xmin);
        positions_tex_quad.push_back(b.ymax);
        positions_tex_quad.push_back(z);
        positions_tex_quad.push_back(1.0);

        //C
        positions_tex_quad.push_back(b.xmax);
        positions_tex_quad.push_back(b.ymin);
        positions_tex_quad.push_back(z);
        positions_tex_quad.push_back(1.0);


        //C
        positions_tex_quad.push_back(b.xmax);
        positions_tex_quad.push_back(b.ymin);
        positions_tex_quad.push_back(z);
        positions_tex_quad.push_back(1.0);

        //B
        positions_tex_quad.push_back(b.xmin);
        positions_tex_quad.push_back(b.ymax);
        positions_tex_quad.push_back(z);
        positions_tex_quad.push_back(1.0);

        //D
        positions_tex_quad.push_back(b.xmax);
        positions_tex_quad.push_back(b.ymax);
        positions_tex_quad.push_back(z);
        positions_tex_quad.push_back(1.0);

        //UV Mapping x2 but I don't know why.
        texture_map.push_back(0.0);
        texture_map.push_back(1.0);

        texture_map.push_back(0.0);
        texture_map.push_back(0.0);

        texture_map.push_back(1.0);
        texture_map.push_back(1.0);

        texture_map.push_back(1.0);
        texture_map.push_back(1.0);

        texture_map.push_back(0.0);
        texture_map.push_back(0.0);

        texture_map.push_back(1.0);
        texture_map.push_back(0.0);

        texture_map.push_back(0.0);
        texture_map.push_back(1.0);

        texture_map.push_back(0.0);
        texture_map.push_back(0.0);

        texture_map.push_back(1.0);
        texture_map.push_back(1.0);

        texture_map.push_back(1.0);
        texture_map.push_back(1.0);

        texture_map.push_back(0.0);
        texture_map.push_back(0.0);

        texture_map.push_back(1.0);
        texture_map.push_back(0.0);
    }
    //The grid
    {

        for(int u = 0; u < 11; u++)
        {

            positions_grid.push_back(b.xmin + x* u);
            positions_grid.push_back(b.ymin);
            positions_grid.push_back(z);
            positions_grid.push_back(1.0);

            positions_grid.push_back(b.xmin + x* u);
            positions_grid.push_back(b.ymax);
            positions_grid.push_back(z);
            positions_grid.push_back(1.0);
        }
        for(int v=0; v<11; v++)
        {

            positions_grid.push_back(b.xmin);
            positions_grid.push_back(b.ymin + v * y);
            positions_grid.push_back(z);
            positions_grid.push_back(1.0);

            positions_grid.push_back(b.xmax);
            positions_grid.push_back(b.ymin + v * y);
            positions_grid.push_back(z);
            positions_grid.push_back(1.0);
        }

    }
    //the Box
    {

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);

        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);
        positions_cube.push_back(1.0);
    }

    //The texture
    for( int i=0 ; i < texture->getWidth() ; i++ )
    {
        for( int j=0 ; j < texture->getHeight() ; j++)
        {
            compute_texture(i,j);
        }
    }

    location[0] = glGetUniformLocation(rendering_program_tex_quad, "mvp_matrix");
    location[2] = glGetUniformLocation(rendering_program_tex_quad, "m_matrix");

    sampler_location = glGetUniformLocation(rendering_program_tex_quad, "s_texture");

    location[1] = glGetUniformLocation(rendering_program_cube, "mvp_matrix");
    location[3] = glGetUniformLocation(rendering_program_grid, "mvp_matrix");
    location[4] = glGetUniformLocation(rendering_program_grid, "m_matrix");

}

Scene_implicit_function_item::
Scene_implicit_function_item(Implicit_function_interface* f)
    : function_(f)
    , frame_(new ManipulatedFrame())
    , need_update_(true)
    , grid_size_(SCENE_IMPLICIT_GRID_SIZE)
    , max_value_(0.)
    , min_value_(0.)
    , blue_color_ramp_()
    , red_color_ramp_()
    , positions_cube(0)
    , positions_grid(0)
    , positions_tex_quad(0)
    , texture_map(0)

{
    texture = new Texture(grid_size_-1,grid_size_-1);
    blue_color_ramp_.build_blue();
    red_color_ramp_.build_red();
    glGenVertexArrays(1, &vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(4, buffer);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &textureId);
    compile_shaders();
    compute_min_max();
    compute_function_grid();
    double offset_x = (bbox().xmin + bbox().xmax) / 2;
    double offset_y = (bbox().ymin + bbox().ymax) / 2;
    double offset_z = (bbox().zmin + bbox().zmax) / 2;
    frame_->setPosition(offset_x, offset_y, offset_z);
    frame_->setOrientation(1., 0, 0, 0);
    connect(frame_, SIGNAL(modified()), this, SLOT(plane_was_moved()));

    changed();
}


Scene_implicit_function_item::~Scene_implicit_function_item()
{
    glDeleteBuffers(4, buffer);
    glDeleteVertexArrays(1, &vao);
    glDeleteProgram(rendering_program_tex_quad);
    glDeleteProgram(rendering_program_cube);
    glDeleteProgram(rendering_program_grid);
    delete frame_;

}


Scene_implicit_function_item::Bbox
Scene_implicit_function_item::bbox() const
{
    return function_->bbox();
}

void
Scene_implicit_function_item::draw(Viewer_interface* viewer) const
{

    //draw_aux(viewer, false);
    glBindVertexArray(vao);
    glUseProgram(rendering_program_tex_quad);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureId);
    uniform_attrib(viewer,0);
    glDrawArrays(GL_TRIANGLES, 0, positions_tex_quad.size()/4);



    // Clean-up
    glUseProgram(0);
    glBindVertexArray(0);

}

void
Scene_implicit_function_item::draw_edges(Viewer_interface* viewer) const
{
    //  draw_aux(viewer, true);
    glBindVertexArray(vao);
    glUseProgram(rendering_program_cube);
    uniform_attrib(viewer,1);
    glDrawArrays(GL_LINES, 0, positions_cube.size()/4);

    glUseProgram(rendering_program_grid);
    uniform_attrib(viewer,2);
    glDrawArrays(GL_LINES, 0, positions_grid.size()/4);

    // Clean-up
    glUseProgram(0);
    glBindVertexArray(0);
}

void
Scene_implicit_function_item::draw_aux(Viewer_interface* viewer, bool edges) const
{
    if(edges) {
        draw_bbox();
        ::glPushMatrix();
        ::glMultMatrixd(frame_->matrix());
        QGLViewer::drawGrid((float)bbox().diagonal_length() * 0.3);
        ::glPopMatrix();
    }

    if(!frame_->isManipulated()) {
        if(need_update_) {
            compute_function_grid();
            need_update_ = false;
        }
        if(!viewer->inFastDrawing()) {
            if(edges)
                Scene_item_with_display_list::draw_edges(viewer);
            else
                Scene_item_with_display_list::draw(viewer);
        }
    }
}

void
Scene_implicit_function_item::direct_draw(Viewer_interface* viewer) const
{
    // draw_function_grid(red_color_ramp_, blue_color_ramp_);

}



QString
Scene_implicit_function_item::toolTip() const
{
    return tr("<p>Function <b>%1</b>")
            .arg(this->name());
}

bool
Scene_implicit_function_item::supportsRenderingMode(RenderingMode m) const
{ 
    switch ( m )
    {
    case Splatting:
    case Gouraud:
        return false;

    case Points:
    case Wireframe:
    case Flat:
    case FlatPlusEdges:
        return true;

    default:
        return false;
    }

    return false;
}

void
Scene_implicit_function_item::
draw_bbox() const
{
    const Bbox& b = bbox();

    ::glDisable(GL_LIGHTING);
    ::glColor3f(0.f,0.f,0.f);
    ::glBegin(GL_LINES);

    ::glVertex3d(b.xmin,b.ymin,b.zmin);
    ::glVertex3d(b.xmin,b.ymin,b.zmax);

    ::glVertex3d(b.xmin,b.ymin,b.zmin);
    ::glVertex3d(b.xmin,b.ymax,b.zmin);

    ::glVertex3d(b.xmin,b.ymin,b.zmin);
    ::glVertex3d(b.xmax,b.ymin,b.zmin);

    ::glVertex3d(b.xmax,b.ymin,b.zmin);
    ::glVertex3d(b.xmax,b.ymax,b.zmin);

    ::glVertex3d(b.xmax,b.ymin,b.zmin);
    ::glVertex3d(b.xmax,b.ymin,b.zmax);

    ::glVertex3d(b.xmin,b.ymax,b.zmin);
    ::glVertex3d(b.xmin,b.ymax,b.zmax);

    ::glVertex3d(b.xmin,b.ymax,b.zmin);
    ::glVertex3d(b.xmax,b.ymax,b.zmin);

    ::glVertex3d(b.xmax,b.ymax,b.zmin);
    ::glVertex3d(b.xmax,b.ymax,b.zmax);

    ::glVertex3d(b.xmin,b.ymin,b.zmax);
    ::glVertex3d(b.xmin,b.ymax,b.zmax);

    ::glVertex3d(b.xmin,b.ymin,b.zmax);
    ::glVertex3d(b.xmax,b.ymin,b.zmax);

    ::glVertex3d(b.xmax,b.ymax,b.zmax);
    ::glVertex3d(b.xmin,b.ymax,b.zmax);

    ::glVertex3d(b.xmax,b.ymax,b.zmax);
    ::glVertex3d(b.xmax,b.ymin,b.zmax);

    ::glEnd();
}

void Scene_implicit_function_item::compute_texture(int i, int j)
{
    const Point& p = (implicit_grid_[i][j]).first;
    double v = (implicit_grid_[i][j]).second;

    if(is_nan(v)) {
        texture->setData(i,j,0.2,0.2,0.2);
    } else
        // determines grey level
        if ( v > 0 )
        {
            v = v/max_value_;
            GLdouble r = red_color_ramp_.r(v), g = red_color_ramp_.g(v), b = red_color_ramp_.b(v);
            texture->setData(i,j,255*r,255*g,255*b);
        }
        else
        {
            v = v/min_value_;
            GLdouble r = blue_color_ramp_.r(v), g = blue_color_ramp_.g(v), b = blue_color_ramp_.b(v);
            texture->setData(i,j,255*r,255*g,255*b);
        }
}

void 
Scene_implicit_function_item::
draw_function_grid(const Color_ramp& ramp_pos,
                   const Color_ramp& ramp_neg) const
{
    ::glDisable(GL_LIGHTING);
    ::glShadeModel(GL_SMOOTH);

    ::glBegin(GL_QUADS);
    const int nb_quads = grid_size_ - 1;
    for( int i=0 ; i < nb_quads ; i++ )
    {
        for( int j=0 ; j < nb_quads ; j++)
        {
            draw_grid_vertex(implicit_grid_[i][j], ramp_pos, ramp_neg);
            draw_grid_vertex(implicit_grid_[i][j+1], ramp_pos, ramp_neg);
            draw_grid_vertex(implicit_grid_[i+1][j+1], ramp_pos, ramp_neg);
            draw_grid_vertex(implicit_grid_[i+1][j], ramp_pos, ramp_neg);
        }
    }
    ::glEnd();
}


void
Scene_implicit_function_item::
draw_grid_vertex(const Point_value& pv,
                 const Color_ramp& ramp_positive,
                 const Color_ramp& ramp_negative) const
{
    const Point& p = pv.first;
    double v = pv.second;

    if(is_nan(v)) {
        ::glColor3f(0.2f, 0.2f, 0.2f);
    } else
        // determines grey level
        if ( v > 0 )
        {
            v = v/max_value_;
            ::glColor3d(ramp_positive.r(v),ramp_positive.g(v),ramp_positive.b(v));
        }
        else
        {
            v = v/min_value_;
            ::glColor3d(ramp_negative.r(v),ramp_negative.g(v),ramp_negative.b(v));
        }

    ::glVertex3d(p.x,p.y,p.z);


}


void
Scene_implicit_function_item::
compute_function_grid() const
{
    typedef CGAL::Simple_cartesian<double>  K;
    typedef K::Aff_transformation_3         Aff_transformation;
    typedef K::Point_3                      Point_3;

    // Get transformation
    const ::GLdouble* m = frame_->matrix();

    // OpenGL matrices are row-major matrices
    Aff_transformation t (m[0], m[4], m[8], m[12],
            m[1], m[5], m[9], m[13],
            m[2], m[6], m[10], m[14]);

    double diag = bbox().diagonal_length() * .6;

    const double dx = diag;
    const double dy = diag;
    const double z (0);

    int nb_quad = grid_size_ - 1;

    for(int i=0 ; i<grid_size_ ; ++i)
    {
        double x = -diag/2. + double(i)/double(nb_quad) * dx;

        for(int j=0 ; j<grid_size_ ; ++j)
        {
            double y = -diag/2. + double(j)/double(nb_quad) * dy;

            Point_3 query = t( Point_3(x, y, z) );
            double v = function_->operator()(query.x(), query.y(), query.z());

            implicit_grid_[i][j] = Point_value(Point(query.x(),query.y(),query.z()),v);
        }
    }

    // Update display list
    const_cast<Scene_implicit_function_item*>(this)->changed();

}

void
Scene_implicit_function_item::
compute_min_max()
{
    if(function_->get_min_max(min_value_, max_value_))
        return;

    double probes_nb = double(grid_size_) / 2;

    // Probe bounding box
    const Bbox& b = bbox();

    for ( int i = 0 ; i <= probes_nb ; ++i )
    {
        double x = b.xmin + double(i) * (b.xmax - b.xmin) / probes_nb;

        for ( int j = 0 ; j <= probes_nb ; ++j )
        {
            double y = b.ymin + double(j) * (b.ymax - b.ymin) / probes_nb;

            for ( int k = 0 ; k <= probes_nb ; ++k )
            {
                double z = b.zmin + double(k) * (b.zmax - b.zmin) / probes_nb;

                double v = (*function_)(x,y,z);
                if(is_nan(v)) continue;
                max_value_ = (std::max)(v,max_value_);
                min_value_ = (std::min)(v,min_value_);
            }
        }
    }
}

void
Scene_implicit_function_item::changed()
{
    Scene_item_with_display_list::changed();
    compute_vertices_and_texmap();
    initialize_buffers();
}

void Scene_implicit_function_item::contextual_changed()
{
    if(!frame_->isManipulated()) {
        if(need_update_) {
            compute_function_grid();
            compute_vertices_and_texmap();
            need_update_ = false;
        }
    }
}

#include "Scene_implicit_function_item.moc"

