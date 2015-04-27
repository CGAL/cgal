#include "Scene_plane_item.h"

#include "Scene_plane_item.moc"


void Scene_plane_item::initialize_buffers()
{
    qFunc.glBindVertexArray(vao[0]);

    qFunc.glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    qFunc.glBufferData(GL_ARRAY_BUFFER,
                 (positions_quad.size())*sizeof(float),
                 positions_quad.data(),
                 GL_STATIC_DRAW);
    qFunc.glVertexAttribPointer(0, //number of the buffer
                          3, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    qFunc.glEnableVertexAttribArray(0);

    qFunc.glBindBuffer(GL_ARRAY_BUFFER, buffer[1]);
    qFunc.glBufferData(GL_ARRAY_BUFFER,
                 (positions_lines.size())*sizeof(float),
                 positions_lines.data(),
                 GL_STATIC_DRAW);
    qFunc.glVertexAttribPointer(1, //number of the buffer
                          3, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    qFunc.glEnableVertexAttribArray(1);
}
void Scene_plane_item::compile_shaders(void)
{
    //fill the vertex shader
    static const GLchar* vertex_shader_source[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec3 positions; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 f_matrix; \n"
        "uniform vec3 color; \n"
        "out highp vec3 fColors; \n"

        "vec4 positions_quad = vec4 (positions, 1.0); \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   fColors = color; \n"
        "   gl_Position =  mvp_matrix * f_matrix * positions_quad; \n"
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

    GLuint vertex_shader = qFunc.glCreateShader(GL_VERTEX_SHADER);
    qFunc.glShaderSource(vertex_shader, 1, vertex_shader_source, NULL);
    qFunc.glCompileShader(vertex_shader);


    GLuint fragment_shader =	qFunc.glCreateShader(GL_FRAGMENT_SHADER);
    qFunc.glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    qFunc.glCompileShader(fragment_shader);


    //creates the program, attaches and links the shaders
    GLuint program= qFunc.glCreateProgram();
    qFunc.glAttachShader(program, vertex_shader);
    qFunc.glAttachShader(program, fragment_shader);
    qFunc.glLinkProgram(program);

    //Clean-up
    qFunc.glDeleteShader(vertex_shader);

    rendering_program_quad = program;


    //For the edges
    static const GLchar* vertex_shader_source_lines[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 1) in vec3 positions; \n"
        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 f_matrix; \n"
        "vec4 positions_lines = vec4(positions,1.0); \n"
        "out highp vec3 fColors; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   fColors = vec3(0.0,0.0,0.0); \n"
        "   gl_Position = mvp_matrix * f_matrix * positions_lines; \n"
        "} \n"
    };

    vertex_shader = qFunc.glCreateShader(GL_VERTEX_SHADER);
    qFunc.glShaderSource(vertex_shader, 1, vertex_shader_source_lines, NULL);
    qFunc.glCompileShader(vertex_shader);


    qFunc.glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    qFunc.glCompileShader(fragment_shader);


    program = qFunc.glCreateProgram();
    qFunc.glAttachShader(program, vertex_shader);
    qFunc.glAttachShader(program, fragment_shader);
    qFunc.glLinkProgram(program);


    //Clean-up
    qFunc.glDeleteShader(vertex_shader);
    qFunc.glDeleteShader(fragment_shader);
    rendering_program_lines = program;

}
void Scene_plane_item::uniform_attrib(Viewer_interface* viewer, int mode) const
{
    GLfloat mvp_mat[16];
    GLdouble d_mat[16];
    GLfloat f_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i){
        mvp_mat[i] = GLfloat(d_mat[i]);
        f_mat[i] = GLfloat(frame->matrix()[i]);
    }

    GLfloat colors[3];

    colors[0] =this->color().redF();
    colors[1] =this->color().greenF();
    colors[2] =this->color().blueF();

    if(mode ==0)
    {
        qFunc.glUseProgram(rendering_program_quad);
        qFunc.glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
        qFunc.glUniform3fv(location[1], 1, colors);
        qFunc.glUniformMatrix4fv(location[2], 1, GL_FALSE, f_mat);

    }
    else if(mode ==1)
    {
        qFunc.glUseProgram(rendering_program_lines);
        qFunc.glUniformMatrix4fv(location[3], 1, GL_FALSE, mvp_mat);
        qFunc.glUniformMatrix4fv(location[4], 1, GL_FALSE, f_mat);
    }

}
void Scene_plane_item::compute_normals_and_vertices(void)
{
    positions_quad.clear();
    positions_lines.clear();

    const double diag = scene_diag();
    //The quad
    {

    positions_quad.push_back(-diag);
    positions_quad.push_back(-diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(-diag);
    positions_quad.push_back(diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(diag);
    positions_quad.push_back(-diag);
    positions_quad.push_back(0.0);

    positions_quad.push_back(-diag);
    positions_quad.push_back(diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(diag);
    positions_quad.push_back(-diag);
    positions_quad.push_back(0.0);
    positions_quad.push_back(diag);
    positions_quad.push_back(diag);
    positions_quad.push_back(0.0);

}
    //The grid
    float x = (2*diag)/10.0;
    float y = (2*diag)/10.0;
    {
        for(int u = 0; u < 11; u++)
        {

            positions_lines.push_back(-diag + x* u);
            positions_lines.push_back(-diag);
            positions_lines.push_back(0.0);

            positions_lines.push_back(-diag + x* u);
            positions_lines.push_back(diag);
            positions_lines.push_back(0.0);
        }
        for(int v=0; v<11; v++)
        {

            positions_lines.push_back(-diag);
            positions_lines.push_back(-diag + v * y);
            positions_lines.push_back(0.0);

            positions_lines.push_back(diag);
            positions_lines.push_back(-diag + v * y);
            positions_lines.push_back(0.0);
        }

    }


    location[0] = qFunc.glGetUniformLocation(rendering_program_quad, "mvp_matrix");
    location[1] = qFunc.glGetUniformLocation(rendering_program_quad, "color");
    location[2] = qFunc.glGetUniformLocation(rendering_program_quad, "f_matrix");
    location[3] = qFunc.glGetUniformLocation(rendering_program_lines, "mvp_matrix");
    location[4] = qFunc.glGetUniformLocation(rendering_program_lines, "f_matrix");
}

void Scene_plane_item::draw(Viewer_interface* viewer)const
{
    qFunc.glBindVertexArray(vao[0]);
    qFunc.glUseProgram(rendering_program_quad);
    uniform_attrib(viewer,0);
    qFunc.glDrawArrays(GL_TRIANGLES, 0, positions_quad.size()/3);\
    qFunc.glUseProgram(0);
    qFunc.glBindVertexArray(0);

}

void Scene_plane_item::draw_edges(Viewer_interface* viewer)const
{
    qFunc.glBindVertexArray(vao[0]);
    qFunc.glUseProgram(rendering_program_lines);
    uniform_attrib(viewer,1);
    qFunc.glDrawArrays(GL_LINES, 0, positions_lines.size()/3);\
    qFunc.glUseProgram(0);
    qFunc.glBindVertexArray(0);
}
