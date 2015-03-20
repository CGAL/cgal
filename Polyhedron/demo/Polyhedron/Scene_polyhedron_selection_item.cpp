#include "Scene_polyhedron_selection_item.h"
#include "Scene_polyhedron_selection_item.moc"
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
    GLfloat spec_power;

};

void Scene_polyhedron_selection_item::initialize_buffers()
{
    glBindVertexArray(vao[0]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_facets.size())*sizeof(float),
                 positions_facets.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0, //number of the buffer
                          3, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[1]);
    glBufferData(GL_ARRAY_BUFFER,
                 (normals.size())*sizeof(float),
                 normals.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(1, //number of the buffer
                          3, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[2]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_lines.size())*sizeof(float),
                 positions_lines.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(2, //number of the buffer
                          3, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[3]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_points.size())*sizeof(float),
                 positions_points.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(3, //number of the buffer
                          3, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(3);

    glBindVertexArray(0);

}
void Scene_polyhedron_selection_item::compile_shaders()
{
    //The facets
    //fill the vertex shader
    static const GLchar* vertex_shader_source_facets[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec3 positions; \n"
        "layout (location = 1) in vec3 vNormals; \n"
        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 mv_matrix; \n"

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
        "vec4 P = mv_matrix * positions_facets; \n"
        "vec3 N = mat3(mv_matrix)* vNormals; \n"
        "vec3 L = light_pos - P.xyz; \n"
        "vec3 V = -P.xyz; \n"

        "N = normalize(N); \n"
        "L = normalize(L); \n"
        "V = normalize(V); \n"

        "vec3 R = reflect(-L, N); \n"
        "  vec3 diffuse; \n"
        "if(is_two_side == 1) \n"
        "   diffuse = abs(dot(N,L)) * light_diff; \n"
        "else \n"
        "   diffuse = max(dot(N,L), 0.0) * light_diff; \n"
        "vec3 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

        "fColors = light_amb + diffuse + specular ; \n"
        "gl_Position = mvp_matrix * positions_facets; \n"
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
    glShaderSource(vertex_shader, 1, vertex_shader_source_facets, NULL);
    glCompileShader(vertex_shader);
    GLuint fragment_shader =	glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    GLuint program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    glDeleteShader(vertex_shader);
    rendering_program_facets = program;

    //The lines
    //fill the vertex shader
    static const GLchar* vertex_shader_source_lines[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 2) in vec3 positions; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform vec3 color; \n"
        "out highp vec3 fColors; \n"
        "vec4 positions_lines = vec4(positions, 1.0); \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color; \n"
        "   gl_Position = mvp_matrix * positions_lines; \n"
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

    glDeleteShader(vertex_shader);
    rendering_program_lines = program;

    //The points
    //fill the vertex shader
    static const GLchar* vertex_shader_source_points[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 3) in vec3 positions; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform vec3 color; \n"
        "out highp vec3 fColors; \n"
        "vec4 positions_points = vec4(positions, 1.0); \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color; \n"
        "   gl_Position = mvp_matrix * positions_points; \n"
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

    glDeleteShader(vertex_shader);
    rendering_program_points = program;


    glDeleteShader(fragment_shader);

}
void Scene_polyhedron_selection_item::uniform_attrib(Viewer_interface* viewer, int mode) const
{
    GLfloat colors[3];
    light_info light;
    GLint is_both_sides = 0;
    GLfloat mvp_mat[16];
    GLfloat mv_mat[16];

    //fills the MVP and MV matrices.

    GLdouble d_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i)
        mvp_mat[i] = GLfloat(d_mat[i]);

    viewer->camera()->getModelViewMatrix(d_mat);
    for (int i=0; i<16; ++i)
        mv_mat[i] = GLfloat(d_mat[i]);


    glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE, &is_both_sides);

    //fills the arraw of colors with the current color
    colors[0] = facet_color.redF();
    colors[1] = facet_color.greenF();
    colors[2] = facet_color.blueF();
    //Gets lighting info :

    //position
    glGetLightfv(GL_LIGHT0, GL_POSITION, light.position);

    //ambient
    glGetLightfv(GL_LIGHT0, GL_AMBIENT, light.ambient);
    light.ambient[0]*=colors[0];
    light.ambient[1]*=colors[1];
    light.ambient[2]*=colors[2];

    //specular
    glGetLightfv(GL_LIGHT0, GL_SPECULAR, light.specular);

    //diffuse
    glGetLightfv(GL_LIGHT0, GL_DIFFUSE, light.diffuse);

    light.diffuse[0]*=colors[0];
    light.diffuse[1]*=colors[1];
    light.diffuse[2]*=colors[2];

    //For the Flat mode
    if(mode ==0)
    {
        glUseProgram(rendering_program_facets);
        glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
        glUniformMatrix4fv(location[1], 1, GL_FALSE, mv_mat);
        glUniform3fv(location[2], 1, light.position);
        glUniform3fv(location[3], 1, light.diffuse);
        glUniform3fv(location[4], 1, light.specular);
        glUniform3fv(location[5], 1, light.ambient);
        glUniform1i(location[6], is_both_sides);


    }
    //For the Wire mode
    else if(mode ==1)
    {
        glUseProgram(rendering_program_lines);
        glUniformMatrix4fv(location[7], 1, GL_FALSE, mvp_mat);
        colors[0] = edge_color.redF();
        colors[1] = edge_color.greenF();
        colors[2] = edge_color.blueF();
        glUniform3fv(location[8],1,colors);
    }
    //For the points mode
    else if(mode ==2)
    {
        glUseProgram(rendering_program_points);
        glUniformMatrix4fv(location[9], 1, GL_FALSE, mvp_mat);
        colors[0] = vertex_color.redF();
        colors[1] = vertex_color.greenF();
        colors[2] = vertex_color.blueF();
        glUniform3fv(location[10],1,colors);
    }
}
void Scene_polyhedron_selection_item::compute_elements()
{
    positions_facets.clear();
    positions_lines.clear();
    positions_points.clear();
    normals.clear();
    //The facets
    {


        for(Selection_set_facet::iterator
            it = selected_facets.begin(),
            end = selected_facets.end();
            it != end; ++it)
        {
            const Kernel::Vector_3 n =
                    compute_facet_normal<Polyhedron::Facet,Kernel>(**it);

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
                    he = (*it)->facet_begin(),
                    cend = he;

            CGAL_For_all(he,cend)
            {
                const Kernel::Point_3& p = he->vertex()->point();
                positions_facets.push_back(p.x());
                positions_facets.push_back(p.y());
                positions_facets.push_back(p.z());
            }
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

    location[0] = glGetUniformLocation(rendering_program_facets, "mvp_matrix");
    location[1] = glGetUniformLocation(rendering_program_facets, "mv_matrix");
    location[2] = glGetUniformLocation(rendering_program_facets, "light_pos");
    location[3] = glGetUniformLocation(rendering_program_facets, "light_diff");
    location[4] = glGetUniformLocation(rendering_program_facets, "light_spec");
    location[5] = glGetUniformLocation(rendering_program_facets, "light_amb");
    location[6] = glGetUniformLocation(rendering_program_facets, "is_two_side");

    location[7] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
    location[8] = glGetUniformLocation(rendering_program_lines, "color");

    location[9] = glGetUniformLocation(rendering_program_points, "mvp_matrix");
    location[10] = glGetUniformLocation(rendering_program_points, "color");
}

void Scene_polyhedron_selection_item::draw(Viewer_interface* viewer) const
{
    draw();

    draw_points(viewer);
    GLfloat offset_factor;
    GLfloat offset_units;
    glGetFloatv( GL_POLYGON_OFFSET_FACTOR, &offset_factor);
    glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
    glPolygonOffset(-1.f, 1.f);
    glBindVertexArray(vao[0]);
    uniform_attrib(viewer,0);
    glUseProgram(rendering_program_facets);
    glDrawArrays(GL_TRIANGLES, 0, positions_facets.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);
    glPolygonOffset(offset_factor, offset_units);
    draw_edges(viewer);


}

void Scene_polyhedron_selection_item::draw_edges(Viewer_interface* viewer) const
{

    glLineWidth(3.f);
    glBindVertexArray(vao[0]);
    uniform_attrib(viewer,1);
    glUseProgram(rendering_program_lines);
    glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);
    glLineWidth(1.f);
}

void Scene_polyhedron_selection_item::draw_points(Viewer_interface* viewer) const
{
    std::cout<<"draw_points"<<std::endl;
    glPointSize(5.f);
    glBindVertexArray(vao[0]);
    uniform_attrib(viewer,2);
    glUseProgram(rendering_program_points);
    glDrawArrays(GL_POINTS, 0, positions_points.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);
    glPointSize(1.f);

}
