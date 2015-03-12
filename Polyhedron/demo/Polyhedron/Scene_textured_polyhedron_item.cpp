#include "Scene_textured_polyhedron_item.h"
#include "Textured_polyhedron_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <QObject>
#include <CGAL/gl_render.h>

typedef EPIC_kernel::Point_3 Point;
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

void Scene_textured_polyhedron_item::initialize_buffers()
{
    glBindVertexArray(vao[0]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_facets.size())*sizeof(float),
                 positions_facets.data(),
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
                 (positions_lines.size())*sizeof(float),
                 positions_lines.data(),
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
                 (normals.size())*sizeof(float),
                 normals.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(2,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[3]);
    glBufferData(GL_ARRAY_BUFFER,
                 (textures_map_facets.size())*sizeof(float),
                 textures_map_facets.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(3,
                          2,
                          GL_FLOAT,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[4]);
    glBufferData(GL_ARRAY_BUFFER,
                 (textures_map_lines.size())*sizeof(float),
                 textures_map_lines.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(4,
                          2,
                          GL_FLOAT,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(4);


    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGB,
                 texture.GetWidth(),
                 texture.GetHeight(),
                 0,
                 GL_RGB,
                 GL_UNSIGNED_BYTE,
                 texture.GetData());
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    // Clean-up
    glBindVertexArray(0);
}

void Scene_textured_polyhedron_item::compile_shaders(void)
{
    //fill the vertex shader
    static const GLchar* vertex_shader_source[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec4 positions_facets; \n"
        "layout (location = 2) in vec3 vNormals; \n"
        "layout (location = 3) in vec2 v_texCoord; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 mv_matrix; \n"
        "uniform int is_two_side; \n"
        "uniform vec3 light_pos;  \n"
        "uniform vec3 light_diff; \n"
        "uniform vec3 light_spec; \n"
        "uniform vec3 light_amb;  \n"
        "uniform vec3 color_facets; \n"
        "float spec_power = 128.0; \n"
        "out highp vec3 fColors; \n"
        "out highp vec2 f_texCoord; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   vec4 P = mv_matrix * positions_facets; \n"
        "   vec3 N = mat3(mv_matrix)* vNormals; \n"
        "   vec3 L = light_pos - P.xyz; \n"
        "   N = normalize(N); \n"
        "   L = normalize(L); \n"
        "   vec3 diffuse; \n"
        "   if(is_two_side == 1) \n"
        "       diffuse = abs(dot(N,L)) * light_diff; \n"
        "   else \n"
        "       diffuse = max(dot(N,L), 0.0) * light_diff; \n"
        "   f_texCoord = v_texCoord; \n"
        "   fColors = color_facets * (light_amb + diffuse); \n"
        "   gl_Position =  mvp_matrix *positions_facets; \n"
        "} \n"
    };
    //fill the fragment shader
    static const GLchar* fragment_shader_source[]=
    {
        "#version 300 es \n"
        " \n"
        "in highp vec3 fColors; \n"
        "in highp vec2 f_texCoord; \n"
        "uniform sampler2D s_texture; \n"
        "out highp vec3 color; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        " color = vec3(texture(s_texture, f_texCoord)) * fColors; \n"
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
    //For the edges
    static const GLchar* vertex_shader_source_lines[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 1) in vec4 positions_lines; \n"
        "layout (location = 4) in vec2 v_texCoord; \n"
        "uniform vec3 color_lines; \n"
        "uniform mat4 mvp_matrix; \n"
        "out highp vec3 fColors; \n"
        "out highp vec2 f_texCoord; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   f_texCoord = v_texCoord; \n"
        "   fColors = color_lines; \n"
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
    rendering_program_lines = program;
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

    //Clean-up
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);


}

void Scene_textured_polyhedron_item::uniform_attrib(Viewer_interface* viewer, int mode) const
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

    //set the colors
    GLfloat colors_facet[3];
    GLfloat colors_lines[3];
    QColor temp = this->color();
    if(is_selected)
    {
        colors_facet[0] = temp.lighter(120).redF();
        colors_facet[1] = temp.lighter(120).greenF();
        colors_facet[2] = temp.lighter(120).blueF();

        colors_lines[0] = 0.0;
        colors_lines[1] = 0.0;
        colors_lines[2] = 0.0;
    }
    else
    {
        colors_facet[0] = temp.redF();
        colors_facet[1] = temp.greenF();
        colors_facet[2] = temp.blueF();

        colors_lines[0] = temp.lighter(50).redF();
        colors_lines[1] = temp.lighter(50).greenF();
        colors_lines[2] = temp.lighter(50).blueF();
    }
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
        glUniform3fv(location[8], 1, colors_facet);
        glUniform1i(sampler_location, 0);

    }
    else if(mode ==1)
    {
        glUseProgram(rendering_program_lines);
        glUniformMatrix4fv(location[7], 1, GL_FALSE, mvp_mat);
        glUniform3fv(location[9], 1, colors_lines);
    }
}
void
Scene_textured_polyhedron_item::compute_normals_and_vertices(void)
{
    positions_facets.clear();
    positions_lines.clear();
    textures_map_facets.clear();
    textures_map_lines.clear();
    normals.clear();

    typedef typename ::EPIC_kernel Kernel;
    typedef typename CGAL::Textured_items Items;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename CGAL::Polyhedron_3<Kernel,Items> Base;

    typedef typename Base::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
    typedef typename Base::Edge_iterator Edge_iterator;
    typedef typename Base::Facet Facet;
    typedef typename Base::Facet_iterator Facet_iterator;

    //Facets
    Facet_iterator f = poly->facets_begin();

    for(f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {

        Halfedge_around_facet_circulator he = f->facet_begin();
        Halfedge_around_facet_circulator end = he;
        CGAL_For_all(he,end)
        {

            // If Flat shading:1 normal per polygon added once per vertex
            if (cur_shading == GL_FLAT)
            {

                Vector n = compute_facet_normal<Facet,Kernel>(*f);
                normals.push_back(n[0]);
                normals.push_back(n[1]);
                normals.push_back(n[2]);
            }

            // If Gouraud shading: 1 normal per vertex
            else if(cur_shading == GL_SMOOTH)
            {

                const typename Facet::Normal_3& n = he->vertex()->normal();
                normals.push_back(n[0]);
                normals.push_back(n[1]);
                normals.push_back(n[2]);
            }

            //position
            const Point& p = he->vertex()->point();
            positions_facets.push_back(p.x());
            positions_facets.push_back(p.y());
            positions_facets.push_back(p.z());
            positions_facets.push_back(1.0);

            const double u = he->vertex()->u();
            const double v = he->vertex()->v();
            textures_map_facets.push_back(u);
            textures_map_facets.push_back(v);
        }


    }
    //Lines
    typedef Kernel::Point_3		Point;
    typedef Base::Edge_iterator	Edge_iterator;

    Edge_iterator he;

    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {

        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        positions_lines.push_back(a.x());
        positions_lines.push_back(a.y());
        positions_lines.push_back(a.z());
        positions_lines.push_back(1.0);

        const double u = he->vertex()->u();
        const double v = he->vertex()->v();
        textures_map_lines.push_back(u);
        textures_map_lines.push_back(v);

        positions_lines.push_back(b.x());
        positions_lines.push_back(b.y());
        positions_lines.push_back(b.z());
        positions_lines.push_back(1.0);

        const double ou = he->opposite()->vertex()->u();
        const double ov = he->opposite()->vertex()->v();
        textures_map_lines.push_back(ou);
        textures_map_lines.push_back(ov);

    }


    location[0] = glGetUniformLocation(rendering_program_facets, "mvp_matrix");
    location[1] = glGetUniformLocation(rendering_program_facets, "mv_matrix");
    location[2] = glGetUniformLocation(rendering_program_facets, "light_pos");
    location[3] = glGetUniformLocation(rendering_program_facets, "light_diff");
    location[4] = glGetUniformLocation(rendering_program_facets, "light_spec");
    location[5] = glGetUniformLocation(rendering_program_facets, "light_amb");
    location[6] = glGetUniformLocation(rendering_program_facets, "is_two_side");
    location[8] = glGetUniformLocation(rendering_program_facets, "color_facets");

    sampler_location = glGetUniformLocation(rendering_program_facets, "s_texture");

    location[7] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
    location[9] = glGetUniformLocation(rendering_program_lines, "color_lines");

}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item()
    : Scene_item(),positions_lines(0),positions_facets(0),normals(0),textures_map_facets(0),
      textures_map_lines(0),poly(new Textured_polyhedron)
{
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    cur_shading=GL_FLAT;
    is_selected=false;
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(5, buffer);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &textureId);

    compile_shaders();
    changed();
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(Textured_polyhedron* const p)
    : Scene_item(),smooth_shading(true),positions_lines(0),positions_facets(0),textures_map_facets(0),
      textures_map_lines(0), poly(p)
{
    cur_shading=GL_FLAT;
    is_selected=false;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(5, buffer);

    compile_shaders();
    changed();
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(const Textured_polyhedron& p)
    : Scene_item(),smooth_shading(true),positions_lines(0),positions_facets(0),textures_map_facets(0),
      textures_map_lines(0), poly(new Textured_polyhedron(p))
{
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    cur_shading=GL_FLAT;
    is_selected=false;
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(5, buffer);

    compile_shaders();
    changed();
}

// Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(const Scene_textured_polyhedron_item& item)
//   : Scene_item_with_display_list(item),
//     poly(new Textured_polyhedron(*item.poly))
// {
// }

Scene_textured_polyhedron_item::~Scene_textured_polyhedron_item()
{
    glDeleteBuffers(5, buffer);
    glDeleteVertexArrays(1, vao);
    glDeleteProgram(rendering_program_lines);
    glDeleteProgram(rendering_program_facets);
    delete poly;
}

Scene_textured_polyhedron_item* 
Scene_textured_polyhedron_item::clone() const {
    return new Scene_textured_polyhedron_item(*poly);
}

// Load textured_polyhedron from .OFF file
bool
Scene_textured_polyhedron_item::load(std::istream& in)
{
    std::cout<<"LOAD"<<std::endl;
    in >> *poly;
    changed();
    return in && !isEmpty();
}

// Write textured_polyhedron to .OFF file
bool 
Scene_textured_polyhedron_item::save(std::ostream& out) const
{
    out << *poly;
    return (bool) out;
}

QString 
Scene_textured_polyhedron_item::toolTip() const
{
    if(!poly)
        return QString();

    return QObject::tr("<p>Textured polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                       "Number of facets: %4</p>")
            .arg(this->name())
            .arg(poly->size_of_vertices())
            .arg(poly->size_of_halfedges()/2)
            .arg(poly->size_of_facets())
            .arg(this->renderingModeName())
            .arg(this->color().name());
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_textured_polyhedron_item::draw(Viewer_interface* viewer) const {

    glBindVertexArray(vao[0]);
    glUseProgram(rendering_program_facets);
    uniform_attrib(viewer,0);
    glDrawArrays(GL_TRIANGLES, 0, positions_facets.size()/4);
    //Clean-up
    glUseProgram(0);
    glBindVertexArray(0);

}
void Scene_textured_polyhedron_item::draw_edges(Viewer_interface* viewer) const {
    glBindVertexArray(vao[0]);

    glUseProgram(rendering_program_lines);
    uniform_attrib(viewer,1);
    glDrawArrays(GL_LINES, 0, positions_lines.size()/4);
    // Clean-up
    glUseProgram(0);
    glBindVertexArray(0);
}

Textured_polyhedron* 
Scene_textured_polyhedron_item::textured_polyhedron()       { return poly; }
const Textured_polyhedron* 
Scene_textured_polyhedron_item::textured_polyhedron() const { return poly; }

bool
Scene_textured_polyhedron_item::isEmpty() const {
    return (poly == 0) || poly->empty();
}

Scene_textured_polyhedron_item::Bbox
Scene_textured_polyhedron_item::bbox() const {
    const Point& p = *(poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Textured_polyhedron::Point_iterator it = poly->points_begin();
        it != poly->points_end();
        ++it) {
        bbox = bbox + it->bbox();
    }
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}
void
Scene_textured_polyhedron_item::changed()
{
    compute_normals_and_vertices();
    initialize_buffers();
}
void
Scene_textured_polyhedron_item::
contextual_changed()
{
    GLint new_shading;
    glGetIntegerv(GL_SHADE_MODEL, &new_shading);
    prev_shading = cur_shading;
    cur_shading = new_shading;
    if(prev_shading != cur_shading)
    {
        changed();
    }
}
void
Scene_textured_polyhedron_item::selection_changed(bool p_is_selected)
{
    if(p_is_selected != is_selected)
    {
        is_selected = p_is_selected;
        initialize_buffers();
    }
    else
        is_selected = p_is_selected;
}
#include "Scene_textured_polyhedron_item.moc"
