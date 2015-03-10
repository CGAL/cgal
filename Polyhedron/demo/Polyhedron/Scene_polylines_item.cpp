#include "Scene_polylines_item.h"

#include <CGAL/bounding_box.h>
#include <CGAL/gl.h>
#include <CGAL/glu.h>
#include <QMenu>
#include <QAction>

#include <QInputDialog>
namespace {
void CGALglcolor(QColor c, int dv = 0)
{
    if ( 0 != dv )
    {
        // workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
        c = c.darker(dv);
#undef darker
    }
    ::glColor4d(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
}
}
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
typedef Scene_polylines_item::K K;
typedef K::Point_3 Point_3;
//Fill the VBO with coordinates of the vertices composing a sphere
void Scene_polylines_item::create_Sphere(double R)
{

    float T, P;
    float x[4],y[4],z[4];


    //Top of the sphere
    for(int t=0; t<360; t+=sectors)
    {

        positions_spheres.push_back(0);
        positions_spheres.push_back(0);
        positions_spheres.push_back(R);
        positions_spheres.push_back(1.0);


        normals_spheres.push_back(0);
        normals_spheres.push_back(0);
        normals_spheres.push_back(1);



        P = rings*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        positions_spheres.push_back(R * x[1]);
        positions_spheres.push_back(R * y[1]);
        positions_spheres.push_back(R * z[1]);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(x[1]);
        normals_spheres.push_back(y[1]);
        normals_spheres.push_back(z[1]);

        //
        P = rings*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        positions_spheres.push_back(R * x[2]);
        positions_spheres.push_back(R * y[2]);
        positions_spheres.push_back(R * z[2]);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(x[2]);
        normals_spheres.push_back(y[2]);
        normals_spheres.push_back(z[2]);

    }

    //Body of the sphere
    for (int p=rings; p<180-rings; p+=rings)
        for(int t=0; t<360; t+=sectors)
        {
            //A
            P = p*M_PI/180.0;
            T = t*M_PI/180.0;
            x[0] = sin(P) * cos(T) ;
            y[0] = sin(P) * sin(T) ;
            z[0] = cos(P);

            positions_spheres.push_back(R * x[0]);
            positions_spheres.push_back(R * y[0]);
            positions_spheres.push_back(R * z[0]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[0]);
            normals_spheres.push_back(y[0]);
            normals_spheres.push_back(z[0]);

            //B
            P = (p+rings)*M_PI/180.0;
            T = t*M_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            positions_spheres.push_back(R * x[1]);
            positions_spheres.push_back(R * y[1]);
            positions_spheres.push_back(R * z[1]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[1]);
            normals_spheres.push_back(y[1]);
            normals_spheres.push_back(z[1]);

            //C
            P = p*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            positions_spheres.push_back(R * x[2]);
            positions_spheres.push_back(R * y[2]);
            positions_spheres.push_back(R * z[2]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[2]);
            normals_spheres.push_back(y[2]);
            normals_spheres.push_back(z[2]);
            //D
            P = (p+rings)*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[3] = sin(P) * cos(T) ;
            y[3] = sin(P) * sin(T) ;
            z[3] = cos(P);
            positions_spheres.push_back(R * x[3]);
            positions_spheres.push_back(R * y[3]);
            positions_spheres.push_back(R * z[3]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[3]);
            normals_spheres.push_back(y[3]);
            normals_spheres.push_back(z[3]);



            positions_spheres.push_back(R * x[1]);
            positions_spheres.push_back(R * y[1]);
            positions_spheres.push_back(R * z[1]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[1]);
            normals_spheres.push_back(y[1]);
            normals_spheres.push_back(z[1]);

            positions_spheres.push_back(R * x[2]);
            positions_spheres.push_back(R * y[2]);
            positions_spheres.push_back(R * z[2]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[2]);
            normals_spheres.push_back(y[2]);
            normals_spheres.push_back(z[2]);

        }
    //Bottom of the sphere
    for(int t=0; t<360; t+=sectors)
    {


        positions_spheres.push_back(0);
        positions_spheres.push_back(0);
        positions_spheres.push_back(-R);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(0);
        normals_spheres.push_back(0);
        normals_spheres.push_back(-1);


        P = (180-rings)*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        positions_spheres.push_back(R * x[1]);
        positions_spheres.push_back(R * y[1]);
        positions_spheres.push_back(R * z[1]);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(x[1]);
        normals_spheres.push_back(y[1]);
        normals_spheres.push_back(z[1]);


        P = (180-rings)*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        positions_spheres.push_back(R * x[2]);
        positions_spheres.push_back(R * y[2]);
        positions_spheres.push_back(R * z[2]);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(x[2]);
        normals_spheres.push_back(y[2]);
        normals_spheres.push_back(z[2]);

    }
}

class Scene_polylines_item_private {
public:
    typedef Scene_polylines_item::K K;
    typedef K::Point_3 Point_3;

    Scene_polylines_item_private() :
        draw_extremities(false),
        spheres_drawn_radius(0),
        sphere_display_list(0),
        quadric(0)
    {}

    ~Scene_polylines_item_private()
    {
        if(quadric != 0)
            gluDeleteQuadric(quadric);
        if(sphere_display_list  != 0)
            glDeleteLists(sphere_display_list, 1);
    }

    void draw_sphere(const K::Point_3&, double) const;
    void draw_spheres(const Scene_polylines_item*) const;

    bool draw_extremities;
    double spheres_drawn_radius;
private:
    mutable GLuint sphere_display_list;
    mutable GLUquadric* quadric;
};

void
Scene_polylines_item::initialize_buffers()
{

    glBindVertexArray(vao[0]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_lines.size())*sizeof(float),
                 positions_lines.data(),
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
                 (positions_spheres.size())*sizeof(float),
                 positions_spheres.data(),
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
                 (normals_spheres.size())*sizeof(float),
                 normals_spheres.data(), GL_STATIC_DRAW);
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
                 (color_spheres.size())*sizeof(float),
                 color_spheres.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(3,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[4]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_center.size())*sizeof(float),
                 positions_center.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(4,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(4);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[5]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_wire_spheres.size())*sizeof(float),
                 positions_wire_spheres.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(5, //number of the buffer
                          4, //number of floats to be taken
                          GL_FLOAT, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(5);

    glVertexAttribDivisor(3, 1);
    glVertexAttribDivisor(4, 1);

    // Clean-up
    glBindVertexArray(0);


}

void
Scene_polylines_item::compile_shaders()
{
    //fill the vertex shader
    static const GLchar* vertex_shader_source[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 1) in vec4 positions_spheres; \n"
        "layout (location = 2) in vec3 vNormals; \n"
        "layout (location = 3) in vec3 color_spheres; \n"
        "layout (location = 4) in vec3 center; \n"
        " \n"
        "uniform mat4 mvp_matrix; \n"
        "uniform mat4 mv_matrix; \n"
        " \n"
        "uniform int is_two_side; \n"
        "uniform vec3 light_pos;  \n"
        "uniform vec3 light_diff; \n"
        "uniform vec3 light_spec; \n"
        "uniform vec3 light_amb;  \n"
        "float spec_power = 128.0; \n"
        " \n"
        "out highp vec3 fColors; \n"
        " \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   vec4 P = mv_matrix * positions_spheres; \n"
        "   vec3 N = mat3(mv_matrix)* vNormals; \n"
        "   vec3 L = light_pos - P.xyz; \n"
        "   vec3 V = -P.xyz; \n"
        " \n"
        "   N = normalize(N); \n"
        "   L = normalize(L); \n"
        "   V = normalize(V); \n"
        " \n"
        "   vec3 R = reflect(-L, N); \n"
        "   vec3 diffuse; \n"
        "   if(is_two_side == 1) \n"
        "       diffuse = abs(dot(N,L)) * light_diff * color_spheres; \n"
        "   else \n"
        "       diffuse = max(dot(N,L), 0.0) * light_diff * color_spheres; \n"
        "   vec3 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"
        " \n"
        "   fColors = light_amb*color_spheres + diffuse + specular ; \n"
        "   gl_Position =  mvp_matrix * vec4(positions_spheres.x + center.x, positions_spheres.y + center.y, positions_spheres.z + center.z, 1.0) ; \n"
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

    //creates and compiles the vertex shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source, NULL);
    glCompileShader(vertex_shader);

    //creates and compiles the fragment shader
    GLuint fragment_shader =	glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    //creates the program, attaches and links the shaders
    GLuint program= glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    //Delete the shaders which are now in the memory
    glDeleteShader(vertex_shader);

    rendering_program_spheres = program;

    //For the lines
    //fill the vertex shader
    static const GLchar* vertex_shader_source_lines[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec4 positions_lines; \n"

        "uniform mat4 mvp_matrix; \n"
        "uniform vec3 color_lines; \n"

        "out highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
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

    glDeleteShader(vertex_shader);
    rendering_program_lines = program;

    //For the wired spheres
    static const GLchar* vertex_shader_source_wire_sphere[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 5) in vec4 positions_WireSpheres; \n"
        "layout (location = 3) in vec3 color_spheres; \n"
        "layout (location = 4) in vec3 center; \n"
        " \n"
        "uniform mat4 mvp_matrix; \n"
        " \n"
        " \n"
        "out highp vec3 fColors; \n"
        " \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        "   fColors = color_spheres; \n"
        "   gl_Position =  mvp_matrix * vec4(positions_WireSpheres.x + center.x, positions_WireSpheres.y + center.y, positions_WireSpheres.z + center.z, 1.0) ; \n"
        "} \n"
    };

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_wire_sphere, NULL);
    glCompileShader(vertex_shader);

    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    rendering_program_WireSpheres = program;

}


void Scene_polylines_item::uniform_attrib(Viewer_interface* viewer, int mode) const
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

    //Program for the Flat mode
    if(mode ==0)
    {
        //Decides if the light is one or both sides
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

        glUseProgram(rendering_program_spheres);
        glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
        glUniformMatrix4fv(location[1], 1, GL_FALSE, mv_mat);
        glUniform3fv(location[2], 1, light.position);
        glUniform3fv(location[3], 1, light.diffuse);
        glUniform3fv(location[4], 1, light.specular);
        glUniform3fv(location[5], 1, light.ambient);
        glUniform1i(location[6], is_both_sides);
    }
    //For the wiremode programs
    else if(mode ==1)
    {
        //Lines
        GLfloat colors[3];
        colors[0] = this->color().redF();
        colors[1] = this->color().greenF();
        colors[2] = this->color().blueF();
        glUseProgram(rendering_program_lines);
        glUniformMatrix4fv(location[7], 1, GL_FALSE, mvp_mat);
        glUniform3fv(location[8], 1, colors);

        glUseProgram(rendering_program_WireSpheres);
        glUniformMatrix4fv(location[9], 1, GL_FALSE, mvp_mat);
        glUniform3fv(location[10], 1, colors);
    }
}

void
Scene_polylines_item::compute_elements()
{
    positions_spheres.clear();
    positions_wire_spheres.clear();
    positions_lines.clear();
    color_spheres.clear();
    normals_spheres.clear();
    positions_center.clear();
    nbSpheres = 0;

    //Fills the VBO with the lines
    for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
        it != polylines.end();
        ++it){
        if(it->empty()) continue;
        for(size_t i = 0, end = it->size()-1;
            i < end; ++i)
        {
            const Point_3& a = (*it)[i];
            const Point_3& b = (*it)[i+1];
            positions_lines.push_back(a.x());
            positions_lines.push_back(a.y());
            positions_lines.push_back(a.z());
            positions_lines.push_back(1.0);

            positions_lines.push_back(b.x());
            positions_lines.push_back(b.y());
            positions_lines.push_back(b.z());
            positions_lines.push_back(1.0);

        }

    }
    //Fills the VBO with the spheres
    if(d->draw_extremities)
    {

        // FIRST, count the number of incident cycles and polylines
        // for all extremities.
        typedef std::map<Point_3, int> Point_to_int_map;
        typedef Point_to_int_map::iterator iterator;
        Point_to_int_map corner_polyline_nb;

        { // scope to fill corner_polyline_nb'
            Point_to_int_map corner_cycles_nb;

            for(std::list<std::vector<Point_3> >::const_iterator
                it = this->polylines.begin(),
                end = this->polylines.end();
                it != end; ++it)
            {
                const K::Point_3& a = *it->begin();
                const K::Point_3& b = *it->rbegin();
                if(a == b) {
                    if ( it->size()>1 )
                        ++corner_cycles_nb[a];
                    else
                        ++corner_polyline_nb[a];
                }
                else {
                    ++corner_polyline_nb[a];
                    ++corner_polyline_nb[b];
                }
            }
            // THEN, ignore points that are incident to one cycle only.
            for(iterator
                c_it = corner_cycles_nb.begin(),
                end = corner_cycles_nb.end();
                c_it != end; ++c_it)
            {
                const Point_3& a = c_it->first;

                iterator p_it = corner_polyline_nb.find(a);

                // If the point 'a'=c_it->first has only incident cycles...
                if(p_it == corner_polyline_nb.end()) {
                    // ...then count it as a corner only if it has two incident cycles
                    // or more.
                    if(c_it->second > 1) {
                        corner_polyline_nb[a] = c_it->second;
                    }
                } else {
                    // else add the number of cycles.
                    p_it->second += c_it->second;
                }
            }
        }
        // At this point, 'corner_polyline_nb' gives the multiplicity of all
        // corners.
        //Finds the centers of the spheres and their color
        for(iterator
            p_it = corner_polyline_nb.begin(),
            end = corner_polyline_nb.end();
            p_it != end; ++p_it)
        {
            nbSpheres++;
            const K::Point_3& centre = p_it->first;
            positions_center.push_back(centre.x());
            positions_center.push_back(centre.y());
            positions_center.push_back(centre.z());

            float colors[3];
            switch(p_it->second) {
            case 1:
                colors[0] = 0.0; // black
                colors[1] = 0.0;
                colors[2] = 0.0;
                break;
            case 2:
                colors[0] = 0.0; // green
                colors[1] = 0.8;
                colors[2] = 0.0;
                break;
            case 3:
                colors[0] = 0.0; // blue
                colors[1] = 0.0;
                colors[2] = 0.8;
                break;
            case 4:
                colors[0] = 0.8; //red
                colors[1] = 0.0;
                colors[2] = 0.0;
                break;
            default:
                colors[0] = 0.8; //fuschia
                colors[1] = 0.0;
                colors[2] = 0.8;
            }

            color_spheres.push_back(colors[0]);
            color_spheres.push_back(colors[1]);
            color_spheres.push_back(colors[2]);
        }
        create_Sphere(d->spheres_drawn_radius);

        //Convert the triangle coordinates to lines coordinates for the
        //Wiremode in the spheres
        for(int i=0; i< positions_spheres.size(); i=i)
        {
            //draw triangles
            if(i< (360/sectors)*12)
            {
                //AB
                for(int j=i; j<i+8; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //BC
                for(int j=i+4; j<i+12; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //CA
                for(int j=i+8; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j%12]);
                }
                i+=12;
            }
            //draw quads
            else if((360/sectors) * 3 * 4 < i < positions_spheres.size() - (360/sectors) * 3 * 4)
            {
                //AB
                for(int j=i; j<i+8; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //BD
                for(int j=i+4; j<i+8; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                for(int j=i+12; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //DC
                for(int j=i+12; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                for(int j=i+8; j<i+12; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //CA
                for(int j=i+8; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j%12]);
                }
                i+=24;
            }
            //draw triangles
            else
            {
                //AB
                for(int j=i; j<i+8; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //BC
                for(int j=i+4; j<i+12; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //CA
                for(int j=i+8; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j%12]);
                }
                i+=12;
            }

        }
    }

    location[0] = glGetUniformLocation(rendering_program_spheres, "mvp_matrix");
    location[1] = glGetUniformLocation(rendering_program_spheres, "mv_matrix");
    location[2] = glGetUniformLocation(rendering_program_spheres, "light_pos");
    location[3] = glGetUniformLocation(rendering_program_spheres, "light_diff");
    location[4] = glGetUniformLocation(rendering_program_spheres, "light_spec");
    location[5] = glGetUniformLocation(rendering_program_spheres, "light_amb");
    location[6] = glGetUniformLocation(rendering_program_spheres, "is_two_side");

    location[7] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
    location[8] = glGetUniformLocation(rendering_program_lines, "color_lines");

    location[9] = glGetUniformLocation( rendering_program_WireSpheres, "mvp_matrix");
    location[10] = glGetUniformLocation(rendering_program_WireSpheres, "color_lines");
}


Scene_polylines_item::Scene_polylines_item() 
    : d(new Scene_polylines_item_private()),positions_lines(0), positions_spheres(0),
      normals_spheres(0), positions_center(0),color_spheres(0), positions_wire_spheres(0),nbSpheres(0),
      rings(18), sectors(36)
{
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(6, buffer);

    compile_shaders();
    changed();

}

Scene_polylines_item::~Scene_polylines_item()
{
    delete d;
}

bool
Scene_polylines_item::isEmpty() const {
    return polylines.empty();
}

Scene_interface::Bbox 
Scene_polylines_item::bbox() const {
    if(isEmpty())
        return Bbox();
    std::list<Point_3> boxes;
    for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
        it != polylines.end();
        ++it){
        if(it->begin() != it->end()) {
            Iso_cuboid_3 cub = CGAL::bounding_box(it->begin(), it->end());
            boxes.push_back((cub.min)());
            boxes.push_back((cub.max)());
        }
    }
    Iso_cuboid_3 bbox =
            boxes.begin() != boxes.end() ?
                CGAL::bounding_box(boxes.begin(), boxes.end()) :
                Iso_cuboid_3();

    return Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax());
}

Scene_polylines_item* 
Scene_polylines_item::clone() const {
    Scene_polylines_item* item = new Scene_polylines_item;
    item->polylines = polylines;
    QVariant metadata_variant = property("polylines metadata");
    if(metadata_variant.type() == QVariant::StringList)
    {
        item->setProperty("polylines metadata", metadata_variant);
    }
    return item;
}

QString
Scene_polylines_item::toolTip() const {
    QString s =
            tr("<p><b>%1</b> (mode: %2, color: %3)<br />"
               "<i>Polylines</i></p>"
               "<p>Number of polylines: %4</p>")
            .arg(this->name())
            .arg(this->renderingModeName())
            .arg(this->color().name())
            .arg(polylines.size());
    if(d->draw_extremities) {
        s += tr("<p>Legende of endpoints colors: <ul>"
                "<li>black: one incident polyline</li>"
                "<li>green: two incident polylines</li>"
                "<li>blue: three incident polylines</li>"
                "<li>red: four incident polylines</li>"
                "<li>fuchsia: five or more incident polylines</li>"
                "</ul></p>");
    }
    return s;
}

bool
Scene_polylines_item::supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe ||
            m == FlatPlusEdges ||
            m == Points);
}

// Shaded OpenGL drawing: only draw spheres
void
Scene_polylines_item::draw(Viewer_interface* viewer) const {

    if(d->draw_extremities)
    {
        glBindVertexArray(vao[0]);
        glUseProgram(rendering_program_spheres);
        uniform_attrib(viewer,0);
        glDrawArraysInstanced(GL_TRIANGLES, 0, positions_spheres.size()/4, nbSpheres);
        glUseProgram(0);
        glBindVertexArray(0);
    }
}

// Wireframe OpenGL drawing
void 
Scene_polylines_item::draw_edges(Viewer_interface* viewer) const {

    glBindVertexArray(vao[0]);
    uniform_attrib(viewer,1);
    glUseProgram(rendering_program_lines);
    glDrawArrays(GL_LINES, 0, positions_lines.size()/4);
    if(d->draw_extremities)
    {
        uniform_attrib(viewer,0);
        glUseProgram(rendering_program_WireSpheres);
        glDrawArraysInstanced(GL_LINES, 0, positions_wire_spheres.size()/4, nbSpheres);
    }
    // Clean-up
    glUseProgram(0);
    glBindVertexArray(0);
}

void 
Scene_polylines_item::draw_points(Viewer_interface* viewer) const {
    glBindVertexArray(vao[0]);
    uniform_attrib(viewer,1);
    glUseProgram(rendering_program_lines);
    glDrawArrays(GL_POINTS, 0, positions_lines.size()/4);
    // Clean-up
    glUseProgram(0);
    glBindVertexArray(0);
}

void
Scene_polylines_item_private::
draw_spheres(const Scene_polylines_item* item) const {
    // FIRST, count the number of incident cycles and polylines
    // for all extremities.
    typedef std::map<Point_3, int> Point_to_int_map;
    typedef Point_to_int_map::iterator iterator;
    Point_to_int_map corner_polyline_nb;

    { // scope to fill corner_polyline_nb'
        Point_to_int_map corner_cycles_nb;

        for(std::list<std::vector<Point_3> >::const_iterator
            it = item->polylines.begin(),
            end = item->polylines.end();
            it != end; ++it)
        {
            const K::Point_3& a = *it->begin();
            const K::Point_3& b = *it->rbegin();
            if(a == b) {
                if ( it->size()>1 )
                    ++corner_cycles_nb[a];
                else
                    ++corner_polyline_nb[a];
            }
            else {
                ++corner_polyline_nb[a];
                ++corner_polyline_nb[b];
            }
        }
        // THEN, ignore points that are incident to one cycle only.
        for(iterator
            c_it = corner_cycles_nb.begin(),
            end = corner_cycles_nb.end();
            c_it != end; ++c_it)
        {
            const Point_3& a = c_it->first;

            iterator p_it = corner_polyline_nb.find(a);

            // If the point 'a'=c_it->first has only incident cycles...
            if(p_it == corner_polyline_nb.end()) {
                // ...then count it as a corner only if it has two incident cycles
                // or more.
                if(c_it->second > 1) {
                    corner_polyline_nb[a] = c_it->second;
                }
            } else {
                // else add the number of cycles.
                p_it->second += c_it->second;
            }
        }
    }
    // At this point, 'corner_polyline_nb' gives the multiplicity of all
    // corners.
    for(iterator
        p_it = corner_polyline_nb.begin(),
        end = corner_polyline_nb.end();
        p_it != end; ++p_it)
    {
        switch(p_it->second) {
        case 1:
            ::glColor3d(0.0, 0.0, 0.0); // black
            break;
        case 2:
            ::glColor3d(0.0, 0.8, 0.0); // green
            break;
        case 3:
            ::glColor3d(0.0, 0.0, 0.8); // blue
            break;
        case 4:
            ::glColor3d(0.8, 0.0, 0.0); //red
            break;
        default:
            ::glColor3d(0.8, 0.0, 0.8); //fuschia
        }
        this->draw_sphere(p_it->first, this->spheres_drawn_radius);
    }
}

void 
Scene_polylines_item_private::draw_sphere(const K::Point_3& p,
                                          double r) const
{
    if(sphere_display_list == 0) {
        sphere_display_list = glGenLists(1);
        if(sphere_display_list == 0)
            std::cerr << "ERROR: Cannot create display list!\n";
        if(quadric == 0)
            quadric = gluNewQuadric();
        if(quadric == 0)
            std::cerr << "ERROR: Cannot create GLU quadric!\n";
        glNewList(sphere_display_list, GL_COMPILE);
        gluSphere(quadric, 1., 10, 10);
        glEndList();
        if(glGetError() != GL_NO_ERROR)
            std::cerr << gluErrorString(glGetError());
    }
    glPushMatrix();
    glTranslated(CGAL::to_double(p.x()),
                 CGAL::to_double(p.y()),
                 CGAL::to_double(p.z()));

    glScaled(r, r, r);
    glCallList(sphere_display_list);
    glPopMatrix();
}

QMenu* Scene_polylines_item::contextMenu() 
{
    const char* prop_name = "Menu modified by Scene_polylines_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // http://doc.trolltech.com/lastest/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
        menu->addSeparator();
        // TODO: add actions to display corners
        QAction* action = menu->addAction(tr("Display corners with radius..."));
        connect(action, SIGNAL(triggered()),
                this, SLOT(change_corner_radii()));

        QAction* actionSmoothPolylines =
                menu->addAction(tr("Smooth polylines"));
        actionSmoothPolylines->setObjectName("actionSmoothPolylines");
        connect(actionSmoothPolylines, SIGNAL(triggered()),this, SLOT(smooth()));
        menu->setProperty(prop_name, true);
    }
    return menu;
}

void Scene_polylines_item::changed()
{
    compute_elements();
    initialize_buffers();


}

void Scene_polylines_item::change_corner_radii() {
    bool ok = true;
    double proposed_radius = d->spheres_drawn_radius;
    if(proposed_radius == 0) {
        Scene_interface::Bbox b = bbox();
        proposed_radius = (std::max)(b.xmax - b.xmin,
                                     proposed_radius);
        proposed_radius = (std::max)(b.ymax - b.ymin,
                                     proposed_radius);
        proposed_radius = (std::max)(b.zmax - b.zmin,
                                     proposed_radius);
        proposed_radius /= 100;
    }
    double r = QInputDialog::getDouble(NULL,
                                       tr("Display corners with new radius..."),
                                       tr("Radius:"),
                                       proposed_radius, // value
                                       0.,          // min
                                       2147483647., // max
                                       10,          // decimals
                                       &ok);
    if(ok) {
        change_corner_radii(r);
    }
}

void Scene_polylines_item::change_corner_radii(double r) {
    if(r >= 0) {
        d->spheres_drawn_radius = r;
        d->draw_extremities = (r > 0);
        this->changed();
        emit itemChanged();
    }
}

void Scene_polylines_item::split_at_sharp_angles()
{
    typedef Polylines_container Bare_polyline_container;
    typedef Polyline Bare_polyline;
    Polylines_container& bare_polylines = polylines;

    int counter = 0;
    for(Bare_polyline_container::iterator
        bare_polyline_it = bare_polylines.begin();
        bare_polyline_it != bare_polylines.end(); // the end changes
        // during the loop
        ++counter /* bare_polyline_it is incremented in the loop */)
    {
        Bare_polyline_container::iterator current_polyline_it =
                bare_polyline_it;
        Bare_polyline& bare_polyline = *bare_polyline_it;
        Bare_polyline::iterator it = boost::next(bare_polyline.begin());

        if(boost::next(bare_polyline.begin()) == bare_polyline.end())
        {
            std::cerr << "WARNING: Isolated point in polylines\n";
            bare_polyline_it = bare_polylines.erase(bare_polyline_it);
            continue;
        }
        else
            ++bare_polyline_it;
        if(it != bare_polyline.end()) {
            for(; it != boost::prior(bare_polyline.end()); ++it) {
                const Point_3 pv = *it;
                const Point_3 pa = *boost::prior(it);
                const Point_3 pb = *boost::next(it);
                const K::Vector_3 av = pv - pa;
                const K::Vector_3 bv = pv - pb;
                const K::FT sc_prod = av * bv;
                if( sc_prod >= 0 ||
                        (sc_prod < 0 &&
                         CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
                {
#ifdef PROTECTION_DEBUG
                    std::cerr << "Split polyline (small angle) "
                              <<  std::acos(sqrt(CGAL::square(sc_prod) /
                                                 ((av*av) * (bv*bv)))) * 180 /CGAL_PI
                               << " degres\n";
#endif
                    Bare_polyline new_polyline;
                    std::copy(it, bare_polyline.end(),
                              std::back_inserter(new_polyline));

                    if(*bare_polyline.begin() == *bare_polyline.rbegin()) {
                        // if the polyline is a cycle, test if its beginning is a sharp
                        // angle...
                        const Point_3 pv = *bare_polyline.begin();
                        const Point_3 pa = *boost::prior(boost::prior(bare_polyline.end()));
                        const Point_3 pb = *boost::next(bare_polyline.begin());
                        const K::Vector_3 av = pv - pa;
                        const K::Vector_3 bv = pv - pb;
                        const K::FT sc_prod = av * bv;
                        if( sc_prod >= 0 ||
                                (sc_prod < 0 &&
                                 CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
                        {
                            // if its beginning is a sharp angle, then split
                            bare_polyline.erase(boost::next(it), bare_polyline.end());
                        }
                        else {
                            // ...if not, modifies its beginning
                            std::copy(boost::next(bare_polyline.begin()),
                                      boost::next(it),
                                      std::back_inserter(new_polyline));
                            bare_polylines.erase(current_polyline_it);
                        }
                    }
                    else {
                        bare_polyline.erase(boost::next(it), bare_polyline.end());
                    }
                    bare_polylines.push_back(new_polyline);
                    break;
                }
            }
        }
    }
    emit itemChanged();
}

void
Scene_polylines_item::merge(Scene_polylines_item* other_item) {
    if(other_item == 0) return;
    std::copy(other_item->polylines.begin(),
              other_item->polylines.end(),
              std::back_inserter(polylines));
    QVariant other_metadata_variant = other_item->property("polylines metadata");
    if(other_metadata_variant.type() == QVariant::StringList)
    {
        QStringList metadata = property("polylines metadata").toStringList();
        metadata.append(other_metadata_variant.toStringList());
        setProperty("polylines metadata", metadata);
    }
    changed();
}

#include "Scene_polylines_item.moc"
