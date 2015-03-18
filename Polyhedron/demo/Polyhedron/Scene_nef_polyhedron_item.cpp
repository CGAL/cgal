#include "Scene_nef_polyhedron_item.h"
#include "Scene_polyhedron_item.h"
#include "Nef_type.h"
#include "Polyhedron_type.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>
// #include <CGAL/OFF_to_nef_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Inverse_index.h>

#include <QObject>
#include "Scene_nef_rendering.h"

#include <CGAL/minkowski_sum_3.h>
#include <CGAL/convex_decomposition_3.h>

#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

typedef Nef_polyhedron::Traits Traits;
typedef Nef_polyhedron::Halffacet Facet;
typedef CGAL::Triangulation_2_filtered_projection_traits_3<Traits>   P_traits;
typedef typename Nef_polyhedron::Halfedge_const_handle Halfedge_handle;
struct Face_info {
    typename Nef_polyhedron::Halfedge_const_handle e[3];
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
struct DPoint {
    DPoint(GLdouble x, GLdouble y, GLdouble z)
    {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    }
    GLdouble coords[3];
};
Scene_nef_polyhedron_item::Scene_nef_polyhedron_item()
    : Scene_item_with_display_list(),
      positions_facets(0),
      positions_lines(0),
      color_lines(0),
      color_facets(0),
      color_points(0),
      normals(0),
      positions_points(0),
      nef_poly(new Nef_polyhedron)
{
    is_selected = true;
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(7, buffer);
    compile_shaders();

}

Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(Nef_polyhedron* const p)
    : Scene_item_with_display_list(),
      positions_facets(0),
      positions_lines(0),
      color_lines(0),
      color_facets(0),
      color_points(0),
      normals(0),
      positions_points(0),
      nef_poly(p)
{
    is_selected = true;
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(7, buffer);
    compile_shaders();

}

Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(const Nef_polyhedron& p)
    : Scene_item_with_display_list(),
      positions_facets(0),
      positions_lines(0),
      normals(0),
      color_lines(0),
      color_facets(0),
      color_points(0),
      positions_points(0),
      nef_poly(new Nef_polyhedron(p))
{
     is_selected = true;
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(7, buffer);
    compile_shaders();

}

// Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(const Scene_nef_polyhedron_item& item)
//   : Scene_item_with_display_list(item),
//     nef_poly(new Nef_polyhedron(*item.nef_poly))
// {
// }

Scene_nef_polyhedron_item::~Scene_nef_polyhedron_item()
{
    glDeleteBuffers(6, buffer);
    glDeleteVertexArrays(1, vao);
    glDeleteProgram(rendering_program_facets);
    glDeleteProgram(rendering_program_lines);
    glDeleteProgram(rendering_program_points);
    delete nef_poly;
}

void Scene_nef_polyhedron_item::initialize_buffers()
{
    glBindVertexArray(vao[0]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_facets.size())*sizeof(double),
                 positions_facets.data(),
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
                 (positions_lines.size())*sizeof(double),
                 positions_lines.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(1, //number of the buffer
                          3, //number of floats to be taken
                          GL_DOUBLE, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[2]);
    glBufferData(GL_ARRAY_BUFFER,
                 (normals.size())*sizeof(double),
                 normals.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(2,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[3]);
    glBufferData(GL_ARRAY_BUFFER,
                 (color_facets.size())*sizeof(double),
                 color_facets.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(3,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[4]);
    glBufferData(GL_ARRAY_BUFFER,
                 (color_lines.size())*sizeof(double),
                 color_lines.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(4,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(4);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[5]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_points.size())*sizeof(double),
                 positions_points.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(5,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(5);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[6]);
        glBufferData(GL_ARRAY_BUFFER,
                     (color_points.size())*sizeof(double),
                     color_points.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(6,
                              3,
                              GL_DOUBLE,
                              GL_FALSE,
                              0,
                              NULL
                              );
        glEnableVertexAttribArray(6);
    // Clean-up
    glBindVertexArray(0);


}

void Scene_nef_polyhedron_item::compile_shaders(void)
{
    //fill the vertex shader
    static const GLchar* vertex_shader_source[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec3 positions; \n"
        "layout (location = 2) in vec3 vNormals; \n"
        "layout (location = 3) in vec3 color_facets; \n"

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
        "       diffuse = abs(dot(N,L)) * light_diff * color_facets; \n"
        "   else \n"
        "       diffuse = max(dot(N,L), 0.0) * light_diff * color_facets; \n"
        "   vec3 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

        "   fColors =light_amb*color_facets+ diffuse + specular ; \n"

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
        "layout (location = 1) in vec3 positions; \n"
        "layout (location = 4) in vec3 color_lines; \n"

        "uniform mat4 mvp_matrix; \n"

        "out highp vec3 fColors; \n"
        "vec4 positions_lines = vec4(positions, 1.0); \n"
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
    //Clean-up
    glDeleteShader(vertex_shader);
    rendering_program_lines = program;

    //For the points
    static const GLchar* vertex_shader_source_points[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 5) in vec3 positions; \n"
        "layout (location = 6) in vec3 color; \n"
        "uniform mat4 mvp_matrix; \n"

        "out highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color; \n"
        "   gl_Position = mvp_matrix * vec4(positions, 1.0); \n"
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

void Scene_nef_polyhedron_item::compute_normals_and_vertices(void)
{
     int count = 0;
    positions_facets.clear();
    positions_points.clear();
    color_lines.clear();
    color_facets.clear();
    color_points.clear();
    normals.clear();
    //The Facets
    {
        for(Nef_polyhedron::Halffacet_const_iterator
            f = nef_poly->halffacets_begin (),
            end = nef_poly->halffacets_end();
            f != end; ++f)
        {
            if(f->is_twin()) continue;
            count++;
            Nef_polyhedron::Vector_3 v = f->plane().orthogonal_vector();
            P_traits cdt_traits(v);
            CDT cdt(cdt_traits);

            for(Nef_polyhedron::Halffacet_cycle_const_iterator
                fc = f->facet_cycles_begin(),
                end = f->facet_cycles_end();
                fc != end; ++fc)
            {
                if ( fc.is_shalfedge() )
                {

                    Nef_polyhedron::SHalfedge_const_handle h = fc;
                    Nef_polyhedron::SHalfedge_around_facet_const_circulator hc(h), he(hc);

                    typename CDT::Vertex_handle previous, first;

                    do {
                        Nef_polyhedron::SVertex_const_handle v = hc->source();
                        const Nef_polyhedron::Point_3& point = v->source()->point();
                        typename CDT::Vertex_handle vh = cdt.insert(point);
                        if(first == 0) {
                            first = vh;
                        }
                        vh->info() = hc->source();
                        if(previous != 0 && previous != vh) {
                            cdt.insert_constraint(previous, vh);
                        }
                        previous = vh;
                    } while( ++hc != he );

                    cdt.insert_constraint(previous, first);

                    // sets mark is_external
                    for(typename CDT::All_faces_iterator
                        fit = cdt.all_faces_begin(),
                        end = cdt.all_faces_end();
                        fit != end; ++fit)
                    {
                        fit->info().is_external = false;

                    }
                    //check if the facet is external or internal
                    std::queue<typename CDT::Face_handle> face_queue;
                    face_queue.push(cdt.infinite_vertex()->face());

                    while(! face_queue.empty() ) {
                        typename CDT::Face_handle fh = face_queue.front();
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

                    for(typename CDT::Finite_faces_iterator
                        ffit = cdt.finite_faces_begin(),
                        end = cdt.finite_faces_end();
                        ffit != end; ++ffit)
                    {


                        if(ffit->info().is_external){ continue;}
                        for(int i = 0; i<3; i++)
                        {
                            positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().x()));
                            positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().y()));
                            positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().z()));

                        }



                        Nef_polyhedron::Vector_3 v = f->plane().orthogonal_vector();
                        GLdouble normal[3];
                        normal[0] = CGAL::to_double(v.x());
                        normal[1] = CGAL::to_double(v.y());
                        normal[2] = CGAL::to_double(v.z());
                        GLdouble norm = normal[0]*normal[0]
                                + normal[1]*normal[1]
                                + normal[2]*normal[2];
                        norm = CGAL::sqrt(norm);
                        normal[0] /= norm;
                        normal[1] /= norm;
                        normal[2] /= norm;

                        normals.push_back(normal[0]);
                        normals.push_back(normal[1]);
                        normals.push_back(normal[2]);

                        normals.push_back(normal[0]);
                        normals.push_back(normal[1]);
                        normals.push_back(normal[2]);

                        normals.push_back(normal[0]);
                        normals.push_back(normal[1]);
                        normals.push_back(normal[2]);

                        if(is_selected)
                        {
                            color_facets.push_back(this->color().lighter(120).redF());
                            color_facets.push_back(this->color().lighter(120).greenF());
                            color_facets.push_back(this->color().lighter(120).blueF());

                            color_facets.push_back(this->color().lighter(120).redF());
                            color_facets.push_back(this->color().lighter(120).greenF());
                            color_facets.push_back(this->color().lighter(120).blueF());

                            color_facets.push_back(this->color().lighter(120).redF());
                            color_facets.push_back(this->color().lighter(120).greenF());
                            color_facets.push_back(this->color().lighter(120).blueF());
                        }
                        else
                        {
                            color_facets.push_back(this->color().redF());
                            color_facets.push_back(this->color().greenF());
                            color_facets.push_back(this->color().blueF());

                            color_facets.push_back(this->color().redF());
                            color_facets.push_back(this->color().greenF());
                            color_facets.push_back(this->color().blueF());

                            color_facets.push_back(this->color().redF());
                            color_facets.push_back(this->color().greenF());
                            color_facets.push_back(this->color().blueF());

                        }

                    }
                }
            }
        }

    } // end facets

    //The Lines
    {
       for(Nef_polyhedron::Halfedge_const_iterator
            e = nef_poly->halfedges_begin(),
            end = nef_poly->halfedges_end();
            e != end; ++e)
        {
            if (e->is_twin()) continue;
            const Nef_polyhedron::Vertex_const_handle& s = e->source();
            const Nef_polyhedron::Vertex_const_handle& t = e->twin()->source();
            const Nef_polyhedron::Point_3& a = s->point();
            const Nef_polyhedron::Point_3& b = t->point();

            positions_lines.push_back(CGAL::to_double(a.x()));
            positions_lines.push_back(CGAL::to_double(a.y()));
            positions_lines.push_back(CGAL::to_double(a.z()));

            positions_lines.push_back(CGAL::to_double(b.x()));
            positions_lines.push_back(CGAL::to_double(b.y()));
            positions_lines.push_back(CGAL::to_double(b.z()));

            if(is_selected)
            {
                color_lines.push_back(this->color().lighter(50).redF());
                color_lines.push_back(this->color().lighter(50).greenF());
                color_lines.push_back(this->color().lighter(50).blueF());

                color_lines.push_back(this->color().lighter(50).redF());
                color_lines.push_back(this->color().lighter(50).greenF());
                color_lines.push_back(this->color().lighter(50).blueF());
            }
            else
            {
                color_lines.push_back(0.0);
                color_lines.push_back(0.0);
                color_lines.push_back(0.0);

                color_lines.push_back(0.0);
                color_lines.push_back(0.0);
                color_lines.push_back(0.0);
            }
        }
    }
    //The points
    {
        for(Nef_polyhedron::Vertex_const_iterator
            v = nef_poly->vertices_begin(),
            end = nef_poly->vertices_end();
            v != end; ++v)
        {
            const Nef_polyhedron::Point_3& p = v->point();
            positions_points.push_back(CGAL::to_double(p.x()));
            positions_points.push_back(CGAL::to_double(p.y()));
            positions_points.push_back(CGAL::to_double(p.z()));

                color_points.push_back(this->color().lighter(50).redF());
                color_points.push_back(this->color().lighter(50).greenF());
                color_points.push_back(this->color().lighter(50).blueF());

                color_points.push_back(this->color().lighter(50).redF());
                color_points.push_back(this->color().lighter(50).greenF());
                color_points.push_back(this->color().lighter(50).blueF());

        }

    } //end points

    location[0] = glGetUniformLocation(rendering_program_facets, "mvp_matrix");
    location[1] = glGetUniformLocation(rendering_program_facets, "mv_matrix");
    location[2] = glGetUniformLocation(rendering_program_facets, "light_pos");
    location[3] = glGetUniformLocation(rendering_program_facets, "light_diff");
    location[4] = glGetUniformLocation(rendering_program_facets, "light_spec");
    location[5] = glGetUniformLocation(rendering_program_facets, "light_amb");
    location[6] = glGetUniformLocation(rendering_program_facets, "is_two_side");

    location[7] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");

    location[8] = glGetUniformLocation(rendering_program_points, "mvp_matrix");


}

void Scene_nef_polyhedron_item::uniform_attrib(Viewer_interface* viewer, int mode) const
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
        glUseProgram(rendering_program_facets);
        glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
        glUniformMatrix4fv(location[1], 1, GL_FALSE, mv_mat);
        glUniform3fv(location[2], 1, light.position);
        glUniform3fv(location[3], 1, light.diffuse);
        glUniform3fv(location[4], 1, light.specular);
        glUniform3fv(location[5], 1, light.ambient);
        glUniform1i(location[6], is_both_sides);
    }
    else if(mode ==1)
    {
        glUseProgram(rendering_program_lines);
        glUniformMatrix4fv(location[7], 1, GL_FALSE, mvp_mat);
    }
    else if(mode ==2)
    {
        glUseProgram(rendering_program_points);
        glUniformMatrix4fv(location[8], 1, GL_FALSE, mvp_mat);
    }
}

Scene_nef_polyhedron_item* 
Scene_nef_polyhedron_item::clone() const {
    return new Scene_nef_polyhedron_item(*nef_poly);
}

bool
Scene_nef_polyhedron_item::load_from_off(std::istream& in)
{
    //   const std::size_t discarded = CGAL::OFF_to_nef_3(in, *nef_poly);
    //   return discarded != 0;

    Exact_polyhedron exact_poly;
    in >> exact_poly;
    *nef_poly = Nef_polyhedron(exact_poly);

    //   Polyhedron poly;
    //   in >> poly;
    //   *nef_poly = Nef_polyhedron(poly);
    changed();
    return (bool) in;
}

QFont 
Scene_nef_polyhedron_item::font() const {
    QFont font;
    font.setItalic(!font.italic());
    return font;
}

bool
Scene_nef_polyhedron_item::load(std::istream& in)
{
    in >> *nef_poly;
    changed();
    return (bool) in;
}

bool
Scene_nef_polyhedron_item::save(std::ostream& in) const
{
    in << *nef_poly;
    return (bool) in;
}

QString 
Scene_nef_polyhedron_item::toolTip() const
{
    if(!nef_poly)
        return QString();

    return QObject::tr("<p><b>%1</b> (mode: %5, color: %6)<br />"
                       "<i>Nef_3 polyhedron</i></p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                       "Number of facets: %4<br />"
                       "number of volumes: %7</p>")
            .arg(this->name())
            .arg(nef_poly->number_of_vertices())
            .arg(nef_poly->number_of_edges())
            .arg(nef_poly->number_of_facets())
            .arg(this->renderingModeName())
            .arg(this->color().name())
            .arg(nef_poly->number_of_volumes());
}

void
Scene_nef_polyhedron_item::direct_draw() const {
    gl_render_nef_facets(nef_poly);

    GLboolean lighting;
    glGetBooleanv(GL_LIGHTING, &lighting);
    glDisable(GL_LIGHTING);

    GLfloat point_size;
    glGetFloatv(GL_POINT_SIZE, &point_size);
    glPointSize(10.f);

    gl_render_nef_vertices(nef_poly);

    if(lighting) {
        glEnable(GL_LIGHTING);
    }
    glPointSize(point_size);
}
void Scene_nef_polyhedron_item::draw(Viewer_interface* viewer) const
{
    glBindVertexArray(vao[0]);

    // tells the GPU to use the program just created
    glUseProgram(rendering_program_facets);
    uniform_attrib(viewer,0);
    //draw the polygons
    // the third argument is the number of vec4 that will be entered
    glDrawArrays(GL_TRIANGLES, 0, positions_facets.size()/3);


    GLfloat point_size;
    glGetFloatv(GL_POINT_SIZE, &point_size);
    glPointSize(10.f);

    draw_points(viewer);
    glPointSize(point_size);

    glUseProgram(0);
    glBindVertexArray(0);

}
void Scene_nef_polyhedron_item::draw_edges(Viewer_interface* viewer) const
{
    glBindVertexArray(vao[0]);
    glUseProgram(rendering_program_lines);
    uniform_attrib(viewer ,1);
    glDrawArrays(GL_LINES,0,positions_lines.size()/3);

    GLfloat point_size;
    glGetFloatv(GL_POINT_SIZE, &point_size);
    glPointSize(10.f);

    draw_points(viewer);
    glPointSize(point_size);

    glUseProgram(0);
    glBindVertexArray(0);
}
void Scene_nef_polyhedron_item::draw_points(Viewer_interface* viewer) const
{
    glBindVertexArray(vao[0]);
    glUseProgram(rendering_program_points);
    uniform_attrib(viewer ,2);
    glDrawArrays(GL_POINTS,0,positions_points.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);

}


Nef_polyhedron* 
Scene_nef_polyhedron_item::nef_polyhedron() {
    return nef_poly;
}

bool
Scene_nef_polyhedron_item::isEmpty() const {
    return (nef_poly == 0) || nef_poly->is_empty();
}

Scene_nef_polyhedron_item::Bbox
Scene_nef_polyhedron_item::bbox() const {
    if(isEmpty())
        return Bbox();
    CGAL::Bbox_3 bbox(nef_poly->vertices_begin()->point().bbox());
    for(Nef_polyhedron::Vertex_const_iterator it = nef_poly->vertices_begin();
        it != nef_poly->vertices_end();
        ++it) {
        bbox = bbox + it->point().bbox();
    }
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}

// quick hacks to convert polyhedra from exact to inexact and vice-versa
template <class Polyhedron_input,
          class Polyhedron_output>
struct Copy_polyhedron_to 
        : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
    Copy_polyhedron_to(const Polyhedron_input& in_poly)
        : in_poly(in_poly) {}

    void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
    {
        typedef typename Polyhedron_output::HalfedgeDS Output_HDS;

        CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

        typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
        typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
        typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

        builder.begin_surface(in_poly.size_of_vertices(),
                              in_poly.size_of_facets(),
                              in_poly.size_of_halfedges());

        for(Vertex_const_iterator
            vi = in_poly.vertices_begin(), end = in_poly.vertices_end();
            vi != end ; ++vi)
        {
            typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
                                                  ::CGAL::to_double( vi->point().y()),
                                                  ::CGAL::to_double( vi->point().z()));
            builder.add_vertex(p);
        }

        typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
        Index index( in_poly.vertices_begin(), in_poly.vertices_end());

        for(Facet_const_iterator
            fi = in_poly.facets_begin(), end = in_poly.facets_end();
            fi != end; ++fi)
        {
            HFCC hc = fi->facet_begin();
            HFCC hc_end = hc;
            //     std::size_t n = circulator_size( hc);
            //     CGAL_assertion( n >= 3);
            builder.begin_facet ();
            do {
                builder.add_vertex_to_facet(index[hc->vertex()]);
                ++hc;
            } while( hc != hc_end);
            builder.end_facet();
        }
        builder.end_surface();
    } // end operator()(..)
private:
    const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_A, class Poly_B>
void copy_to(const Poly_A& poly_a, Poly_B& poly_b)
{
    Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
    poly_b.delegate(modifier);
}

void from_exact(Exact_polyhedron& in,
                Polyhedron& out)
{
    copy_to(in, out);
    CGAL_assertion(out.is_valid());
}

void to_exact(Polyhedron& in,
              Exact_polyhedron& out)
{
    copy_to(in, out);
    CGAL_assertion(out.is_valid());
}

bool Scene_nef_polyhedron_item::is_simple() const
{
    return nef_poly->is_simple();
}

// [static]
Scene_nef_polyhedron_item* 
Scene_nef_polyhedron_item::from_polyhedron(Scene_polyhedron_item* item)
{
    Polyhedron* poly = item->polyhedron();
    if(!poly) return 0;

    Exact_polyhedron exact_poly;
    to_exact(*poly, exact_poly);
    Nef_polyhedron* nef_poly = new Nef_polyhedron(exact_poly);
    exact_poly.clear();

    return new Scene_nef_polyhedron_item(nef_poly);
}

Scene_polyhedron_item*
Scene_nef_polyhedron_item::convert_to_polyhedron() const {
    Exact_polyhedron exact_poly;
    nef_poly->convert_to_Polyhedron(exact_poly);
    Polyhedron* poly = new Polyhedron;
    from_exact(exact_poly, *poly);
    exact_poly.clear();
    return new Scene_polyhedron_item(poly);
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator+=(const Scene_nef_polyhedron_item& other)
{
    (*nef_poly) += (*other.nef_poly);
    return *this;
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator*=(const Scene_nef_polyhedron_item& other)
{
    (*nef_poly) *= (*other.nef_poly);
    return *this;
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator-=(const Scene_nef_polyhedron_item& other)
{
    (*nef_poly) -= (*other.nef_poly);
    return *this;
}

Scene_nef_polyhedron_item*
Scene_nef_polyhedron_item::
sum(const Scene_nef_polyhedron_item& a, 
    const Scene_nef_polyhedron_item& b)
{
    return new Scene_nef_polyhedron_item(CGAL::minkowski_sum_3(*a.nef_poly,
                                                               *b.nef_poly));
}

void
Scene_nef_polyhedron_item::
convex_decomposition(std::list< Scene_polyhedron_item*>& convex_parts)
{
    // copy the Nef polyhedron, as markers are added
    Nef_polyhedron N(*nef_poly);
    CGAL::convex_decomposition_3(N);

    typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;

    Volume_const_iterator ci = ++N.volumes_begin();
    for( ; ci != N.volumes_end(); ++ci) {
        if(ci->mark()) {
            Exact_polyhedron P;
            N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
            Polyhedron* poly = new Polyhedron;
            from_exact(P, *poly);
            Scene_polyhedron_item *spoly = new Scene_polyhedron_item(poly);
            convex_parts.push_back(spoly);
            spoly->changed();
        }

    }
}

void
Scene_nef_polyhedron_item::
changed()
{
    //  init();
    Base::changed();
    compute_normals_and_vertices();
    initialize_buffers();
}
#include "Scene_nef_polyhedron_item.moc"
void
Scene_nef_polyhedron_item::selection_changed(bool p_is_selected)
{

    if(p_is_selected != is_selected)
    {
        is_selected = p_is_selected;
        changed();
    }

}
