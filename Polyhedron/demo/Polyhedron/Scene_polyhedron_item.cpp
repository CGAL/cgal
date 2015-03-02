#include <CGAL/AABB_intersections.h>
#include "Scene_polyhedron_item.h"
#include "Kernel_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>



#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>

#include <CGAL/internal/Operations_on_polyhedra/compute_normal.h>

#include <QVariant>
#include <list>
#include <queue>
#include <iostream>

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Input_facets_AABB_tree;
const char* aabb_property_name = "Scene_polyhedron_item aabb tree";

Input_facets_AABB_tree* get_aabb_tree(Scene_polyhedron_item* item)
{
    QVariant aabb_tree_property = item->property(aabb_property_name);
    if(aabb_tree_property.isValid()) {
        void* ptr = aabb_tree_property.value<void*>();
        return static_cast<Input_facets_AABB_tree*>(ptr);
    }
    else {
        Polyhedron* poly = item->polyhedron();
        if(poly) {
            Input_facets_AABB_tree* tree =
                    new Input_facets_AABB_tree(faces(*poly).first,
                                               faces(*poly).second,
                                               *poly);
            item->setProperty(aabb_property_name,
                              QVariant::fromValue<void*>(tree));
            return tree;
        }
        else return 0;
    }
}

void delete_aabb_tree(Scene_polyhedron_item* item)
{
    QVariant aabb_tree_property = item->property(aabb_property_name);
    if(aabb_tree_property.isValid()) {
        void* ptr = aabb_tree_property.value<void*>();
        Input_facets_AABB_tree* tree = static_cast<Input_facets_AABB_tree*>(ptr);
        if(tree) {
            delete tree;
            tree = 0;
        }
        item->setProperty(aabb_property_name, QVariant());
    }
}

typedef typename Polyhedron::Traits Traits;
typedef typename Polyhedron::Facet Facet;
typedef CGAL::Triangulation_2_filtered_projection_traits_3<Traits>   P_traits;
typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
struct Face_info {
    typename Polyhedron::Halfedge_handle e[3];
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

//Make sure all the facets are triangles
void
Scene_polyhedron_item::triangulate_facet(Facet_iterator fit)
{
    //Computes the normal of the facet
    typename Traits::Vector_3 normal =
            compute_facet_normal<Facet,Traits>(*fit);

    P_traits cdt_traits(normal);
    CDT cdt(cdt_traits);

    typename Facet::Halfedge_around_facet_circulator
            he_circ = fit->facet_begin(),
            he_circ_end(he_circ);

    // Iterates on the vector of facet handles
    typename CDT::Vertex_handle previous, first;
    do {
        typename CDT::Vertex_handle vh = cdt.insert(he_circ->vertex()->point());
        if(first == 0) {
            first = vh;
        }
        vh->info() = he_circ;
        if(previous != 0 && previous != vh) {
            cdt.insert_constraint(previous, vh);
        }
        previous = vh;
    } while( ++he_circ != he_circ_end );
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
        if(ffit->info().is_external)
            continue;

        positions_facets.push_back(ffit->vertex(0)->point().x());
        positions_facets.push_back(ffit->vertex(0)->point().y());
        positions_facets.push_back(ffit->vertex(0)->point().z());
        positions_facets.push_back(1.0);

        positions_facets.push_back(ffit->vertex(1)->point().x());
        positions_facets.push_back(ffit->vertex(1)->point().y());
        positions_facets.push_back(ffit->vertex(1)->point().z());
        positions_facets.push_back(1.0);

        positions_facets.push_back(ffit->vertex(2)->point().x());
        positions_facets.push_back(ffit->vertex(2)->point().y());
        positions_facets.push_back(ffit->vertex(2)->point().z());
        positions_facets.push_back(1.0);

        if (cur_shading == GL_FLAT || cur_shading == GL_SMOOTH)
        {

            typedef typename Kernel::Vector_3	    Vector;
            Vector n = compute_facet_normal<Facet,Traits>(*fit);
            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());

            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());

            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());
        }




    }
}

void
Scene_polyhedron_item::triangulate_facet_color(Facet_iterator fit)
{
    typename Traits::Vector_3 normal =
            compute_facet_normal<Facet,Traits>(*fit);

    P_traits cdt_traits(normal);
    CDT cdt(cdt_traits);

    typename Facet::Halfedge_around_facet_circulator
            he_circ = fit->facet_begin(),
            he_circ_end(he_circ);

    // Iterates on the vector of facet handles
    typename CDT::Vertex_handle previous, first;
    do {
        typename CDT::Vertex_handle vh = cdt.insert(he_circ->vertex()->point());
        if(first == 0) {
            first = vh;
        }
        vh->info() = he_circ;
        if(previous != 0 && previous != vh) {
            cdt.insert_constraint(previous, vh);
        }
        previous = vh;
    } while( ++he_circ != he_circ_end );
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

    //iterates on the internal faces to add the vertices to the positions vector
    for(typename CDT::Finite_faces_iterator
        ffit = cdt.finite_faces_begin(),
        end = cdt.finite_faces_end();
        ffit != end; ++ffit)
    {
        if(ffit->info().is_external)
            continue;
        //Add Colors
        for(int i = 0; i<3; ++i)
        {
            const int this_patch_id = fit->patch_id();
            if(is_selected)
            {
                color_facets.push_back(colors_[this_patch_id].lighter(120).redF());
                color_facets.push_back(colors_[this_patch_id].lighter(120).greenF());
                color_facets.push_back(colors_[this_patch_id].lighter(120).blueF());

                color_facets.push_back(colors_[this_patch_id].lighter(120).redF());
                color_facets.push_back(colors_[this_patch_id].lighter(120).greenF());
                color_facets.push_back(colors_[this_patch_id].lighter(120).blueF());
            }
            else
            {
                color_facets.push_back(colors_[this_patch_id].redF());
                color_facets.push_back(colors_[this_patch_id].greenF());
                color_facets.push_back(colors_[this_patch_id].blueF());

                color_facets.push_back(colors_[this_patch_id].redF());
                color_facets.push_back(colors_[this_patch_id].greenF());
                color_facets.push_back(colors_[this_patch_id].blueF());

            }
        }
    }
}

#include <QObject>
#include <QMenu>
#include <QAction>
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

void
Scene_polyhedron_item::initialize_buffers()
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
                 (color_facets.size())*sizeof(float),
                 color_facets.data(), GL_STATIC_DRAW);
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
                 (color_lines.size())*sizeof(float),
                 color_lines.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(4,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(4);

    // Clean-up
    glBindVertexArray(0);

}

void
Scene_polyhedron_item::compile_shaders(void)
{
    //fill the vertex shader
    static const GLchar* vertex_shader_source[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec4 positions_facets; \n"
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

        "   fColors = light_amb*color_facets + diffuse + specular ; \n"

        "   gl_Position =  mvp_matrix *positions_facets; \n"
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
        " color = fColors; \n" //polygons de couleur rouge
        "} \n"
    };

    //creates and compiles the vertex shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source, NULL);
    glCompileShader(vertex_shader);

    /*  GLint result;
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
    }*/
    GLuint fragment_shader =	glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    /*   glGetShaderiv(fragment_shader,GL_COMPILE_STATUS,&result);
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
*/

    //creates the program, attaches and links the shaders
    GLuint program= glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);


    /*   glGetProgramiv(program,GL_LINK_STATUS,&result);
    if(result == GL_TRUE){
        std::cout<<"Link OK"<<std::endl;
    } else {
        int maxLength;
        int length;
        glGetProgramiv(program,GL_INFO_LOG_LENGTH,&maxLength);
        char* log = new char[maxLength];
        glGetProgramInfoLog(program,maxLength,&length,log);
        std::cout<<"link error : Length = "<<length<<", log ="<<log<<std::endl;
    }*/
    //Delete the shaders which are now in the memory
    glDeleteShader(vertex_shader);

    rendering_program_facets = program;

    //For the edges
    //fill the vertex shader
    static const GLchar* vertex_shader_source_lines[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 1) in vec4 positions_lines; \n"
        "layout (location = 4) in vec3 color_lines; \n"

        "uniform mat4 mvp_matrix; \n"

        "out highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color_lines; \n"
        "   gl_Position = mvp_matrix * positions_lines; \n"
        "} \n"
    };

    //creates and compiles the vertex shader
    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_lines, NULL);
    glCompileShader(vertex_shader);


    /*  glGetShaderiv(vertex_shader,GL_COMPILE_STATUS,&result);
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
*/
    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    /*   glGetShaderiv(fragment_shader,GL_COMPILE_STATUS,&result);
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
*/

    //creates the program, attaches and links the shaders
    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);


    /* glGetProgramiv(program,GL_LINK_STATUS,&result);
    if(result == GL_TRUE){
        std::cout<<"Link OK"<<std::endl;
    } else {
        int maxLength;
        int length;
        glGetProgramiv(program,GL_INFO_LOG_LENGTH,&maxLength);
        char* log = new char[maxLength];
        glGetProgramInfoLog(program,maxLength,&length,log);
        std::cout<<"link error : Length = "<<length<<", log ="<<log<<std::endl;
    }*/
    //Delete the shaders which are now in the memory
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    rendering_program_lines = program;

}

void
Scene_polyhedron_item::uniform_attrib(Viewer_interface* viewer, int mode) const
{
    light_info light;
    GLint is_both_sides = 0;
    GLfloat mvp_mat[16];
    GLfloat mv_mat[16];

    //fills the MVP and MV matrices.

    GLdouble d_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    //Convert the GLdoubles matrix in GLfloats
    for (int i=0; i<16; ++i){
        mvp_mat[i] = GLfloat(d_mat[i]);
    }

    viewer->camera()->getModelViewMatrix(d_mat);
    for (int i=0; i<16; ++i)
        mv_mat[i] = GLfloat(d_mat[i]);

    glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE, &is_both_sides);

    //fills the arraw of colors with the current color


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
        //Set the light infos
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

}
void
Scene_polyhedron_item::compute_normals_and_vertices(void)
{
    positions_facets.clear();
    positions_lines.clear();
    normals.clear();


    //Facets


    typedef typename Polyhedron::Traits	    Kernel;
    typedef typename Kernel::Point_3	    Point;
    typedef typename Kernel::Vector_3	    Vector;
    typedef typename Polyhedron::Facet	    Facet;
    typedef typename Polyhedron::Facet_iterator Facet_iterator;
    typedef typename Polyhedron::Halfedge_around_facet_circulator HF_circulator;



    Facet_iterator f = poly->facets_begin();


    for(f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {
        triangulate_facet(f);
        /* HF_circulator he = f->facet_begin();
         HF_circulator end = he;
       CGAL_For_all(he,end)
        {

            // If Flat shading:1 normal per polygon added once per vertex
            if (cur_shading == GL_FLAT)
            {

                Vector n = compute_facet_normal<Facet,Kernel>(*f);
                normals.push_back(n.x());
                normals.push_back(n.y());
                normals.push_back(n.z());
            }


            // revolve around current face to get vertices


            // If Gouraud shading: 1 normal per vertex
            else if (cur_shading == GL_SMOOTH)
            {

                Vector n = compute_vertex_normal<typename Polyhedron::Vertex,Kernel>(*he->vertex());
                normals.push_back(n.x());
                normals.push_back(n.y());
                normals.push_back(n.z());
            }
            const int this_patch_id = f->patch_id();

            //position
            const Point& p = he->vertex()->point();
            positions_facets.push_back(p.x());
            positions_facets.push_back(p.y());
            positions_facets.push_back(p.z());
            positions_facets.push_back(1.0);
        }
        */
        std::cout<<"Nombre de vertices = "<<positions_facets.size()/4<<std::endl;
    }
    //Lines
    typedef Kernel::Point_3		Point;
    typedef Polyhedron::Edge_iterator	Edge_iterator;

    Edge_iterator he;
    if(!show_only_feature_edges_m) {
        for(he = poly->edges_begin();
            he != poly->edges_end();
            he++)
        {
            if(he->is_feature_edge()) continue;
            const Point& a = he->vertex()->point();
            const Point& b = he->opposite()->vertex()->point();
            positions_lines.push_back(a.x());
            positions_lines.push_back(a.y());
            positions_lines.push_back(a.z());
            positions_lines.push_back(1.0);

            positions_lines.push_back(b.x());
            positions_lines.push_back(b.y());
            positions_lines.push_back(b.z());
            positions_lines.push_back(1.0);

            /* color_lines.push_back(colors[0]);
            color_lines.push_back(colors[1]);
            color_lines.push_back(colors[2]);

            color_lines.push_back(colors[0]);
            color_lines.push_back(colors[1]);
            color_lines.push_back(colors[2]);*/

        }
    }
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
        if(!he->is_feature_edge()) continue;
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();

        positions_lines.push_back(a.x());
        positions_lines.push_back(a.y());
        positions_lines.push_back(a.z());
        positions_lines.push_back(1.0);

        positions_lines.push_back(b.x());
        positions_lines.push_back(b.y());
        positions_lines.push_back(b.z());
        positions_lines.push_back(1.0);

        /*  color_lines.push_back(1.0);
        color_lines.push_back(0.0);
        color_lines.push_back(0.0);

        color_lines.push_back(1.0);
        color_lines.push_back(0.0);
        color_lines.push_back(0.0);*/
    }


    //set the colors
    compute_colors();
    //Allocates a uniform location for the MVP and MV matrices
    location[0] = glGetUniformLocation(rendering_program_facets, "mvp_matrix");
    location[1] = glGetUniformLocation(rendering_program_facets, "mv_matrix");

    //Allocates a uniform location for the light values
    location[2] = glGetUniformLocation(rendering_program_facets, "light_pos");
    location[3] = glGetUniformLocation(rendering_program_facets, "light_diff");
    location[4] = glGetUniformLocation(rendering_program_facets, "light_spec");
    location[5] = glGetUniformLocation(rendering_program_facets, "light_amb");
    location[6] = glGetUniformLocation(rendering_program_facets, "is_two_side");

    location[7] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");

}

void
Scene_polyhedron_item::compute_colors()
{
    GLfloat colors[4];
    color_lines.clear();
    color_facets.clear();
    //Facets


    typedef typename Polyhedron::Traits	    Kernel;
    typedef typename Kernel::Point_3	    Point;
    typedef typename Kernel::Vector_3	    Vector;
    typedef typename Polyhedron::Facet	    Facet;
    typedef typename Polyhedron::Facet_iterator Facet_iterator;
    typedef typename Polyhedron::Halfedge_around_facet_circulator HF_circulator;



    // int patch_id = -1;

    Facet_iterator f = poly->facets_begin();
    QColor temp = colors_[f->patch_id()];
    if(is_selected)
    {
        colors[0]=0.0;
        colors[1]=0.0;
        colors[2]=0.0;
    }
    else
    {
        colors[0]=temp.lighter(50).redF();
        colors[1]=temp.lighter(50).greenF();
        colors[2]=temp.lighter(50).blueF();
    }
    for(f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {
        triangulate_facet_color(f);
        /*HF_circulator he = f->facet_begin();
        HF_circulator end = he;
        CGAL_For_all(he,end)
        {

            const int this_patch_id = f->patch_id();
            if(is_selected)
            {
                color_facets.push_back(colors_[this_patch_id].lighter(120).redF());
                color_facets.push_back(colors_[this_patch_id].lighter(120).greenF());
                color_facets.push_back(colors_[this_patch_id].lighter(120).blueF());
            }
            else
            {
                color_facets.push_back(colors_[this_patch_id].redF());
                color_facets.push_back(colors_[this_patch_id].greenF());
                color_facets.push_back(colors_[this_patch_id].blueF());
            }
        }*/
    }
    //Lines
    typedef Polyhedron::Edge_iterator	Edge_iterator;

    Edge_iterator he;
    if(!show_only_feature_edges_m) {
        for(he = poly->edges_begin();
            he != poly->edges_end();
            he++)
        {
            if(he->is_feature_edge()) continue;
            color_lines.push_back(colors[0]);
            color_lines.push_back(colors[1]);
            color_lines.push_back(colors[2]);

            color_lines.push_back(colors[0]);
            color_lines.push_back(colors[1]);
            color_lines.push_back(colors[2]);

        }
    }
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
        if(!he->is_feature_edge()) continue;
        color_lines.push_back(1.0);
        color_lines.push_back(0.0);
        color_lines.push_back(0.0);

        color_lines.push_back(1.0);
        color_lines.push_back(0.0);
        color_lines.push_back(0.0);
    }
}

Scene_polyhedron_item::Scene_polyhedron_item()
    : Scene_item(),
      positions_facets(0),
      positions_lines(0),
      normals(0),
      poly(new Polyhedron),
      show_only_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    cur_shading=GL_FLAT;
    is_selected = false;
    //init();
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(5, buffer);
    compile_shaders();
}

Scene_polyhedron_item::Scene_polyhedron_item(Polyhedron* const p)
    : Scene_item(),
      positions_facets(0),
      positions_lines(0),
      normals(0),
      poly(p),
      show_only_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    cur_shading=GL_FLAT;
    is_selected = false;
    init();
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(5, buffer);

    compile_shaders();
}

Scene_polyhedron_item::Scene_polyhedron_item(const Polyhedron& p)
    : Scene_item(),
      positions_facets(0),
      positions_lines(0),
      normals(0),
      poly(new Polyhedron(p)),
      show_only_feature_edges_m(false),
      facet_picking_m(false),
      erase_next_picked_facet_m(false),
      plugin_has_set_color_vector_m(false)
{
    cur_shading=GL_FLAT;
    is_selected=false;
    init();
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(5, buffer);

    compile_shaders();
}

// Scene_polyhedron_item::Scene_polyhedron_item(const Scene_polyhedron_item& item)
//   : Scene_item(item),
//     poly(new Polyhedron(*item.poly)),
//     show_only_feature_edges_m(false)
// {
// }

Scene_polyhedron_item::~Scene_polyhedron_item()
{
    glDeleteBuffers(5, buffer);
    glDeleteVertexArrays(1, vao);
    glDeleteProgram(rendering_program_facets);
    glDeleteProgram(rendering_program_lines);

    delete_aabb_tree(this);
    delete poly;
}

#include "Color_map.h"

void
Scene_polyhedron_item::
init()
{
    typedef Polyhedron::Facet_iterator Facet_iterator;

    if ( !plugin_has_set_color_vector_m )
    {
        // Fill indices map and get max subdomain value
        int max = 0;
        for(Facet_iterator fit = poly->facets_begin(), end = poly->facets_end() ;
            fit != end; ++fit)
        {
            max = (std::max)(max, fit->patch_id());
        }

        colors_.clear();
        compute_color_map(this->color(), max + 1,
                          std::back_inserter(colors_));
    }
}


Scene_polyhedron_item*
Scene_polyhedron_item::clone() const {
    return new Scene_polyhedron_item(*poly);}

// Load polyhedron from .OFF file
bool
Scene_polyhedron_item::load(std::istream& in)
{


    in >> *poly;

    if ( in && !isEmpty() )
    {
        changed();
        return true;
    }
    return false;
}

// Write polyhedron to .OFF file
bool
Scene_polyhedron_item::save(std::ostream& out) const
{
    out.precision(13);
    out << *poly;
    return (bool) out;
}

QString
Scene_polyhedron_item::toolTip() const
{
    if(!poly)
        return QString();

    return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
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

QMenu* Scene_polyhedron_item::contextMenu()
{
    const char* prop_name = "Menu modified by Scene_polyhedron_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // http://doc.trolltech.com/lastest/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {

        QAction* actionShowOnlyFeatureEdges =
                menu->addAction(tr("Show only &feature edges"));
        actionShowOnlyFeatureEdges->setCheckable(true);
        actionShowOnlyFeatureEdges->setObjectName("actionShowOnlyFeatureEdges");
        connect(actionShowOnlyFeatureEdges, SIGNAL(toggled(bool)),
                this, SLOT(show_only_feature_edges(bool)));

        QAction* actionPickFacets =
                menu->addAction(tr("Facets picking"));
        actionPickFacets->setCheckable(true);
        actionPickFacets->setObjectName("actionPickFacets");
        connect(actionPickFacets, SIGNAL(toggled(bool)),
                this, SLOT(enable_facets_picking(bool)));

        QAction* actionEraseNextFacet =
                menu->addAction(tr("Erase next picked facet"));
        actionEraseNextFacet->setCheckable(true);
        actionEraseNextFacet->setObjectName("actionEraseNextFacet");
        connect(actionEraseNextFacet, SIGNAL(toggled(bool)),
                this, SLOT(set_erase_next_picked_facet(bool)));

        menu->setProperty(prop_name, true);
    }
    QAction* action = menu->findChild<QAction*>("actionPickFacets");
    if(action) action->setChecked(facet_picking_m);
    action = menu->findChild<QAction*>("actionEraseNextFacet");
    if(action) action->setChecked(erase_next_picked_facet_m);
    return menu;
}

void Scene_polyhedron_item::show_only_feature_edges(bool b)
{
    show_only_feature_edges_m = b;
    emit itemChanged();
}

void Scene_polyhedron_item::enable_facets_picking(bool b)
{
    facet_picking_m = b;
}

void Scene_polyhedron_item::set_erase_next_picked_facet(bool b)
{
    if(b) { facet_picking_m = true; } // automatically activate facet_picking
    erase_next_picked_facet_m = b;
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::draw(Viewer_interface* viewer) const {
    glBindVertexArray(vao[0]);

    // tells the GPU to use the program just created
    //  glUseProgram(rendering_program_facets);
    uniform_attrib(viewer,0);
    //draw the polygons
    // the third argument is the number of vec4 that will be entered
    //viewer->displayMessage(tr("TEST"),5000);
    glDrawArrays(GL_TRIANGLES, 0, positions_facets.size()/4);
    glUseProgram(0);
    glBindVertexArray(0);

}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::draw_edges(Viewer_interface* viewer) const {
    glBindVertexArray(vao[0]);

    // tells the GPU to use the program just created
    // glUseProgram(rendering_program_lines);

    uniform_attrib(viewer,1);

    //draw the edges
    glDrawArrays(GL_LINES, 0, positions_lines.size()/4);
    // Clean-up
    glUseProgram(0);
    glBindVertexArray(0);
}

void
Scene_polyhedron_item::draw_points(Viewer_interface* viewer) const {

    glBindVertexArray(vao[0]);

    // tells the GPU to use the program just created
    //glUseProgram(rendering_program_lines);

    uniform_attrib(viewer,1);

    //draw the points
    glDrawArrays(GL_POINTS, 0, positions_facets.size());

    // Clean-up
    glBindVertexArray(0);
}

Polyhedron*
Scene_polyhedron_item::polyhedron()       { return poly; }
const Polyhedron*
Scene_polyhedron_item::polyhedron() const { return poly; }

bool
Scene_polyhedron_item::isEmpty() const {
    return (poly == 0) || poly->empty();
}

Scene_polyhedron_item::Bbox
Scene_polyhedron_item::bbox() const {
    const Kernel::Point_3& p = *(poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedron::Point_iterator it = poly->points_begin();
        it != poly->points_end();
        ++it) {
        bbox = bbox + it->bbox();
    }
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}


void
Scene_polyhedron_item::
changed()
{
    emit item_is_about_to_be_changed();
    delete_aabb_tree(this);
    init();
    Base::changed();
    compute_normals_and_vertices();
    initialize_buffers();

}
void
Scene_polyhedron_item::
shading_mode_changed()
{
    GLint new_shading;
    glGetIntegerv(GL_SHADE_MODEL, &new_shading);
    prev_shading = cur_shading;
    cur_shading = new_shading;
    if(prev_shading != cur_shading)
        if(cur_shading == GL_SMOOTH || cur_shading == GL_FLAT && prev_shading == GL_SMOOTH )
        {
            //Change the normals
            changed();
        }
}
void
Scene_polyhedron_item::selection_changed(bool p_is_selected)
{
    if(p_is_selected != is_selected)
    {is_selected = p_is_selected;
        compute_colors();
        initialize_buffers();
    }
    else
        is_selected = p_is_selected;
}

void
Scene_polyhedron_item::select(double orig_x,
                              double orig_y,
                              double orig_z,
                              double dir_x,
                              double dir_y,
                              double dir_z)
{
    if(facet_picking_m) {
        typedef Input_facets_AABB_tree Tree;
        typedef Tree::Object_and_primitive_id Object_and_primitive_id;

        Tree* aabb_tree = get_aabb_tree(this);
        if(aabb_tree) {
            const Kernel::Point_3 ray_origin(orig_x, orig_y, orig_z);
            const Kernel::Vector_3 ray_dir(dir_x, dir_y, dir_z);
            const Kernel::Ray_3 ray(ray_origin, ray_dir);

            typedef std::list<Object_and_primitive_id> Intersections;
            Intersections intersections;

            aabb_tree->all_intersections(ray, std::back_inserter(intersections));

            Intersections::iterator closest = intersections.begin();
            if(closest != intersections.end()) {
                const Kernel::Point_3* closest_point =
                        CGAL::object_cast<Kernel::Point_3>(&closest->first);

                for(Intersections::iterator
                    it = boost::next(intersections.begin()),
                    end = intersections.end();
                    it != end; ++it)
                {
                    if(! closest_point) {
                        closest = it;
                    }
                    else {
                        const Kernel::Point_3* it_point =
                                CGAL::object_cast<Kernel::Point_3>(&it->first);
                        if(it_point &&
                                (ray_dir * (*it_point - *closest_point)) < 0)
                        {
                            closest = it;
                            closest_point = it_point;
                        }
                    }
                }
                if(closest_point) {
                    Polyhedron::Facet_handle selected_fh = closest->second;

                    // The computation of the nearest vertex may be costly.  Only
                    // do it if some objects are connected to the signal
                    // 'selected_vertex'.
                    if(QObject::receivers(SIGNAL(selected_vertex(void*))) > 0)
                    {
                        Polyhedron::Halfedge_around_facet_circulator
                                he_it = selected_fh->facet_begin(),
                                around_end = he_it;

                        Polyhedron::Vertex_handle v = he_it->vertex(), nearest_v = v;

                        Kernel::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                                    v->point());

                        while(++he_it != around_end) {
                            v = he_it->vertex();
                            Kernel::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                                            v->point());
                            if(new_sq_dist < sq_dist) {
                                sq_dist = new_sq_dist;
                                nearest_v = v;
                            }
                        }

                        emit selected_vertex((void*)(&*nearest_v));
                    }

                    if(QObject::receivers(SIGNAL(selected_edge(void*))) > 0
                            || QObject::receivers(SIGNAL(selected_halfedge(void*))) > 0)
                    {
                        Polyhedron::Halfedge_around_facet_circulator
                                he_it = selected_fh->facet_begin(),
                                around_end = he_it;

                        Polyhedron::Halfedge_handle nearest_h = he_it;
                        Kernel::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                                    Kernel::Segment_3(he_it->vertex()->point(), he_it->opposite()->vertex()->point()));

                        while(++he_it != around_end) {
                            Kernel::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                                            Kernel::Segment_3(he_it->vertex()->point(), he_it->opposite()->vertex()->point()));
                            if(new_sq_dist < sq_dist) {
                                sq_dist = new_sq_dist;
                                nearest_h = he_it;
                            }
                        }

                        emit selected_halfedge((void*)(&*nearest_h));
                        emit selected_edge((void*)(std::min)(&*nearest_h, &*nearest_h->opposite()));
                    }

                    emit selected_facet((void*)(&*selected_fh));
                    if(erase_next_picked_facet_m) {
                        polyhedron()->erase_facet(selected_fh->halfedge());
                        polyhedron()->normalize_border();
                        //set_erase_next_picked_facet(false);
                        changed();
                        emit itemChanged();
                    }
                }
            }
        }
    }
    Base::select(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z);
}

void Scene_polyhedron_item::update_vertex_indices()
{
    std::size_t id=0;
    for (Polyhedron::Vertex_iterator vit = polyhedron()->vertices_begin(),
         vit_end = polyhedron()->vertices_end(); vit != vit_end; ++vit)
    {
        vit->id()=id++;
    }
}
void Scene_polyhedron_item::update_facet_indices()
{
    std::size_t id=0;
    for (Polyhedron::Facet_iterator  fit = polyhedron()->facets_begin(),
         fit_end = polyhedron()->facets_end(); fit != fit_end; ++fit)
    {
        fit->id()=id++;
    }
}
void Scene_polyhedron_item::update_halfedge_indices()
{
    std::size_t id=0;
    for (Polyhedron::Halfedge_iterator hit = polyhedron()->halfedges_begin(),
         hit_end = polyhedron()->halfedges_end(); hit != hit_end; ++hit)
    {
        hit->id()=id++;
    }
}

#include "Scene_polyhedron_item.moc"

