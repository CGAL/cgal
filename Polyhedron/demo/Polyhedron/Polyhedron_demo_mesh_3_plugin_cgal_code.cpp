#include <CGAL/AABB_intersections.h>
#include <CGAL/AABB_tree.h>

#include "Polyhedron_type.h"

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <Scene_polyhedron_item.h>
#include <Scene_polygon_soup_item.h>
#include <fstream>
#include <sstream>

#include <CGAL/Timer.h>

#include <QMenu>

typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel,
Polyhedron> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

// 3D complex
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Edge_criteria Edge_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

typedef Tr::Point Point_3;

#include "Scene_item.h"
#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

namespace {
void CGALglcolor(QColor c)
{
    ::glColor4d(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
}
}

class Q_DECL_EXPORT Scene_c3t3_item : public Scene_item
{
    Q_OBJECT
public:
    typedef qglviewer::ManipulatedFrame ManipulatedFrame;

    Scene_c3t3_item(const C3t3& c3t3)
        : c3t3_(c3t3), frame(new ManipulatedFrame()), last_known_scene(NULL)
    {
        positions_lines.resize(0);
        positions_poly.resize(0);
        color_lines.resize(0);
        color_poly.resize(0);
        color_grid.resize(0);
        normals.resize(0);
        glGenVertexArrays(2, vao);
        //Generates an integer which will be used as ID for each buffer
        glGenBuffers(7, buffer);
        compile_shaders();
    }

    ~Scene_c3t3_item()
    {
        glDeleteBuffers(7, buffer);
        glDeleteVertexArrays(2, vao);
        glDeleteProgram(rendering_program_lines);
        glDeleteProgram(rendering_program_poly);
        delete frame;
    }

    void changed()
    {
        compute_elements();
        initialize_buffers();
    }

    void contextual_changed()
    {
        if(frame->isInMouseGrabberPool())
            changed();
    }
    const C3t3& c3t3() const {
        return c3t3_;
    }

    bool manipulatable() const {
        return true;
    }
    ManipulatedFrame* manipulatedFrame() {
        return frame;
    }

    void setPosition(float x, float y, float z) {
        frame->setPosition(x, y, z);
    }

    void setNormal(float x, float y, float z) {
        frame->setOrientation(x, y, z, 0.f);
    }

    Kernel::Plane_3 plane() const {
        const qglviewer::Vec& pos = frame->position();
        const qglviewer::Vec& n =
                frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
        return Kernel::Plane_3(n[0], n[1],  n[2], - n * pos);
    }

    bool isFinite() const { return true; }
    bool isEmpty() const {
        return c3t3().triangulation().number_of_vertices() == 0;
    }

    Bbox bbox() const {
        if(isEmpty())
            return Bbox();
        else {
            CGAL::Bbox_3 result = c3t3().triangulation().finite_vertices_begin()->point().bbox();
            for(Tr::Finite_vertices_iterator
                vit = ++c3t3().triangulation().finite_vertices_begin(),
                end = c3t3().triangulation().finite_vertices_end();
                vit != end; ++vit)
            {
                result = result + vit->point().bbox();
            }
            return Bbox(result.xmin(), result.ymin(), result.zmin(),
                        result.xmax(), result.ymax(), result.zmax());
        }
    }

    Scene_c3t3_item* clone() const {
        return 0;
    }

    QString toolTip() const {
        int number_of_tets = 0;
        for(Tr::Finite_cells_iterator
            cit = c3t3().triangulation().finite_cells_begin(),
            end = c3t3().triangulation().finite_cells_end();
            cit != end; ++cit)
        {
            if( c3t3().is_in_complex(cit) )
                ++number_of_tets;
        }
        return tr("<p><b>3D complex in a 3D triangulation</b></p>"
                  "<p>Number of vertices: %1<br />"
                  "Number of surface facets: %2<br />"
                  "Number of volume tetrahedra: %3</p>")
                .arg(c3t3().triangulation().number_of_vertices())
                .arg(c3t3().number_of_facets())
                .arg(number_of_tets);
    }

    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const {
        return (m != Gouraud && m!=PointsPlusNormals && m!=Splatting); // CHECK THIS!
    }

    void draw() const {
        ::glPushMatrix();
        ::glMultMatrixd(frame->matrix());
        QGLViewer::drawGrid((float)complex_diag());
        ::glPopMatrix();

        if(isEmpty())
            return;

        GLboolean two_side;
        ::glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &two_side);
        ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

        const Kernel::Plane_3& plane = this->plane();
        GLdouble clip_plane[4];
        clip_plane[0] = -plane.a();
        clip_plane[1] = -plane.b();
        clip_plane[2] = -plane.c();
        clip_plane[3] = -plane.d();

        ::glClipPlane(GL_CLIP_PLANE0, clip_plane);
        ::glEnable(GL_CLIP_PLANE0);
        ::glBegin(GL_TRIANGLES);
        for(C3t3::Facet_iterator
            fit = c3t3().facets_begin(),
            end = c3t3().facets_end();
            fit != end; ++fit)
        {
            const Tr::Cell_handle& cell = fit->first;
            const int& index = fit->second;
            const Kernel::Point_3& pa = cell->vertex((index+1)&3)->point();
            const Kernel::Point_3& pb = cell->vertex((index+2)&3)->point();
            const Kernel::Point_3& pc = cell->vertex((index+3)&3)->point();
            typedef Kernel::Oriented_side Side;
            using CGAL::ON_ORIENTED_BOUNDARY;
            const Side sa = plane.oriented_side(pa);
            const Side sb = plane.oriented_side(pb);
            const Side sc = plane.oriented_side(pc);
            if( sa != ON_ORIENTED_BOUNDARY &&
                    sb != ON_ORIENTED_BOUNDARY &&
                    sc != ON_ORIENTED_BOUNDARY &&
                    sb == sa && sc == sa )
            {
                // draw_triangle(pa, pb, pc);
            }
        }
        ::glEnd();
        ::glDisable(GL_CLIP_PLANE0);

        ::glBegin(GL_TRIANGLES);
        // workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
        CGALglcolor(this->color().darker(150));
#undef darker
        for(Tr::Finite_cells_iterator
            cit = c3t3().triangulation().finite_cells_begin(),
            end = c3t3().triangulation().finite_cells_end();
            cit != end; ++cit)
        {
            if(! c3t3().is_in_complex(cit) )
                continue;

            const Kernel::Point_3& pa = cit->vertex(0)->point();
            const Kernel::Point_3& pb = cit->vertex(1)->point();
            const Kernel::Point_3& pc = cit->vertex(2)->point();
            const Kernel::Point_3& pd = cit->vertex(3)->point();
            typedef Kernel::Oriented_side Side;
            using CGAL::ON_ORIENTED_BOUNDARY;
            const Side sa = plane.oriented_side(pa);
            const Side sb = plane.oriented_side(pb);
            const Side sc = plane.oriented_side(pc);
            const Side sd = plane.oriented_side(pd);

            if( sa == ON_ORIENTED_BOUNDARY ||
                    sb == ON_ORIENTED_BOUNDARY ||
                    sc == ON_ORIENTED_BOUNDARY ||
                    sd == ON_ORIENTED_BOUNDARY ||
                    sb != sa || sc != sa || sd != sa)
            {
                // draw_triangle(pa, pb, pc);
                // draw_triangle(pa, pb, pd);
                // draw_triangle(pa, pc, pd);
                // draw_triangle(pb, pc, pd);
            }

            //       for(int i = 0; i < 4; ++i) {
            //         if(c3t3().is_in_complex(cit, i)) continue;
            //         const Point_3& pa = cit->vertex((i+1)&3)->point();
            //         const Point_3& pb = cit->vertex((i+2)&3)->point();
            //         const Point_3& pc= cit->vertex((i+3)&3)->point();
            //         typedef Kernel::Oriented_side Side;
            //         using CGAL::ON_ORIENTED_BOUNDARY;
            //         const Side sa = plane.oriented_side(pa);
            //         const Side sb = plane.oriented_side(pb);
            //         const Side sc = plane.oriented_side(pc);

            //         if( sa == ON_ORIENTED_BOUNDARY ||
            //             sb == ON_ORIENTED_BOUNDARY ||
            //             sc == ON_ORIENTED_BOUNDARY ||
            //             sb != sa || sc != sa )
            //         {
            //           draw_triangle(pa, pb, pc);
            //         }
            //       }
        }
        ::glEnd();
        if(!two_side)
            ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    };
    void draw(Viewer_interface* viewer) const {

        glBindVertexArray(vao[0]);
        glUseProgram(rendering_program_poly);
        uniform_attrib(viewer,0);
        glDrawArrays(GL_TRIANGLES, 0, positions_poly.size()/3);
        glUseProgram(0);
        glBindVertexArray(0);


    }
    void draw_edges(Viewer_interface* viewer) const {
        glBindVertexArray(vao[1]);
        glUseProgram(rendering_program_lines);
        uniform_attrib(viewer,2);
        glDrawArrays(GL_LINES, 0, positions_grid.size()/3);
        glUseProgram(0);
        glBindVertexArray(0);

        glBindVertexArray(vao[0]);
        glUseProgram(rendering_program_lines);
        uniform_attrib(viewer,1);
        glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
        glUseProgram(0);
        glBindVertexArray(0);

    }
    void draw_points(Viewer_interface * viewer) const
    {
        glBindVertexArray(vao[0]);
        glUseProgram(rendering_program_lines);
        uniform_attrib(viewer,1);
        glDrawArrays(GL_POINTS, 0, positions_lines.size()/3);
        glUseProgram(0);
        glBindVertexArray(0);

        glBindVertexArray(vao[1]);
        glUseProgram(rendering_program_lines);
        uniform_attrib(viewer,2);
        glDrawArrays(GL_LINES, 0, positions_grid.size()/3);
        glUseProgram(0);
        glBindVertexArray(0);
    }
private:
    void draw_triangle(const Kernel::Point_3& pa,
                       const Kernel::Point_3& pb,
                       const Kernel::Point_3& pc, bool is_cut) {
        // workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
#undef darker
        Kernel::Vector_3 n = cross_product(pb - pa, pc - pa);
        n = n / CGAL::sqrt(n*n);

if(!is_cut)
{
    for(int i=0; i<3; i++)
    {

        color_poly.push_back(this->color().redF());
        color_poly.push_back(this->color().greenF());
        color_poly.push_back(this->color().blueF());
    }
}
else
{
    for(int i=0; i<3; i++)
    {

        color_poly.push_back(this->color().darker(150).redF());
        color_poly.push_back(this->color().darker(150).greenF());
        color_poly.push_back(this->color().darker(150).blueF());
    }
}
        for(int i=0; i<3; i++)
        {
            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());
        }
        positions_poly.push_back(pa.x());
        positions_poly.push_back(pa.y());
        positions_poly.push_back(pa.z());

        positions_poly.push_back(pb.x());
        positions_poly.push_back(pb.y());
        positions_poly.push_back(pb.z());

        positions_poly.push_back(pc.x());
        positions_poly.push_back(pc.y());
        positions_poly.push_back(pc.z());

    }

    void draw_triangle_edges(const Kernel::Point_3& pa,
                       const Kernel::Point_3& pb,
                       const Kernel::Point_3& pc) {
        // workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
#undef darker
        Kernel::Vector_3 n = cross_product(pb - pa, pc - pa);
        n = n / CGAL::sqrt(n*n);
        for(int i=0; i<6; i++)
        {
            color_lines.push_back(0.0);
            color_lines.push_back(0.0);
            color_lines.push_back(0.0);
        }
        positions_lines.push_back(pa.x());
        positions_lines.push_back(pa.y());
        positions_lines.push_back(pa.z());

        positions_lines.push_back(pb.x());
        positions_lines.push_back(pb.y());
        positions_lines.push_back(pb.z());

        positions_lines.push_back(pb.x());
        positions_lines.push_back(pb.y());
        positions_lines.push_back(pb.z());

        positions_lines.push_back(pc.x());
        positions_lines.push_back(pc.y());
        positions_lines.push_back(pc.z());

        positions_lines.push_back(pc.x());
        positions_lines.push_back(pc.y());
        positions_lines.push_back(pc.z());

        positions_lines.push_back(pa.x());
        positions_lines.push_back(pa.y());
        positions_lines.push_back(pa.z());

    }



    double complex_diag() const {
        const Bbox& bbox = this->bbox();
        const double& xdelta = bbox.xmax-bbox.xmin;
        const double& ydelta = bbox.ymax-bbox.ymin;
        const double& zdelta = bbox.zmax-bbox.zmin;
        const double diag = std::sqrt(xdelta*xdelta +
                                      ydelta*ydelta +
                                      zdelta*zdelta);
        return diag * 0.7;
    }

public slots:
    void export_facets_in_complex()
    {
        std::stringstream off_sstream;
        c3t3().output_facets_in_complex_to_off(off_sstream);
        std::string backup = off_sstream.str();
        // Try to read .off in a polyhedron
        Scene_polyhedron_item* item = new Scene_polyhedron_item();
        if(!item->load(off_sstream))
        {
            delete item;
            off_sstream.str(backup);

            // Try to read .off in a polygon soup
            Scene_polygon_soup_item* soup_item = new Scene_polygon_soup_item;

            if(!soup_item->load(off_sstream)) {
                delete soup_item;
                return;
            }

            soup_item->setName(QString("%1_%2").arg(this->name()).arg("facets"));
            last_known_scene->addItem(soup_item);
        }
        else{
            item->setName(QString("%1_%2").arg(this->name()).arg("facets"));
            last_known_scene->addItem(item);
        }
    }

public:

    QMenu* contextMenu()
    {
        const char* prop_name = "Menu modified by Scene_c3t3_item.";

        QMenu* menu = Scene_item::contextMenu();

        // Use dynamic properties:
        // http://doc.trolltech.com/lastest/qobject.html#property
        bool menuChanged = menu->property(prop_name).toBool();

        if(!menuChanged) {
            QAction* actionExportFacetsInComplex =
                    menu->addAction(tr("Export facets in complex"));
            actionExportFacetsInComplex->setObjectName("actionExportFacetsInComplex");
            connect(actionExportFacetsInComplex,
                    SIGNAL(triggered()),this,
                    SLOT(export_facets_in_complex()));
        }
        return menu;
    }

    void set_scene(Scene_interface* scene){ last_known_scene=scene; }

private:

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
    C3t3 c3t3_;
    qglviewer::ManipulatedFrame* frame;
    Scene_interface* last_known_scene;


    std::vector<float> positions_lines;
    std::vector<float> positions_grid;
    std::vector<float> positions_poly;
    std::vector<float> normals;
    std::vector<float> color_lines;
    std::vector<float> color_poly;
    std::vector<float> color_grid;

    GLint location[10];
    GLuint vao[2];
    GLuint buffer[7];
    GLuint rendering_program_lines;
    GLuint rendering_program_poly;

    void initialize_buffers()
    {
        glBindVertexArray(vao[0]);

        glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
        glBufferData(GL_ARRAY_BUFFER,
                     (positions_poly.size())*sizeof(float),
                     positions_poly.data(),
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
                     (positions_lines.size())*sizeof(float),
                     positions_lines.data(),
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
                     (normals.size())*sizeof(float),
                     normals.data(),
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
                     (color_poly.size())*sizeof(float),
                     color_poly.data(),
                     GL_STATIC_DRAW);
        glVertexAttribPointer(3, //number of the buffer
                              3, //number of floats to be taken
                              GL_FLOAT, // type of data
                              GL_FALSE, //not normalized
                              0, //compact data (not in a struct)
                              NULL //no offset (seperated in several buffers)
                              );
        glEnableVertexAttribArray(3);

        glBindBuffer(GL_ARRAY_BUFFER, buffer[4]);
        glBufferData(GL_ARRAY_BUFFER,
                     (color_lines.size())*sizeof(float),
                     color_lines.data(),
                     GL_STATIC_DRAW);
        glVertexAttribPointer(4, //number of the buffer
                              3, //number of floats to be taken
                              GL_FLOAT, // type of data
                              GL_FALSE, //not normalized
                              0, //compact data (not in a struct)
                              NULL //no offset (seperated in several buffers)
                              );
        glEnableVertexAttribArray(4);

        glBindVertexArray(vao[1]);

        glBindBuffer(GL_ARRAY_BUFFER, buffer[5]);
        glBufferData(GL_ARRAY_BUFFER,
                     (positions_grid.size())*sizeof(float),
                     positions_grid.data(),
                     GL_STATIC_DRAW);
        glVertexAttribPointer(1, //number of the buffer
                              3, //number of floats to be taken
                              GL_FLOAT, // type of data
                              GL_FALSE, //not normalized
                              0, //compact data (not in a struct)
                              NULL //no offset (seperated in several buffers)
                              );
        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, buffer[6]);
        glBufferData(GL_ARRAY_BUFFER,
                     (color_grid.size())*sizeof(float),
                     color_grid.data(),
                     GL_STATIC_DRAW);
        glVertexAttribPointer(4, //number of the buffer
                              3, //number of floats to be taken
                              GL_FLOAT, // type of data
                              GL_FALSE, //not normalized
                              0, //compact data (not in a struct)
                              NULL //no offset (seperated in several buffers)
                              );
        glEnableVertexAttribArray(4);
        glBindVertexArray(0);
    }
    void compile_shaders()
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
            "   diffuse = abs(dot(N,L)) * light_diff * color_facets; \n"
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

        rendering_program_poly = program;

        //For the edges
        static const GLchar* vertex_shader_source_lines[] =
        {
            "#version 300 es \n"
            " \n"
            "layout (location = 1) in vec3 positions; \n"
            "layout (location = 4) in vec3 color_lines; \n"

            "uniform mat4 mvp_matrix; \n"
            "uniform mat4 f_matrix; \n"
            "vec4 positions_lines = vec4(positions, 1.0); \n"
            "out highp vec3 fColors; \n"
            " \n"

            "void main(void) \n"
            "{ \n"
            "   fColors = color_lines; \n"
            "   gl_Position = mvp_matrix * f_matrix * positions_lines; \n"
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
        glDeleteShader(fragment_shader);
        rendering_program_lines = program;

    }
    void uniform_attrib(Viewer_interface* viewer, int mode) const
    {
        light_info light;
        GLint is_both_sides = 0;
        GLfloat mvp_mat[16];
        GLfloat mv_mat[16];
        GLfloat f_mat[16];

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
        //Polys
        if(mode ==0)
        {

            glUseProgram(rendering_program_poly);
            glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
            glUniformMatrix4fv(location[1], 1, GL_FALSE, mv_mat);
            glUniform3fv(location[2], 1, light.position);
            glUniform3fv(location[3], 1, light.diffuse);
            glUniform3fv(location[4], 1, light.specular);
            glUniform3fv(location[5], 1, light.ambient);
        }
        //Edges
        else if(mode ==1)
        {
            for(int i=0; i<16; i++)
            {
                if(i%5 == 0)
                    f_mat[i] = 1.0;
                else
                    f_mat[i] = 0.0;
            }
            glUseProgram(rendering_program_lines);
            glUniformMatrix4fv(location[6], 1, GL_FALSE, mvp_mat);
            glUniformMatrix4fv(location[7], 1, GL_FALSE, f_mat);

        }
        //Grid
        else if(mode ==2)
        {
            for(int i=0; i<16; i++)
                f_mat[i] = frame->matrix()[i];
            glUseProgram(rendering_program_lines);
            glUniformMatrix4fv(location[6], 1, GL_FALSE, mvp_mat);
            glUniformMatrix4fv(location[7], 1, GL_FALSE, f_mat);

        }

    }
    void compute_elements()
    {
        positions_lines.clear();
        positions_poly.clear();
        color_lines.clear();
        color_grid.clear();
        color_poly.clear();
        normals.clear();

        //The grid
        {
            float x = (2*(float)complex_diag())/10.0;
            float y = (2*(float)complex_diag())/10.0;
            for(int u = 0; u < 11; u++)
            {

                positions_grid.push_back(-(float)complex_diag() + x* u);
                positions_grid.push_back(-(float)complex_diag());
                positions_grid.push_back(0.0);

                positions_grid.push_back(-(float)complex_diag() + x* u);
                positions_grid.push_back((float)complex_diag());
                positions_grid.push_back(0.0);
            }
            for(int v=0; v<11; v++)
            {

                positions_grid.push_back(-(float)complex_diag());
                positions_grid.push_back(-(float)complex_diag() + v * y);
                positions_grid.push_back(0.0);

                positions_grid.push_back((float)complex_diag());
                positions_grid.push_back(-(float)complex_diag() + v * y);
                positions_grid.push_back(0.0);
            }
            float colors[3];
            colors[0] = this->color().redF();
            colors[1] = this->color().greenF();
            colors[2] = this->color().blueF();

            for(int i=0; i< 132; i++)
            {
                color_grid.push_back(colors[i%3]);
            }
        }

        //The facets
        {
            if(isEmpty())
                return;


            const Kernel::Plane_3& plane = this->plane();
            GLdouble clip_plane[4];
            clip_plane[0] = -plane.a();
            clip_plane[1] = -plane.b();
            clip_plane[2] = -plane.c();
            clip_plane[3] = -plane.d();



            for(C3t3::Facet_iterator
                fit = c3t3().facets_begin(),
                end = c3t3().facets_end();
                fit != end; ++fit)
            {
                const Tr::Cell_handle& cell = fit->first;
                const int& index = fit->second;
                const Kernel::Point_3& pa = cell->vertex((index+1)&3)->point();
                const Kernel::Point_3& pb = cell->vertex((index+2)&3)->point();
                const Kernel::Point_3& pc = cell->vertex((index+3)&3)->point();
                typedef Kernel::Oriented_side Side;
                using CGAL::ON_ORIENTED_BOUNDARY;
                const Side sa = plane.oriented_side(pa);
                const Side sb = plane.oriented_side(pb);
                const Side sc = plane.oriented_side(pc);
                bool is_showned = false;
                if(pa.x() * clip_plane[0]  + pa.y() * clip_plane[1]  + pa.z() * clip_plane[2] + clip_plane[3]  > 0
                        && pb.x() * clip_plane[0]  + pb.y() * clip_plane[1]  + pb.z() * clip_plane[2] + clip_plane[3]  > 0
                        && pc.x() * clip_plane[0]  + pc.y() * clip_plane[1]  + pc.z() * clip_plane[2] + clip_plane[3]  > 0)
                    is_showned = true;

                if(is_showned && sa != ON_ORIENTED_BOUNDARY &&
                        sb != ON_ORIENTED_BOUNDARY &&
                        sc != ON_ORIENTED_BOUNDARY &&
                        sb == sa && sc == sa )
                {
                    draw_triangle(pa, pb, pc, false);
                    draw_triangle_edges(pa, pb, pc);
                }

            }


            for(Tr::Finite_cells_iterator
                cit = c3t3().triangulation().finite_cells_begin(),
                end = c3t3().triangulation().finite_cells_end();
                cit != end; ++cit)
            {
                if(! c3t3().is_in_complex(cit) )
                    continue;

                const Kernel::Point_3& pa = cit->vertex(0)->point();
                const Kernel::Point_3& pb = cit->vertex(1)->point();
                const Kernel::Point_3& pc = cit->vertex(2)->point();
                const Kernel::Point_3& pd = cit->vertex(3)->point();
                typedef Kernel::Oriented_side Side;
                using CGAL::ON_ORIENTED_BOUNDARY;
                const Side sa = plane.oriented_side(pa);
                const Side sb = plane.oriented_side(pb);
                const Side sc = plane.oriented_side(pc);
                const Side sd = plane.oriented_side(pd);

                if( sa == ON_ORIENTED_BOUNDARY ||
                        sb == ON_ORIENTED_BOUNDARY ||
                        sc == ON_ORIENTED_BOUNDARY ||
                        sd == ON_ORIENTED_BOUNDARY ||
                        sb != sa || sc != sa || sd != sa)
                {
                    draw_triangle(pa,pb,pc, true);
                    draw_triangle(pa,pb,pd, true);
                    draw_triangle(pa,pc,pd, true);
                    draw_triangle(pb,pc,pd, true);

                    draw_triangle_edges(pa,pb,pc);
                    draw_triangle_edges(pa,pb,pd);
                    draw_triangle_edges(pa,pc,pd);
                    draw_triangle_edges(pb,pc,pd);


                }


            }

        }



        location[0] = glGetUniformLocation(rendering_program_poly, "mvp_matrix");
        location[1] = glGetUniformLocation(rendering_program_poly, "mv_matrix");
        location[2] = glGetUniformLocation(rendering_program_poly, "light_pos");
        location[3] = glGetUniformLocation(rendering_program_poly, "light_diff");
        location[4] = glGetUniformLocation(rendering_program_poly, "light_spec");
        location[5] = glGetUniformLocation(rendering_program_poly, "light_amb");
        location[6] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
        location[7] = glGetUniformLocation(rendering_program_lines, "f_matrix");
    }

};

Scene_item* cgal_code_mesh_3(const Polyhedron* pMesh,
                             QString filename,
                             const double angle,
                             const double facet_sizing,
                             const double approx,
                             const double tet_sizing,
                             const double tet_shape,
                             const bool protect_features,
                             Scene_interface* scene)
{
    if(!pMesh) return 0;

    // remesh

    // Set mesh criteria
    Edge_criteria edge_criteria(facet_sizing);
    Facet_criteria facet_criteria(angle, facet_sizing, approx); // angle, size, approximation
    Cell_criteria cell_criteria(tet_shape, tet_sizing); // radius-edge ratio, size
    Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

    CGAL::Timer timer;
    timer.start();
    std::cerr << "Meshing file \"" << qPrintable(filename) << "\"\n";
    std::cerr << "  angle: " << angle << std::endl
              << "  facets size bound: " << facet_sizing << std::endl
              << "  approximation bound: " << approx << std::endl
              << "  tetrahedra size bound: " << tet_sizing << std::endl;
    std::cerr << "Build AABB tree...";
    // Create domain
    Mesh_domain domain(*pMesh);
    if(protect_features) {
        domain.detect_features();
    }
    std::cerr << "done (" << timer.time() << " ms)" << std::endl;

    // Meshing
    std::cerr << "Mesh...";
    CGAL::parameters::internal::Features_options features =
            protect_features ?
                CGAL::parameters::features(domain) :
                CGAL::parameters::no_features();

    Scene_c3t3_item* new_item =
            new Scene_c3t3_item(CGAL::make_mesh_3<C3t3>(domain, criteria, features));
    new_item->set_scene(scene);
    std::cerr << "done (" << timer.time() << " ms, " << new_item->c3t3().triangulation().number_of_vertices() << " vertices)" << std::endl;

    if(new_item->c3t3().triangulation().number_of_vertices() > 0)
    {
        std::ofstream medit_out("out.mesh");
        new_item->c3t3().output_to_medit(medit_out);

        const Scene_item::Bbox& bbox = new_item->bbox();
        new_item->setPosition((float)(bbox.xmin + bbox.xmax)/2.f,
                              (float)(bbox.ymin + bbox.ymax)/2.f,
                              (float)(bbox.zmin + bbox.zmax)/2.f);
        return new_item;
    }
    else {
        delete new_item;
        return 0;
    }
}

#include "Scene_c3t3_item.moc"
