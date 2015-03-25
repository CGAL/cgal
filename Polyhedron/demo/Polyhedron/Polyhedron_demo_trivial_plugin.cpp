#include <GL/glew.h>
#include <QtCore/qglobal.h>

#include "Scene_item.h"
#include "Scene_interface.h"
#include <CGAL/gl.h>

#include "Viewer_interface.h"
#include <QAction>
#include <QMainWindow>

class Q_DECL_EXPORT Scene_bbox_item : public Scene_item
{
    Q_OBJECT
public:
    Scene_bbox_item(const Scene_interface* scene_interface)
        : scene(scene_interface)
    {

        positions_lines.resize(0);
        glGenVertexArrays(1, vao);
        //Generates an integer which will be used as ID for each buffer
        glGenBuffers(1, buffer);
        compile_shaders();
    }
    ~Scene_bbox_item()
    {
        glDeleteBuffers(1, buffer);
        glDeleteVertexArrays(1, vao);
        glDeleteProgram(rendering_program_lines);
    }
    bool isFinite() const { return true; }
    bool isEmpty() const { return true; }
    Bbox bbox() const { return Bbox(); }

    Scene_bbox_item* clone() const {
        return 0;
    }

    QString toolTip() const {
        const Bbox& bb = scene->bbox();
        return QString("<p><b>Scene bounding box</b></p>"
                       "<p>x range: (%1, %2)<br />"
                       "y range: (%3, %4)<br />"
                       "z range: (%5, %6)</p>")
                .arg(bb.xmin).arg(bb.xmax)
                .arg(bb.ymin).arg(bb.ymax)
                .arg(bb.zmin).arg(bb.zmax);
    }

    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const {
        return (m == Wireframe);
    }

    // Flat/Gouraud OpenGL drawing
    void draw() const {}

    // Wireframe OpenGL drawing
    void draw_edges() const {
        const Bbox& bb = scene->bbox();
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
    void draw_edges(Viewer_interface* viewer) const
    {
        glBindVertexArray(vao[0]);
        glUseProgram(rendering_program_lines);
        uniform_attrib(viewer);
        glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
        glUseProgram(0);
        glBindVertexArray(0);

    }

    void changed()
    {
        compute_elements();
        initialize_buffers();
    }

private:
    static void gl_draw_edge(double px, double py, double pz,
                             double qx, double qy, double qz)
    {
        ::glVertex3d(px,py,pz);
        ::glVertex3d(qx,qy,qz);
    }

    std::vector<float> positions_lines;
    GLint location[2];
    GLuint vao[1];
    GLuint buffer[1];
    GLuint rendering_program_lines;

    void initialize_buffers()
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
        glShaderSource(vertex_shader, 1, vertex_shader_source, NULL);
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
        rendering_program_lines = program;
    }
    void uniform_attrib(Viewer_interface* viewer) const
    {
        GLfloat colors[3];
        GLfloat mvp_mat[16];

        //fills the MVP and MV matrices.

        GLdouble d_mat[16];
        viewer->camera()->getModelViewProjectionMatrix(d_mat);
        for (int i=0; i<16; ++i)
            mvp_mat[i] = GLfloat(d_mat[i]);


        //fills the arraw of colors with the current color
        colors[0] = this->color().redF();
        colors[1] = this->color().greenF();
        colors[2] = this->color().blueF();


        glUseProgram(rendering_program_lines);
        glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
        glUniform3fv(location[1],1,colors);
    }
    void compute_elements()
    {
        positions_lines.clear();
        const Bbox& bb = scene->bbox();
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);

        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);

        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);

        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);

        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);

        location[0] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
        location[1] = glGetUniformLocation(rendering_program_lines, "color");
    }

    const Scene_interface* scene;
};

#include "Polyhedron_demo_plugin_interface.h"

class Polyhedron_demo_trivial_plugin : 
        public QObject,
        public Polyhedron_demo_plugin_interface
{
    Q_OBJECT
    Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
    void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
    QList<QAction*> actions() const {
        return QList<QAction*>() << actionBbox;
    }

    bool applicable(QAction*) const {
        return true;
    }
public slots:

    void bbox();
    void enableAction();

private:
    Scene_interface* scene;
    QAction* actionBbox;


}; // end Polyhedron_demo_trivial_plugin

void Polyhedron_demo_trivial_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface)
{
    scene = scene_interface;
    actionBbox = new QAction(tr("Create bbox"), mainWindow);
    connect(actionBbox, SIGNAL(triggered()),
            this, SLOT(bbox()));
}

void Polyhedron_demo_trivial_plugin::bbox()
{
    for(int i = 0, end = scene->numberOfEntries();
        i < end; ++i)
    {
        if(qobject_cast<Scene_bbox_item*>(scene->item(i)))
            return;
    }
    Scene_item* item = new Scene_bbox_item(scene);
    connect(item, SIGNAL(destroyed()),
            this, SLOT(enableAction()));
    item->setName("Scene bbox");
    item->setColor(Qt::black);
    item->setRenderingMode(Wireframe);
    scene->addItem(item);
    actionBbox->setEnabled(false);
}

void Polyhedron_demo_trivial_plugin::enableAction() {
    actionBbox->setEnabled(true);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_trivial_plugin, Polyhedron_demo_trivial_plugin)

#include "Polyhedron_demo_trivial_plugin.moc"
