//! \file Polyhedron_demo_example_plugin.cpp
#include <QApplication>
#include <QMainWindow>
#include <QAction>

#include "Scene_item.h"
#include <CGAL/Three/Viewer_interface.h>

class Q_DECL_EXPORT Scene_triangle_item : public Scene_item
{

    Q_OBJECT
public :
    Scene_triangle_item()
        :  Scene_item(1,1)
    {

        vertices.resize(0);
        changed();

    }
    ~Scene_triangle_item()
    {
    }
    bool isFinite() const { return true; }
    bool isEmpty() const { return true; }
    Bbox bbox() const { return Bbox(); }

    Scene_triangle_item* clone() const {
        return 0;
    }

    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const {
        return (m == Flat);
    }

    QString toolTip() const {    
      QString str =
             QObject::tr( "<p>Number of vertices: %1<br />"
                           "Number of edges: %2<br />"
                         "Number of facets: %3")
                .arg(this->name())
                .arg(3)
                .arg(3)
                .arg(1)
                .arg(this->renderingModeName())
                .arg(this->color().name());
        str+=QObject::tr("<br />Volume: %1").arg(0);
        str+=QObject::tr("<br />Area: %1").arg(0.5);
      str+="</p>";
      str += QString("Bounding box: min (%1,%2,%3), max(%4,%5,%6)")
           .arg(bbox().xmin)
           .arg(bbox().ymin)
           .arg(bbox().zmin)
           .arg(bbox().xmax)
           .arg(bbox().ymax)
           .arg(bbox().zmax);
      str += QString("<br />Number of isolated vertices : 0<br />");

      return str;
    }

    void draw(CGAL::Three::Viewer_interface* viewer) const
    {
        if(!are_buffers_filled)
            initialize_buffers(viewer);
        vaos[0]->bind();
        program = getShaderProgram(PROGRAM_WITH_LIGHT);
        attrib_buffers(viewer, PROGRAM_WITH_LIGHT);
        program->bind();
        program->setAttributeValue("colors", this->color());
        viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()/3));
        vaos[0]->release();
        program->release();

    }

    void changed()
    {
        compute_elements();
        are_buffers_filled = false;
    }

private:

    std::vector<float> vertices;
    mutable QOpenGLShaderProgram *program;
    using Scene_item::initialize_buffers;
    void initialize_buffers(CGAL::Three::Viewer_interface *viewer)const
    {

        //vao containing the data for the lines
        {
            program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
            program->bind();

            vaos[0]->bind();
            buffers[0].bind();
            buffers[0].allocate(vertices.data(),
                                static_cast<GLsizei>(vertices.size()*sizeof(float)));
            program->enableAttributeArray("vertex");
            program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
            buffers[0].release();

            vaos[0]->release();
            program->release();

        }
        are_buffers_filled = true;
    }

    void compute_elements()
    {
        vertices.resize(9);
        vertices[0] = 0.0; vertices[1] = 0.0; vertices[2] = 0.0;
        vertices[3] = 0.5; vertices[4] = 1.0; vertices[5] = 0.0;
        vertices[6] = 1.0; vertices[7] = 0.0; vertices[8] = 0.0;
    }

}; //end of class Scene_triangle_item

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
using namespace CGAL::Three;
class Polyhedron_demo_example_plugin :
        public QObject,
        public Polyhedron_demo_plugin_helper
{
    //Configures CMake to use MOC correctly
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")


public :
    // To silent a warning -Woverloaded-virtual
    // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
    using Polyhedron_demo_plugin_helper::init;

    void init(QMainWindow* mainWindow,
              CGAL::Three::Scene_interface* scene_interface) {
      //get the references
      this->scene = scene_interface;
      this->mw = mainWindow;
      //creates and link the actions
      actionDrawTriangle= new QAction("Draw Triangle", mw);
      if(actionDrawTriangle) {
        connect(actionDrawTriangle, SIGNAL(triggered()),
                this, SLOT(draw_triangle()));
      }
    }

    bool applicable(QAction*) const
    {
        return true;
    }
    QList<QAction*> actions() const {
      return QList<QAction*>() << actionDrawTriangle;
    }

  public Q_SLOTS:

  void draw_triangle() {
    triangle = new Scene_triangle_item();
    scene->addItem(triangle);
  }

private:
  Scene_item* triangle;
  QAction* actionDrawTriangle;

}; //end of class Polyhedron_demo_example_plugin

#include "Polyhedron_demo_example_plugin.moc"
