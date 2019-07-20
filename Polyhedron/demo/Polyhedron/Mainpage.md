Main Page                        {#mainpage}
============

These pages are not documenting the whole Polyhedron demo but only the API that can be useful to create and add a new plugin.

Understanding the Polyhedron demo
============

There are several levels in this demo. 

- The MainWindow, which contains the UI elements. 

- Among these elements is the Viewer, which is the drawable surface that handles all the drawing and all the keyboard and mouse events. 

- The Viewer has a reference to the Scene, which contains the Scene_item list, which is a list of the drawn elements. 

A plugin usually defines an object that inherits from Scene_item or uses some of them to demonstrate a CGAL feature, so it might have to deal with the above elements.

Creating a simple Plugin
============
A basic plugin will inherit from Polyhedron_demo_plugin_interface. It can also inherits from the Polyhedron_demo_plugin_helper instead, for a more detailed model of plugin.
Its name must be of the form Polyhedron_demo_xxxx_yyyy_plugin. \n
<b>In the CMakeList.txt file, in the section Plugins, add the following lines :</b>

    polyhedron_demo_plugin(xxxx_yyyy_plugin Polyhedron_demo_xxxx_yyyy_plugin)
    target_link_libraries(xxxx_yyyy_plugin scene_polyhedron_item) 
  
  [init]: @ref Polyhedron_demo_plugin_helper#init(QMainWindow *, Scene_interface *)
  The class must contain the following lines :\n
  
    Q_OBJECT\n
    Q_INTERFACES(Polyhedron_demo_plugin_interface)\n
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")\n
    
In the function [init], get a reference to the Scene  and to the MainWindow. Then, create and link the actions of the plugin.\n
Create a list of QActions containing the actions of the plugin.\n
Add the following line:

    actionName->setProperty("submenuName", "Name_you_want_for_your_submenu");
    
to place your action in a submenu in the Operation Menu.\n
If the plugin implements a new Scene_item, please notice that a Scene_itam have a number of functions that will need a reference to the Viewer through the Viewer_interface type.

A plugin must always contain
~~~~~~~~~~~~~{.cpp}
#include "Polyhedron_demo_xxxx_yyyy_plugin.moc"
~~~~~~~~~~~~~

List of useful classes :
========================
- MainWindow
- Viewer_interface
- Scene_interface
- Scene_item
- Polyhedron_demo_plugin_helper
- Polyhedron_demo_plugin_interface


Example : 
============
The following code will create a plugin that adds an action to the MainWindow. This action is called "Draw Triangle" and adds a triangle to the scene.


~~~~~~~~~~~~~{.cpp}
#include <QApplication>
#include <QMainWindow>
#include <QAction>

#include  <CGAL/Three/Scene_item.h>
#include "Viewer_interface.h"
class Q_DECL_EXPORT Scene_triangle_item : public CGAL::Three::Scene_item
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
             QObject::tr( "<p>Number of vertices: %3<br />"
                           "Number of edges: %3<br />"
                         "Number of facets: %1")
                .arg(this->name())
                .arg(poly->size_of_vertices())
                .arg(poly->size_of_halfedges()/2)
                .arg(poly->size_of_facets())
                .arg(this->renderingModeName())
                .arg(this->color().name());
      if (volume!=-std::numeric_limits<double>::infinity())
        str+=QObject::tr("<br />Volume: %1").arg(volume);
      if (area!=-std::numeric_limits<double>::infinity())
        str+=QObject::tr("<br />Area: %1").arg(area);
      str+="</p>";
      item_text += QString("Bounding box: min (%1,%2,%3), max(%4,%5,%6)")
           .arg(item->bbox().xmin)
           .arg(item->bbox().ymin)
           .arg(item->bbox().zmin)
           .arg(item->bbox().xmax)
           .arg(item->bbox().ymax)
           .arg(item->bbox().zmax);
      m_text += QString("<br />Number of isolated vertices: 0<br />");

      return str;
    }

    void draw(Viewer_interface* viewer) const
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
    using CGAL::Three::Scene_item::initialize_buffers;
    void initialize_buffers(Viewer_interface *viewer)const
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
class Polyhedron_demo_example_plugin :
        public QObject,
        public Polyhedron_demo_plugin_helper
{
    //Configures CMake to use MOC correctly
    Q_OBJECT
    Q_INTERFACES(Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")


public :
    // To silent a warning -Woverloaded-virtual
    // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
    using Polyhedron_demo_plugin_helper::init;

    void init(QMainWindow* mainWindow,
              Scene_interface* scene_interface) {
      //get the references
      this->scene = scene_interface;
      this->mw = mainWindow;
      //creates and link the actions
      actionDrawTriangle= new QAction("Draw Triangle", mw);
      actionDrawTriangle->setProperty("subMenuName", "Object creation");
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

~~~~~~~~~~~~~
