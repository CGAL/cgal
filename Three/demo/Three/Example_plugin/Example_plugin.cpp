//! \file Polyhedron_demo_example_plugin.cpp

#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_group_item.h>

//! [itemdeclaration]
// The special Scene_item only for triangles

//this is used by the Qt's MOC system to manage the metadata.
#ifdef scene_triangle_item_EXPORTS
#  define SCENE_TRIANGLE_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_TRIANGLE_ITEM_EXPORT Q_DECL_IMPORT
#endif

class Scene_triangle_item : public CGAL::Three::Scene_item
{

  Q_OBJECT
public :
  Scene_triangle_item(double ax,double ay, double az,
                      double bx,double by, double bz,
                      double cx,double cy, double cz);

  // Indicates if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Flat);
  }

  //Displays the item
  void draw(CGAL::Three::Viewer_interface* viewer) const;

  //Specifies that the buffers need to be initialized again.
  //Is mostly called after a change of geometry in the data.
  void invalidateOpenGLBuffers();

  //fills the std::vector
  void computeElements(double ax,double ay, double az,
                        double bx,double by, double bz,
                        double cx,double cy, double cz) const;

  Scene_item* clone() const {return 0;}
  QString toolTip() const {return QString();}


private:
  //contains the data
  mutable std::vector<float> vertices;
  mutable int nb_pos;
  mutable QOpenGLShaderProgram *program;
  using CGAL::Three::Scene_item::initializeBuffers;
  //Fills the buffers with data. The buffers allow us to give data to the shaders.
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const;
}; //end of class Scene_triangle_item
//! [itemdeclaration]
Scene_triangle_item::Scene_triangle_item(double ax,double ay, double az,
                                         double bx,double by, double bz,
                                         double cx,double cy, double cz)
  :  CGAL::Three::Scene_item(1,1)
{

  //Color is uniform, no need for a buffer. Changing the color will not re-compute the data
  is_monochrome = true;
  nb_pos = 0;
  are_buffers_filled = false;
  computeElements(ax, ay, az,
                   bx, by, bz,
                   cx, cy, cz);
  invalidateOpenGLBuffers();
}

//! [computeelements]
//Fills the position vector with data.
void Scene_triangle_item::computeElements(double ax, double ay, double az,
                                           double bx, double by, double bz,
                                           double cx, double cy, double cz)const
{
  vertices.resize(9);
  vertices[0] = ax; vertices[1] = ay; vertices[2] = az;
  vertices[3] = bx; vertices[4] = by; vertices[5] = bz;
  vertices[6] = cx; vertices[7] = cy; vertices[8] = cz;
}

//! [computeelements]
//! [draw]

void Scene_triangle_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  //The filling of the buffers should be performed in this function, because it needs a valid openGL context, and we are certain to have one in this function.
  if(!are_buffers_filled)
  {
    computeElements(0, 0, 0,
                     1, 0, 0,
                     0.5, 0.5, 0);
    initializeBuffers(viewer);
  }
  //Binds the vao corresponding to the type of data we are drawing.
  vaos[0]->bind();
  //Gets the program corresponding to the type of data we are drawing.
  //Here we want triangles with light effects.
  program = getShaderProgram(PROGRAM_WITH_LIGHT);
  //Gives most of the uniform values to the shaders.
  attribBuffers(viewer, PROGRAM_WITH_LIGHT);
  //Binds the program chosen before to use the right shaders.
  program->bind();
  //Gives the wanted color to the fragment shader as uniform value.
  program->setAttributeValue("colors", this->color());
  //Draws the items
  viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_pos/3));
  //clean up
  vaos[0]->release();
  program->release();

}
//! [draw]
//Specifies that the buffers need to be initialized again.
//Is mostly called after a change of geometry in the data.
void Scene_triangle_item::invalidateOpenGLBuffers()
{
  are_buffers_filled = false;
}


//! [fillbuffers]
void Scene_triangle_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
{

  //vao containing the data for the facets
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

//once the buffers are fill, we can empty the vectors to optimize memory consumption
  nb_pos = vertices.size();
  vertices.resize(0);
  //"Swap trick" insures that the memory is indeed freed and not kept available
  std::vector<float>(vertices).swap(vertices);
  are_buffers_filled = true;
}
//! [fillbuffers]
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
//The actual plugin
using namespace CGAL::Three;
class Q_DECL_EXPORT Polyhedron_demo_example_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  //Configures CMake to use MOC correctly
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")


public :
  // Adds an action to the menu and configures the widget
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    //get the references
    this->scene = scene_interface;
    this->mw = mainWindow;

    //creates and link the actions
    QAction* actionDrawTriangle= new QAction("Draw Triangle", mw);
    if(actionDrawTriangle) {
      connect(actionDrawTriangle, SIGNAL(triggered()),
              this, SLOT(draw_triangle()));
      _actions << actionDrawTriangle;
    }
  }
  bool applicable(QAction*) const
  {
    return true;
  }
  QList<QAction*> actions() const {
    return _actions;
  }

public Q_SLOTS:


  void draw_triangle() {

    triangle = new Scene_triangle_item(0, 0, 0,
                                       1, 0, 0,
                                       0.5, 0.5, 0);
    scene->addItem(triangle);
  }
private:
  CGAL::Three::Scene_item* triangle;
  QList<QAction*> _actions;

}; //end of class Polyhedron_demo_example_plugin
#include "Example_plugin.moc"

