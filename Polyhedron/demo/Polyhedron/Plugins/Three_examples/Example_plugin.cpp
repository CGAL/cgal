//! \file Polyhedron_demo_example_plugin.cpp

#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_group_item.h>

#include <CGAL/Three/Triangle_container.h>

//! [itemdeclaration]
// The special Scene_item only for triangles

using namespace CGAL::Three;
typedef Triangle_container Tc;
typedef Viewer_interface VI;


class Scene_triangle_item
    : public CGAL::Three::Scene_item_rendering_helper
{

  Q_OBJECT
public :
  Scene_triangle_item(double ax,double ay, double az,
                      double bx,double by, double bz,
                      double cx,double cy, double cz);

  // Indicates if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE {
    return (m == Flat);
  }

  //Displays the item
  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;

  //Specifies that the buffers need to be initialized again.
  //Is mostly called after a change of geometry in the data.
  void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;

  //fills the std::vector
  void computeElements() const Q_DECL_OVERRIDE;

  Scene_item* clone() const Q_DECL_OVERRIDE {return 0;}
  QString toolTip() const Q_DECL_OVERRIDE {return QString();}
  void compute_bbox() const Q_DECL_OVERRIDE { _bbox = Bbox(); }

private:
  //contains the data
  mutable std::vector<float> vertices;
  mutable std::size_t nb_pos;
  //Fills the buffers with data. The buffers allow us to give data to the shaders.
  void initializeBuffers(Viewer_interface *viewer) const Q_DECL_OVERRIDE;
}; //end of class Scene_triangle_item
//! [itemdeclaration]
Scene_triangle_item::Scene_triangle_item(double ax,double ay, double az,
                                         double bx,double by, double bz,
                                         double cx,double cy, double cz)
{
  //! [creation]
  //Prepare a single TriangleContainer, as we will only use one program.
  setTriangleContainer(0, new Tc(VI::PROGRAM_NO_SELECTION,
                                 false));
  //! [creation]
//! [computeelements]
  //Fills the position vector with data.
  nb_pos = 0;
  vertices.resize(9);
  vertices[0] = ax; vertices[1] = ay; vertices[2] = az;
  vertices[3] = bx; vertices[4] = by; vertices[5] = bz;
  vertices[6] = cx; vertices[7] = cy; vertices[8] = cz;
  nb_pos=vertices.size();
//! [computeelements]
  //be sure the data will be computed next draw call
  invalidateOpenGLBuffers();
}


//prepare the TriangleContainer with the computed data.
void Scene_triangle_item::computeElements()const
{
//! [allocateelements]
getTriangleContainer(0)->allocate(Tc::Flat_vertices, vertices.data(),
                                  static_cast<int>(vertices.size()
                                                   * sizeof(float)));
setBuffersFilled(true);
//! [allocateelements]
}

//! [draw]

void Scene_triangle_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  //Initializes the OpenGL context for `viewer` if needed.
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }

  //set the uniform properties for the TriangleContainer.
  //Uniform values are setted at each draw call and are defined for the whole item.
  //Values per simplex are computed as buffers in ComputeElements() and bound in initializeBuffers().
  getTriangleContainer(0)->setColor(this->color());
  getTriangleContainer(0)->draw(viewer, true);

}
//! [draw]
//Specifies that the buffers need to be initialized again.
//Is mostly called after a change of geometry in the data.
void Scene_triangle_item::invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getTriangleContainer(0)->reset_vbos(ALL);
}


//! [fillbuffers]
void Scene_triangle_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
{
//Bind the buffers filled in ComputeElements() for the TriangleContainer.
  getTriangleContainer(0)->initializeBuffers(viewer);
  getTriangleContainer(0)->setFlatDataSize(nb_pos);

//once the buffers are bound, we can clear the vectors to optimize memory consumption
  vertices.clear();
  vertices.shrink_to_fit();
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
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) Q_DECL_OVERRIDE{
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
  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    return true;
  }
  QList<QAction*> actions() const Q_DECL_OVERRIDE{
    return _actions;
  }

public Q_SLOTS:


  void draw_triangle() {

    triangle = new Scene_triangle_item(0, 0, 0,
                                       1, 0, 0,
                                       0.5, 0.5, 0);
    triangle->setName(QString("Basic triangle"));
    scene->addItem(triangle);
  }
private:
  CGAL::Three::Scene_item* triangle;
  QList<QAction*> _actions;

}; //end of class Polyhedron_demo_example_plugin
#include "Example_plugin.moc"

