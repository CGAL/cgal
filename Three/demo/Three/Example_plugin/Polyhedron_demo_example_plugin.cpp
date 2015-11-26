//! \file Polyhedron_demo_example_plugin.cpp
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_group_item.h>
#include "ui_Polyhedron_demo_example_plugin.h"
#include "ui_dock_example.h"

#ifdef scene_triangle_item_EXPORTS
#  define SCENE_TRIANGLE_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_TRIANGLE_ITEM_EXPORT Q_DECL_IMPORT
#endif

//! The dialog class to get the coordinates of the triangle
class Polyhedron_demo_example_plugin_dialog : public QDialog, private Ui::Polyhedron_demo_example_plugin_dialog
{
  Q_OBJECT
  public:
    Polyhedron_demo_example_plugin_dialog(QWidget* /*parent*/ = 0)
    {
      setupUi(this);

    }

    double getAx() const { return doubleSpinBox_Ax->value(); }
    double getAy() const { return doubleSpinBox_Ay->value(); }
    double getAz() const { return doubleSpinBox_Az->value(); }

    double getBx() const { return doubleSpinBox_Bx->value(); }
    double getBy() const { return doubleSpinBox_By->value(); }
    double getBz() const { return doubleSpinBox_Bz->value(); }

    double getCx() const { return doubleSpinBox_Cx->value(); }
    double getCy() const { return doubleSpinBox_Cy->value(); }
    double getCz() const { return doubleSpinBox_Cz->value(); }
    void setTriangleName(QString name) {triangle_name_label->setText(name); }
};

//! The special Scene_item only for triangles
class SCENE_TRIANGLE_ITEM_EXPORT Scene_triangle_item : public CGAL::Three::Scene_item
{

    Q_OBJECT
public :
    Scene_triangle_item(double ax,double ay, double az,
                        double bx,double by, double bz,
                        double cx,double cy, double cz)
        :  CGAL::Three::Scene_item(1,1)
    {
        vertices.resize(9);
        vertices[0] = ax; vertices[1] = ay; vertices[2] = az;
        vertices[3] = bx; vertices[4] = by; vertices[5] = bz;
        vertices[6] = cx; vertices[7] = cy; vertices[8] = cz;
        changed();

    }
    Scene_triangle_item()
        :  CGAL::Three::Scene_item(1,1)
    {

        Scene_triangle_item(
                    0,0,0,
                    0.5,1,0,
                    1,0,0);

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
        are_buffers_filled = false;
    }

private:

    std::vector<float> vertices;
    mutable QOpenGLShaderProgram *program;
    using CGAL::Three::Scene_item::initialize_buffers;
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

}; //end of class Scene_triangle_item

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
//!The actual plugin
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
    //! Adds an action to the menu and configures the widget
    void init(QMainWindow* mainWindow,
              CGAL::Three::Scene_interface* scene_interface) {
      //get the references
      this->scene = scene_interface;
      this->mw = mainWindow;
        triangle_array.resize(0);

      //creates and link the actions
      actionDrawTriangle= new QAction("Draw Triangle", mw);
      if(actionDrawTriangle) {
        connect(actionDrawTriangle, SIGNAL(triggered()),
                this, SLOT(populate_dialog()));
      }

      //Dock widget initialization
      dock_widget = new QDockWidget("Triangle Creation", mw);
      dock_widget->setVisible(false); // do not show at the beginning

      ui_widget.setupUi(dock_widget);
      mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);
      connect(ui_widget.buttonBox, SIGNAL(accepted()), this, SLOT(on_dock_OK_clicked()));
    }

    bool applicable(QAction*) const
    {
        return true;
    }
    QList<QAction*> actions() const {
      return QList<QAction*>() << actionDrawTriangle;
    }

  public Q_SLOTS:

    void populate_dialog()
    {
        // dock widget should be constructed in init()
        if(dock_widget->isVisible()) { dock_widget->hide(); }
        else                         { dock_widget->show(); }
        dialog = new Polyhedron_demo_example_plugin_dialog();
        dialog->setTriangleName(QString("Triangle %1").arg(triangle_array.size()));
        if(!dialog->exec())
            return;
        ax = dialog->getAx();
        ay = dialog->getAy();
        az = dialog->getAz();
        bx = dialog->getBx();
        by = dialog->getBy();
        bz = dialog->getBz();
        cx = dialog->getCx();
        cy = dialog->getCy();
        cz = dialog->getCz();
        triangle_array.push_back(new Scene_triangle_item(ax, ay, az,
                                                               bx, by, bz,
                                                               cx, cy, cz));
        if(triangle_array.size()%3 > 0)
        {
            scene->addItem(triangle_array[triangle_array.size()-1]);
            populate_dialog();
        }
        else
        {
            Scene_group_item *group = new Scene_group_item("triangle group");
            for(int i =0; i<3; i++)
                scene->changeGroup(triangle_array[triangle_array.size()-i-1], group);
            scene->addItem(group);
        }
    }

  void draw_triangle() {

    triangle = new Scene_triangle_item(ax, ay, az,
                                       bx, by, bz,
                                       cx, cy, cz);
    scene->addItem(triangle);
  }

  void on_dock_OK_clicked()
  {
      ax = ui_widget.doubleSpinBox_Ax->value();
      ay = ui_widget.doubleSpinBox_Ay->value();
      az = ui_widget.doubleSpinBox_Az->value();

      bx = ui_widget.doubleSpinBox_Bx->value();
      by = ui_widget.doubleSpinBox_By->value();
      bz = ui_widget.doubleSpinBox_Bz->value();

      cx = ui_widget.doubleSpinBox_Cx->value();
      cy = ui_widget.doubleSpinBox_Cy->value();
      cz = ui_widget.doubleSpinBox_Cz->value();
      draw_triangle();

  }

private:
  CGAL::Three::Scene_item* triangle;
  QAction* actionDrawTriangle;
  Polyhedron_demo_example_plugin_dialog *dialog;
  double ax;
  double ay;
  double az;
  double bx;
  double by;
  double bz;
  double cx;
  double cy;
  double cz;
  QVector<Scene_triangle_item*> triangle_array;
Ui::Dock_example ui_widget;
QDockWidget* dock_widget;

}; //end of class Polyhedron_demo_example_plugin
#include "Polyhedron_demo_example_plugin.moc"
