#include "Messages_interface.h"
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <QByteArray>
#include <QThread>
#include <QMessageBox>
#include <QFileDialog>
#include <QFileInfo>
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include "ui_Animate_widget.h"
#include <fstream>

using namespace CGAL::Three;
typedef Viewer_interface Vi;


class DockWidget :
    public QDockWidget,
    public Ui::AnimateWidget
{
public:
  DockWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};



class Q_DECL_EXPORT Animate_mesh_plugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public :
  void init(QMainWindow* mw,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* mi) override {
    //get the references
    this->scene = scene_interface;
    this->messages = mi;
    this->mw = mw;
    actionAnimate = new QAction("Animate Surface Mesh", mw);
    actionAnimate->setProperty("subMenuName","Operations on Polyhedra");

    dock_widget = new DockWidget("Mesh Animation", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    addDockWidget(dock_widget);
    connect(actionAnimate, &QAction::triggered,
            this, &Animate_mesh_plugin::init_animation);

    connect(dock_widget->startButton, &QPushButton::clicked,
            this, &Animate_mesh_plugin::start_animation);

    connect(dock_widget->resetButton, &QPushButton::clicked,
            this, &Animate_mesh_plugin::reset_animation);

    connect(dock_widget->nextButton, &QPushButton::clicked,
            this, &Animate_mesh_plugin::next_frame);
    connect(dock_widget->prevButton, &QPushButton::clicked,
            this, &Animate_mesh_plugin::prev_frame);
    connect(dock_widget->stopButton, &QPushButton::clicked,
            this, [this](){
      timer.stop();
    });
    connect(&timer, SIGNAL (timeout()), this, SLOT (next_frame()));
    frame = -1;
    connect(dock_widget->helpButton, &QPushButton::clicked,this, [this](){
      QMessageBox::information(dock_widget, QString("Animation"),
                             QString("The TRJS format contains informations for a succession of modifications on a Surface Mesh. "
                                     "Such a modification is called a frame, and every frame is composed with a ligne for the "
                                     "number of points modified, and one ligne per modified point and its index.\n\n"
                                     "Example:\n\n"
                                     "n \n"
                                     "id1 x1 y1 z1 \n"
                                     "...           \n"
                                     "idn xn yn zn  \n"
                                     "m             \n"
                                     "id1 x1 y1 z1  \n"
                                     "...           \n"
                                     "idm xm ym zm"));
  });

  }

  bool applicable(QAction*) const override
  {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    if(!item)
      return false;
    return qobject_cast<Scene_surface_mesh_item*>(item);
  }

  QList<QAction*> actions() const override{
    return QList<QAction*>() << actionAnimate;
  }

public Q_SLOTS:

  void read_frame()
  {
    typedef boost::property_map<SMesh,CGAL::vertex_point_t>::type VPmap;
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();

    std::ifstream is(filepath.toUtf8());
    is.seekg(position);
    if(!is.good())
    {
      timer.stop();
      return;
    }
    int nb_verts;
    is >> nb_verts;
    if(!is.good())
    {
      timer.stop();
      return;
    }

    if(nb_verts != 0)
    {
      EPICK::FT x,y,z;
      for(int i = 0; i< nb_verts; ++i)
      {
        int id;
        is >> id;
        if(!is.good())
        {
          timer.stop();
          return;
        }
        is >> x;
        if(!is.good())
        {
          timer.stop();
          return;
        }
        is >> y;
        if(!is.good())
        {
          timer.stop();
          return;
        }
        is >> z;
        if(!is.good())
        {
          timer.stop();
          return;
        }
        if(id <= static_cast<int>(sm_item->face_graph()->number_of_vertices()))
        {
          SMesh::Vertex_index vh(id);
          VPmap vpm = get(CGAL::vertex_point, *sm_item->face_graph());
          put(vpm, vh, Point_3(x-offset.x,
                               y-offset.y,
                               z-offset.z));
          sm_item->updateVertex(vh);
        }
      }
    }

    sm_item->redraw();
    if(!is.good())
    {
      timer.stop();
      return;
    }
  }

  void next_frame()
  {
    if(frame == static_cast<int>(frame_pos.size()-1))
    {
      timer.stop();
      return;
    }
    position = frame_pos[++frame];
    dock_widget->frameSlider->setValue(frame);
    dock_widget->frameLabel->setText(QString("%1/%2").arg(frame).arg(frame_pos.size()-1));
    read_frame();
  }

  void prev_frame()
  {
    if(frame <= 0)
      return;
    position = frame_pos[--frame];
    dock_widget->frameSlider->setValue(frame);
    dock_widget->frameLabel->setText(QString("%1/%2").arg(frame).arg(frame_pos.size()-1));
    read_frame();
  }

  void start_animation()
  {
    timer.start(dock_widget->frame_time->value());
  }

  void reset_animation()
  {
    timer.stop();
    frame=-1;
    dock_widget->frameSlider->setValue(frame);
    dock_widget->frameLabel->setText(QString("%1/%2").arg(frame).arg(frame_pos.size()-1));
    for(std::size_t id = 0; id<initial_points.size();++id)
    {
      sm_item->face_graph()->points()[SMesh::Vertex_index(
            static_cast<CGAL::SM_Index<CGAL::SM_Vertex_index>::size_type>(id))]//not size_t on windows.
          =initial_points[id];
    }
    sm_item->invalidateOpenGLBuffers();
    sm_item->redraw();
  }

  void init_animation()
  {
    clean_up();
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    if(!item)
      return;
    sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    connect(sm_item, &Scene_surface_mesh_item::aboutToBeDestroyed, this,
            [this](){
      clean_up();
    });
    initial_points.reserve(sm_item->face_graph()->number_of_vertices());
    for(const auto& p : sm_item->face_graph()->points())
      initial_points.push_back(p);

    //ask the file (*.trjs)
    QFileDialog dialog(mw);
    QStringList filters;
    filters << "Trajectories file(*.trjs)";
    dialog.setNameFilters(filters);
    dialog.setFileMode(QFileDialog::ExistingFiles);
    if(dialog.exec() != QDialog::Accepted) { return; }
    QFileInfo info(dialog.selectedFiles().first());
    filepath = info.filePath();
    position = 0;
    dock_widget->resetButton->setEnabled(true);
    dock_widget->startButton->setEnabled(true);
    dock_widget->prevButton->setEnabled(true);
    dock_widget->nextButton->setEnabled(true);
    dock_widget->stopButton->setEnabled(true);
    dock_widget->frameSlider->setEnabled(true);
    dock_widget->setWindowTitle(QString("Animate %1").arg(sm_item->name()));
    if(!info.exists())
    {
      QMessageBox::warning(mw, "Error","File does not exist.");
      return;
    }

    if(!dock_widget->isVisible()) { dock_widget->show(); }

    //pre-process to count frames
    QApplication::setOverrideCursor(Qt::WaitCursor);
    std::ifstream is(filepath.toUtf8());
    frame_pos.clear();
    while(is.good())
    {
      std::string line;
      std::streampos pos = is.tellg();
      std::getline(is, line);
      if(line.length() > 0 && line.find(" ") == std::string::npos)
        frame_pos.push_back(pos);
    }
    is.close();
    dock_widget->frameSlider->setMaximum(static_cast<int>(frame_pos.size())-1);
    dock_widget->frameLabel->setText(QString("%1/%2").arg(frame).arg(frame_pos.size()-1));

    connect(dock_widget->frameSlider, &QSlider::sliderMoved, [this](int i){
      frame = i;
      position = frame_pos[frame];
      dock_widget->frameLabel->setText(QString("%1/%2").arg(frame).arg(frame_pos.size()-1));
      read_frame();
    });

    QApplication::restoreOverrideCursor();
  }

  void closure() override
  {
    dock_widget->hide();
  }

  void clean_up()
  {
    sm_item=nullptr;
    position=0;
    frame_pos.clear();
    filepath="";
    frame = -1;
    dock_widget->resetButton->setEnabled(false);
    dock_widget->startButton->setEnabled(false);
    dock_widget->prevButton->setEnabled(false);
    dock_widget->nextButton->setEnabled(false);
    dock_widget->stopButton->setEnabled(false);
    dock_widget->frameSlider->setEnabled(false);
    dock_widget->setWindowTitle(QString("Animate"));
    initial_points.clear();
    initial_points.shrink_to_fit();
  }

private:
  QAction* actionAnimate;
  Messages_interface* messages;
  Scene_interface* scene;
  DockWidget* dock_widget;
  QString filepath;
  std::streampos position;
  Scene_surface_mesh_item* sm_item;
  std::vector<std::streampos> frame_pos;
  QTimer timer;
  int frame;
  std::vector<SMesh::Point> initial_points;
}; //end of plugin class
#include "Animate_mesh_plugin.moc"
