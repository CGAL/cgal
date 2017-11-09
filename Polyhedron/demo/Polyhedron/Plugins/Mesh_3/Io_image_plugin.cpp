#ifdef _MSC_VER
#  pragma warning(disable:4244) // conversion with loss of data
#endif

#include "Volume_plane.h"
#include "Volume_plane_thread.h"
#include "Volume_plane_intersection.h"
#include "Messages_interface.h"
#include "config.h"
#include "Scene_image_item.h"
#include "Image_type.h"
#include "ui_Image_res_dialog.h"

#include <CGAL/Image_3.h>
#include <CGAL/ImageIO.h>
#include <CGAL/read_sep_image_data.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/config.h>
#include <CGAL/use.h>
#include <QAction>
#include <QMenu>
#include <QMouseEvent>
#include <QList>
#include <QInputDialog>
#include <QSlider>
#include <QLabel>
#include <QLineEdit>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QDockWidget>
#include <QMainWindow>
#include <QMessageBox>
#include <QString>
#include <QFontMetrics>
#include <QFileDialog>
#include <QPushButton>

#include <cassert>
#include <iostream>

#include <boost/type_traits.hpp>
#include <boost/optional.hpp>
#include "Scene.h"

#include <QSettings>
#include <QUrl>
#include "Raw_image_dialog.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <fstream>
#include <cstdlib>
#ifdef CGAL_USE_VTK
#include <CGAL/read_vtk_image_data.h>

#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkDemandDrivenPipeline.h>
#endif

// Covariant return types don't work for scalar types and we cannot
// have templates here, hence this unfortunate hack.

// The input float value we are reading is always in
// 0..1 and min_max is the range it came from.
struct IntConverter {
  std::pair<int, int> min_max;

  int operator()(float f) {
    float s = f * (min_max.second - min_max.first);
    return s + min_max.first;
  }
};

struct DoubleConverter {
  std::pair<float, float> min_max;

  float operator()(float f) {
    float s = f * (min_max.second - min_max.first);
    return s + min_max.first;
  }
};

class PixelReader : public QObject
{
Q_OBJECT
public Q_SLOTS:
  void update(const QMouseEvent *e) {
    getPixel(e->pos());
  }
Q_SIGNALS:
  void x(QString);

public:
  void setIC(const IntConverter& x) { ic = x; fc = boost::optional<DoubleConverter>(); }
  void setFC(const DoubleConverter& x) { fc = x; ic = boost::optional<IntConverter>(); }
  void setViewer(Viewer_interface* viewer) { this->viewer = viewer; }

private:
  boost::optional<IntConverter> ic;
  boost::optional<DoubleConverter> fc;
  Viewer_interface* viewer;
  void getPixel(const QPoint& e) {
    float data[3];
    int vp[4];
    viewer->glGetIntegerv(GL_VIEWPORT, vp);
    viewer->glReadPixels(e.x(), vp[3] - e.y(), 1, 1, GL_RGB, GL_FLOAT, data);

    if(fc) {
      Q_EMIT x(QString::number((*fc)(data[0]), 'f', 6 ));
    } else if(ic) {
      Q_EMIT x( QString::number((*ic)(data[0]) ));
    }
  }
};


class Plane_slider : public QSlider
{
  Q_OBJECT
public:
  Plane_slider(const qglviewer::Vec& v, int id, Scene_interface* scene,
               qglviewer::ManipulatedFrame* frame, Qt::Orientation ori, QWidget* widget)
    : QSlider(ori, widget), v(v), id(id), scene(scene), frame(frame) {
    this->setTracking(true);
    connect(frame,  SIGNAL(manipulated()), this, SLOT(updateCutPlane()));
  }

public Q_SLOTS:
  void updateCutPlane()
  {
     ready_to_cut = true;
     QTimer::singleShot(0,this,SLOT(updateValue()));
  }

  void setFramePosition()
  {
    if(!ready_to_move)
      return;
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    qglviewer::Vec v2 = v * (this->value() / scale);
    v2+=offset;
    frame->setTranslationWithConstraint(v2);
    scene->itemChanged(id);

    Q_EMIT realChange(this->value() / scale);
    ready_to_move = false;
  }
  void updateValue() {
    if(!ready_to_cut)
      return;
#if QGLVIEWER_VERSION >= 0x020600
    typedef qreal qglviewer_real;
#else // QGLViewer < 2.6.0
    typedef float qglviewer_real;
#endif // QGLViewer < 2.6.0
    qglviewer_real a, b, c;
    frame->getPosition(a, b, c);
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    a-=offset.x;
    b-=offset.y;
    c-=offset.z;
    float sum1 = float(a + b + c);
    float sum2 = float(v.x + v.y + v.z);
    sum1 /= sum2;
    setValue(sum1 * scale);
    ready_to_cut = false;
  }

  void updateFramePosition()
  {
    ready_to_move = true;
    QTimer::singleShot(0,this,SLOT(setFramePosition()));
  }
  unsigned int getScale() const { return scale; }
Q_SIGNALS:
  void realChange(int);

private:
  static const unsigned int scale;
  bool ready_to_cut;
  bool ready_to_move;
  qglviewer::Vec v;
  int id;
  Scene_interface* scene;
  qglviewer::ManipulatedFrame* frame;
};

const unsigned int Plane_slider::scale = 100;

class Io_image_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_helper,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:

  bool applicable(QAction*) const {
    return qobject_cast<Scene_image_item*>(scene->item(scene->mainSelectionIndex()));
  }


  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface *mi) {
    this->message_interface = mi;
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->is_gray = false;
    x_control = NULL;
    y_control = NULL;
    z_control = NULL;
    current_control = NULL;
    planeSwitch = new QAction("Add Volume Planes", mw);
    QAction *actionLoadDCM = new QAction("Open directory", mw);
    connect(actionLoadDCM, SIGNAL(triggered()), this, SLOT(on_actionLoadDCM_triggered()));
    if(planeSwitch) {
      planeSwitch->setProperty("subMenuName", "3D Mesh Generation");
      connect(planeSwitch, SIGNAL(triggered()),
              this, SLOT(selectPlanes()));
      Scene* true_scene = dynamic_cast<Scene*>(scene);
      connect(true_scene,SIGNAL(itemIndexSelected(int)),
              this, SLOT(connect_controls(int)));
    }
    Viewer_interface* v = mw->findChild<Viewer_interface*>("viewer");
    CGAL_assertion(v != 0);
    pxr_.setViewer(v);
    connect(v, SIGNAL(pointSelected(const QMouseEvent *)), &pxr_, SLOT(update(const QMouseEvent *)));
    createOrGetDockLayout();

    QMenu* menuFile = mw->findChild<QMenu*>("menuFile");
    if ( NULL != menuFile )
    {
      QList<QAction*> menuFileActions = menuFile->actions();

      // Look for action just after "Load..." action
      QAction* actionAfterLoad = NULL;
      for ( QList<QAction*>::iterator it_action = menuFileActions.begin(),
           end = menuFileActions.end() ; it_action != end ; ++ it_action ) //Q_FOREACH( QAction* action, menuFileActions)
      {
        if ( NULL != *it_action && (*it_action)->text().contains("Load Plugin") )
        {
          ++it_action;
          if ( it_action != end && NULL != *it_action )
          {
            actionAfterLoad = *it_action;
          }
        }
      }

      // Insert "Load implicit function" action
      if ( NULL != actionAfterLoad )
      {
        menuFile->insertAction(actionAfterLoad,actionLoadDCM);
      }
    }
  }
  QList<QAction*> actions() const {
    return QList<QAction*>() << planeSwitch;
  }
  virtual void closure()
  {
      QDockWidget* controlDockWidget = mw->findChild<QDockWidget*>("volumePlanesControl");
      if(controlDockWidget)
          controlDockWidget->hide();
  }
  Io_image_plugin() : planeSwitch(NULL) {}

  QString nameFilters() const;
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item* item, QFileInfo fi) {
    const Scene_image_item* im_item = qobject_cast<const Scene_image_item*>(item);

    point_image p_im = *im_item->image()->image();
    return _writeImage(&p_im, fi.filePath().toUtf8()) == 0;
  }
  QString name() const { return "segmented images"; }


public Q_SLOTS:

  void setXNum(int i)
  {
   x_cubeLabel->setText(QString("%1").arg(i));
  }
  void setYNum(int i)
  {
   y_cubeLabel->setText(QString("%1").arg(i));
  }
  void setZNum(int i)
  {
   z_cubeLabel->setText(QString("%1").arg(i));
  }

  void XTextEdited(QString s)
  {
    int i = s.toInt();
    x_slider->setValue(i*qobject_cast<Plane_slider*>(x_slider)->getScale());
    x_slider->sliderMoved(i);
  }
  void YTextEdited(QString s)
  {
    int i = s.toInt();
    y_slider->setValue(i*qobject_cast<Plane_slider*>(y_slider)->getScale());
    y_slider->sliderMoved(i);
  }
  void ZTextEdited(QString s)
  {
    int i = s.toInt();
    z_slider->setValue(i*qobject_cast<Plane_slider*>(z_slider)->getScale());
    z_slider->sliderMoved(i);

  }

  void on_imageType_changed(int index)
  {
    if(index == 0)
      ui.groupBox->setVisible(true);
    else
      ui.groupBox->setVisible(false);
  }
  void selectPlanes() {

    std::vector< Scene_image_item* > seg_items;
    Scene_image_item* seg_img;
    seg_img = NULL;
    for(int i = 0; i < scene->numberOfEntries(); ++i) {
      Scene_image_item* tmp = qobject_cast<Scene_image_item*>(scene->item(i));
      if(tmp != NULL){
        seg_items.push_back(tmp);
      }
    }
    if(seg_items.empty()) {
      QMessageBox::warning(mw, tr("No suitable item found"), tr("Load an inrimage or hdr file to enable Volume Planes."));
      return;
    } else {
      QList<QString> items;
      for(std::vector< Scene_image_item* >::const_iterator it = seg_items.begin();
          it != seg_items.end(); ++it) {
        items << (*it)->name();
      }
      bool ok;
      QString selected = QInputDialog::getItem(mw, tr("Select a dataset:"), tr("Items"), items, 0, false, &ok);
      if(!ok || selected.isEmpty())
        return;
      for(std::vector< Scene_image_item*>::const_iterator it = seg_items.begin();
          it != seg_items.end(); ++it) {
        if(selected == (*it)->name())
          seg_img = *it;
      }
    }
    if(group_map.keys().contains(seg_img))
      this->message_interface->warning("This item already has volume planes.");
    else
    {
      // Opens a modal Dialog to prevent the user from manipulating things that could mess with the planes creation and cause a segfault.
      msgBox.setText("Planes created : 0/3");
      msgBox.setStandardButtons(QMessageBox::NoButton);
      msgBox.show();
      createPlanes(seg_img);
    }
  }

  void addVP(Volume_plane_thread* thread) {
    Volume_plane_interface* plane = thread->getItem();
    plane->init();
    // add the interface for this Volume_plane
    int id = scene->addItem(plane);
    scene->changeGroup(plane, group);
    group->lockChild(plane);
    switch(thread->type())
    {
    case 'x':
      delete x_slider;
      x_slider = new Plane_slider(plane->translationVector(), id, scene, plane->manipulatedFrame(),
                                         Qt::Horizontal, x_control);
      x_slider->setRange(0, (plane->cDim() - 1) * 100);
      connect(x_slider, SIGNAL(realChange(int)), this, SLOT(setXNum(int)));
      connect(x_cubeLabel, SIGNAL(textEdited(QString)), this, SLOT(XTextEdited(QString)));
      connect(x_slider, SIGNAL(realChange(int)), this, SLOT(set_value()));
      connect(x_slider, SIGNAL(sliderMoved(int)), x_slider, SLOT(updateFramePosition()));
      connect(plane, SIGNAL(manipulated(int)), this, SLOT(setXNum(int)));
      connect(plane, SIGNAL(aboutToBeDestroyed()), this, SLOT(destroy_x_item()));
      x_box->addWidget(x_slider);
      x_box->addWidget(x_cubeLabel);
      x_control->setVisible(true);
      break;
    case 'y':
      delete y_slider;
      y_slider = new Plane_slider(plane->translationVector(), id, scene, plane->manipulatedFrame(),
                                         Qt::Horizontal, y_control);
      y_slider->setRange(0, (plane->cDim() - 1) * 100);
      connect(y_slider, SIGNAL(realChange(int)), this, SLOT(setYNum(int)));
      connect(y_cubeLabel, SIGNAL(textEdited(QString)), this, SLOT(YTextEdited(QString)));
      connect(y_slider, SIGNAL(realChange(int)), this, SLOT(set_value()));
      connect(y_slider, SIGNAL(sliderMoved(int)), y_slider, SLOT(updateFramePosition()));
      connect(plane, SIGNAL(manipulated(int)), this, SLOT(setYNum(int)));
      connect(plane, SIGNAL(aboutToBeDestroyed()), this, SLOT(destroy_y_item()));
      y_box->addWidget(y_slider);
      y_box->addWidget(y_cubeLabel);
      y_control->setVisible(true);
      break;
    case 'z':
      delete z_slider;
      z_slider = new Plane_slider(plane->translationVector(), id, scene, plane->manipulatedFrame(),
                                         Qt::Horizontal, z_control);
      z_slider->setRange(0, (plane->cDim() - 1) * 100);
      connect(z_slider, SIGNAL(realChange(int)), this, SLOT(setZNum(int)));
      connect(z_cubeLabel, SIGNAL(textEdited(QString)), this, SLOT(ZTextEdited(QString)));
      connect(z_slider, SIGNAL(realChange(int)), this, SLOT(set_value()));
      connect(z_slider, SIGNAL(sliderMoved(int)), z_slider, SLOT(updateFramePosition()));
      connect(plane, SIGNAL(manipulated(int)), this, SLOT(setZNum(int)));
      connect(plane, SIGNAL(aboutToBeDestroyed()), this, SLOT(destroy_z_item()));
      z_box->addWidget(z_slider);
      z_box->addWidget(z_cubeLabel);
      z_control->setVisible(true);
      break;
    default:
      break;
    }
    std::vector<Volume_plane_thread*>::iterator it = std::find(threads.begin(), threads.end(), thread);

    // this slot has been connected to a thread that hasn't been
    // registered here.
    assert(it != threads.end());
    delete *it;
    threads.erase(it);

    update_msgBox();
    Volume_plane_intersection* intersection = dynamic_cast<Volume_plane_intersection*>(scene->item(intersection_id));
    if(!intersection) {
      // the intersection is gone before it was initialized
      return;
    }
    // FIXME downcasting mode
    // FIXME this will bug if two volume planes are generated simultaneously by the plugin
    if(Volume_plane<x_tag>* p = dynamic_cast< Volume_plane<x_tag>* >(plane)) {
      intersection->setX(p);
    } else if(Volume_plane<y_tag>* p = dynamic_cast< Volume_plane<y_tag>* >(plane)) {
      intersection->setY(p);
    } else if(Volume_plane<z_tag>* p = dynamic_cast< Volume_plane<z_tag>* >(plane)) {
      intersection->setZ(p);
    }
    connect(plane, SIGNAL(planeDestructionIncoming(Volume_plane_interface*)),
            intersection, SLOT(planeRemoved(Volume_plane_interface*)));

  }

  void on_actionLoadDCM_triggered()
  {
    QSettings settings;
    QString start_dir = settings.value("Open directory",
                                       QDir::current().dirName()).toString();
    QString dir =
      QFileDialog::getExistingDirectory(mw,
                                        tr("Open directory"),
                                        start_dir,
                                        QFileDialog::ShowDirsOnly
                                        | QFileDialog::DontResolveSymlinks);

    if (!dir.isEmpty()) {
      QFileInfo fileinfo(dir);
      if (fileinfo.isDir() && fileinfo.isReadable())
      {
        settings.setValue("Open directory",
          fileinfo.absoluteDir().absolutePath());
        QApplication::setOverrideCursor(Qt::WaitCursor);
        QApplication::processEvents();
        loadDCM(dir);
        QApplication::restoreOverrideCursor();
      }
    }

  }
private:
  qglviewer::Vec first_offset;
#ifdef CGAL_USE_VTK
  vtkImageData* vtk_image;
  vtkDICOMImageReader* dicom_reader;
  vtkDemandDrivenPipeline* executive;
  vtkImageGaussianSmooth* smoother;
#endif // CGAL_USE_VTK
  bool is_gray;
  Messages_interface* message_interface;
  QMessageBox msgBox;
  QAction* planeSwitch;
  QWidget *x_control, *y_control, *z_control;
  QSlider *x_slider, *y_slider, *z_slider;
  QLineEdit*x_cubeLabel, *y_cubeLabel, *z_cubeLabel;
  QHBoxLayout *x_box, *y_box, *z_box;
  PixelReader pxr_;
  Ui::ImagePrecisionDialog ui;

  CGAL::Three::Scene_group_item* group;
  std::vector<Volume_plane_thread*> threads;
  struct Controls{
    CGAL::Three::Scene_item* group;
    CGAL::Three::Scene_item* x_item;
    CGAL::Three::Scene_item* y_item;
    CGAL::Three::Scene_item* z_item;
    int x_value;
    int y_value;
    int z_value;
  };
  Controls *current_control;
  QMap<CGAL::Three::Scene_item*, Controls> group_map;
  unsigned int intersection_id;
  bool loadDCM(QString filename);
  Image* createDCMImage(QString dirname);
  QLayout* createOrGetDockLayout() {
    QLayout* layout = NULL;
    QDockWidget* controlDockWidget = mw->findChild<QDockWidget*>("volumePlanesControl");;

    if(!controlDockWidget) {
      controlDockWidget = new QDockWidget(mw);
      controlDockWidget->setObjectName("volumePlanesControl");
      QWidget* content = new QWidget(controlDockWidget);
      layout = new QVBoxLayout(content);
      layout->setObjectName("vpSliderLayout");
      controlDockWidget->setWindowTitle("Volume Planes Control");
      mw->addDockWidget(Qt::LeftDockWidgetArea, controlDockWidget);

      QWidget* vlabels = new QWidget(content);
      layout->addWidget(vlabels);
      layout->setAlignment(Qt::AlignLeft);
      QHBoxLayout* vbox = new QHBoxLayout(vlabels);
      vbox->setAlignment(Qt::AlignLeft);
      QLabel* text = new QLabel(vlabels);
      text->setText("Value of that pixel:");
      QLabel* help = new QLabel(vlabels);
      help->setText("Cut planes for the selected image:");
      QLabel* x = new QLabel(vlabels);

      connect(&pxr_, SIGNAL(x(QString)), x, SLOT(setText(QString)));

      layout->addWidget(help); vbox->addWidget(text); vbox->addWidget(x);
      controlDockWidget->setWidget(content);
      controlDockWidget->hide();

    } else {
      layout = controlDockWidget->findChild<QLayout*>("vpSliderLayout");
      controlDockWidget->show();
    }

    return layout;
  }

  void createPlanes(Scene_image_item* seg_img) {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    //Control widgets creation
    QLayout* layout = createOrGetDockLayout();
    QRegExpValidator* validator = new QRegExpValidator(QRegExp("\\d*"), this);
    if(x_control == NULL)
    {
      x_control = new QWidget;
      x_box = new QHBoxLayout(x_control);
      layout->addWidget(x_control);

      QLabel* label = new QLabel(x_control);
      label->setStyleSheet("QLabel { color : red; }");
      label->setText("X Slice");

      x_cubeLabel = new QLineEdit(x_control);

      // Find the right width for the label to accommodate at least 9999
      QFontMetrics metric = x_cubeLabel->fontMetrics();
      x_cubeLabel->setFixedWidth(metric.width(QString("9999")));
      x_cubeLabel->setText("0");
      x_cubeLabel->setValidator(validator);

      x_slider = new QSlider(mw);

      x_box->addWidget(label);
      x_box->addWidget(x_slider);
      x_box->addWidget(x_cubeLabel);
    }

    if(y_control == NULL)
    {
      y_control = new QWidget;
      y_box = new QHBoxLayout(y_control);
      layout->addWidget(y_control);

      QLabel* label = new QLabel(y_control);
      label->setStyleSheet("QLabel { color : green; }");
      label->setText("Y Slice");

      y_cubeLabel = new QLineEdit(y_control);

      // Find the right width for the label to accommodate at least 9999
      QFontMetrics metric = y_cubeLabel->fontMetrics();
      y_cubeLabel->setFixedWidth(metric.width(QString("9999")));
      y_cubeLabel->setText("0");
      y_cubeLabel->setValidator(validator);
      y_slider = new QSlider(mw);

      y_box->addWidget(label);
      y_box->addWidget(y_slider);
      y_box->addWidget(y_cubeLabel);
    }

    if(z_control == NULL)
    {
      z_control = new QWidget;
      z_box = new QHBoxLayout(z_control);
      layout->addWidget(z_control);

      QLabel* label = new QLabel(z_control);
      label->setStyleSheet("QLabel { color : blue; }");
      label->setText("Z Slice");

      z_cubeLabel = new QLineEdit(z_control);

      // Find the right width for the label to accommodate at least 9999
      QFontMetrics metric = z_cubeLabel->fontMetrics();
      z_cubeLabel->setFixedWidth(metric.width(QString("9999")));
      z_cubeLabel->setText("0");
      z_cubeLabel->setValidator(validator);
      z_slider = new QSlider(mw);

      z_box->addWidget(label);
      z_box->addWidget(z_slider);
      z_box->addWidget(z_cubeLabel);
    }
    if(!(seg_img == NULL)) {
      const CGAL::Image_3* img = seg_img->image();
      CGAL_IMAGE_IO_CASE(img->image(), this->launchAdders<Word>(seg_img, seg_img->name()))

          Volume_plane_intersection* i
          = new Volume_plane_intersection(img->xdim() * img->vx()-1,
                                          img->ydim() * img->vy()-1,
                                          img->zdim() * img->vz()-1);
      this->intersection_id = scene->addItem(i);
      scene->changeGroup(i, group);
      group->lockChild(i);
    } else {
      QMessageBox::warning(mw, tr("Something went wrong"), tr("Selected a suitable Object but couldn't get an image pointer."));
      return;
    }
  }


  template<typename Word>
  void launchAdders(Scene_image_item* seg_img, const QString& name) {
    const CGAL::Image_3* img = seg_img->image();
    const Word* begin = (const Word*)img->data();
    const Word* end = (const Word*)img->data() + img->size();

    std::pair<const Word, const Word> minmax(*std::min_element(begin, end), *std::max_element(begin, end));

    Clamp_to_one_zero_range clamper = { minmax };

    switchReaderConverter< Word >(minmax);

    Volume_plane<x_tag> *x_item = new Volume_plane<x_tag>();
    Volume_plane<y_tag> *y_item = new Volume_plane<y_tag>();
    Volume_plane<z_tag> *z_item = new Volume_plane<z_tag>();

    x_item->setProperty("img",qVariantFromValue((void*)seg_img));
    y_item->setProperty("img",qVariantFromValue((void*)seg_img));
    z_item->setProperty("img",qVariantFromValue((void*)seg_img));

    x_item->setColor(QColor("red"));
    y_item->setColor(QColor("green"));
    z_item->setColor(QColor("blue"));
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first());
    viewer->installEventFilter(x_item);
    viewer->installEventFilter(y_item);
    viewer->installEventFilter(z_item);
    connect(x_item, SIGNAL(selected(CGAL::Three::Scene_item*)), this, SLOT(select_plane(CGAL::Three::Scene_item*)));
    connect(y_item, SIGNAL(selected(CGAL::Three::Scene_item*)), this, SLOT(select_plane(CGAL::Three::Scene_item*)));
    connect(z_item, SIGNAL(selected(CGAL::Three::Scene_item*)), this, SLOT(select_plane(CGAL::Three::Scene_item*)));
    scene->setSelectedItem(-1);
    group = new Scene_group_item(QString("Planes for %1").arg(seg_img->name()));
    connect(group, SIGNAL(aboutToBeDestroyed()),
            this, SLOT(erase_group()));
    scene->addItem(group);
    Controls c;
    c.group = group;
    c.x_item = x_item;
    c.y_item = y_item;
    c.z_item = z_item;
    c.x_value = 0;
    c.y_value = 0;
    c.z_value = 0;
    group_map[seg_img] = c;
    current_control = &group_map[seg_img];
    connect(seg_img, SIGNAL(aboutToBeDestroyed()),
            this, SLOT(on_img_detroyed()));

    threads.push_back(new X_plane_thread<Word>(x_item, img, clamper, name));

    connect(threads.back(), SIGNAL(finished(Volume_plane_thread*)), this, SLOT(addVP(Volume_plane_thread*)));
    threads.back()->start();
    threads.push_back(new Y_plane_thread<Word>(y_item,img, clamper, name));
    connect(threads.back(), SIGNAL(finished(Volume_plane_thread*)), this, SLOT(addVP(Volume_plane_thread*)));
    threads.back()->start();

    threads.push_back(new Z_plane_thread<Word>(z_item,img, clamper, name));
    connect(threads.back(), SIGNAL(finished(Volume_plane_thread*)), this, SLOT(addVP(Volume_plane_thread*)));
    threads.back()->start();

    first_offset = viewer->offset();

  }
  template<typename T>
  void switchReaderConverter(std::pair<T, T> minmax) {
    switchReaderConverter(minmax, typename boost::is_integral<T>::type());
  }

  template<typename T>
  void switchReaderConverter(std::pair<T, T> minmax, boost::true_type) {
    // IntConverter
    IntConverter x = { minmax }; pxr_.setIC(x);
  }

  template<typename T>
  void switchReaderConverter(std::pair<T, T> minmax, boost::false_type) {
    // IntConverter
    DoubleConverter x = { minmax }; pxr_.setFC(x);
  }

private Q_SLOTS:
  void select_plane(CGAL::Three::Scene_item* item)
  {
    //Scene* true_scene = dynamic_cast<Scene*>(scene);
    Scene_image_item* img = (Scene_image_item*)item->property("img").value<void*>();
    if(img)
      scene->setSelectedItem(scene->item_id(img));
  }
  //updates the msgBox
  void update_msgBox()
  {
    static int nbPlanes =0;
    nbPlanes ++;
    msgBox.setText(QString("Planes created : %1/3").arg(nbPlanes));
    if(nbPlanes == 3)
    {
      const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
      if(offset != first_offset)
      {
        for(int i=0; i<scene->numberOfEntries(); ++i)
        {
          scene->item(i)->invalidateOpenGLBuffers();
        }
      }
      msgBox.hide();
      nbPlanes = 0;
      QApplication::restoreOverrideCursor();
    }
  }
  // Avoids the segfault after the deletion of an item
  void erase_group()
  {

    CGAL::Three::Scene_group_item* group_item = qobject_cast<CGAL::Three::Scene_group_item*>(sender());
    if(group_item)
    {

      Q_FOREACH(CGAL::Three::Scene_item* key, group_map.keys())
      {
        if(group_map[key].group == group_item)
        {
          group_map[key].x_item = NULL;
          group_map[key].y_item = NULL;
          group_map[key].z_item = NULL;
          group_map.remove(key);
          break;
        }
      }
    }
      //try to re-connect to another group
      if(!group_map.isEmpty())
      {
        int id = scene->item_id(group_map.keys().first());
        connect_controls(id);
      }
  }
    //destroy planes on image deletion
    void on_img_detroyed()
    {
      Scene_image_item* img_item = qobject_cast<Scene_image_item*>(sender());
      if(img_item)
      {
        Scene_group_item* group = qobject_cast<Scene_group_item*>(group_map[img_item].group);
        if(!group)
          return;
        group_map[img_item].x_item = NULL;
        group_map[img_item].y_item = NULL;
        group_map[img_item].z_item = NULL;
        disconnect(group_map[img_item].group, SIGNAL(aboutToBeDestroyed()),
                           this, SLOT(erase_group()));
        group_map.remove(img_item);
        QList<int> deletion;
        Q_FOREACH(Scene_item* child, group->getChildren())
        {
          group->unlockChild(child);
          deletion.append(scene->item_id(child));
        }
        deletion.append(scene->item_id(group));
        scene->erase(deletion);
      }
    //try to re-connect to another group
    if(!group_map.isEmpty())
    {
      int id = scene->item_id(group_map.keys().first());
      connect_controls(id);
    }
  }
  void connect_controls(int id)
  {
    CGAL::Three::Scene_item* sel_itm = scene->item(id);
    if(!sel_itm)
      return;
    if(!group_map.contains(sel_itm)) //the planes are not yet created or the selected item is not a segmented_image
    {
      Scene_image_item* img = (Scene_image_item*)sel_itm->property("img").value<void*>();
      if(img)
        sel_itm = img;
      else
        return;
    }
    Controls c = group_map[sel_itm];
    current_control = &group_map[sel_itm];
    // x line
    if(c.x_item != NULL)
    {
      Volume_plane_interface* x_plane = qobject_cast<Volume_plane_interface*>(c.x_item);
      if(x_slider)
        delete x_slider;
      x_slider = new Plane_slider(x_plane->translationVector(), scene->item_id(x_plane), scene, x_plane->manipulatedFrame(),
                                  Qt::Horizontal, x_control);
      x_slider->setRange(0, (x_plane->cDim() - 1) * 100);
      connect(x_slider, SIGNAL(realChange(int)), this, SLOT(setXNum(int)));
      connect(x_plane, SIGNAL(manipulated(int)), this, SLOT(setXNum(int)));
      connect(x_cubeLabel, SIGNAL(textEdited(QString)), this, SLOT(XTextEdited(QString)));
      connect(x_plane, SIGNAL(aboutToBeDestroyed()), this, SLOT(destroy_x_item()));
      connect(x_slider, SIGNAL(realChange(int)), this, SLOT(set_value()));
      connect(x_slider, SIGNAL(sliderMoved(int)), x_slider, SLOT(updateFramePosition()));
      x_slider->setValue(c.x_value);

      x_box->addWidget(x_slider);
      x_box->addWidget(x_cubeLabel);
    }
    //y line
    if(c.y_item != NULL)
    {
      Volume_plane_interface* y_plane = qobject_cast<Volume_plane_interface*>(c.y_item);
      if(y_slider)
        delete y_slider;
      y_slider = new Plane_slider(y_plane->translationVector(), scene->item_id(y_plane), scene, y_plane->manipulatedFrame(),
                                  Qt::Horizontal, z_control);
      y_slider->setRange(0, (y_plane->cDim() - 1) * 100);
      connect(y_slider, SIGNAL(realChange(int)), this, SLOT(setYNum(int)));
      connect(y_plane, SIGNAL(manipulated(int)), this, SLOT(setYNum(int)));
      connect(y_cubeLabel, SIGNAL(textEdited(QString)), this, SLOT(YTextEdited(QString)));
      connect(y_plane, SIGNAL(aboutToBeDestroyed()), this, SLOT(destroy_y_item()));
      connect(y_slider, SIGNAL(realChange(int)), this, SLOT(set_value()));
      connect(y_slider, SIGNAL(sliderMoved(int)), y_slider, SLOT(updateFramePosition()));
      y_slider->setValue(c.y_value);
      y_box->addWidget(y_slider);
      y_box->addWidget(y_cubeLabel);
    }
    // z line
    if(c.z_item != NULL)
    {
      Volume_plane_interface* z_plane = qobject_cast<Volume_plane_interface*>(c.z_item);
      if(z_slider)
        delete z_slider;
      z_slider = new Plane_slider(z_plane->translationVector(), scene->item_id(z_plane), scene, z_plane->manipulatedFrame(),
                                  Qt::Horizontal, z_control);
      z_slider->setRange(0, (z_plane->cDim() - 1) * 100);
      connect(z_slider, SIGNAL(realChange(int)), this, SLOT(setZNum(int)));
      connect(z_plane, SIGNAL(manipulated(int)), this, SLOT(setZNum(int)));
      connect(z_cubeLabel, SIGNAL(textEdited(QString)), this, SLOT(ZTextEdited(QString)));
      connect(z_plane, SIGNAL(aboutToBeDestroyed()), this, SLOT(destroy_z_item()));
      connect(z_slider, SIGNAL(sliderMoved(int)), z_slider, SLOT(updateFramePosition()));
      connect(z_slider, SIGNAL(realChange(int)), this, SLOT(set_value()));
      z_slider->setValue(c.z_value);
      z_box->addWidget(z_slider);
      z_box->addWidget(z_cubeLabel);
    }
  }
//Keeps the position of the planes for the next time
  void set_value()
  {
    current_control->x_value = x_slider->value();
    current_control->y_value = y_slider->value();
    current_control->z_value = z_slider->value();
  }

  void destroy_x_item()
  {
    current_control->x_item = NULL;
    if(group_map.isEmpty())
      x_control->hide();
  }
  void destroy_y_item()
  {
    current_control->y_item = NULL;
    if(group_map.isEmpty())
      y_control->hide();
  }
  void destroy_z_item()
  {
    current_control->z_item = NULL;
    if(group_map.isEmpty())
      z_control->hide();
  }

};


QString Io_image_plugin::nameFilters() const {
  return QString("Inrimage files (*.inr *.inr.gz) ;; "
                 "Analyze files (*.hdr *.img *img.gz) ;; "
                 "Stanford Exploration Project files (*.H *.HH)");
}


bool Io_image_plugin::canLoad() const {
  return true;
}

template<typename Word>
void convert(Image* image)
{
  float *f_data = (float*)ImageIO_alloc(image->xdim()*image->ydim()*image->zdim()*sizeof(float));
  Word* d_data = (Word*)(image->data());
  //convert image from double to float
  for(std::size_t x = 0; x<image->xdim(); ++x)
    for(std::size_t y = 0; y<image->ydim(); ++y)
      for(std::size_t z = 0; z<image->zdim(); ++z)
      {
        std::size_t i =(z * image->ydim() + y) * image->xdim() + x;
        f_data[i] =(float)d_data[i];
      }
  ImageIO_free(d_data);
  image->image()->data = (void*)f_data;
  image->image()->wdim = 4;
  image->image()->wordKind = WK_FLOAT;
}
CGAL::Three::Scene_item*
Io_image_plugin::load(QFileInfo fileinfo) {
  QApplication::restoreOverrideCursor();
  Image* image = new Image;
  if(fileinfo.suffix() != "H" && fileinfo.suffix() != "HH" &&
     !image->read(fileinfo.filePath().toUtf8()))
    {
      QMessageBox qmb(QMessageBox::NoIcon,
                      "Raw Dialog",
                      tr("Error with file <tt>%1</tt>:\n"
                         "unknown file format!\n"
                         "\n"
                         "Open it as a raw image?").arg(fileinfo.fileName()),
                      QMessageBox::Yes|QMessageBox::No);

      bool success = true;
      if(qmb.exec() == QMessageBox::Yes) {
        Raw_image_dialog raw_dialog;
        raw_dialog.label_file_size->setText(QString("%1 B").arg(fileinfo.size()));
        raw_dialog.buttonBox->button(QDialogButtonBox::Open)->setEnabled(false);
        if( raw_dialog.exec() ){

          QApplication::setOverrideCursor(Qt::WaitCursor);
          QApplication::processEvents();

          if(image->read_raw(fileinfo.filePath().toUtf8(),
                             raw_dialog.dim_x->value(),
                             raw_dialog.dim_y->value(),
                             raw_dialog.dim_z->value(),
                             raw_dialog.spacing_x->value(),
                             raw_dialog.spacing_y->value(),
                             raw_dialog.spacing_z->value(),
                             raw_dialog.offset->value(),
                             raw_dialog.image_word_size(),
                             raw_dialog.image_word_kind(),
                             raw_dialog.image_sign())
             ){
            switch(raw_dialog.image_word_kind())
            {
            case WK_FLOAT:
              is_gray = true;
              convert<double>(image);
              break;
            case WK_FIXED:
            {
              switch(raw_dialog.image_word_size())
              {
              case 2:
                is_gray = true;
                convert<short>(image);
                break;
              case 4:
                is_gray = true;
                convert<int>(image);
                break;
              default:
                is_gray = false;
                break;
              }
              break;
            }
            default:
              break;
            }
            QSettings settings;
            settings.beginGroup(QUrl::toPercentEncoding(fileinfo.absoluteFilePath()));
            settings.setValue("is_raw", true);
            settings.setValue("dim_x", raw_dialog.dim_x->value());
            settings.setValue("dim_y", raw_dialog.dim_y->value());
            settings.setValue("dim_z", raw_dialog.dim_z->value());
            settings.setValue("spacing_x", raw_dialog.spacing_x->value());
            settings.setValue("spacing_y", raw_dialog.spacing_y->value());
            settings.setValue("spacing_z", raw_dialog.spacing_z->value());
            settings.setValue("offset", raw_dialog.offset->value());
            settings.setValue("wdim", QVariant::fromValue(raw_dialog.image_word_size()));
            settings.setValue("wk", raw_dialog.image_word_kind());
            settings.setValue("sign", raw_dialog.image_sign());
            settings.endGroup();
          }else {
            success = false;
          }
        }else {
          success = false;
        }
      }else {
        success = false;
      }
      if(!success){
        delete image;
        return NULL;
      }
    }
  //read a sep file
  else if(fileinfo.suffix() == "H" || fileinfo.suffix() == "HH")
  {
    Sep_reader<float> reader(fileinfo.filePath().toUtf8().data());
    *image = *reader.cgal_image();
    is_gray = true;
  }
  // Get display precision
  QDialog dialog;
  ui.setupUi(&dialog);

  connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));
  connect(ui.imageType, SIGNAL(currentIndexChanged(int)),
          this, SLOT(on_imageType_changed(int)));
  dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);

  // Add precision values to the dialog
  for ( int i=1 ; i<9 ; ++i )
  {
    QString s = tr("1:%1").arg(i*i*i);
    ui.precisionList->addItem(s);
  }

  //Adds Image type
  ui.imageType->addItem(QString("Segmented image"));
  ui.imageType->addItem(QString("Gray-level image"));

  QString type;
  int voxel_scale = 0;
  // Open window
  QApplication::restoreOverrideCursor();
  if(!is_gray)
  {
    int return_code = dialog.exec();
    if(return_code != QDialog::Accepted)
    {
      delete image;
      return NULL;
    }

    // Get selected precision
    voxel_scale = ui.precisionList->currentIndex() + 1;

    //Get the image type
    type = ui.imageType->currentText();
  }
  else
    type = "Gray-level image";
  QApplication::setOverrideCursor(Qt::WaitCursor);
  QApplication::processEvents();
  Scene_image_item* image_item;
  if(type == "Gray-level image")
  {
    //Create planes
    image_item = new Scene_image_item(image,0, true);
    image_item->setName(fileinfo.baseName());
    msgBox.setText("Planes created : 0/3");
    msgBox.setStandardButtons(QMessageBox::NoButton);
    msgBox.show();

    createPlanes(image_item);
  }
  else
    image_item = new Scene_image_item(image,voxel_scale, false);
  image_item->setName(fileinfo.baseName());
  return image_item;
}

bool Io_image_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_image_item*>(item);
}

bool Io_image_plugin::loadDCM(QString dirname)
{
  QApplication::restoreOverrideCursor();
#ifdef CGAL_USE_VTK
  QFileInfo fileinfo;
  fileinfo.setFile(dirname);
  bool result = true;
  if(!fileinfo.isReadable())
  {
    QMessageBox::warning(mw, mw->windowTitle(),
                         tr("Cannot read directory <tt>%1</tt>!").arg(dirname));
    message_interface->warning(tr("Opening of directory %1 failed!").arg(dirname));
    result = false;
  }
  else
  {
    // Get display precision
    QDialog dialog;
    ui.setupUi(&dialog);

    connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));
    connect(ui.imageType, SIGNAL(currentIndexChanged(int)),
            this, SLOT(on_imageType_changed(int)));


    // Add precision values to the dialog
    for ( int i=1 ; i<9 ; ++i )
    {
      QString s = tr("1:%1").arg(i*i*i);
      ui.precisionList->addItem(s);
    }

    //Adds Image type
    ui.imageType->addItem(QString("Segmented image"));
    ui.imageType->addItem(QString("Gray-level image"));


    // Open window
    QApplication::restoreOverrideCursor();
    int return_code = dialog.exec();
    if(return_code != QDialog::Accepted)
    {
      return false;
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();

    // Get selected precision
    int voxel_scale = ui.precisionList->currentIndex() + 1;

    //Get the image type
    QString type = ui.imageType->currentText();
    Scene_image_item* image_item;
    if(type == "Gray-level image")
    {

      Image *image = createDCMImage(dirname);
      if(image->image() == 0)
      {
        QMessageBox::warning(mw, mw->windowTitle(),
                             tr("Error with file <tt>%1/</tt>:\nunknown file format!").arg(dirname));
        message_interface->warning(tr("Opening of file %1/ failed!").arg(dirname));
        result = false;
      }
      else
      {
        message_interface->information(tr("File %1/ successfully opened.").arg(dirname));
      }
      if(result)
      {
      //Create planes
      image_item = new Scene_image_item(image,125, true);
      msgBox.setText("Planes created : 0/3");
      msgBox.setStandardButtons(QMessageBox::NoButton);
      msgBox.show();
      createPlanes(image_item);
      image_item->setName(fileinfo.baseName());
      scene->addItem(image_item);
      }
    }
    else
    {
      Image *image = createDCMImage(dirname);
      if(image->image() == 0)
      {
        QMessageBox::warning(mw, mw->windowTitle(),
                             tr("Error with file <tt>%1/</tt>:\nunknown file format!").arg(dirname));
        message_interface->warning(tr("Opening of file %1/ failed!").arg(dirname));
        result = false;
      }
      else
      {
        message_interface->information(tr("File %1/ successfully opened.").arg(dirname));
      }
      if(result)
      {
        image_item = new Scene_image_item(image,voxel_scale, false);
        image_item->setName(fileinfo.baseName());
        scene->addItem(image_item);
      }
    }
  }
  return result;
#else
  message_interface->warning("You need VTK to read a DCM file");
  CGAL_USE(dirname);
  return false;
#endif
}
Image* Io_image_plugin::createDCMImage(QString dirname)
{
  Image* image = NULL;
#ifdef CGAL_USE_VTK
  dicom_reader = vtkDICOMImageReader::New();
  dicom_reader->SetDirectoryName(dirname.toUtf8());

  executive =
    vtkDemandDrivenPipeline::SafeDownCast(dicom_reader->GetExecutive());
  if (executive)
  {
    executive->SetReleaseDataFlag(0, 0); // where 0 is the port index
  }

  smoother = vtkImageGaussianSmooth::New();
  smoother->SetStandardDeviations(1., 1., 1.);
  smoother->SetInputConnection(dicom_reader->GetOutputPort());
  smoother->Update();
  vtk_image = smoother->GetOutput();
  vtk_image->Print(std::cerr);
  image = new Image;
  *image = CGAL::read_vtk_image_data(vtk_image);
#else
  message_interface->warning("You need VTK to read a DCM file");
  CGAL_USE(dirname);
#endif
  return image;

}

#include "Io_image_plugin.moc"
