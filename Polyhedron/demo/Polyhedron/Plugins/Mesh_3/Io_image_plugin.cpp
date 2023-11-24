#ifdef _MSC_VER
#  pragma warning(disable:4244) // conversion with loss of data
#  pragma warning(disable:4996) // boost_1_65_1\boost/iostreams/positioning.hpp(96):
                                // warning C4996: 'std::fpos<_Mbstatet>::seekpos': warning STL4019:
                                // The member std::fpos::seekpos() is non-Standard
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
#include <CGAL/SEP_to_ImageIO.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/config.h>
#include <CGAL/use.h>
#include <CGAL/IO/helpers.h>
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

#include <QSettings>
#include <QUrl>
#include "Raw_image_dialog.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>

#ifdef CGAL_USE_VTK
#include <CGAL/IO/read_vtk_image_data.h>

#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDICOMImageReader.h>
#include <vtkBMPReader.h>
#include <vtkNIFTIImageReader.h>
#include <vtkImageReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkNrrdReader.h>
#endif

#include <CGAL/Three/Three.h>
#include <CGAL/use.h>

#include <boost/type_traits.hpp>
#include <optional>
#include <boost/filesystem.hpp>

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <locale>

template <class vtkReader>
bool
load_vtk_file(QFileInfo fileinfo, Image* image)
{
#ifdef CGAL_USE_VTK
    vtkNew<vtkReader> reader;
    reader->SetFileName(fileinfo.filePath().toUtf8());
    reader->Update();
    auto vtk_image = reader->GetOutput();
    vtk_image->Print(std::cerr);
    *image = CGAL::IO::read_vtk_image_data(vtk_image); // copy the image data
    return true;
#else
    CGAL_USE(fileinfo);
    CGAL_USE(image);
    return false;
#endif
}

// Covariant return types don't work for scalar types and we cannot
// have templates here, hence this unfortunate hack.

// The input float value we are reading is always in
// 0..1 and min_max is the range it came from.
struct IntConverter
{
  std::pair<int, int> min_max;

  int operator()(float f)
  {
    float s = f * float((min_max.second - min_max.first));
    //approximate instead of just floor.
    if (s - floor(s) >= 0.5)
      return int(s)+1 + min_max.first;
    else
      return s + float(min_max.first);
  }
};

struct DoubleConverter
{
  std::pair<float, float> min_max;

  float operator()(float f)
  {
    float s = f * (min_max.second - min_max.first);
    return s + min_max.first;
  }
};

class PixelReader
  : public QObject
{
Q_OBJECT

public Q_SLOTS:
  void update(const QMouseEvent *e) { getPixel(e->pos()); }

Q_SIGNALS:
  void x(QString);

public:
  void setIC(const IntConverter& x) { ic = x; fc = std::optional<DoubleConverter>(); }
  void setFC(const DoubleConverter& x) { fc = x; ic = std::optional<IntConverter>(); }
  void setViewer(Viewer_interface* viewer) { this->viewer = viewer; }

private:
  std::optional<IntConverter> ic;
  std::optional<DoubleConverter> fc;
  Viewer_interface* viewer;
  void getPixel(const QPoint& e)
  {
    const auto data = read_pixel_as_float_rgb(e, viewer, viewer->camera());
    if(fc) {
      Q_EMIT x(QString::number((*fc)(data[0]), 'f', 6 ));
    } else if(ic) {
      Q_EMIT x( QString::number((*ic)(data[0]) ));
    }
  }
};

class Plane_slider
  : public QSlider
{
  Q_OBJECT

public:
  Plane_slider(const CGAL::qglviewer::Vec& v,
               int id,
               Scene_interface* scene,
               CGAL::qglviewer::ManipulatedFrame* frame,
               Qt::Orientation ori,
               QWidget* widget)
    : QSlider(ori, widget), v(v), id(id), scene(scene), frame(frame)
  {
    this->setTracking(true);
    connect(frame, SIGNAL(manipulated()), this, SLOT(updateCutPlane()));
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

    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    CGAL::qglviewer::Vec v2 = v * (this->value() / scale);
    v2 += offset;
    frame->setTranslationWithConstraint(v2);
    scene->itemChanged(id);

    Q_EMIT realChange(this->value() / scale);
    ready_to_move = false;
  }

  void updateValue()
  {
    if(!ready_to_cut)
      return;

    typedef qreal qglviewer_real;
    qglviewer_real a, b, c;

    frame->getPosition(a, b, c);
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    a -= offset.x;
    b -= offset.y;
    c -= offset.z;
    float sum1 = float(a + b + c);
    float sum2 = float(v.x + v.y + v.z);
    sum1 /= sum2;
    setValue(sum1 * float(scale));
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
  CGAL::qglviewer::Vec v;
  int id;

  Scene_interface* scene;
  CGAL::qglviewer::ManipulatedFrame* frame;
};

const unsigned int Plane_slider::scale = 100;

enum class Directory_extension_type
{
  DCM = 0,
  BMP
};

class Io_image_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_helper,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "io_image_plugin.json")

public:
  bool applicable(QAction*) const override
  {
    return qobject_cast<Scene_image_item*>(scene->item(scene->mainSelectionIndex()));
  }

  using Polyhedron_demo_io_plugin_interface::init;
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface *mi) override
  {
    this->message_interface = mi;
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->is_gray = false;

    x_control = nullptr;
    y_control = nullptr;
    z_control = nullptr;
    current_control = nullptr;

    planeSwitch = new QAction("Add Volume Planes", mw);

    QAction *actionLoadDCM = new QAction("Open Directory (DCM files)", mw);
    connect(actionLoadDCM, SIGNAL(triggered()), this, SLOT(on_actionLoadDCM_triggered()));

    QAction *actionLoadBMP = new QAction("Open Directory (BMP files)", mw);
    connect(actionLoadBMP, SIGNAL(triggered()), this, SLOT(on_actionLoadBMP_triggered()));

    if(planeSwitch)
    {
      planeSwitch->setProperty("subMenuName", "3D Mesh Generation");
      connect(planeSwitch, SIGNAL(triggered()),
              this, SLOT(selectPlanes()));
      connect(CGAL::Three::Three::connectableScene(),SIGNAL(itemIndexSelected(int)),
              this, SLOT(connect_controls(int)));
    }

    Viewer_interface* v = CGAL::Three::Three::mainViewer();
    CGAL_assertion(v != nullptr);

    pxr_.setViewer(v);
    connect(v, SIGNAL(pointSelected(const QMouseEvent *)), &pxr_, SLOT(update(const QMouseEvent *)));
    createOrGetDockLayout();
    connect(mw, SIGNAL(newViewerCreated(QObject*)),
            this, SLOT(connectNewViewer(QObject*)));

    QMenu* menuFile = mw->findChild<QMenu*>("menuFile");
    if(nullptr != menuFile )
    {
      QList<QAction*> menuFileActions = menuFile->actions();

      // Look for action just after "Load..." action
      QAction* actionAfterLoad = nullptr;
      for(QList<QAction*>::iterator it_action = menuFileActions.begin(),
          end = menuFileActions.end() ; it_action != end ; ++ it_action ) //for( QAction* action : menuFileActions)
      {
        if(NULL != *it_action && (*it_action)->text().contains("Load..."))
        {
          ++it_action;
          if(it_action != end && NULL != *it_action)
          {
            actionAfterLoad = *it_action;
          }
        }
      }

      // Insert "Load implicit function" action
      if(nullptr != actionAfterLoad)
      {
        menuFile->insertAction(actionAfterLoad, actionLoadDCM);
        menuFile->insertAction(actionAfterLoad, actionLoadBMP);
      }
    }
  }

  QList<QAction*> actions() const override
  {
    return QList<QAction*>() << planeSwitch;
  }

  virtual void closure() override
  {
    QDockWidget* controlDockWidget = mw->findChild<QDockWidget*>("volumePlanesControl");
    if(controlDockWidget)
      controlDockWidget->hide();
  }

  Io_image_plugin() : planeSwitch(nullptr) { }

  QString nameFilters() const override;
  bool canLoad(QFileInfo) const override;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene = true) override;

  bool canSave(const CGAL::Three::Scene_item*) override;
  bool save(QFileInfo fileinfo, QList<CGAL::Three::Scene_item*>& items ) override
  {
    Scene_item* item = items.front();
    const Scene_image_item* im_item = qobject_cast<const Scene_image_item*>(item);

    point_image p_im = *im_item->image()->image();
    bool ok = _writeImage(&p_im, fileinfo.filePath().toUtf8()) == 0;
    if(ok)
      items.pop_front();
    return ok;
  }

  bool isDefaultLoader(const Scene_item* item) const override
  {
    if(qobject_cast<const Scene_image_item*>(item))
      return true;
    return false;
  }

  QString name() const override{ return "segmented images"; }

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

  void selectPlanes()
  {
    std::vector< Scene_image_item* > seg_items;
    Scene_image_item* seg_img;
    seg_img = nullptr;
    for(int i = 0; i < scene->numberOfEntries(); ++i)
    {
      Scene_image_item* tmp = qobject_cast<Scene_image_item*>(scene->item(i));
      if(tmp != nullptr)
        seg_items.push_back(tmp);
    }

    if(seg_items.empty())
    {
      QMessageBox::warning(mw, tr("No suitable item found"), tr("Load an inrimage or hdr file to enable Volume Planes."));
      return;
    }
    else
    {
      QList<QString> items;
      for(std::vector< Scene_image_item* >::const_iterator it = seg_items.begin();
          it != seg_items.end(); ++it)
      {
        items << (*it)->name();
      }

      bool ok;
      QString selected = QInputDialog::getItem(mw, tr("Select a dataset:"), tr("Items"), items, 0, false, &ok);
      if(!ok || selected.isEmpty())
        return;

      for(std::vector< Scene_image_item*>::const_iterator it = seg_items.begin();
          it != seg_items.end(); ++it)
      {
        if(selected == (*it)->name())
          seg_img = *it;
      }
    }

    if(group_map.keys().contains(seg_img))
    {
      CGAL::Three::Three::warning("This item already has volume planes.");
    }
    else
    {
      // Opens a modal Dialog to prevent the user from manipulating things that could mess with the planes creation and cause a segfault.
      msgBox.setText("Planes created : 0/3");
      msgBox.setStandardButtons(QMessageBox::NoButton);
      msgBox.show();
      createPlanes(seg_img);
    }
  }

  void addVP(Volume_plane_thread* thread)
  {
    Volume_plane_interface* plane = thread->getItem();
    plane->init(Three::mainViewer());

    // add the interface for this Volume_plane
    int id = scene->addItem(plane);
    scene->changeGroup(plane, group);
    group->lockChild(plane);
    //connect(plane->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::manipulated,
    //        plane, &Volume_plane_interface::redraw);

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

    // this slot has been connected to a thread that hasn't been registered here.
    assert(it != threads.end());
    delete *it;
    threads.erase(it);

    update_msgBox();
    Volume_plane_intersection* intersection = dynamic_cast<Volume_plane_intersection*>(scene->item(intersection_id));
    if(!intersection)
    {
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

  void loadDirectory(const Directory_extension_type ext)
  {
    QSettings settings;
    QString start_dir = settings.value("Open directory", QDir::current().dirName()).toString();
    QString dir = QFileDialog::getExistingDirectory(mw, tr("Open directory"),
                                                    start_dir, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    if(!dir.isEmpty())
    {
      QFileInfo fileinfo(dir);
      if(fileinfo.isDir() && fileinfo.isReadable())
      {
        settings.setValue("Open directory", fileinfo.absoluteDir().absolutePath());
        QApplication::setOverrideCursor(Qt::WaitCursor);
        QApplication::processEvents();

        loadDirectory(dir, ext);
      }
    }
  }

  void on_actionLoadDCM_triggered()
  {
    return loadDirectory(Directory_extension_type::DCM);
  }

  void on_actionLoadBMP_triggered()
  {
    return loadDirectory(Directory_extension_type::BMP);
  }

  void connectNewViewer(QObject* o)
  {
    for(Controls c : group_map.values())
    {
      o->installEventFilter(c.x_item);
      o->installEventFilter(c.y_item);
      o->installEventFilter(c.z_item);
    }
  }

private:
  CGAL::qglviewer::Vec first_offset;
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

  struct Controls
  {
    CGAL::Three::Scene_item* group;
    CGAL::Three::Scene_item* x_item;
    CGAL::Three::Scene_item* y_item;
    CGAL::Three::Scene_item* z_item;
    int x_value;
    int y_value;
    int z_value;
  };

  Controls* current_control;
  QMap<CGAL::Three::Scene_item*, Controls> group_map;
  unsigned int intersection_id;

  bool loadDirectory(const QString& filename, const Directory_extension_type ext);
  Image* createDirectoryImage(const QString& dirname, const Directory_extension_type ext, const bool smooth);

  QLayout* createOrGetDockLayout()
  {
    QLayout* layout = nullptr;
    QDockWidget* controlDockWidget = mw->findChild<QDockWidget*>("volumePlanesControl");

    if(!controlDockWidget)
    {
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

    }
    else
    {
      layout = controlDockWidget->findChild<QLayout*>("vpSliderLayout");
      controlDockWidget->show();
      controlDockWidget->raise();
    }

    return layout;
  }

  void createPlanes(Scene_image_item* seg_img)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    //Control widgets creation
    QLayout* layout = createOrGetDockLayout();
    QRegularExpressionValidator* validator = new QRegularExpressionValidator(QRegularExpression("\\d*"), this);
    bool show_sliders = true;
    if(x_control == nullptr)
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
#if QT_VERSION >= QT_VERSION_CHECK(5, 11, 0)
      x_cubeLabel->setFixedWidth(metric.horizontalAdvance(QString(".9999.")));
#else
      x_cubeLabel->setFixedWidth(metric.width(QString(".9999.")));
#endif
      x_cubeLabel->setText("0");
      x_cubeLabel->setValidator(validator);

      x_slider = new QSlider(mw);

      x_box->addWidget(label);
      x_box->addWidget(x_slider);
      x_box->addWidget(x_cubeLabel);
      show_sliders &= seg_img->image()->xdim() > 1;
    }

    if(y_control == nullptr)
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
#if QT_VERSION >= QT_VERSION_CHECK(5, 11, 0)
      y_cubeLabel->setFixedWidth(metric.horizontalAdvance(QString(".9999.")));
#else
      y_cubeLabel->setFixedWidth(metric.width(QString(".9999.")));
#endif
      y_cubeLabel->setText("0");
      y_cubeLabel->setValidator(validator);
      y_slider = new QSlider(mw);

      y_box->addWidget(label);
      y_box->addWidget(y_slider);
      y_box->addWidget(y_cubeLabel);
      show_sliders &= seg_img->image()->ydim() > 1;
    }

    if(z_control == nullptr)
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
#if QT_VERSION >= QT_VERSION_CHECK(5, 11, 0)
      z_cubeLabel->setFixedWidth(metric.horizontalAdvance(QString(".9999.")));
#else
      z_cubeLabel->setFixedWidth(metric.width(QString(".9999.")));
#endif
      z_cubeLabel->setText("0");
      z_cubeLabel->setValidator(validator);
      z_slider = new QSlider(mw);

      z_box->addWidget(label);
      z_box->addWidget(z_slider);
      z_box->addWidget(z_cubeLabel);
      show_sliders &= seg_img->image()->zdim() > 1;
    }

    x_control->setEnabled(show_sliders);
    y_control->setEnabled(show_sliders);
    z_control->setEnabled(show_sliders);

    if(!(seg_img == nullptr))
    {
      const CGAL::Image_3* img = seg_img->image();
      CGAL_IMAGE_IO_CASE(img->image(), this->launchAdders<Word>(seg_img, seg_img->name()))

          Volume_plane_intersection* i
          = new Volume_plane_intersection(img->xdim() * img->vx()-1,
                                          img->ydim() * img->vy()-1,
                                          img->zdim() * img->vz()-1,
                                          img->image()->tx,
                                          img->image()->ty,
                                          img->image()->tz);
      this->intersection_id = scene->addItem(i);
      scene->changeGroup(i, group);
      group->lockChild(i);
    }
    else
    {
      QMessageBox::warning(mw, tr("Something went wrong"), tr("Selected a suitable Object but couldn't get an image pointer."));
      return;
    }
  }

  template<typename Word>
  void launchAdders(Scene_image_item* seg_img, const QString& name)
  {
    const CGAL::Image_3* img = seg_img->image();
    const Word* begin = (const Word*)img->data();
    const Word* end = (const Word*)img->data() + img->size();

    std::pair<const Word, const Word> minmax(*std::min_element(begin, end), *std::max_element(begin, end));

    Clamp_to_one_zero_range clamper = { minmax };

    switchReaderConverter< Word >(minmax);

    Volume_plane<x_tag> *x_item = new Volume_plane<x_tag>(img->image()->tx,img->image()->ty, img->image()->tz);
    Volume_plane<y_tag> *y_item = new Volume_plane<y_tag>(img->image()->tx,img->image()->ty, img->image()->tz);
    Volume_plane<z_tag> *z_item = new Volume_plane<z_tag>(img->image()->tx,img->image()->ty, img->image()->tz);

    x_item->setProperty("img", QVariant::fromValue((void*)seg_img));
    y_item->setProperty("img", QVariant::fromValue((void*)seg_img));
    z_item->setProperty("img", QVariant::fromValue((void*)seg_img));

    x_item->setColor(QColor("red"));
    y_item->setColor(QColor("green"));
    z_item->setColor(QColor("blue"));
    for(CGAL::QGLViewer* viewer : CGAL::QGLViewer::QGLViewerPool())
    {
      viewer->installEventFilter(x_item);
      viewer->installEventFilter(y_item);
      viewer->installEventFilter(z_item);
    }

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

    first_offset = Three::mainViewer()->offset();

  }

  template<typename T>
  void switchReaderConverter(std::pair<T, T> minmax)
  {
    switchReaderConverter(minmax, typename std::is_integral<T>::type());
  }

  template<typename T>
  void switchReaderConverter(std::pair<T, T> minmax, std::true_type)
  {
    // IntConverter
    IntConverter x = { minmax }; pxr_.setIC(x);
  }

  template<typename T>
  void switchReaderConverter(std::pair<T, T> minmax, std::false_type)
  {
    // IntConverter
    DoubleConverter x = { minmax }; pxr_.setFC(x);
  }

private Q_SLOTS:
  void select_plane(CGAL::Three::Scene_item* item)
  {
    Scene_image_item* img = (Scene_image_item*)item->property("img").value<void*>();
    if(img)
      scene->setSelectedItem(scene->item_id(img));
  }
  //updates the msgBox
  void update_msgBox()
  {
    static int nbPlanes =0;
    ++nbPlanes;
    msgBox.setText(QString("Planes created : %1/3").arg(nbPlanes));
    if(nbPlanes == 3)
    {
      const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
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
      for(CGAL::Three::Scene_item* key : group_map.keys())
      {
        if(group_map[key].group == group_item)
        {
          group_map[key].x_item = nullptr;
          group_map[key].y_item = nullptr;
          group_map[key].z_item = nullptr;
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

  // destroy planes on image deletion
  void on_img_detroyed()
  {
    Scene_image_item* img_item = qobject_cast<Scene_image_item*>(sender());
    if(img_item)
    {
      Scene_group_item* group = qobject_cast<Scene_group_item*>(group_map[img_item].group);
      if(!group)
        return;

      group_map[img_item].x_item = nullptr;
      group_map[img_item].y_item = nullptr;
      group_map[img_item].z_item = nullptr;
      disconnect(group_map[img_item].group, SIGNAL(aboutToBeDestroyed()),
                          this, SLOT(erase_group()));
      group_map.remove(img_item);

      QList<int> deletion;
      for(Scene_interface::Item_id id : group->getChildren())
      {
        Scene_item* child = group->getChild(id);
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
    bool show_sliders = true;

    // x line
    if(c.x_item != nullptr)
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
      show_sliders &= qobject_cast<Scene_image_item*>(sel_itm)->image()->xdim() > 1;
    }

    //y line
    if(c.y_item != nullptr)
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
      show_sliders &= qobject_cast<Scene_image_item*>(sel_itm)->image()->ydim() > 1;
    }

    // z line
    if(c.z_item != nullptr)
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
      show_sliders &= qobject_cast<Scene_image_item*>(sel_itm)->image()->zdim() > 1;
    }

      x_control->setEnabled(show_sliders);
      y_control->setEnabled(show_sliders);
      z_control->setEnabled(show_sliders);
  }

  // Keeps the position of the planes for the next time
  void set_value()
  {
    current_control->x_value = x_slider->value();
    current_control->y_value = y_slider->value();
    current_control->z_value = z_slider->value();
  }

  void destroy_x_item()
  {
    current_control->x_item = nullptr;
    if(group_map.isEmpty())
      x_control->hide();
  }
  void destroy_y_item()
  {
    current_control->y_item = nullptr;
    if(group_map.isEmpty())
      y_control->hide();
  }
  void destroy_z_item()
  {
    current_control->z_item = nullptr;
    if(group_map.isEmpty())
      z_control->hide();
  }
};

QString Io_image_plugin::nameFilters() const
{
  return QString("Inrimage files (*.inr *.inr.gz) ;; "
                 "Analyze files (*.hdr *.img *.img.gz) ;; "
                 "Stanford Exploration Project files (*.H *.HH) ;; "
                 "VTK image files (*.vti) ;; "
                 "NRRD image files (*.nrrd) ;; "
                 "NIFTI image files (*.nii *.nii.gz)");
}

bool Io_image_plugin::canLoad(QFileInfo) const
{
  return true;
}

template<typename Word>
void convert(Image* image)
{
  float *f_data = (float*)ImageIO_alloc(image->xdim()*image->ydim()*image->zdim()*sizeof(float));
  Word* d_data = (Word*)(image->data());

  // convert image from double to float
  for(std::size_t x = 0; x<image->xdim(); ++x) {
    for(std::size_t y = 0; y<image->ydim(); ++y) {
      for(std::size_t z = 0; z<image->zdim(); ++z)
      {
        std::size_t i =(z * image->ydim() + y) * image->xdim() + x;
        f_data[i] =(float)d_data[i];
      }
    }
  }

  ImageIO_free(d_data);
  image->image()->data = (void*)f_data;
  image->image()->wdim = 4;
  image->image()->wordKind = WK_FLOAT;
}


QList<Scene_item*>
Io_image_plugin::load(QFileInfo fileinfo, bool& ok, bool add_to_scene)
{
  ok = true;
  QApplication::restoreOverrideCursor();
  Image* image = new Image;

  QString warningMessage;
  // read a vti file
  if(fileinfo.suffix() == "vti")
  {
#ifdef CGAL_USE_VTK
    ok = load_vtk_file<vtkXMLImageDataReader>(fileinfo, image);
#else
    ok = false;
    warningMessage = "VTK is required to read VTI files";
#endif
  }

  // read a nrrd file
  else if(fileinfo.suffix() == "nrrd")
  {
#ifdef CGAL_USE_VTK
    ok = load_vtk_file<vtkNrrdReader>(fileinfo, image);
#else
    ok = false;
    warningMessage = "VTK is required to read NRRD files";
#endif
  }

  // read a NIFTI file
  else if(fileinfo.suffix() == "nii"
    || (   fileinfo.suffix() == "gz"
        && fileinfo.fileName().endsWith(QString(".nii.gz"), Qt::CaseInsensitive)))
  {
#ifdef CGAL_USE_VTK
    ok = load_vtk_file<vtkNIFTIImageReader>(fileinfo, image);
#else
    ok = false;
    warningMessage = "VTK is required to read NifTI files";
#endif
  }

  // read a sep file
  else if (fileinfo.suffix() == "H" || fileinfo.suffix() == "HH")
  {
    CGAL::SEP_to_ImageIO<float> reader(fileinfo.filePath().toUtf8().data());
    *image = *reader.cgal_image();
    is_gray = true;
  }

  else if(fileinfo.suffix() != "H" && fileinfo.suffix() != "HH" &&
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
    if(qmb.exec() == QMessageBox::Yes)
    {
      Raw_image_dialog raw_dialog;
      raw_dialog.label_file_size->setText(QString("%1 B").arg(fileinfo.size()));
      raw_dialog.buttonBox->button(QDialogButtonBox::Open)->setEnabled(false);
      if(raw_dialog.exec())
      {
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
                            raw_dialog.image_sign()))
        {
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
        }
        else
        {
          success = false;
        }
      }
      else
      {
        success = false;
      }
    }
    else
    {
      success = false;
    }

    if(!success)
    {
      ok = false;
    }
  }

  if (!ok) {
    if (warningMessage.length() > 0)
        CGAL::Three::Three::warning(warningMessage);
    delete image;
    return QList<Scene_item*>();
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
  for(int i=1 ; i<9; ++i)
  {
    QString s = tr("1:%1").arg(i*i*i);
    ui.precisionList->addItem(s);
  }

  // Adds Image type
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
      ok = false;
      return QList<Scene_item*>();
    }

    // Get selected precision
    voxel_scale = ui.precisionList->currentIndex() + 1;

    //Get the image type
    type = ui.imageType->currentText();
  }
  else
  {
    type = "Gray-level image";
  }

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
  {
    image_item = new Scene_image_item(image,voxel_scale, false);
  }
  image_item->setName(fileinfo.baseName());

  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(image_item);

  return QList<Scene_item*>() << image_item;
}

bool Io_image_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_image_item*>(item);
}

bool Io_image_plugin::loadDirectory(const QString& dirname,
                                    const Directory_extension_type ext)
{
#ifndef CGAL_USE_VTK
  QApplication::restoreOverrideCursor();
  CGAL::Three::Three::warning("VTK is required to read DCM and BMP files");
  CGAL_USE(dirname);
  CGAL_USE(ext);
  return false;
#else
  QFileInfo fileinfo;
  fileinfo.setFile(dirname);
  bool result = true;
  if(!fileinfo.isReadable())
  {
    QMessageBox::warning(mw, mw->windowTitle(),
                         tr("Cannot read directory <tt>%1</tt>!").arg(dirname));
    CGAL::Three::Three::warning(tr("Opening of directory %1 failed!").arg(dirname));
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

    // Adds Image type
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

    bool smooth = ui.smoothImage->isChecked();

    Image *image = createDirectoryImage(dirname, ext, smooth);
    if(image->image() == nullptr)
    {
      QMessageBox::warning(mw, mw->windowTitle(), tr("Error opening directory <tt>%1/</tt>!").arg(dirname));
      CGAL::Three::Three::warning(tr("Opening of directory %1/ failed!").arg(dirname));
      result = false;
    }
    else
    {
      CGAL::Three::Three::information(tr("Directory %1/ successfully opened.").arg(dirname));
    }

    if(result)
    {
      // Get the image type
      QString type = ui.imageType->currentText();
      Scene_image_item* image_item;
      if(type == "Gray-level image")
      {
        // Create planes
        image_item = new Scene_image_item(image,125, true);
        image_item->setName(fileinfo.baseName());
        msgBox.setText("Planes created : 0/3");
        msgBox.setStandardButtons(QMessageBox::NoButton);
        msgBox.show();

        scene->addItem(image_item);
        createPlanes(image_item);
      }
      else
      {
        // Get selected precision
        int voxel_scale = ui.precisionList->currentIndex() + 1;

        image_item = new Scene_image_item(image,voxel_scale, false);
        image_item->setName(fileinfo.baseName());
        scene->addItem(image_item);
      }
    }
  }

  QApplication::restoreOverrideCursor();
  return result;
#endif
}

Image* Io_image_plugin::createDirectoryImage(const QString& dirname,
                                             const Directory_extension_type ext,
                                             const bool smooth)
{
  Image* image = nullptr;

#ifndef CGAL_USE_VTK
  CGAL::Three::Three::warning("VTK is required to read DCM and BMP files");
  CGAL_USE(dirname);
  CGAL_USE(ext);
  CGAL_USE(smooth);
#else
  auto create_image = [&](auto&& reader) -> void
  {
    auto executive = vtkDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
    if(executive)
      executive->SetReleaseDataFlag(0, 0); // where 0 is the port index

    vtkImageData* vtk_image = nullptr;
    vtkNew<vtkImageGaussianSmooth> smoother; // must be here because it will own the vtk image

    if(smooth)
    {
      smoother->SetStandardDeviations(1., 1., 1.);
      smoother->SetInputConnection(reader->GetOutputPort());
      smoother->Update();
      vtk_image = smoother->GetOutput();
    }
    else
    {
      reader->Update();
      vtk_image = reader->GetOutput();
    }

    vtk_image->Print(std::cerr);

    *image = CGAL::IO::read_vtk_image_data(vtk_image); // copy the image data
  };

  image = new Image;
  if(ext == Directory_extension_type::DCM)
  {
    vtkNew<vtkDICOMImageReader> dicom_reader;
    dicom_reader->SetDirectoryName(dirname.toUtf8());
    create_image(dicom_reader);
  }
  else
  {
    CGAL_assertion(ext == Directory_extension_type::BMP);

    // vtkBMPReader does not provide SetDirectoryName()...
    std::vector<boost::filesystem::path> paths;
    vtkStringArray* files = vtkStringArray::New();
    boost::filesystem::path p(dirname.toUtf8().data());
    for(boost::filesystem::directory_entry& x : boost::filesystem::directory_iterator(p))
    {
      std::string s = x.path().string();
      if(CGAL::IO::internal::get_file_extension(s) != "bmp")
        continue;

      paths.push_back(x.path());
    }

    // boost::filesystem::directory_iterator does not guarantee a sorted order
    std::sort(std::begin(paths), std::end(paths));

    for(const boost::filesystem::path& p : paths)
      files->InsertNextValue(p.string());

    if(files->GetSize() == 0)
      return image;

    vtkNew<vtkBMPReader> bmp_reader;
    bmp_reader->SetFileNames(files);
    create_image(bmp_reader);
  }
#endif

  return image;
}

#include "Io_image_plugin.moc"
