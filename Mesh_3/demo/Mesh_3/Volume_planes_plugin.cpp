#include <CGAL/Image_3.h>

#include "Volume_plane.h"
#include "Volume_plane_thread.h"
#include "Volume_plane_intersection.h"

#include "Scene_segmented_image_item.h"
#include <CGAL_demo/Plugin_helper.h>
#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Messages_interface.h>
#include <CGAL_demo/Scene_interface.h>
#include <CGAL_demo/Scene_item.h>

#include <QAction>
#include <QMenu>
#include <QList>
#include <QInputDialog>
#include <QSlider>
#include <QLabel>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QDockWidget>
#include <QMainWindow>
#include <QMessageBox>
#include <QString>

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <map>
#include <string>

class Plane_slider : public QSlider
{
  Q_OBJECT
public:
  Plane_slider(const qglviewer::Vec& v, int id, Scene_interface* scene,
               qglviewer::ManipulatedFrame* frame, Qt::Orientation ori, QWidget* widget) 
    : QSlider(ori, widget), v(v), id(id), scene(scene), frame(frame) { 
    this->setTracking(true);
    connect(frame,  SIGNAL(manipulated()), this, SLOT(updateValue()));
  }

public:
  void sliderChange(SliderChange c) {
    QSlider::sliderChange(c);
    if(c == SliderValueChange) {
      qglviewer::Vec v2 = v * (this->value() / scale);
      frame->setTranslationWithConstraint(v2);
      scene->itemChanged(id);
    }

    emit realChange(this->value() / scale);
  }

public slots:
  void updateValue() {
    float a, b, c;
    frame->getPosition(a, b, c);
    float sum1 = a + b + c;
    float sum2 = v.x + v.y + v.z;
    sum1 /= sum2;
    setValue(sum1 * scale);
  }

signals:
  void realChange(int);

private:
  const static unsigned int scale;

  qglviewer::Vec v;
  int id;
  Scene_interface* scene;
  qglviewer::ManipulatedFrame* frame;
};

const unsigned int Plane_slider::scale = 100;

class Volume_plane_plugin :
  public QObject,
  public Plugin_interface 
{
  Q_OBJECT
  Q_INTERFACES(Plugin_interface);
public:
  Volume_plane_plugin() : planeSwitch(NULL), sc(NULL), mw(NULL) {
  }
  
  virtual ~Volume_plane_plugin() { }

  virtual void init(QMainWindow* mw, Scene_interface* sc) {
    assert(mw != NULL);
    assert(sc != NULL);
    this->sc = sc;
    this->mw = mw;

    QList<QMenu*> menus = mw->findChildren<QMenu*>();

    planeSwitch = new QAction(mw);
    planeSwitch->setText("Add Volume Planes");

    connect(planeSwitch, SIGNAL(triggered()), this, SLOT(selectPlanes()));
  }

  virtual QList<QAction*> actions() const {
    return QList<QAction*>() << planeSwitch;
  }

public slots:
  void selectPlanes() {
    std::vector< Scene_segmented_image_item* > seg_items;
    Scene_segmented_image_item* seg_img = NULL;

    for(unsigned int i = 0; i < sc->numberOfEntries(); ++i) {
      Scene_segmented_image_item* tmp = qobject_cast<Scene_segmented_image_item*>(sc->item(i));
      if(tmp != NULL)
        seg_items.push_back(tmp);
    }
    
    if(seg_items.empty()) {
      QMessageBox::warning(mw, tr("No suitable item found"), tr("Load an inrimage or hdr file to enable Volume Planes."));
      return;
    } else {
      QList<QString> items;
      for(std::vector< Scene_segmented_image_item* >::const_iterator it = seg_items.begin(); 
          it != seg_items.end(); ++it) { 
        items << (*it)->name();
      }

      bool ok;
      QString selected = QInputDialog::getItem(mw, tr("Select a dataset:"), tr("Items"), items, 0, false, &ok);
      
      if(!ok || selected.isEmpty())
        return;

      for(std::vector< Scene_segmented_image_item*>::const_iterator it = seg_items.begin(); 
          it != seg_items.end(); ++it) { 
        if(selected == (*it)->name())
          seg_img = *it;
      }
    }
    
    if(!(seg_img == NULL)) {
      const CGAL::Image_3* img = seg_img->image();
      CGAL_IMAGE_IO_CASE(img->image(), this->launchAdders<Word>(img, seg_img->name()))

      Volume_plane_intersection* i = new Volume_plane_intersection(img->xdim() * img->vx(), 
                                                                   img->ydim() * img->vy(), 
                                                                   img->zdim() * img->vz());
      this->intersectionId = sc->addItem(i);
    } else {
      QMessageBox::warning(mw, tr("Something went wrong"), tr("Selected a suitable Object but couldn't get an image pointer."));
      return;
    }

  }

  void addVP(Volume_plane_thread* thread) {
    Volume_plane_interface* plane = thread->getItem();
    plane->init();

    // add the interface for this Volume_plane
    int id = sc->addItem(plane);
    
    QLayout* layout = createOrGetDockLayout();
    
    QWidget* controls = new QWidget;
    QHBoxLayout* box = new QHBoxLayout(controls);
    layout->addWidget(controls);

    QLabel* label = new QLabel(controls);
    label->setText(plane->name());

    QLabel* cubeLabel = new QLabel(controls);
    cubeLabel->setNum(static_cast<int>(plane->getCurrentCube()));

    QSlider* slider = new Plane_slider(plane->translationVector(), id, sc, plane->manipulatedFrame(), Qt::Horizontal, controls);
    slider->setRange(0, (plane->cDim() - 1) * 100);

    connect(slider, SIGNAL(realChange(int)), cubeLabel, SLOT(setNum(int)));
    connect(plane, SIGNAL(manipulated(int)), cubeLabel, SLOT(setNum(int)));
    
    box->addWidget(cubeLabel);
    box->addWidget(label);
    box->addWidget(slider);
    
    connect(plane, SIGNAL(aboutToBeDestroyed()), controls, SLOT(deleteLater()));

    std::vector<Volume_plane_thread*>::iterator it = std::find(threads.begin(), threads.end(), thread);

    // this slot has been connected to a thread that hasn't been
    // registered here.
    assert(it != threads.end());
    delete *it;
    threads.erase(it);

    Volume_plane_intersection* intersection = dynamic_cast<Volume_plane_intersection*>(sc->item(intersectionId));
    if(!intersection) {
      // the intersection is gone before it was initialized
      return;
    }
    // FIXME downcasting mode
    // FIXME this will bug if two volume planes are generated simultaneously by the plug
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
 
private:
  QAction* planeSwitch;
  Scene_interface* sc;
  QMainWindow* mw;
  std::vector<Volume_plane_thread*> threads;
  unsigned int intersectionId;

  QLayout* createOrGetDockLayout() {
    QLayout* layout = NULL;
    QDockWidget* controlDockWidget = mw->findChild<QDockWidget*>("volumePlanesControl");
    if(!controlDockWidget) {
      controlDockWidget = new QDockWidget(mw);
      controlDockWidget->setObjectName("volumePlanesControl");
      QWidget* content = new QWidget(controlDockWidget);
      layout = new QVBoxLayout(content);
      layout->setObjectName("vpSliderLayout");
      controlDockWidget->setWindowTitle("Control Widget");
      controlDockWidget->setWidget(content);
      mw->addDockWidget(Qt::LeftDockWidgetArea, controlDockWidget);
    } else {
      layout = controlDockWidget->findChild<QLayout*>("vpSliderLayout");
    }

    return layout;
  }

  template<typename Word>
  void launchAdders(const CGAL::Image_3* img, const QString& name) {
    const Word* begin = (const Word*)img->data();
    const Word* end = (const Word*)img->data() + img->size();

    Clamp_to_one_zero_range clamper = { std::make_pair(*std::max_element(begin, end),
                                                       *std::min_element(begin, end)) };

    threads.push_back(new X_plane_thread<Word>(img, clamper, name));
    connect(threads.back(), SIGNAL(finished(Volume_plane_thread*)), this, SLOT(addVP(Volume_plane_thread*)));
    threads.back()->start();
    
    threads.push_back(new Y_plane_thread<Word>(img, clamper, name));
    connect(threads.back(), SIGNAL(finished(Volume_plane_thread*)), this, SLOT(addVP(Volume_plane_thread*)));
    threads.back()->start();

    threads.push_back(new Z_plane_thread<Word>(img, clamper, name));
    connect(threads.back(), SIGNAL(finished(Volume_plane_thread*)), this, SLOT(addVP(Volume_plane_thread*)));
    threads.back()->start();

  }
};


Q_EXPORT_PLUGIN2(Volume_plane_plugin, Volume_plane_plugin)

#include "Volume_planes_plugin.moc"
