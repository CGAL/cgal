#ifndef CGAL_VOLUME_PLANE_INTERFACE_H_
#define CGAL_VOLUME_PLANE_INTERFACE_H_

#include <QObject>
#include <QGLViewer/qglviewer.h>
#include <CGAL/Three/Scene_item.h>
#include <iostream>
#include <QGLViewer/manipulatedFrame.h>
using namespace CGAL::Three;
class Volume_plane_interface : public Scene_item {
Q_OBJECT
public:
  Volume_plane_interface(qglviewer::ManipulatedFrame* f) : mFrame_(f) { 
    connect(mFrame_, SIGNAL(manipulated()), this, SLOT(propagateManipulation()));
  }

  virtual ~Volume_plane_interface() {
    delete mFrame_;
  }

  virtual void init() = 0;

  virtual void setData(unsigned int adim, unsigned int bdim, unsigned int cdim, float xscale, float yscale, float zscale, std::vector<float> &colors) =0;
  virtual unsigned int aDim() const = 0;
  virtual unsigned int bDim() const = 0;
  virtual unsigned int cDim() const = 0;
  virtual qglviewer::Vec translationVector() const = 0;
  void itemAboutToBeDestroyed(CGAL::Three::Scene_item* item) {
    if(this == item) {
      Q_EMIT planeDestructionIncoming(this);
      Q_EMIT aboutToBeDestroyed();
    }
  }

  virtual unsigned int getCurrentCube() const = 0;

  virtual qglviewer::ManipulatedFrame* manipulatedFrame() { return mFrame_; }

Q_SIGNALS:
  void planeDestructionIncoming(Volume_plane_interface*);
  void manipulated(int);
private Q_SLOTS:
  void propagateManipulation() {
    Q_EMIT manipulated(getCurrentCube());
  }
protected:
  qglviewer::ManipulatedFrame* mFrame_;
};


#endif /* CGAL_VOLUME_PLANE_INTERFACE_H_ */

