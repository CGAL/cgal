#ifndef CGAL_VOLUME_PLANE_INTERFACE_H_
#define CGAL_VOLUME_PLANE_INTERFACE_H_

#include <QObject>
#include <QGLViewer/qglviewer.h>
#include <CGAL_demo/Scene_item.h>
#include <iostream>

class Volume_plane_interface : public Scene_item {
Q_OBJECT
public:
  virtual void init() = 0;

  virtual unsigned int aDim() const = 0;
  virtual unsigned int bDim() const = 0;
  virtual unsigned int cDim() const = 0;
  virtual qglviewer::Vec translationVector() const = 0;
  void itemAboutToBeDestroyed(Scene_item* item) {
    if(this == item) {
      emit planeDestructionIncoming(this);
      emit aboutToBeDestroyed();
    }
  }

signals:
  void planeDestructionIncoming(Volume_plane_interface*);
};


#endif /* CGAL_VOLUME_PLANE_INTERFACE_H_ */

