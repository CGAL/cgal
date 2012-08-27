#ifndef CGAL_VOLUME_PLANE_INTERSECTION_H_
#define CGAL_VOLUME_PLANE_INTERSECTION_H_

#include <CGAL_demo/Scene_item.h>

#include <QColor>
#include <QString>

class Volume_plane_interface;

class Volume_plane_intersection
  : public Scene_item {
  typedef std::pair<Volume_plane_interface*, Volume_plane_interface*> Interface_pair;
Q_OBJECT
public:
  Volume_plane_intersection(float x, float y, float z)
    : a(NULL), b(NULL), c(NULL), x(x), y(y), z(z) {
    setColor(QColor(255, 0, 0));
    setName("Volume plane intersection");
  }

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  bool manipulatable() const { return false; }
  Volume_plane_intersection* clone() const { return 0; }
  bool supportsRenderingMode(RenderingMode) const { return true; }
  QString toolTip() const { return "Tooling"; }

  void draw() const;

  void setX(Volume_plane_interface* x) { a = x; }
  void setY(Volume_plane_interface* x) { b = x; }
  void setZ(Volume_plane_interface* x) { c = x; }

public slots:
  void planeRemoved(Volume_plane_interface* i) {
    if(a == i) {
      a = NULL;
    } else if(b == i) {
      b = NULL;
    } else if(c == i) {
      c = NULL;
    }
  }

private:
  Volume_plane_interface *a, *b, *c;
  float x, y, z;
};

#endif /* CGAL_VOLUME_PLANE_INTERSECTION_H_ */

