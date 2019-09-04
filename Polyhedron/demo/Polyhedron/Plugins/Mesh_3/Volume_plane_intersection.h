#ifndef CGAL_VOLUME_PLANE_INTERSECTION_H_
#define CGAL_VOLUME_PLANE_INTERSECTION_H_

#include <CGAL/Three/Scene_item_rendering_helper.h>

#include <QColor>
#include <QString>
#include <CGAL/Qt/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <CGAL/Three/Viewer_interface.h>
using namespace CGAL::Three;
class Volume_plane_interface;
struct Volume_plane_intersection_priv;
class Volume_plane_intersection
  : public Scene_item_rendering_helper {
  typedef std::pair<Volume_plane_interface*, Volume_plane_interface*> Interface_pair;
Q_OBJECT
public:
  Volume_plane_intersection(float x, float y, float z, float tx, float ty, float tz);
  ~Volume_plane_intersection();

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  bool manipulatable() const { return false; }
  Volume_plane_intersection* clone() const { return 0; }
  bool supportsRenderingMode(RenderingMode) const { return true; }
  QString toolTip() const { return "Tooling"; }
  void compute_bbox() const{}

  void draw(Viewer_interface*)const;

  void setX(Volume_plane_interface* x);
  void setY(Volume_plane_interface* x);
  void setZ(Volume_plane_interface* x);
  void initializeBuffers(Viewer_interface *) const;
  void computeElements() const;
  void invalidateOpenGLBuffers();
public Q_SLOTS:
  void planeRemoved(Volume_plane_interface* i);
protected:
  friend struct Volume_plane_intersection_priv;
  Volume_plane_intersection_priv *d;
};

#endif /* CGAL_VOLUME_PLANE_INTERSECTION_H_ */

