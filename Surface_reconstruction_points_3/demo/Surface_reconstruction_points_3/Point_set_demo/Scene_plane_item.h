#ifndef SCENE_PLANE_ITEM_H
#define SCENE_PLANE_ITEM_H

#include "Scene_item.h"
#include "Scene_interface.h"

#include "Scene_basic_objects_config.h"

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

#include <cmath>

#include "Kernel_type.h"

class SCENE_BASIC_OBJECTS_EXPORT Scene_plane_item 
  : public Scene_item
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_plane_item(const Scene_interface* scene_interface) 
    : scene(scene_interface),
      manipulable(false),
      can_clone(true),
      frame(new ManipulatedFrame())
  {
    setNormal(0., 0., 1.);
  }

  ~Scene_plane_item() {
    delete frame;
  }

  bool isFinite() const { return false; }
  bool isEmpty() const { return false; }
  Bbox bbox() const { return Bbox(); }

  bool manipulatable() const {
    return manipulable;
  }
  ManipulatedFrame* manipulatedFrame() {
    return frame;
  }

  Scene_plane_item* clone() const {
    if(can_clone)
    {
      Scene_plane_item* item = new Scene_plane_item(scene);
      item->manipulable = manipulable;
      item->can_clone = true;
      item->frame = new ManipulatedFrame;
      item->frame->setPosition(frame->position());
      item->frame->setOrientation(frame->orientation());
      return item;
    }
    else
      return 0;
  }

  QString toolTip() const {
    const qglviewer::Vec& pos = frame->position();
    const qglviewer::Vec& n = frame->orientation().axis(); 
    return
      tr("<p><b>%1</b> (mode: %2, color: %3)<br />")
      .arg(this->name())
      .arg(this->renderingModeName())
      .arg(this->color().name())

      +
      tr("<i>Plane</i></p>"
         "<p>Equation: %1*x + %2*y + %3*z + %4 = 0<br />"
         "Normal vector: (%1, %2, %3)<br />"
         "Point: (%5, %6, %7)</p>")
      .arg(n[0]).arg(n[1]).arg(n[2])
      .arg( - pos * n)
      .arg(pos[0]).arg(pos[1]).arg(pos[2])
      +
      tr("<p>Can clone: %1<br />"
         "Manipulatable: %2</p>")
      .arg(can_clone?tr("true"):tr("false"))
      .arg(manipulable?tr("true"):tr("false"));
  }

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe || m == Flat); 
  }

  // Flat OpenGL drawing
  void draw() const {
    const double diag = scene_diag();
    ::glPushMatrix();
    ::glMultMatrixd(frame->matrix());
    GLboolean lighting;
    ::glGetBooleanv(GL_LIGHTING, &lighting);
    ::glDisable(GL_LIGHTING);
    ::glBegin(GL_POLYGON);
    ::glVertex3d(-diag, -diag, 0.);
    ::glVertex3d(-diag,  diag, 0.);
    ::glVertex3d( diag,  diag, 0.);
    ::glVertex3d( diag, -diag, 0.);
    ::glEnd();
    if(lighting)
      ::glEnable(GL_LIGHTING);
    ::glPopMatrix();
  };

  // Wireframe OpenGL drawing
  void draw_edges() const {
    ::glPushMatrix();
    ::glMultMatrixd(frame->matrix());
    QGLViewer::drawGrid((float)scene_diag());
    ::glPopMatrix();
  }

  Plane_3 plane() const {
    const qglviewer::Vec& pos = frame->position();
    const qglviewer::Vec& n = 
      frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
    return Plane_3(n[0], n[1],  n[2], - n * pos);
  }

private:
  double scene_diag() const {
    const Bbox& bbox = scene->bbox();
    const double& xdelta = bbox.xmax-bbox.xmin;
    const double& ydelta = bbox.ymax-bbox.ymin;
    const double& zdelta = bbox.zmax-bbox.zmin;
    const double diag = std::sqrt(xdelta*xdelta + 
                            ydelta*ydelta +
                            zdelta*zdelta);
    return diag * 0.7;
  }

public slots:
  void setPosition(float x, float y, float z) {
    frame->setPosition(x, y, z);
  }
  
  void setPosition(double x, double y, double z) {
    frame->setPosition((float)x, (float)y, (float)z);
  }
  
  void setNormal(float x, float y, float z) {
    frame->setOrientation(x, y, z, 0.f);
  }

  void setNormal(double x, double y, double z) {
    frame->setOrientation((float)x, (float)y, (float)z, 0.f);
  }

  void setClonable(bool b = true) {
    can_clone = b;
  }

  void setManipulatable(bool b = true) {
    manipulable = b;
  }
private:
  const Scene_interface* scene;
  bool manipulable;
  bool can_clone;
  qglviewer::ManipulatedFrame* frame;
};

#endif // SCENE_PLANE_ITEM_H
