#include "config.h"

#include "Scene_c3t3_item.h"

#include <QVector>
#include <QColor>
#include <QPixmap>
#include <QPainter>
#include <QtCore/qglobal.h>

#include <map>
#include <vector>
#include <CGAL/gl.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>

#include <CGAL/Three/Scene_interface.h>

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

struct Scene_c3t3_item_priv {
  Scene_c3t3_item_priv() : c3t3() {}
  Scene_c3t3_item_priv(const C3t3& c3t3_) : c3t3(c3t3_) {}

  C3t3 c3t3;
  QVector<QColor> colors;
};
double complex_diag(const Scene_item* item) {
  const Scene_item::Bbox& bbox = item->bbox();
  const double& xdelta = bbox.xmax-bbox.xmin;
  const double& ydelta = bbox.ymax-bbox.ymin;
  const double& zdelta = bbox.zmax-bbox.zmin;
  const double diag = std::sqrt(xdelta*xdelta +
                                ydelta*ydelta +
                                zdelta*zdelta);
  return diag * 0.7;
}




