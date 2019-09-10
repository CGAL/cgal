#ifndef SCENE_H
#define SCENE_H

#include "typedefs.h"
#include <QMap>
#include <CGAL/Qt/qglviewer.h>

class Scene {

  friend class Viewer;

public:
  Scene() {}
  ~Scene() {  eraseOldData();  }

public:
  inline void setViewer(CGAL::QGLViewer* v) {  m_viewer = v;  }
  inline void showError(const QString & msg) {
    if(!m_viewer) m_viewer->displayMessage( msg );
  }
  inline bool isDTEmpty() {  return m_dt.number_of_vertices()==0;  }
  inline void eraseOldData() {  m_dt.clear();  m_vhArray.clear();  }

public:
  void generatePoints(int);
  void loadPointsOFF(const char*);
  void loadPointsXYZ(const char*);
  void savePointsOFF(const char*);
  void savePointsXYZ(const char*);

  void readOFFPointsandFacets(const char*, std::list<Point_3> &);

private:
  //added for T3 demo
  DT3 m_dt;
  QList<Vertex_handle> m_vhArray;

  CGAL::QGLViewer* m_viewer;
};

#endif
