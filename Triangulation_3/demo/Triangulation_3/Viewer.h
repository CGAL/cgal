#ifndef VIEWER_H
#define VIEWER_H

#include "Scene.h"
#include <QGLViewer/qglviewer.h>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QSettings>
#include "PreferenceDlg.h"

#include <QOpenGLFunctions_2_1>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <CGAL/Qt/CreateOpenGLContext.h>
#include <iostream>
using namespace qglviewer;

class MainWindow;

class Viewer : public QGLViewer, QOpenGLFunctions_2_1 {

  Q_OBJECT

public:
  Viewer(QWidget* parent)
    : QGLViewer(CGAL::Qt::createOpenGLContext(),parent)
    , m_showAxis(false)
    , m_showVertex(true)
    , m_showDEdge(true)
    , m_showVEdge(false)
    , m_showFacet(false)
    , m_isFlat(false)
    , m_fRadius(1.)
    , m_curMode(NONE)
    , m_isPress(false)
    , m_isMoving(false)
    , m_hasNewPt(false)
    , m_selMode(NORMAL)
    , m_nearestNb(NULL)
    , m_hasEmptyS(false)
    , m_showTrackball(true)
    , m_pDlgPrefer(NULL)
    {
      pos_emptyFacet = new std::vector<float>();
      pos_emptySphere= new std::vector<float>();
      points_emptySphere = new std::vector<float>();
      pos_points = new std::vector<float>();
      pos_newPoint = new std::vector<float>();
      pos_selectedVertex = new std::vector<float>();
      pos_movingPoint = new std::vector<float>();
      pos_queryPoint = new std::vector<float>();
      pos_trackBall  = new std::vector<float>();
      pos_voronoi = new std::vector<float>();
      pos_delaunay = new std::vector<float>();
      pos_facets = new std::vector<float>();
      pos_newFacet = new std::vector<float>();
      pos_nearest_neighbor = new std::vector<float>();
      points_locationSphere = new std::vector<float>();
      points_cylinder = new std::vector<float>();
      normals_cylinder = new std::vector<float>();
      points_sphere = new std::vector<float>();
      points_trackBall = new std::vector<float>();
      normals_sphere = new std::vector<float>();
      normals_emptySphere = new std::vector<float>();
      normals_trackBall= new std::vector<float>();
      transfo1_voronoi = new std::vector<float>();
      transfo2_voronoi = new std::vector<float>();
      transfo3_voronoi = new std::vector<float>();
      transfo4_voronoi = new std::vector<float>();
      transfo1_delaunay = new std::vector<float>();
      transfo2_delaunay = new std::vector<float>();
      transfo3_delaunay = new std::vector<float>();
      transfo4_delaunay = new std::vector<float>();
      incremental_points = new std::vector<float>();
      incremental_next_point = new std::vector<float>();
      incremental_facet = new std::vector<float>();
      incremental_conflict = new std::vector<float>();
    }

  ~Viewer()
  {
      for(int i=0; i< vboSize; i++)
          buffers[i].destroy();

      for(int i=0; i< vaoSize; i++)
          vao[i].destroy();


  }
  enum Mode { NONE, INSERT_V, INSERT_PT, MOVE, SELECT, FINDNB, EMPTYSPH };

public:
  inline void setScene(Scene* pScene) { m_pScene = pScene; }

  // set current mode
  inline void setMode(Mode m) {
    m_curMode = m;
    m_isMoving = false;
    m_hasEmptyS = false;
    m_nearestNb = NULL;
    updateGL();
  }

  // set selectBuffer size (if necessary)
  inline void setSelBuffSize() {
    // Default selectBuffer size is 4000
    // (i.e. 1000 objects in selection region, since each object pushes 4 values).
    if( m_pScene->m_vhArray.size() > 900 )
      // The previous selectBuffer is deleted and a new one is created.
      setSelectBufferSize( 4*(m_pScene->m_vhArray.size() + 100) );
  }

  void readSettings() {
    // read from an .ini file
    QSettings settings("settings.ini", QSettings::IniFormat);
    // QVariant value ( const QString & key, const QVariant & defaultValue = QVariant() )
    // Because QVariant is part of the QtCore library,
    //   it cannot provide conversion functions to data types such as QColor and QImage,
    //   which are part of QtGui.
    // In other words, there is no toColor(), toImage(), or toPixmap() functions in QVariant.
    // Instead, use the QVariant::value() or the qVariantValue() template function
    m_colorVertex = settings.value( "Show/vertexcolor", QColor(255, 150, 0) ).value<QColor>();
    m_fSizeVertex = settings.value( "Show/vertexsize", 0.04f ).toFloat();
    m_colorDEdge = settings.value( "Show/dedgecolor", QColor(0, 255, 0) ).value<QColor>();
    m_fSizeDEdge = settings.value( "Show/dedgesize", 0.01f ).toFloat();
    m_colorVEdge = settings.value( "Show/vedgecolor", QColor(0, 0, 255) ).value<QColor>();
    m_fSizeVEdge = settings.value( "Show/vedgesize", 0.01f ).toFloat();
    m_colorFacet = settings.value( "Show/facetcolor",
                             QColor(255, 255, 0, 96) ).value<QColor>();
    m_colorTrackball = settings.value( "Show/ballcolor",
                                 QColor(150, 150, 150, 128) ).value<QColor>();
    m_iStep = settings.value( "Show/ballstep", 4000 ).toInt();
    m_colorEmptySphere = settings.value( "Show/spherecolor",
                                   QColor(180, 50, 180, 64) ).value<QColor>();
  }

  void writeSettings() {
    // write to an .ini file
    QSettings settings("settings.ini", QSettings::IniFormat);
    // The inverse conversion (e.g., from QColor to QVariant) is automatic
    //   for all data types supported by QVariant, including GUI-related types
    settings.setValue("Show/vertexcolor", m_colorVertex);
    settings.setValue("Show/vertexsize", m_fSizeVertex);
    settings.setValue("Show/dedgecolor", m_colorDEdge);
    settings.setValue("Show/dedgesize", m_fSizeDEdge);
    settings.setValue("Show/vedgecolor", m_colorVEdge);
    settings.setValue("Show/vedgesize", m_fSizeVEdge);
    settings.setValue("Show/facetcolor", m_colorFacet);
    settings.setValue("Show/ballcolor", m_colorTrackball);
    settings.setValue("Show/ballstep", m_iStep);
    settings.setValue("Show/spherecolor", m_colorEmptySphere);
  }

public Q_SLOTS:
  // clear scene
  void changed()
  {
      compute_elements();
      are_buffers_initialized = false;
  }
  void clear() {
    m_pScene->eraseOldData();
    m_hasNewPt = false;
    m_boundaryFacets.clear();
    m_conflictCells.clear();
    m_vidSeled.clear();
    m_isMoving = false;
    m_nearestNb = NULL;
    m_hasEmptyS = false;
    if( !m_incrementalPts.isEmpty() ) {
      Q_EMIT( stopIncAnimation() );
      m_incrementalPts.clear();
    }
  }

  // play/pause incremental construction
  void toggleIncremental(bool on);
  // clean up old data and information
  void initClean() {
    m_pScene->eraseOldData();
    m_hasNewPt = false;
    m_boundaryFacets.clear();
    m_conflictCells.clear();
    m_vidSeled.clear();
    m_isMoving = false;
    m_nearestNb = NULL;
    m_hasEmptyS = false;
  }
  // stop incremental construction
  void stopIncremental();
  // incremental insert a vertex (invoked by Timer)
  void incremental_insert();

  // show options
  inline void toggleShowAxis(bool flag)  { m_showAxis = flag; updateGL(); }
  inline void toggleShowVertex(bool flag)  { m_showVertex = flag; updateGL(); }
  inline void toggleShowDEdge(bool flag)  { m_showDEdge = flag; updateGL(); }
  inline void toggleShowVEdge(bool flag)  { m_showVEdge = flag; updateGL(); }
  inline void toggleShowFacet(bool flag)  { m_showFacet = flag; updateGL(); }
  inline void toggleFlat(bool flag)  { m_isFlat = flag; updateGL(); }

  // set preferences
  void setPreferences() {
    if (!m_pDlgPrefer) {
      m_pDlgPrefer = new PreferenceDlg(this);
      m_pDlgPrefer->init( m_colorVertex, m_fSizeVertex, m_colorDEdge, m_fSizeDEdge, m_colorVEdge, m_fSizeVEdge,
              m_colorFacet, m_colorTrackball, m_colorEmptySphere, m_iStep/40 ); // 5*8, 5 degrees of wheel
      connect(m_pDlgPrefer, SIGNAL(applyChanges()), this, SLOT(acceptChanges()) );
    }
    m_pDlgPrefer->show();
    m_pDlgPrefer->raise();
    m_pDlgPrefer->activateWindow();
  }

  void acceptChanges() {
    m_colorVertex = m_pDlgPrefer->m_colorVertex;
    m_fSizeVertex = m_pDlgPrefer->m_fSizeVertex;
    m_colorDEdge = m_pDlgPrefer->m_colorDEdge;
    m_fSizeDEdge = m_pDlgPrefer->m_fSizeDEdge;
    m_colorVEdge = m_pDlgPrefer->m_colorVEdge;
    m_fSizeVEdge = m_pDlgPrefer->m_fSizeVEdge;
    m_colorFacet = m_pDlgPrefer->m_colorFacet;
    m_colorTrackball = m_pDlgPrefer->m_colorTrackball;
    m_iStep = m_pDlgPrefer->m_iStep*40;
    m_colorEmptySphere = m_pDlgPrefer->m_colorEmptySphere;
    // redraw
    updateGL();
  }

  Q_SIGNALS:
  void stopIncAnimation();

// overloading QGLViewer virtual functions
protected:
  // initialize Viewer OpenGL context
  // Note: the default implement is empty and this is overloading.
  void init();
  // draw points, segments, and polygons
  void draw();

  // customize selection process
  void drawWithNames();
  void endSelection(const QPoint& point);
  // customize mouse events
  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *event);
  // customize key events
  void keyPressEvent(QKeyEvent *event);

  // customize help message
  QString helpString() const;

private:
  // draw a 3d effect vertex
  void drawVertex(const Point_3& p, std::vector<float> *vertices);
  // draw a 3d effect edge
  void drawEdge(const Point_3& from, const Point_3 &to, std::vector<float> *vertices);
  // draw a facet
  void drawFacet(const Triangle_3& t, std::vector<float> *vertices);
  // test whether the give 3D point is on the sphere
  inline bool isOnSphere( const Point_3 & pt ) {
    return ( (pt.x()*pt.x() + pt.y()*pt.y() + pt.z()*pt.z()) == (m_fRadius*m_fRadius) );
  }

  // compute the intersection point with the sphere
  bool computeIntersect( const QPoint & pos, Vec & pt );
  // compute the conflict region
  void computeConflict( Point_3 pt );


private:
  Scene* m_pScene;
  // play timer
  QTimer* m_pTimer;
  Point_3 m_curIncPt;
  QList<Point_3> m_incrementalPts;
  enum Step { INIT, NEWPT, CELL, CONFLICT };
  Step m_curStep;
  Cell_handle m_cellContain;
  QList<Facet> m_boundaryFacets;
  QList<Cell_handle> m_conflictCells;
  // show options
  bool m_showAxis;
  bool m_showVertex;
  bool m_showDEdge;
  bool m_showVEdge;
  bool m_showFacet;
  bool m_isFlat;
  // trackball
  float m_fRadius;
  // mode
  Mode m_curMode;
  bool m_isPress;
  bool m_isMoving;
  // insert point
  bool m_hasNewPt;
  Point_3 m_newPt;
  // select vertex
  enum selectionMode { NORMAL, ADD };
  selectionMode m_selMode;
  QRect m_rectSel;
  QList<int> m_vidSeled;
  // move vertex/point
  int m_vidMoving;
  // nearest neighbor
  Point_3 m_queryPt;
  Vertex_handle m_nearestNb;
  // empty sphere
  bool m_hasEmptyS;
  bool m_showTrackball;
  Point_3 m_centerPt;
  float m_fREmptyS;
  // change colors
  PreferenceDlg* m_pDlgPrefer;
  float m_fSizeVertex;
  float m_fSizeDEdge;
  float m_fSizeVEdge;
  QColor m_colorVertex;
  QColor m_colorDEdge;
  QColor m_colorVEdge;
  QColor m_colorFacet;
  QColor m_colorTrackball;
  QColor m_colorEmptySphere;
  // trackball resizing fineness
  int m_iStep;


  QColor color;
  static const int vaoSize = 29;
  static const int vboSize = 33;
  // define material
   QVector4D	ambient;
   QVector4D	diffuse;
   QVector4D	specular;
   GLfloat      shininess ;
      int poly_vertexLocation[3];
      int normalsLocation[3];
      int mvpLocation[3];
      int mvLocation[2];
      int centerLocation[5];
      int colorLocation[3];
      int lightLocation[5*2];

      bool are_buffers_initialized;
      bool extension_is_found;

      std::vector<float> *pos_emptyFacet;
      std::vector<float> *pos_emptySphere;
      std::vector<float> *points_emptySphere;
      std::vector<float> *pos_points;
      std::vector<float> *pos_newPoint;
      std::vector<float> *pos_selectedVertex;
      std::vector<float> *pos_movingPoint;
      std::vector<float> *pos_queryPoint;
      std::vector<float> *pos_trackBall;
      std::vector<float> *points_trackBall;
      std::vector<float> *pos_voronoi;
      std::vector<float> *pos_delaunay;
      std::vector<float> *pos_facets;
      std::vector<float> *pos_newFacet;
      std::vector<float> *pos_nearest_neighbor;
      std::vector<float> *points_locationSphere;
      std::vector<float> *points_cylinder;
      std::vector<float> *normals_cylinder;
      std::vector<float> *points_sphere;
      std::vector<float> *normals_sphere;
      std::vector<float> *normals_emptySphere;
      std::vector<float> *normals_trackBall;
      std::vector<float> *transfo1_voronoi;
      std::vector<float> *transfo2_voronoi;
      std::vector<float> *transfo3_voronoi;
      std::vector<float> *transfo4_voronoi;
      std::vector<float> *transfo1_delaunay;
      std::vector<float> *transfo2_delaunay;
      std::vector<float> *transfo3_delaunay;
      std::vector<float> *transfo4_delaunay;
      std::vector<float> *incremental_points;
      std::vector<float> *incremental_next_point;
      std::vector<float> *incremental_facet;
      std::vector<float> *incremental_conflict;

      QOpenGLBuffer buffers[vboSize];
      QOpenGLVertexArrayObject vao[vaoSize];
      QOpenGLShaderProgram rendering_program;
      QOpenGLShaderProgram rendering_program_spheres;
      QOpenGLShaderProgram rendering_program_cylinders;
      typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
      typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
      PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
      PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;
      void initialize_buffers();
      void compute_elements();
      void attrib_buffers(QGLViewer*);
      void compile_shaders();
      void draw_cylinder(float R, int prec, std::vector<float> *vertices, std::vector<float> *normals);
      void draw_sphere(float R, int prec, std::vector<float> *vertices, std::vector<float> *normals);
};

#endif
