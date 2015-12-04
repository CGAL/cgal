//! \file Viewer.h

#ifndef VIEWER_H
#define VIEWER_H

#include <CGAL/Three/Viewer_config.h>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLShaderProgram>
#include <CGAL/Three/Viewer_interface.h>

#include <QGLViewer/qglviewer.h>
#include <QPoint>
#include <QFont>
#include <QOpenGLFramebufferObject>

// forward declarations
class QWidget;
namespace CGAL{
namespace Three{
class Scene_draw_interface;
}
}
class QMouseEvent;
class QKeyEvent;
class QContextMenuEvent;

class Viewer_impl;
//!The viewer class. Deals with all the openGL rendering and the mouse/keyboard events.
class VIEWER_EXPORT Viewer : public CGAL::Three::Viewer_interface {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);
  ~Viewer();

  // overload several QGLViewer virtual functions
  //! Deprecated and does nothing.
  void draw();
  //!This step happens after draw(). It is here that the axis system is
  //!displayed.
  void drawVisualHints();
  //! Deprecated. Does the same as draw().
  void fastDraw();
  //! Initializes the OpenGL functions and sets the backGround color.
  void initializeGL();
  //! Deprecated and does nothing.
  void drawWithNames();
  /*! Uses the parameter pixel's coordinates to get the corresponding point
   * in the World frame. If this point is found, emits selectedPoint, selected,
   * and selectionRay signals.
   */
  void postSelection(const QPoint&);
  //! Sets the picking matrix to allow the picking.
  void beginSelection(const QPoint &point);
  //! Sets the pick matrix to Identity once the picking is done.
  void endSelection(const QPoint &point);
  //! Sets the scene for the viewer.
  void setScene(CGAL::Three::Scene_draw_interface* scene);
  //! @returns the antialiasing state.
  bool antiAliasing() const;
  //! @returns the fastDrawing state.
  bool inFastDrawing() const;
  //! Implementation of `Viewer_interface::inDrawWithNames()`
  bool inDrawWithNames() const;
  //! Implementation of `Viewer_interface::attrib_buffers()`
  void attrib_buffers(int program_name) const;
  //! Implementation of `Viewer_interface::getShaderProgram()`
  QOpenGLShaderProgram* getShaderProgram(int name) const;

public Q_SLOTS:
  //! Sets the antialiasing to true or false.
  void setAntiAliasing(bool b);
  //! If b is true, facets will be ligted from both internal and external sides.
  //! If b is false, only the side that is exposed to the light source will be lighted.
  void setTwoSides(bool b);
  //! If b is true, some items are displayed in a simplified version when moving the camera.
  //! If b is false, items display is never altered, even when moving.
  void setFastDrawing(bool b);
  //! Make the camera turn around.
  void turnCameraBy180Degres();
  //! @returns a QString containing the position and orientation of the camera.
  QString dumpCameraCoordinates();
  //!Moves the camera to the new coordinates (position and orientation) through an animation.
  bool moveCameraToCoordinates(QString, 
                               float animation_duration = 0.5f);

protected:
  //! Holds useful data to draw the axis system
  struct AxisData
  {
      std::vector<float> *vertices;
      std::vector<float> *normals;
      std::vector<float> *colors;
  };

  //! The buffers used to draw the axis system
  QOpenGLBuffer buffers[3];
  //! The VAO used to draw the axis system
  QOpenGLVertexArrayObject vao[1];
  //! The rendering program used to draw the axis system
  QOpenGLShaderProgram rendering_program;
  //! Holds the vertices data for the axis system
  std::vector<float> v_Axis;
  //! Holds the normals data for the axis system
  std::vector<float> n_Axis;
  //! Holds the color data for the axis system
  std::vector<float> c_Axis;
  //! Decides if the axis system must be drawn or not
  bool axis_are_displayed;
  //!Defines the behaviour for the mouse press events
  void mousePressEvent(QMouseEvent*);
  void wheelEvent(QWheelEvent *);
  //!Defines the behaviour for the key press events
  void keyPressEvent(QKeyEvent*);
  //!Deal with context menu events
  void contextMenuEvent(QContextMenuEvent*);
  /*!
   * \brief makeArrow creates an arrow and stores it in a struct of vectors.
   * \param R the radius of the arrow.
   * \param prec the precision of the quadric. The lower this value is, the higher precision you get.
   * It can be any int between 1 and 360.
   * \param from the starting point of the arrow.
   * \param to the destination point of the arrow (the pointed extremity).
   * \param color the RGB color of the arrow.
   * \param data the struct of std::vector that will contain the results.
   */

  void makeArrow(double R, int prec, qglviewer::Vec from, qglviewer::Vec to, qglviewer::Vec color, AxisData &data);
  void resizeGL(int w, int h);


protected:
  Viewer_impl* d;
  double prev_radius;
}; // end class Viewer

//!This class holds the properties of each line of text to be rendered.
class TextItem{
public :
    TextItem(float p_x, float p_y, float p_z, QString p_text, double id, QFont font = QFont(), QColor p_color = Qt::black)
        :x(p_x), y(p_y), z(p_z), m_text(p_text), m_id(id), m_font(font), m_color(p_color)
    {
       QFontMetrics fm(m_font);
       _width = fm.width(m_text);
       _height = fm.height();
    }
    QString text()const {return m_text;}
    //!Returns the position of the center of the text, in world coordinates.
    QVector3D *position(){return new QVector3D(x,y,z);}
    float width(){return _width;}
    float height(){return _height;}
    QFont font(){return m_font;}
    QColor color() {return m_color;}
    double id(){return m_id;}
private:
    float x;
    float y;
    float z;
    float _width;
    float _height;
    QString m_text;
    double m_id;
    QFont m_font;
    QColor m_color;
};//end class TextItem

//!This class draws all the textItems.
class TextRenderer{
public:
    TextRenderer()
    {
        m_lastId = 0.0;
    }
    /*!
      * Projects each textItem from the world coordinates to the Screen coordinates
      * and draws it.
     */
    void draw(CGAL::Three::Viewer_interface* viewer);
    void draw(CGAL::Three::Viewer_interface* viewer, TextItem* item);
    void addText(TextItem* ti);
    void addText(float p_x, float p_y, float p_z, QString p_text, double p_id, QFont font = QFont(), QColor p_color = Qt::black);
    void removeText(double id);
    double lastId()const{return m_lastId;}
    QMap<double,TextItem*> items() const{return textItems;}
    void printFacetId(QPoint pt, CGAL::Three::Viewer_interface*);
private:
    QList<double> displayList;
    QMap<double,TextItem*> textItems;
    QOpenGLFramebufferObject *fbo;
    double m_lastId;

};//end class TextRenderer
#endif // VIEWER_H
