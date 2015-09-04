#ifndef SCENE_ITEM_H
#define SCENE_ITEM_H
#include "Scene_item_config.h"
#include "Scene_interface.h"
#include <QString>
#include <QPixmap>
#include <QFont>
#include <QOpenGLBuffer>
#include <QOpenGLShader>
#include <QOpenGLVertexArrayObject>
#include <vector>
#include <QMap>
#define PROGRAM_WITH_LIGHT 0
#define PROGRAM_WITHOUT_LIGHT 1
#define PROGRAM_WITH_TEXTURE 2
#define PROGRAM_WITH_TEXTURED_EDGES 3
#define PROGRAM_INSTANCED 4
#define PROGRAM_INSTANCED_WIRE 5


namespace qglviewer {
  class ManipulatedFrame;
}

class QMenu;
class QKeyEvent;
class Viewer_interface;

//! This class represents an object in the OpenGL scene
class SCENE_ITEM_EXPORT Scene_item : public QObject {
  Q_OBJECT
  Q_PROPERTY(QColor color READ color WRITE setColor)
  Q_PROPERTY(QString name READ name WRITE setName)
  Q_PROPERTY(bool visible READ visible WRITE setVisible)
  Q_ENUMS(RenderingMode)
  Q_PROPERTY(RenderingMode renderingMode READ renderingMode WRITE setRenderingMode)
public:
  typedef Scene_interface::Bbox Bbox;
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;
  //! The default color of a scene_item.
  static const QColor defaultColor; // defined in Scene_item.cpp

  //!The default Constructor.
  /*!
   * Initializes the number of VBOs to 20 and VAOs to 10 and creates them.
   */
  Scene_item()
    : name_("unamed"),
      color_(defaultColor),
      visible_(true),
      are_buffers_filled(false),
      rendering_mode(FlatPlusEdges),
      defaultContextMenu(0),
      buffersSize(20),
      vaosSize(10)
  {

      nbVaos = 0;
      for(int i=0; i<vaosSize; i++)
      {
          addVaos(i);
          vaos[i]->create();
      }

      for(int i=0; i<buffersSize; i++)
      {
          QOpenGLBuffer n_buf;
          buffers.push_back(n_buf);
          buffers[i].create();
      }
  }
  //!The Constructor.
  /*!
   * Initializes the number of VBOs and VAOs and creates them.
   */
  Scene_item(int buffers_size, int vaos_size)
    : name_("unamed"),
      color_(defaultColor),
      visible_(true),
      are_buffers_filled(false),
      rendering_mode(FlatPlusEdges),
      defaultContextMenu(0),
      buffersSize(buffers_size),
      vaosSize(vaos_size)
  {
      nbVaos = 0;
      for(int i=0; i<vaosSize; i++)
      {
          addVaos(i);
          vaos[i]->create();
      }

      for(int i=0; i<buffersSize; i++)
      {
          QOpenGLBuffer n_buf;
          buffers.push_back(n_buf);
          buffers[i].create();
      }
  }
  //!The destructor. It is virtual as the item is virtual.
  virtual ~Scene_item();
  //! Creates a new item as a copy of the sender. Must be overloaded.
  virtual Scene_item* clone() const = 0;

  //! Indicates if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const = 0;
  //! Deprecated. Does nothing.
  virtual void draw() const {}
  /*! \brief The drawing function.
   * Draws the facets of the item in the viewer using OpenGL functions. The data
   * for the drawing is gathered in compute_elements(), and is sent
   * to buffers in initialize_buffers().
   * @see compute_elements()
   * @see initialize_buffers()
   */
  virtual void draw(Viewer_interface*) const  { draw(); }
  //! Deprecated. Does nothing.
  virtual void draw_edges() const { draw(); }
  /*! \brief The drawing function.
   * Draws the edges and lines of the item in the viewer using OpenGL functions. The data
   * for the drawing is gathered in compute_elements(), and is sent
   * to buffers in initialize_buffers().
   * @see compute_elements()
   * @see initialize_buffers()
   */
  virtual void draw_edges(Viewer_interface* viewer) const { draw(viewer); }
  //! Deprecated. Does nothing.
  virtual void draw_points() const { draw(); }
  /*! \brief The drawing function.
   * Draws the points of the item in the viewer using OpenGL functions. The data
   * for the drawing is gathered in compute_elements(), and is sent
   * to buffers in initialize_buffers().
   * @see compute_elements()
   * @see initialize_buffers()
   */
  virtual void draw_points(Viewer_interface*) const { draw_points(); }

  //! Specifies which data must be updated when selection has changed.
  //! Must be overloaded.
  virtual void selection_changed(bool);

  // Functions for displaying meta-data of the item
  //! @returns a QString containing meta-data about the item.
  //!//! Must be overloaded.
  virtual QString toolTip() const = 0;
  //! @returns a QPixmap containing graphical meta-data about the item.
  virtual QPixmap graphicalToolTip() const { return QPixmap(); }
  //! @returns a QFont containing the font used for the data of the item.
  virtual QFont font() const { return QFont(); }

  // Functions that help the Scene to compute its bbox
  //! If isFinite() returns false, the BBox is not computed.
  virtual bool isFinite() const { return true; }
  //! Specifies if the item is empty or null.
  virtual bool isEmpty() const { return true; }
  //!@returns the item's bounding box.
  virtual Bbox bbox() const { return Bbox(); }

  // Function about manipulation
  //! Decides if the item can have a ManipulatedFrame.
  virtual bool manipulatable() const { return false; }
  //!@returns the manipulatedFrame of the item.
  virtual ManipulatedFrame* manipulatedFrame() { return 0; }

  // Getters for the four basic properties
  //! @returns the current color of the item.
  virtual QColor color() const { return color_; }
  //! @returns the current name of the item.
  virtual QString name() const { return name_; }
  //! @returns the current visibility of the item.
  virtual bool visible() const { return visible_; }
  //! @returns the current rendering mode of the item.
  virtual RenderingMode renderingMode() const { return rendering_mode; }
  //! @returns the current rendering mode of the item as a human readable string.
  virtual QString renderingModeName() const;

  //! Context menu
  virtual QMenu* contextMenu();

  //!Handles key press events.
  virtual bool keyPressEvent(QKeyEvent*){return false;}
public Q_SLOTS:
  // Call that once you have finished changing something in the item
  // (either the properties or internal data)
  //!Specifies what to do when the item has changed. Typically calls
  //! compute_elements() and initialize_buffers(). Must be overloaded.
  virtual void changed();
  //!When changed() is not enough.
  virtual void contextual_changed(){}

  // Setters for the four basic properties
  virtual void setColor(QColor c) { color_ = c; changed(); }
  void setRbgColor(int r, int g, int b) { setColor(QColor(r, g, b)); }
  virtual void setName(QString n) { name_ = n; }
  virtual void setVisible(bool b) { visible_ = b; }
  virtual void setRenderingMode(RenderingMode m) { 
    if (supportsRenderingMode(m))
      rendering_mode = m; 
  }
  void setPointsMode() {
    setRenderingMode(Points);
  }

  void setWireframeMode() {
    setRenderingMode(Wireframe);
  }
  void setWireframe() {
    setRenderingMode(Wireframe);
  }

  void setFlat() {
    setRenderingMode(Flat);
  }
  void setFlatMode() {
    setRenderingMode(Flat);
  }

  void setFlatPlusEdgesMode() {
    setRenderingMode(FlatPlusEdges);
  }

  void setGouraudMode() {
    setRenderingMode(Gouraud);
  }

  void setPointsPlusNormalsMode(){
    setRenderingMode(PointsPlusNormals);
  }
  
  virtual void itemAboutToBeDestroyed(Scene_item*);

  virtual void select(double orig_x,
                      double orig_y,
                      double orig_z,
                      double dir_x,
                      double dir_y,
                      double dir_z);

Q_SIGNALS:
  void itemChanged();
  void aboutToBeDestroyed();

protected:
  // The four basic properties
  QString name_;
  QColor color_;
  bool visible_;
  bool is_selected;
  mutable bool are_buffers_filled;
  RenderingMode rendering_mode;
  QMenu* defaultContextMenu;

  RenderingMode prev_shading;
  RenderingMode cur_shading;

  int buffersSize;
  int vaosSize;
  mutable std::vector<QOpenGLBuffer> buffers;
  //not allowed to use vectors of VAO for some reason
  //mutable QOpenGLVertexArrayObject vaos[10];
  QMap<int,QOpenGLVertexArrayObject*> vaos;
  int nbVaos;
  void addVaos(int i)
  {
      QOpenGLVertexArrayObject* n_vao = new QOpenGLVertexArrayObject();
      vaos[i] = n_vao;
      nbVaos ++;
  }


  mutable QMap<int, QOpenGLShaderProgram*> shader_programs;
  QOpenGLShaderProgram* getShaderProgram(int , Viewer_interface *viewer = 0) const;

  int vertexLoc;
  int normalLoc;
  int colorLoc;

  virtual void initialize_buffers(){}
  virtual void compute_elements(){}
  virtual void attrib_buffers(Viewer_interface*, int program_name) const;



}; // end class Scene_item


#include <QMetaType>
Q_DECLARE_METATYPE(Scene_item*)

#endif // SCENE_ITEM_H
