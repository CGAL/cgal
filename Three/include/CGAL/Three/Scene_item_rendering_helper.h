#ifndef SCENE_ITEM_RENDERING_HELPER_H
#define SCENE_ITEM_RENDERING_HELPER_H

#include <CGAL/Three/Scene_item.h>
#include <QThread>


#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif

struct D;

namespace CGAL {
namespace Three{

//!
//! \brief The Scene_item_rendering_helper class is a convenience class for  constructing an item.
//! It is more elaborated than a `Scene_item` and facilitates the process of creating an item that can
//! be rendered.
//!
class DEMO_FRAMEWORK_EXPORT Scene_item_rendering_helper
    :public Scene_item
{
  Q_OBJECT
public:
  Scene_item_rendering_helper();
  ~Scene_item_rendering_helper();

  //!
  //! \brief newViewer adds Vaos for `viewer`.
  //! \param viewer the new viewer.
  //!
  void newViewer(Viewer_interface *viewer) Q_DECL_OVERRIDE;

  //! \brief removeViewer removes the Vaos for `viewer`.
  //! \param viewer the viewer to be removed.
  void removeViewer(Viewer_interface *viewer) Q_DECL_OVERRIDE;
  QMenu* contextMenu() Q_DECL_OVERRIDE;
  /*!
    * \brief processData calls `computeElements()` in a dedicated thread so the
    * application does not get stuck while the processing is performed.
    * Emits `dataProcessed()`.
    */
  void processData(Gl_data_names name)const;
  //!
  //! \brief setAlpha sets the integer value of the alpha channel of this item.
  //! Also updates the slider value.
  //! It must be between 0 and 255.
  //! \param alpha the integer value for the alpha channel.
  //!
  void setAlpha(int alpha) Q_DECL_OVERRIDE;
  //! \brief The item's bounding box.
  //!
  //! If the Bbox has never been computed, computes it and
  //! saves the result for further calls.
  //! @returns the item's bounding box.
  Scene_item::Bbox bbox()const Q_DECL_OVERRIDE;
  //!
  //! \brief getTriangleContainer returns the `id`th `Triangle_container`.
  //!
  CGAL::Three::Triangle_container* getTriangleContainer(std::size_t id) const;
  //!
  //! \brief getEdgeContainer returns the `id`th `Edge_container`.
  //!
  CGAL::Three::Edge_container* getEdgeContainer(std::size_t id)const;

  //!
  //! \brief setTriangleContainer sets the `id`th `Triangle_container` to `tc`.
  //!
  //! If `id` is bigger than the current size of the container vector, this vector is
  //! resized accordingly. This means that for optimisation reasons, containers should be created
  //! decreasingly.
  //!
  void setTriangleContainer(std::size_t id,
                            Triangle_container* tc);

  //!
  //! \brief setEdgeContainer sets the `id`th `Edge_container` to `tc`.
  //!
  //! If `id` is bigger than the current size of the container vector, this vector is
  //! resized accordingly. This means that for optimisation reasons, containers should be created
  //! decreasingly.
  //!
  void setEdgeContainer(std::size_t id,
                        Edge_container* tc);

  //!
  //! \brief setBuffersFilled specifies if the data should be re-computed.
  //!
  //! If called with `false`, the item rendering data will be re-computed at the next `draw()`.
  //! If called with `true`, the item rendering data is considered ready and will not be computed
  //! until `setBuffersFilled()` is called with `false` again.
  //!
  void setBuffersFilled(bool b) const;

  //!
  //! \brief getBuffersFilled returns `false` if the item rendering data needs to be re-computed.,
  //! `true` otherwise.
  //! \see `setBuffersFilled()`
  bool getBuffersFilled()const;

  //!
  //! \brief getBuffersInit returns true if the `Vao`s of `viewer` are ready
  //! for rendering.
  //!
  bool getBuffersInit(Viewer_interface *viewer)const;

  //!
  //! \brief setBuffersInit specifies if the `Vbo`s need to be initialized.
  //!
  //! If called with `false`, the item `Vbo`s will be refilled at the next `draw()`.
  //! If called with `true`, the item `Vbo`s are considered ready and will not be refilled
  //! until `setBuffersInit()` is called with `false` again.
  //!
  //! This function should be called in the drawing functions, when `getBuffersFilled()` is `true`.
  //!
  void setBuffersInit(Viewer_interface *viewer, bool val) const;
protected:

  /*!
   * \return a pointer to the slider initialized in initGL();
   */
  QSlider* alphaSlider();

  //!
  //! \return `true` if `initGL()` was called.
  //!
  bool isInit()const;
  //! \brief the item's bounding box's diagonal length.
  //!
  //! If the diagonal's length has never been computed, computes it and
  //! saves the result for further calls.
  //! @returns the item's bounding box's diagonal length.
  virtual double diagonalBbox() const;

  //!Returns the float alpha value of an item.
  //! This value is between 0.0f and 1.0f.
  float alpha() const Q_DECL_OVERRIDE;

  /*! Fills the `Vbo`s with data. Must be called after each call to #computeElements().
     * @see computeElements()
     */
  virtual void initializeBuffers(Viewer_interface*)const{}

  //!Creates the VAOs and VBOs for each existing viewer.
  virtual void initGL() const;
  //!
  //! Computes the items Bbox and stores the result. Must be overridden.
  //!
  virtual void compute_bbox() const= 0;
  //!
  //! \brief setBbox allows to set the Bbox in compute_bbox();
  //! \param b
  //!
  void setBbox(Bbox b) const;
private:
  friend struct D;
  mutable D* d;
};//end Scene_item_rendering_helper

}}


//!
//! \brief The WorkerThread class computes the data of this item in a separated
//! thread. It allows to keep the hand on the GUI and to manage several items at the same time.
//!
class WorkerThread : public QThread
{
  Q_OBJECT
  CGAL::Three::Scene_item* item;
  CGAL::Three::Scene_item::Gl_data_names name;
public:
  //!
  //! \brief The `WorkerThread` constructor.
  //! \param item the `Scene_item` with the data that needs computation.
  //! \param name specifies which type of data must be re-computed.
  //!
  WorkerThread(CGAL::Three::Scene_item* item, CGAL::Three::Scene_item::Gl_data_names name):
    item(item),
    name(name){}
private:
  void run() Q_DECL_OVERRIDE{
    item->computeElements(name);
  }
};

#endif // SCENE_ITEM_RENDERING_HELPER_H
