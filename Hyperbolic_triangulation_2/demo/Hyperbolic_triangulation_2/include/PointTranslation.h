
#ifndef CGAL_QT_POINT_TRANSLATION
#define CGAL_QT_POINT_TRANSLATION

#include <CGAL/Qt/GraphicsViewInput.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>

#include <PointGraphicsItem.h>

namespace CGAL {
namespace Qt {

template <typename DT, typename Translation>
class PointTranslation : public GraphicsViewInput
{
public:
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Point Point;

  PointTranslation(const Translation& translation_, QGraphicsScene& scene_, DT  * dt_, QObject* parent);
  
  void show();
  void hide();
  
protected:
  Point do_point_translation() const;
  virtual Vertex_handle insert_point(const Point& p);

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  Translation translation;
  Vertex_handle chosen_vertex;

private:
  void highlight_chosen_vertex();

  DT * dt;

  QGraphicsScene& scene;
  QGraphicsEllipseItem* circle;

  // highlight the chosen vertex
  PointGraphicsItem<DT>* highlight_point;
};


template <typename DT, typename Translation>
PointTranslation<DT, Translation>::PointTranslation(const Translation& translation_,
                                                    QGraphicsScene& scene_, 
                                                    DT * dt_,
                                                    QObject* parent)
  :  GraphicsViewInput(parent), translation(translation_), dt(dt_), scene(scene_)
{
  highlight_point = new PointGraphicsItem<DT>();
  highlight_point->hide();
  
  // scene is the holder of the added item
  scene.addItem(highlight_point);
}


template <typename DT, typename Translation>
void
PointTranslation<DT, Translation>::show()
{
  highlight_point->show();
}
  
  
template <typename DT, typename Translation>
void
PointTranslation<DT, Translation>::hide()
{
  highlight_point->hide();
} 
  
template <typename DT, typename Translation>
typename PointTranslation<DT, Translation>::Point 
PointTranslation<DT, Translation>::do_point_translation() const
{
  return translation.DoAction(chosen_vertex->point());
}

  
template <typename DT, typename Translation>
typename PointTranslation<DT, Translation>::Vertex_handle 
PointTranslation<DT, Translation>::insert_point(const Point& p)
{
  Vertex_handle v = dt->insert(p, chosen_vertex->face());
  
  emit(modelChanged());
  
  return v;
}
  

template <typename DT, typename Translation>
void 
PointTranslation<DT, Translation>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if(chosen_vertex == Vertex_handle() ||
     event->modifiers() != 0 ||
     event->button() != ::Qt::LeftButton) {
    return;
  }
  
  // translate chosen point
  Point translated_point = do_point_translation();
  // insert to the triangulation
  insert_point(translated_point);
}



template <typename DT, typename Translation>
void 
PointTranslation<DT, Translation>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  Point scene_pos = Point(event->scenePos().x(), event->scenePos().y());
  
  Face_handle start_face = Face_handle();
  if(chosen_vertex != Vertex_handle()) {
    start_face = chosen_vertex->face();
  }
  
  Vertex_handle new_chosen_vertex = dt->nearest_vertex(scene_pos, start_face);
  
  // highlight the vertex if only it's been changed
  if(new_chosen_vertex != chosen_vertex) {
    chosen_vertex = new_chosen_vertex;
    highlight_chosen_vertex();
  }
}


template <typename DT, typename Translation>
bool 
PointTranslation<DT, Translation>::eventFilter(QObject *obj, QEvent *event)
{  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return true;
  } else if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return false; // do not eat move event!
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
}
  
template <typename DT, typename Translation>
void 
PointTranslation<DT, Translation>::highlight_chosen_vertex()
{
  assert(chosen_vertex != Vertex_handle());
  
  highlight_point->setPoint(chosen_vertex->point());
  highlight_point->show();
  emit(modelChanged());
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POINT_TRANSLATION
