
#ifndef CGAL_QT_POINT_TRANSLATION_WITH_INFO
#define CGAL_QT_POINT_TRANSLATION_WITH_INFO

#include <PointTranslation.h>

namespace CGAL {
namespace Qt {

template <typename DT, typename Translation>
class PointTranslationWithInfo : public PointTranslation<DT, Translation>
{
  typedef PointTranslation<DT, Translation> Base;
public:
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Point Point;
  
  typedef typename DT::Vertex::Info Info;

  PointTranslationWithInfo(const Translation& translation_, QGraphicsScene& scene_, DT  * dt_, QObject* parent);
  
  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
  
protected:
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  virtual Vertex_handle insert_point(const Point& p);
  
  Info compute_info() const;
  void add_info(Vertex_handle v) const;
  
  void show_info(const QPointF& pos, Vertex_handle v) const;
  
private:
  Info _info;
};


template <typename DT, typename Translation>
PointTranslationWithInfo<DT, Translation>::PointTranslationWithInfo(const Translation& translation_,
                                                    QGraphicsScene& scene_, 
                                                    DT * dt_,
                                                    QObject* parent)
  :  Base(translation_, scene_, dt_, parent)
{
}

  
template <typename DT, typename Translation>
typename PointTranslationWithInfo<DT, Translation>::Info 
PointTranslationWithInfo<DT, Translation>::compute_info() const
{
  return this->chosen_vertex->info() + info();
}

  
template <typename DT, typename Translation>
void 
PointTranslationWithInfo<DT, Translation>::add_info(Vertex_handle v) const
{
  Info new_info = compute_info();
  v->info().setString(new_info.toString());
}
  
template <typename DT, typename Translation>
void 
PointTranslationWithInfo<DT, Translation>::show_info(const QPointF& pos, Vertex_handle v) const
{
  if(v == Vertex_handle()) {
    return;
  }
  
  const Info& vertex_info = this->chosen_vertex->info();
  
  QToolTip::showText(pos.toPoint(), QString::fromStdWString(vertex_info.toString()));
}
  
  
  
template <typename DT, typename Translation>
typename PointTranslationWithInfo<DT, Translation>::Vertex_handle 
PointTranslationWithInfo<DT, Translation>::insert_point(const Point& p)
{
  Vertex_handle v = Base::insert_point(p);
  
  add_info(v);
  
  return v;
}
  

template <typename DT, typename Translation>
void 
PointTranslationWithInfo<DT, Translation>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  
  Base::mouseMoveEvent(event);
  
  show_info(event->scenePos(), this->chosen_vertex);
}
  
  
template <typename DT, typename Translation>
bool 
PointTranslationWithInfo<DT, Translation>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return false; // do not eat move event!
  }
  return Base::eventFilter(obj, event);
}
  

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POINT_TRANSLATION_WITH_INFO

