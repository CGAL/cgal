#include <QDebug>
#include <vector>

#include <QWidget>
#include <QObject>
#include <QPainter>

#include <QMouseEvent>
#include <QEvent>

#include <QPoint>
#include <QMatrix3x3>
#include <QVector2D>
#include <QVector3D>

#include <cmath>

struct State{
 bool is_left_pressed;
 bool is_right_pressed;
};

class UVProjector:public QWidget
{
public:
  UVProjector(QWidget* parent = 0, Qt::WindowFlags flags = Qt::WindowType(0))
    :QWidget(parent,flags)
  {
    setMouseTracking(true);
    state = {false,false};
    translation = QVector3D(0,0,1);
    rotation = 0;
    prev_pos = QPoint(0,0);
    rot_center = QPoint(width()/2, height()/2);

  }
  void setPoints(const std::vector<QPointF>& p){points = p;}
protected:
  void paintEvent(QPaintEvent*)
  {
    QPainter painter(this);
    QFont font;
    font.setPointSize(10);
    painter.setFont(font);
    painter.setBrush(QBrush(Qt::white));
    painter.setPen(Qt::white);
    painter.drawRect(QRect(QPoint(0,0), QPoint(this->width(),this->height())));

    painter.setPen(QPen(Qt::black));
    Q_FOREACH(QPointF p, points)
    {
     /*Translation(-w/2, -h/2) to recenter the scene, then
      * Scaling then Rotation and finaly the Translation
      * + Translation(w/2+h/2) to put it back.    */
     //scaled values
     qreal sx(translation.z()* (p.x()-width() /2.0)), sy(translation.z()* (p.y()-height()/2.0)) ;

     painter.drawPoint(
        translation.x() + width() /2.0 + cos(rotation)*sx + sin(rotation)*sy,
        translation.y() + height()/2.0 + -sin(rotation)*sx + cos(rotation)*sy
        );
    }
    painter.end();
  }
  void mouseMoveEvent(QMouseEvent* me)
  {
    if(state.is_left_pressed)
    {
     QVector2D prev_dir(prev_pos.x() - width()/2,
                   prev_pos.y() - height()/2);
     QVector2D dir(me->pos().x() - width()/2,
                   me->pos().y() - height()/2);
     prev_dir.normalize();
     dir.normalize();
     qreal a = acos(QVector2D::dotProduct(dir,prev_dir)/(dir.length()*prev_dir.length()));
     qreal det = dir.x()*prev_dir.y()-prev_dir.x()*dir.y();
     if(det > 0)
      rotation += a;
     else
      rotation -= a;
     update();
    }
    else if(state.is_right_pressed)
    {
      QVector2D dir(me->pos().x() - prev_pos.x(),
                    me->pos().y() - prev_pos.y());

      translation[0] += dir.x();
      translation[1] += dir.y();
      update();
    }
    prev_pos = me->pos();
  }

  void mousePressEvent(QMouseEvent* me)
  {
    if(me->button() == Qt::LeftButton)
      state.is_left_pressed = true;
    else if (me->button() == Qt::RightButton)
      state.is_right_pressed = true;
  }
  void mouseReleaseEvent(QMouseEvent* me)
  {
    if(me->button() == Qt::LeftButton)
      state.is_left_pressed = false;
    else if (me->button() == Qt::RightButton)
      state.is_right_pressed = false;
  }
  void wheelEvent(QWheelEvent *event)
  {
   if(event->angleDelta().y() >0)
     translation[2] *= 1.2;
   else
    translation[2] /= 1.2;
    update();
  }

  void mouseDoubleClickEvent(QMouseEvent * me)
  {
   if(state.is_right_pressed)
   {
     rot_center = me->pos();
   }
   else{
     translation = QVector3D(0,0,1);
     rotation = 0;
     rot_center = QPoint(width()/2, height()/2);
   }
    update();
  }
private:
  QPoint prev_pos;
  QPoint rot_center;
  std::vector<QPointF> points;
  qreal rotation;
  QVector3D translation;
  State state;

};
