/*
#include <QApplication>
#include <QtGui/QGraphicsView>
#include <QtGui/QGraphicsScene>
#include <QKeyEvent>
#include <iostream>

struct  Up : public QObject {

  Q_OBJECT

public:
  Up(QGraphicsView* v_)
    : v(v_)
  {}
  
  bool eventFilter(QObject *obj, QEvent *event)
  {
    if (event->type() == QEvent::KeyPress) {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      if (keyEvent->key() == Qt::Key_Up){
	std::cout << "up" << std::endl;
	v->translate(0, -10);
	return true;
      }  
    }
    return QObject::eventFilter(obj, event);
  }

  QGraphicsView* v;

};

int main(int argc, char **argv)
{
  QApplication app(argc, argv);


  QGraphicsScene scene;
  scene.setSceneRect(0,0, 100, 100);
  scene.addRect(20, 20, 30, 5);

  QGraphicsView* view = new QGraphicsView(&scene); 
  view->installEventFilter(new Up(view));
  view->setTransformationAnchor(QGraphicsView::NoAnchor);

  view->show();
  return app.exec();
}
*/

#include <QtGui>

struct  Up : public QObject {

    Q_OBJECT

public:
    Up(QGraphicsView* v_)
        : v(v_)
    {}

    bool eventFilter(QObject *obj, QEvent *event)
    {
        if (event->type() == QEvent::KeyPress) {
            QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
            switch (keyEvent->key()) {
            case Qt::Key_Up:
                v->translate(0, -10);
                break;
            case Qt::Key_Down:
                v->translate(0, 10);
                break;
            case Qt::Key_Left:
                v->translate(-10, 0);
                break;
            case Qt::Key_Right:
                v->translate(10, 0);
                break;
            case Qt::Key_PageUp:
                v->rotate(-6);
                break;
            case Qt::Key_PageDown:
                v->rotate(6);
                break;
            default:
                return false;
            }
            return true;
        }
        return false;
    }

    QGraphicsView* v;
};

int main(int argc, char **argv)
{
    QApplication app(argc, argv);


    QGraphicsScene scene;
    scene.setSceneRect(0,0, 100, 100);
    //scene.setBackgroundBrush(Qt::red);
    scene.addRect(20, 20, 30, 5);
    scene.addRect(80, 80, 135, 5);

    QGraphicsView* view = new QGraphicsView(&scene);
    view->installEventFilter(new Up(view));
    view->setTransformationAnchor(QGraphicsView::NoAnchor);

    view->show();
    return app.exec();
}

#include "min.moc"
