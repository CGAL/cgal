#ifndef CGAL_QT_CALLBACK_HPP
#define CGAL_QT_CALLBACK_HPP
#include <QObject>

class QEvent;
class QKeyEvent;
class QGraphicsSceneMouseEvent;

namespace CGAL {
namespace Qt {

class Callback : public QObject
{
Q_OBJECT

public:
    Callback( QObject* parent );

signals:
    void modelChanged( );

protected:
    virtual bool eventFilter( QObject* object, QEvent* event );
    virtual void mousePressEvent( QGraphicsSceneMouseEvent* event );
    virtual void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
    virtual void mouseReleaseEvent( QGraphicsSceneMouseEvent* event );
    virtual void keyPressEvent( QKeyEvent* event );
};

} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_CALLBACK_HPP
