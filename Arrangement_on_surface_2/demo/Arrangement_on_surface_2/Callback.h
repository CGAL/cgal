#ifndef CGAL_QT_CALLBACK_H
#define CGAL_QT_CALLBACK_H
#include <QObject>

class QRectF;
class QEvent;
class QKeyEvent;
class QGraphicsScene;
class QGraphicsSceneMouseEvent;

namespace CGAL {
namespace Qt {

class Callback : public QObject
{
Q_OBJECT

public:
    Callback( QObject* parent );
    virtual void setScene( QGraphicsScene* scene_ );
    virtual QGraphicsScene* getScene( ) const;
    virtual void reset( );

public slots:
    virtual void slotModelChanged( );

signals:
    void modelChanged( );

protected:
    virtual bool eventFilter( QObject* object, QEvent* event );
    virtual void mousePressEvent( QGraphicsSceneMouseEvent* event );
    virtual void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
    virtual void mouseReleaseEvent( QGraphicsSceneMouseEvent* event );
    virtual void keyPressEvent( QKeyEvent* event );

    /**
    Return the bounding box of the visible scene.
    */
    QRectF viewportRect( ) const;

    QGraphicsScene* scene;
};

} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_CALLBACK_H
