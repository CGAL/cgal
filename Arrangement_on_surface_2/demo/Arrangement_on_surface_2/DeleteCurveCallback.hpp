#ifndef DELETE_CURVE_CALLBACK_HPP
#define DELETE_CURVE_CALLBACK_HPP
#include "Callback.hpp"
#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

/**
Handles deletion of arrangement curves selected from the scene.

The template parameter is the Arrangement type.

TODO: Create a base callback class.
*/
template < class TArr >
class DeleteCurveCallback : public CGAL::Qt::Callback
{

public:
    DeleteCurveCallback( TArr* arr_, QObject* parent_ ):
        CGAL::Qt::Callback( parent_ ),
        arr( arr_ )
    { }

    void setScene( QGraphicsScene* scene_ )
    {
        this->scene = scene_;
    }

    QGraphicsScene* getScene( ) const
    {
        return this->scene;
    }

protected:
    void mousePressEvent( QGraphicsSceneMouseEvent *event )
    {
        std::cout << "mouse press event stub" << std::endl;
        // TODO: delete the curve-to-delete here
    }
    
    void mouseMoveEvent( QGraphicsSceneMouseEvent *event )
    {
        std::cout << "mouse move event stub" << std::endl;
        // TODO: highlight the curve-to-delete here
    }

    QGraphicsScene* scene;
    TArr* arr;
};
#endif // DELETE_CURVE_CALLBACK_HPP
