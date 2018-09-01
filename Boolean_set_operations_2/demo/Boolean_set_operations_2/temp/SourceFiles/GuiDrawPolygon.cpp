#include "GuiDrawPolygon.h"
 
GuiDrawPolygon::GuiDrawPolygon(QObject *parent) : QGraphicsScene(parent)// , QDialog(parent)
{
  //setMouseTracking(true);
}
 
GuiDrawPolygon::~GuiDrawPolygon()
{

}
 
void GuiDrawPolygon::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    cout<<"click detected"<<endl;
    if(event->buttons() & Qt::LeftButton)
    {
        GUI_Polygon.push_back(event->scenePos());
        addLine(previousPoint.x(),previousPoint.y(),event->scenePos().x(),
            event->scenePos().y(),QPen(Qt::green,2,Qt::SolidLine,Qt::RoundCap));
    }
    else
    {
        addPolygon(GUI_Polygon,QPen(Qt::yellow,2,Qt::SolidLine,Qt::RoundCap),QBrush(Qt::red));
    }
    CGAL_Polygon.push_back(Point(event->scenePos().x(),event->scenePos().y()));
    // Save the coordinates of the point of pressing
    previousPoint=event->scenePos();
}
/* 
void GuiDrawPolygon::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
    //cout<<"mouse move detected"<<endl;
    //addLine(previousPoint.x(),previousPoint.y(),event->scenePos().x(),
      //      event->scenePos().y(),QPen(Qt::green,5,Qt::SolidLine,Qt::RoundCap));
    //QTimer::singleShot(5000, this, SLOT(quit()));
    //addLine(previousPoint.x(),previousPoint.y(),event->scenePos().x(),
      //      event->scenePos().y(),QPen(Qt::white,5,Qt::SolidLine,Qt::RoundCap));
    //~addLine();
    //previousPoint=event->scenePos();
}
*/
Polygon GuiDrawPolygon::getPolygon()
{
    return CGAL_Polygon;
}
