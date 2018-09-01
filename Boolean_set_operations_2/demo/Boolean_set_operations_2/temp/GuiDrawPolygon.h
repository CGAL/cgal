#ifndef CGAL_GUIDRAWPOLYGON_H
#define CGAL_GUIDRAWPOLYGON_H

#include <QDialog>
#include <QPolygon>
#include <QPainterPath>
#include <QPainter>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QDebug>//may be needed *******Remove at time of final submission******
#include "Typedefs.h"
#include <iostream>
#include <QTimer>
using namespace std;

class GuiDrawPolygon : public QGraphicsScene// ,public QDialog
{
 
    Q_OBJECT 

public:
    explicit GuiDrawPolygon(QObject *parent = 0);
    ~GuiDrawPolygon();
    Polygon getPolygon();
 
private:
	QPointF previousPoint;// The coordinates of the previous point
	Polygon CGAL_Polygon;
	QPolygonF GUI_Polygon;
    //QGraphicsLineItem* itemToDraw;
 
private:
    void mousePressEvent(QGraphicsSceneMouseEvent * event);
    //void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
 
};



#endif//CGAL_GUIDRAWPOLYGON_H
