#ifndef ARRANGEMENT_DEMO_GRAPHICS_VIEW_H
#define ARRANGEMENT_DEMO_GRAPHICS_VIEW_H

#include <QGraphicsView>
#include <QColor>

class ArrangementDemoGraphicsView : public QGraphicsView
{
public:
  ArrangementDemoGraphicsView( QWidget* parent = 0 );

  void setShowGrid( bool b );
  bool getShowGrid( ) const;
  void setGridSize( int size );
  int getGridSize( ) const;
  void setGridColor( QColor color );
  QColor getGridColor( ) const;
  void setBackgroundColor( QColor color );
  QColor getBackgroundColor( ) const;


protected:
  void drawForeground( QPainter* painter, const QRectF& rect );
  QRectF getViewportRect( ) const;

  bool showGrid;
  int gridSize;
  QColor gridColor;
  QColor backgroundColor;
};

#endif // ARRANGEMENT_DEMO_GRAPHICS_VIEW_H
