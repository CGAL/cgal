#include <CGAL/IO/Qt_widget_layer.h>

#include "cgal_types1.h"

class QTabWidget;

/////////////////////////////////////////////////////////////////////////////////////////////
class Qt_layer : public CGAL::Qt_widget_layer
{
public:

    Qt_layer( QTabWidget * );
	void draw();
 	
private:
	QTabWidget *myBar;
}; //end class 
