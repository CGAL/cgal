/*! class Qt_layer is the main layer in the program.
 *  all the tab widget are attached to it.
 */
#include <CGAL/IO/Qt_widget_layer.h>

#include "cgal_types.h"

class QTabWidget;

class Qt_layer : public CGAL::Qt_widget_layer
{
public:
    Qt_layer( QTabWidget * );
	void draw();
 	
private:
	QTabWidget *myBar;
}; 
