7th tutorial
-----------------

This tutorial is very similar to the previous one. Insert new points
in a Delaunay triangulation when you click the mouse on the widget,
using a tool, and also let you use the standard toolbar.

The difference is that this example is using a generic tool, developed
by CGAL. This tool creates new CGAL points every time you click the
mouse. The coordinates of the point are the coordinates of the real
world. This means that the mouse coordinates are transformed using the
current scales to the real world coordinates.

The generic tools are documented in the manual. They are templetized
by a kernel of CGAL. In this tutorial the generic tool get_point it is
templetized by Cartesian<double>.

First comes the include statement:

	#include <CGAL/IO/Qt_widget_get_point.h>

In the class My_Window, it is declared as a private member:

	CGAL::Qt_widget_get_point<Rep> get_point;

In the constructor of My_Window also we attached this generic tool:
	    widget->attach(&get_point);

The connect is still there in the constructor, the good news is that
no matter how many tools you use, you will have to connect only once
the SIGNAL with your SLOT.
