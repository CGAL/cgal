5th tutorial
-----------------

The fifth tutorial includes the standard toolbar in the
application. The application has the same functionality but as you
see you can now zoom in, zoom out and translate, using the tools from
the standard toolbar.

The layer is still there, doing nothing but drawing the triangulation.

The class My_Widget is like in the previous example an intance of
Qt_widget, that you need in the class My_Window.

The only difference between this example and the previous one is the
use of standard toolbar. It is declared as private in My_Window class:

	CGAL::Qt_widget_standard_toolbar *stoolbar;

To use it, in the constructor of My_Window, it is added:

	stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this);
	this->addToolBar(stoolbar->toolbar(), Top, FALSE);

In this tutorial you can play a little bit with the standard toolbar
but you will see probably something that is not quite pleasant. If you
select the hand, and try to use it, you'll see that when you click to
select the first point, also a new CGAL::Point_2 is inserted in the
triangulation. This is due to the fact that the mousePressEvent
implemented in My_Window it is called every time you click on the
widget, even if a tool from the standard toolbar is active.

You will see later a solution to this problem.
