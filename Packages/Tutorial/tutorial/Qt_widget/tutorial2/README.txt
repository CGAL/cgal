Second tutorial
-----------------

In this tutorial you can see a different method to draw on the window,
and you'll pass over the limitations of the previous example.

This tutorials insert a point in a Delaunay triangulation every time
you press the mouse button over the widget, and calls redraw(). This
means that you'll see the results of your insertions immediately. The
advantage of this approach is that you can resize the window, the
painting will not disappear.

As you see the main entry point is the same as in the previous
tutorial, but instead of using an instance of Qt_widget class, you'll
use an instance of a class derived from Qt_widget.
In this sense, we create a new class My_Window as a child for
Qt_widget. The resize function has been moved into the constructor of
My_Window.

We define a triangulation as a global variable:
typedef CGAL::Cartesian<double>		    Rep;
typedef CGAL::Point_2<Rep>		    Point;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

The private member redraw() is used to output the triangulation on the screen. 
private:
  void redraw()
  {
    Qt_widget::redraw();
    *this << dt;
  }

There is another line in the redraw():
	Qt_widget::redraw();

This line is meant to call the redraw() from the Qt_widget base class,
used to draw the attached layers as you'll see later. This function
also clears the screen. If you remove this line it's necessary to add
clear() instead.

Another private member catches the mouse press event from the Window
system (X11/Windows). The code put in this function will be executed when
you press the mouse on the widget.

private:
  void mousePressEvent(QMouseEvent *e)
  {
    Qt_widget::mousePressEvent(e);
    dt.insert(Point(x_real(e->x()), y_real(e->y())));
    redraw();
  }

As you see the code inserts a point in the triangulation, and calls
redraw(). You can comment redraw() to see what happens: You'll see
the results only when you'll resize the window.

Of course if you want to trigger the mouse press event in the base
class you have to put this in the code of the function:

    Qt_widget::mousePressEvent(e);


