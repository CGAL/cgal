Third tutorial
-----------------

This tutorial guides you through the first contact with Views.

This tutorial is doing exactly what the previous tutorial does, but
instead of putting the code for drawing in the redraw() function, it use a
view.

In this example the view is:

class My_View : public CGAL::Qt_widget_view{
  void draw(CGAL::Qt_widget& win){
    win << dt;
  }
};

As you see you have to provide a draw() function in your view, in
order to create output on the screen. The Qt_widget will call this
draw() function for every attached and active view. A view is active
in the moment of attaching by default.

In the code of the constructor of My_Window the view is attached to
the widget:

	attach(&v);

Because the view is attached, the triangulation will appear on the
screen every time you call redraw(), as the widget redraws all the
views. This means that in this tutorial every time you press the mouse
button, the triangulation will be redrawn.

Try to detach the view to see what happens. Or try to deactivate
the view: The triangulation will not be shown anymore.
