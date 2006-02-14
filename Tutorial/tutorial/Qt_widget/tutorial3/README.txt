Third tutorial
-----------------

This tutorial guides you through the first contact with Layers.

This tutorial is doing exactly what the previous tutorial does, but
instead of putting the code for drawing in the redraw() function, it use a
layer.

In this example the layer is:

class My_Layer : public CGAL::Qt_widget_layer{
  void draw(CGAL::Qt_widget& win){
    win << dt;
  }
};

As you see you have to provide a draw() function in your layer, in
order to create output on the screen. The Qt_widget will call this
draw() function for every attached and active layer. A layer is active
in the moment of attaching by default.

In the code of the constructor of My_Window the layer is attached to
the widget:

	attach(&v);

Because the layer is attached, the triangulation will appear on the
screen every time you call redraw(), as the widget redraws all the
layers. This means that in this tutorial every time you press the mouse
button, the triangulation will be redrawn.

Try to detach the layer to see what happens. Or try to deactivate
the layer: The triangulation will not be shown anymore.
