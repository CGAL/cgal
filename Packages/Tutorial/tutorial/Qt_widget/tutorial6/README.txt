6th tutorial
-----------------

The sixth tutorial uses for the first time a tool. It is declared the
class My_Tool derived from Qt_widget_tool, that is used to create a
CGAL point every time you click on the widget.

class My_Tool : public CGAL::Qt_widget_tool{
public:
  My_Tool(){};
private:
  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton)
    {
      double
	x=static_cast<double>(widget->x_real(e->x())),
	y=static_cast<double>(widget->y_real(e->y()));
      widget->new_object(CGAL::make_object(Point(x, y)));
    }
  }
};

This class has a member function mousePressEvent(QMouseEvent *e) that
is called by Qt_widget if the tool is attached.

In My_Window class an instance of My_Tool is created :

	My_Tool	t;

In the constructor of My_Window we attach the tool:

	widget->attach(&t);

To receive the object that the tool creates, in the constructor of My_Window, 
we connect:

    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
	    this, SLOT(get_object(CGAL::Object)));

In My_Window class it is declared the private slot get_object(CGAL::Object obj). 
This is what you have to do to receive CGAL objects from a tool. The
tool calls the new_object() member function from Qt_widget, that emits
a signal new_cgal_object(CGAL::Object). To receive this signal, in the
constructor of My_Window a connect is declared:

	connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
		this, SLOT(get_object(CGAL::Object)));

This connect tells the application that every time Qt_widget emits the
signal new_cgal_object, the function get_object is called, having as
parameter a CGAL::Object. It could be anything, and to verify what
object have been received you have to use CGAL::assign. If you are
already familiar with CGAL this will be a piece of cake to you.

You can insert new points in the triangulation by clicking on the
widget. When you use the Standard toolbar, your already attached tool will be put on waiting state, till you finish working with it. For example if you click on hand tool, and you use it, then you click once again for detaching it, your attached tool will receive the focus being brought to life.

It is written in the documentation that when a new tool is attached to the
Qt_widget, the last tool attached, is detached. Two tools will
never be attached in the same time on the widget.

