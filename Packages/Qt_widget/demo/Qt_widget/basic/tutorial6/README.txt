6th tutorial
-----------------

The sixth tutorial uses for the first time a tool. It is declared the class My_Tool derived from Qt_widget_tool, that is used to create a Cgal point every time you click on the widget.

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

This class has a member function mousePressEvent(QMouseEvent *e) that is called by Qt_widget if the tool is attached.

In My_Window class it is created an instance of My_Tool:
	My_Tool	t;
In the constructor of My_Window we attach My_Tool:
	win.attach(t);
To receive the object that the tool creates, in the constructor of My_Window, it is the following line:
    connect(&win, SIGNAL(new_cgal_object(CGAL::Object)), 
	this, SLOT(get_object(CGAL::Object)));

In My_Window class it is declared the private slot get_object(CGAL::Object obj). This is what you have to do to receive Cgal objects from a tool. The tool calls the new_object() member function from Qt_widget that emit a signal new_cgal_object(CGAL::Object). To receive this signal, it is declared in the constructor of My_Window a connect, it is the SIGNAL-SLOT mechanism of Qt:
	connect(&win, SIGNAL(new_cgal_object(CGAL::Object)), 
		this, SLOT(get_object(CGAL::Object)));
This connect tells the application that every time Qt_widget emits the signal new_cgal_object, the function get_object is called, having the parameter an CGAL::Object. It could be anything, and to verify what object have been received you have to use CGAL::assign. If you are already familiar with CGAL this will be a peace of cake to you.

You can insert new points in the triangulation by clicking on the widget. You can notice that after the first use of the standard toolbar, you cannot insert more points in the triangulation. It's because some tools from standard toolbar detach the current tool before attaching. Due to the conflicts that may appear, it is preferable this way of doing it.

It is written in documantation that when a new tool is attached to the Qt_widget, the last tool attached, is detached. Never two tools will be attached in the same time on the widget.

To attach again the tool created in this tutorial you have to create some button that let you attach it or detach it at run time.