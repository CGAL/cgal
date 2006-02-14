Fourth tutorial
-----------------

The fourth tutorial shows how to create a more complex application
using the Qt class QMainWindow. With this class you can create a
MDI(Multiple document interface) application.

For drawing, the tutorial use the same layer as for the previous one
but this time we use another class QMainWindow, as the main frame of
the application. Further down is described how the flow is:

The entry point is the same:

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using QMainWindow QT class");
    return app.exec();
}

This time you see that W is no longer an instance of Qt_widget but is
an instance of QMainWindow. QMainWindow is also a widget but provides
functionalities like you can use a toolbar, a status bar ... in other
words you can create a complex MDI application. To use an instance of
Qt_widget you have to say in the constructor:

	setCentralWidget(widget);

where widget is an instance of Qt_widget. As you see it is also declared
in My_window.

There is one more thing, at the constructor you see:

	widget = new My_Widget(this);

Also the constructor of My_widget is adapted:

	My_widget(QMainWindow* c) : CGAL::Qt_widget(c) {};

This lines tells the application that My_window is a parent for
My_widget. Try to comment this lines to see what happens. Two
windows will appear, one for My_window and one for My_widget.

In the constructor it is widget->attach(&v); This way, the layer is
attached by My_widget. The rest of the code does the same thing as the
previous tutorials: insert a new point in a Delaunay triangulation and
draw the triangulation every time you click on the window.
