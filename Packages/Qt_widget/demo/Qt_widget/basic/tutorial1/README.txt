First tutorial
-----------------

In this tutorial you can see how can you use Qt_widget like a stream,
for outputing Cgal objects.  Of course I recomend to read the tutorial
from Trolltech, that is the original Qt tutorial, but I think that you
can pass this tutorials without having strong skills of Qt
programming. Anyway, the code that belongs to Qt it is explained in
this tutorials.

The following is a typical template of how to create a window using Qt
and Qt_widget.

#include <CGAL/IO/Qt_widget.h>
#include <qapplication.h>

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    CGAL::Qt_widget * W = new CGAL::Qt_widget();
    app.setMainWidget( W );
    W.resize(600, 600);
    W.set_window(0, 600, 0, 600);
    W.show();

    return app.exec();
}

You'll allways need to include the header:

#include <qapplication.h>

The entry point for a typical Qt application is the function main. In
this function you should define an application object of Qt:

    QApplication app( argc, argv );

You will run the Qt application with the line:

    return app.exec();

To use Qt_widget you need an instance and tell the application to use
that instance:

    CGAL::Qt_widget *W = new CGAL::Qt_widget();
    app.setMainWidget( W );

To resize and set the scales of the window you'll use:

    W->resize(600, 600);
    W->set_window(0, 600, 0, 600);

At the end you need to show the window when the initialization have been done:

    W->show();

All the drawing code should be put betwen Qt_Widget's lock() and
unlock() functions. See the manual reference pages of Qt_widget. Doing
like this, the window will be updated only once, when Qt_widget find
the last unlock(). This way you can avoid the window flickering.

As you'll notice, this tutorial has some limitations. If you try to
resize the window you'll see that what you have been painted will
disappear. This is not a very pleasant thing but you'll see in the
next tutorial how you can solve this problem.

Applications following this approach are only usefull when you quickly
want to see how the output of a computation looks like.
