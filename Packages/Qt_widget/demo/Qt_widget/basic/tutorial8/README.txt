8th tutorial
-----------------

In this tutorial you learn how to use a button to attach and
detach a Qt_widget_tool. As in the previous tutorial we use the
generic tool Qt_widget_get_point and a toolbar button that will
control the process of attaching and detaching.

To use a toolbar, here is what you have to do in the constructor of My_Window:

    QToolBar  *tools_toolbar;
    tools_toolbar = new QToolBar("Tools", this, QMainWindow::Top, TRUE, "Tools");
    addToolBar(tools_toolbar, Top, FALSE);

To add a button in the toolbar you have to:
-declare the button in My_Window:

	QToolButton *get_point_but;	//the toolbar button

-add the button in the toolbar:

    get_point_but =  new QToolButton(QPixmap( (const char**)point_xpm ),
				  "Point Tool", 
				  0, 
				  this, 
				  SLOT(pointtool()), 
				  tools_toolbar, 
				  "Point Tool");

The constructor of QToolButton needs a pointer to QMainWindow, that is
"this", and a pointer to a toolbar, that is "tools_toolbar".

To make the button a toggle button:

    get_point_but->setToggleButton(TRUE);

Also in the constructor of QToolButton is needed to declare a function
that will be called every time the button is pressed. That is
pointtool(). The first parameter is the icon found in <CGAL/IO/pixmaps/point.xpm>.


