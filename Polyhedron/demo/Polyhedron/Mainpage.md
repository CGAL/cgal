Main Page                        {#mainpage}
============

These pages are not documenting the whole Polyhedron demo but only the API that can be useful to create and add a new plugin.

Understanding the Polyhedron demo
============

There are several levels in this demo. 

- The MainWindow, which contains the UI elements. 

- Among these elements is the Viewer, which is the drawable surface that handles all the drawing and all the keyboard and mouse events. 

- The Viewer has a reference to the Scene, which contains the Scene_item list, which is a list of the drawn elements. 

A plugin usually defines an object that inherits from Scene_item or uses some of them to demonstrate a CGAL feature, so it might have to deal with the above elements.

Creating a simple Plugin
============
A basic plugin will inherit from Polyhedron_demo_plugin_interface. It can also inherits from the Polyhedron_demo_plugin_helper instead, for a more detailed model of plugin.
Its name must be of the form Polyhedron_demo_xxxx_yyyy_plugin. The next steps will assume the plugin's name is Polyhedron_demo_example_plugin
In the CMakeList.txt file, in the section Plugins, add the following lines :
  polyhedron_demo_plugin(example_plugin Polyhedron_demo_example_plugin)
  target_link_libraries(example_plugin scene_polyhedron_item)
It can get a reference to the Scene through the Scene_interface type, and to the MainWindow through the basic QMainWindow type. 
If the plugin implements a new Scene_item, please notice that a Scene_itam have a number of functions that will need a reference to the Viewer through the Viewer_interface type.

Example : 
============