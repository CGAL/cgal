
User
====

The CMakeLists.txt of this application allows some customization. Let's review the key Cache variables :

- `LINK_WITH_TBB`: set this to `true` to allow the application to use multi-threading.
- `SINGLE_PLUGIN_*`: those variables represent the plugins of the application. Check one of them to set the corresponding plugin individually. 
- `CATEGORY_PLUGINS_IO`: set this variable to `true` to set all the IO plugins. 
- `CATEGORY_PLUGINS_all`: set this variable to `true` to set all the plugins. 
- `CATEGORY_PLUGINS_mesh_processing`: set this variable to `true` to set all the plugins that deal with polygon meshes. 
- `CATEGORY_PLUGINS_meshing`: set this variable to `true` to set all the plugins that deal with volumic and surfacic meshing.
- `CATEGORY_PLUGINS_point_set_processing`: set this variable to `true` to set all the plugins that deal with point sets.

Developper
==========

When you want to add a plugin to the demo, you will have to use the `polyhedron_demo_plugin()` macro. Its first argument is the name of the target, 
its second argument is the corresponding source file, its third argument is one of the `CATEGORY_PLUGINS_*` prefixes. Add the UI files after that.
example : 
    qt5_wrap_ui( exampleUI_FILES Example_widget.ui)
    polyhedron_demo_plugin(example_plugin Example_plugin.cpp CATEGORY-PLUGINS_example ${exampleUI_FILES})
This will add your plugin to the list of plugins that will be set when the `CATEGRY_PLUGINS_example` is `true`.
