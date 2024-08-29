# GLFW
- Auto detection of keyboard layout
- Set pivot point position (Shift+Right) 

# Qt 
- Fix drop of performance when normals are displayed
- Remove CGAL_USE_BASIC_VIEWER target (Ask Laurent Rineau)

# Qt & GLFW viewer 
- Mesh selection 
- Screenshot with Wayland (?) 
- Add compatibility shaders in `Basic_shaders.h` for : 
  - `VERTEX_SOURCE_SHAPE`   
  - `GEOMETRY_SOURCE_SPHERE`
  - `GEOMETRY_SOURCE_CYLINDER`
  - `VERTEX_SOURCE_CLIPPING_PLANE`
  - `FRAGMENT_SOURCE_CLIPPING_PLANE`
  - `VERTEX_SOURCE_LINE`
  - `FRAGMENT_SOURCE_LINE`
  - `GEOMETRY_SOURCE_LINE`
  - `GEOMETRY_SOURCE_ARROW`
  - `VERTEX_SOURCE_NORMAL`
  - `GEOMETRY_SOURCE_NORMAL`
  - `VERTEX_SOURCE_TRIANGLE`
  - `GEOMETRY_SOURCE_TRIANGLE`
  - `VERTEX_SOURCE_LINE_WIDTH`
  - `GEOMETRY_SOURCE_LINE_WIDTH`

# Graphics scene
- Add a default size component for points, segments (, rays and lines ?) 
- Add methods to set size of points, segments (, rays and lines ?) 
