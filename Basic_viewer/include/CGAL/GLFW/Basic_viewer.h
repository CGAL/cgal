// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Théo Benard <benard320@gmail.com>

//                 Théo Grillon <theogrillon6f9@gmail.com>

#ifndef CGAL_BASIC_VIEWER_GLFW_H
#define CGAL_BASIC_VIEWER_GLFW_H
#include <iostream>
#include <stdlib.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Basic_shaders.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "bv_settings.h"
#include "internal/Shader.h"
#include "internal/Input.h"
#include "internal/utils.h"
#include "internal/Camera.h"
#include "internal/Line_renderer.h"
#include "internal/Clipping_plane.h"
#include "internal/Animation_controller.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace CGAL 
{
namespace GLFW 
{
  const int WINDOW_SAMPLES = CGAL_WINDOW_SAMPLES;

  class Basic_viewer : public Input
  {
  public:
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;

    enum class RenderingMode 
    {                  
      DRAW_ALL=-1,      // draw all
      DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
      DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
    };

    enum class DisplayMode 
    { 
      CLIPPING_PLANE_OFF=0, 
      CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF,
      CLIPPING_PLANE_SOLID_HALF_WIRE_HALF,
      CLIPPING_PLANE_SOLID_HALF_ONLY,  
    };

  public:
    Basic_viewer(
      const Graphics_scene* graphicScene,
      const char* title = "",
      bool drawVertices = false,
      bool drawEdges = true,
      bool drawFaces = true,
      bool drawRays = true,
      bool drawLines = true,
      bool useMonoColor = false,
      bool inverseNormal = false,
      bool flatShading = true
    );

    void show();
    void initialize(bool screenshotOnly=false);
    void make_screenshot(const std::string& filePath);

    void clear_application();

    /***** Getter & Setter ****/

    // Setter Section
    inline void set_scene(const Graphics_scene* scene)
    {
      m_scene = scene;
      m_areBuffersInitialized = false;
    }

    inline void window_size(const vec2f &size) { window_size_callback(m_window, size.x(), size.y()); }

    inline void vertices_mono_color(const CGAL::IO::Color& c) { m_verticeMonoColor = c; }
    inline void edges_mono_color(const CGAL::IO::Color& c) { m_edgesMonoColor = c; }
    inline void rays_mono_color(const CGAL::IO::Color& c) { m_raysMonoColor = c; }
    inline void lines_mono_color(const CGAL::IO::Color& c) { m_linesMonoColor = c; }
    inline void faces_mono_color(const CGAL::IO::Color& c) { m_facesMonoColor = c; }

    inline void size_points(const float size) { m_sizeVertices = size; }
    inline void size_edges(const float size) { m_sizeEdges = size; }
    inline void size_rays(const float size) { m_sizeRays = size; }
    inline void size_lines(const float size) { m_sizeLines = size; }

    inline void light_position(const vec4f& pos) { m_lightPosition = pos; }
    inline void light_ambient(const vec4f& color) { m_ambientColor = color; }
    inline void light_diffuse(const vec4f& color) { m_diffuseColor = color; }
    inline void light_specular(const vec4f& color) { m_specularColor = color; }
    inline void light_shininess(const float shininess) { m_shininess = shininess; }

    inline void draw_vertices(bool b) { m_drawVertices = b; }
    inline void draw_edges(bool b) { m_drawEdges = b; }
    inline void draw_rays(bool b) { m_drawRays = b; }
    inline void draw_lines(bool b) { m_drawLines = b; }
    inline void draw_faces(bool b) { m_drawFaces = b; }
    inline void use_mono_color(bool b) { m_useMonoColor = b; }
    inline void inverse_normal(bool b) { m_inverseNormal = b; }
    inline void flat_shading(bool b) { m_flatShading = b; }

    inline void set_azerty_layout() { set_keyboard_layout(KeyboardLayout::AZERTY); }
    inline void animation_duration(std::chrono::milliseconds duration) { m_animationController.set_duration(duration); }

    inline void center(const vec3f& center) { m_camera.set_center(center); }
    inline void center(float x, float y, float z) { m_camera.set_center({x, y, z}); }
    inline void position(const vec3f& position) { m_camera.set_position(position); m_camera.set_default_position(position); }
    inline void position(float x, float y, float z) { m_camera.set_position({x, y, z}); m_camera.set_default_position({x, y, z}); }
    inline void direction(const vec3f& direction) { m_camera.set_orientation(direction); m_camera.set_default_orientation(direction); }
    inline void direction(float x, float y, float z) { m_camera.set_orientation({x, y, z}); m_camera.set_default_orientation({x, y, z});}

    // Getter section
    inline vec3f position() const { return m_camera.get_position(); }
    inline vec3f forward() const { return m_camera.get_forward(); }

    inline CGAL::IO::Color vertices_mono_color() const { return m_verticeMonoColor; }
    inline CGAL::IO::Color edges_mono_color() const { return m_edgesMonoColor; }
    inline CGAL::IO::Color rays_mono_color() const { return m_raysMonoColor; }
    inline CGAL::IO::Color lines_mono_color() const { return m_linesMonoColor; }
    inline CGAL::IO::Color faces_mono_color() const { return m_facesMonoColor; }

    inline float size_points() const { return m_sizeVertices; }
    inline float size_edges() const { return m_sizeEdges; }
    inline float size_rays() const { return m_sizeRays; }
    inline float size_lines() const { return m_sizeLines; }

    inline vec4f light_position() const { return m_lightPosition; }
    inline vec4f light_ambient() const { return m_ambientColor; }
    inline vec4f light_diffuse() const { return m_diffuseColor; }
    inline vec4f light_specular() const { return m_specularColor; }
    inline float light_shininess() const { return m_shininess; }

    inline bool draw_vertices() const { return m_drawVertices; }
    inline bool draw_edges() const { return m_drawEdges; }
    inline bool draw_rays() const { return m_drawRays; }
    inline bool draw_lines() const { return m_drawLines; }
    inline bool draw_faces() const { return m_drawFaces; }
    inline bool use_mono_color() const { return m_useMonoColor; }
    inline bool inverse_normal() const { return m_inverseNormal; }
    inline bool flat_shading() const { return m_flatShading; }

    inline bool clipping_plane_enable() const { return m_displayMode != DisplayMode::CLIPPING_PLANE_OFF; }
    inline bool is_orthograpic() const { return m_camera.is_orthographic(); }

    CGAL::Plane_3<Local_kernel> clipping_plane() const;

    void exit_app();

  private:
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void cursor_callback(GLFWwindow* window, double xpos, double ypo);
    static void mouse_btn_callback(GLFWwindow* window, int button, int action, int mods);
    static void window_size_callback(GLFWwindow* window, int width, int height);
    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

    static GLFWwindow* create_window(int width, int height, const char* title, bool hidden = false);
    static void error_callback(int error, const char* description);

    void compile_shaders();
    void load_buffer(int i, int location, int gsEnum, int dataCount);
    void load_buffer(int i, int location, const std::vector<float>& vector, int dataCount);
    void initialize_buffers();
    void initialize_camera();
    void initialize_and_load_world_axis();
    void load_scene();

    void print_application_state(float& elapsedTime, const float deltaTime);

    void compute_model_view_projection_matrix(const float deltaTime);
    void update_uniforms(const float deltaTime);

    void update_face_uniforms();
    void update_pl_uniforms();
    void update_clipping_uniforms();
    void set_world_axis_uniforms();
    void set_XY_grid_uniforms();

    void render_scene(const float deltaTime=0.0);
    void draw_faces();
    void draw_rays();
    void draw_lines();

    void draw_faces_bis(RenderingMode mode);
    void draw_vertices(RenderingMode mode);
    void draw_edges(RenderingMode mode);

    void initialize_and_load_clipping_plane();
    void render_clipping_plane();

    void initialize_keys_actions();

    void start_action(int action, const float deltaTime) override;
    void action_event(int action, const float deltaTime) override;
    void end_action(int action, const float deltaTime) override;

    void double_click_event(int btn) override;
    void scroll_event(const float deltaTime) override;

    void rotate_camera();
    void rotate_clipping_plane();
    void translate_camera(const float deltaTime);
    void translate_clipping_plane(const float deltaTime, bool useCameraForward=false);

    void save_key_frame(); 
    void run_or_stop_animation();

    void increase_light_all(const float deltaTime);
    void increase_red_component(const float deltaTime);
    void increase_green_component(const float deltaTime);
    void increase_blue_component(const float deltaTime);

    void switch_display_mode();

    void fullscreen();
    void capture_screenshot(const std::string& filePath);

    vec4f color_to_vec4(const CGAL::IO::Color& c) const;

    void draw_world_axis();
    void draw_xy_grid();

    void reset_camera_and_clipping_plane();

  private:
    GLFWwindow* m_window;

    const Graphics_scene* m_scene;
    const char* m_title;
    bool m_drawVertices { false };
    bool m_drawEdges { true };
    bool m_drawFaces { true };
    bool m_drawRays { true };
    bool m_drawLines { true };
    bool m_useMonoColor { false };
    bool m_inverseNormal { false };
    bool m_flatShading { true };

    bool m_printApplicationState { true };

    bool m_areBuffersInitialized { false };

    bool m_drawWorldAxis { true };
    bool m_drawXYGrid { false };

    bool m_isOpengl4_3 { false };

    bool m_isFullscreen { false };
    
    Line_renderer m_worldAxisRenderer; 
    Line_renderer m_XYGridRenderer; 
    Line_renderer m_XYAxisRenderer; 
    Line_renderer m_clippingPlaneRenderer;

    float m_sizeVertices { CGAL_SIZE_POINTS };
    float m_sizeEdges { CGAL_SIZE_EDGES };
    float m_sizeRays { CGAL_SIZE_RAYS };
    float m_sizeLines { CGAL_SIZE_LINES };

    CGAL::IO::Color m_facesMonoColor = CGAL_FACES_MONO_COLOR;
    CGAL::IO::Color m_verticeMonoColor = CGAL_VERTICES_MONO_COLOR;
    CGAL::IO::Color m_edgesMonoColor = CGAL_EDGES_MONO_COLOR;
    CGAL::IO::Color m_raysMonoColor = CGAL_RAYS_MONO_COLOR;
    CGAL::IO::Color m_linesMonoColor = CGAL_LINES_MONO_COLOR;

    vec4f m_lightPosition { CGAL_LIGHT_POSITION };
    vec4f m_ambientColor { CGAL_AMBIENT_COLOR };
    vec4f m_diffuseColor { CGAL_DIFFUSE_COLOR };
    vec4f m_specularColor { CGAL_SPECULAR_COLOR };

    float m_shininess { CGAL_SHININESS };

    vec4f m_clipPlane { 0, 0, 1, 0 };
    vec4f m_pointPlane { 0, 0, 0, 1 };

    mat4f m_modelViewMatrix;
    mat4f m_modelViewProjectionMatrix;

    Shader m_plShader;
    Shader m_faceShader; 
    Shader m_planeShader;
    Shader m_lineShader;

    vec2i m_windowSize { CGAL_WINDOW_WIDTH_INIT, CGAL_WINDOW_HEIGHT_INIT };
    vec2i m_oldWindowSize;
    vec2i m_oldWindowPosition;

    Camera m_camera;

    /***************CLIPPING PLANE****************/

    Clipping_plane m_clippingPlane;

    DisplayMode m_displayMode { DisplayMode::CLIPPING_PLANE_OFF };
    
    bool m_drawClippingPlane { true };  // will be toggled when alt+c is pressed, which is used for indicating whether or not to render the clipping plane ;

    /*********************/

    Animation_controller m_animationController;

    enum VAOEnum 
    {
      VAO_MONO_POINTS = 0,
      VAO_COLORED_POINTS,
      VAO_MONO_SEGMENTS,
      VAO_COLORED_SEGMENTS,
      VAO_MONO_RAYS,
      VAO_COLORED_RAYS,
      VAO_MONO_LINES,
      VAO_COLORED_LINES,
      VAO_MONO_FACES,
      VAO_COLORED_FACES,
      NB_VAO_BUFFERS
    };

    GLuint m_vao[NB_VAO_BUFFERS];

    static const unsigned int NB_GL_BUFFERS = (Graphics_scene::END_POS - Graphics_scene::BEGIN_POS) +
                                              (Graphics_scene::END_COLOR - Graphics_scene::BEGIN_COLOR) + 2; // +2 for normals (mono and color)

    GLuint m_vbo[NB_GL_BUFFERS]; // +1 for the vbo buffer of clipping plane

    enum ActionEnum
    {
      /*APPLICATION*/
      EXIT,
      PRINT_APPLICATION_STATE,

      /*WINDOW*/
      FULLSCREEN,
      SCREENSHOT,

      /*SCENE*/
      VERTICES_DISPLAY,
      FACES_DISPLAY,
      EDGES_DISPLAY,
      INVERSE_NORMAL,
      SHADING_MODE,
      MONO_COLOR,
      INC_LIGHT_ALL,
      INC_LIGHT_R,
      INC_LIGHT_G,
      INC_LIGHT_B,
      DEC_LIGHT_ALL,
      DEC_LIGHT_R,
      DEC_LIGHT_G,
      DEC_LIGHT_B,
      INC_POINTS_SIZE,
      DEC_POINTS_SIZE,
      INC_EDGES_SIZE,
      DEC_EDGES_SIZE,

      DISPLAY_WORLD_AXIS, 
      DISPLAY_XY_GRID,

      /*CAMERA*/
      ROTATE_CAMERA,
      TRANSLATE_CAMERA,
      
      UP,
      LEFT,
      RIGHT,
      DOWN,
      FORWARD,
      BACKWARDS,

      SWITCH_CAM_MODE,
      SWITCH_CAM_TYPE,   
      RESET_CAM,

      INC_CAMERA_TRANSLATION_SPEED,
      DEC_CAMERA_TRANSLATION_SPEED,
      INC_CAMERA_ROTATION_SPEED,
      DEC_CAMERA_ROTATION_SPEED,

      INC_CAMERA_ROTATION_SMOOTHNESS, 
      DEC_CAMERA_ROTATION_SMOOTHNESS, 
      INC_CAMERA_TRANSLATION_SMOOTHNESS, 
      DEC_CAMERA_TRANSLATION_SMOOTHNESS, 

      /*CLIPPING PLANE*/
      ROTATE_CLIPPING_PLANE,
      TRANSLATE_CLIPPING_PLANE,
      TRANSLATE_CP_IN_CAMERA_DIRECTION,
      TRANSLATE_CP_IN_NORMAL_DIRECTION,
      TOGGLE_CLIPPING_PLANE_MODE,
      TOGGLE_CLIPPING_PLANE_DISPLAY,
      TOGGLE_CP_CONSTRAINT_AXIS,
      INC_CP_TRANSLATION_SPEED,
      INC_CP_ROTATION_SPEED,
      DEC_CP_TRANSLATION_SPEED,
      DEC_CP_ROTATION_SPEED,

      /*ANIMATION*/
      SAVE_KEY_FRAME,
      RUN_OR_STOP_ANIMATION,
      CLEAR_ANIMATION, 
    };
  };
} // end namespace GLFW
  using GLFW::Basic_viewer;

  inline void draw_graphics_scene(const Graphics_scene& graphics_scene, const char* title="CGAL Basic Viewer (GLFW)")
  {
    Basic_viewer basic_viewer(&graphics_scene, title);
    basic_viewer.show();
  } 
} // end namespace CGAL 

#include "Basic_viewer_impl.h"
  
#endif // CGAL_BASIC_VIEWER_GLFW_H 