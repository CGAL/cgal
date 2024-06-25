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

#include <CGAL/Graphics_scene.h>
#include <CGAL/Basic_shaders.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>

#include "Bv_Settings.h"
#include "internal/Shader.h"
#include "internal/Input.h"
#include "internal/utils.h"
#include "internal/Camera.h"
#include "internal/Line_renderer.h"

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

    enum class Rendering_mode : int 
    {                  
      DRAW_ALL=-1,      // draw all
      DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
      DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
    };

    enum class Display_mode : int 
    { 
      CLIPPING_PLANE_OFF=0, 
      CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF,
      CLIPPING_PLANE_SOLID_HALF_WIRE_HALF,
      CLIPPING_PLANE_SOLID_HALF_ONLY,  
    };

    enum class Constraint_axis : int 
    {
      NO_CONSTRAINT=0,
      X_AXIS,
      Y_AXIS,
      Z_AXIS, 
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
    void make_screenshot(const std::string& filePath);

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

    inline void size_points(const float size) { m_sizePoints = size; }
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

    // Getter section
    inline vec3f position() const { return m_camera.position(); }
    inline vec3f forward() const { return m_camera.forward(); }

    inline CGAL::IO::Color vertices_mono_color() const { return m_verticeMonoColor; }
    inline CGAL::IO::Color edges_mono_color() const { return m_edgesMonoColor; }
    inline CGAL::IO::Color rays_mono_color() const { return m_raysMonoColor; }
    inline CGAL::IO::Color lines_mono_color() const { return m_linesMonoColor; }
    inline CGAL::IO::Color faces_mono_color() const { return m_facesMonoColor; }

    inline float size_points() const { return m_sizePoints; }
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

    inline bool clipping_plane_enable() const { return m_displayModeEnum != Display_mode::CLIPPING_PLANE_OFF; }
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
    void init_and_load_renderers();
    void load_scene();

    void update_uniforms(const double deltaTime=0.00);

    void set_face_uniforms();
    void set_pl_uniforms();
    void set_clipping_uniforms();
    void set_world_axis_uniforms();
    void set_XY_grid_uniforms();

    void render_scene(const double deltaTime=0.00);
    void draw_faces();
    void draw_rays();
    void draw_lines();

    void draw_faces_bis(Rendering_mode mode);
    void draw_vertices(Rendering_mode mode);
    void draw_edges(Rendering_mode mode);

    void init_and_load_clipping_plane();
    void render_clipping_plane();

    void init_keys_actions();

    void start_action(int action, const double deltaTime) override;
    void action_event(int action, const double deltaTime) override;
    void end_action(int action, const double deltaTime) override;

    void double_click_event(int btn) override;
    void scroll_event(const double deltaTime) override;

    void translate(const vec3f dir);
    void mouse_rotate();
    void mouse_translate(const double deltaTime);

    vec2f to_ndc(double, double);
    vec3f mapping_cursor_toHemisphere(double x, double y);
    mat4f get_rotation(const vec3f& start, const vec3f& end);
    void rotate_clipping_plane();

    void translate_clipping_plane();
    void translate_clipping_plane_cam_dir();

    void switch_constraint_axis();
    void switch_display_mode();

    void fullscreen();
    void screenshot(const std::string& filePath);

    void print_help();

    vec4f color_to_vec4(const CGAL::IO::Color& c) const;

    void draw_world_axis();
    void draw_xy_grid();

  private:
    enum ActionEnum
    {
      MOUSE_ROTATE,
      MOUSE_TRANSLATE,
      UP,
      LEFT,
      RIGHT,
      DOWN,
      FORWARD,
      BACKWARDS,
      SWITCH_CAM_MODE,
      SWITCH_CAM_ROTATION,

      FULLSCREEN,
      SCREENSHOT,
      INC_ZOOM,
      DEC_ZOOM,
      INC_MOVE_SPEED_1,
      DEC_MOVE_SPEED_1,
      INC_ROT_SPEED_1,
      DEC_ROT_SPEED_1,
      RESET_CAM,

      CLIPPING_PLANE_MODE,
      CLIPPING_PLANE_DISPLAY,
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

      CP_ROTATION,
      CP_TRANSLATION,
      CP_TRANS_CAM_DIR,
      CP_TRANS_N_DIR,
      CONSTRAINT_AXIS,

      INC_ROT_SMOOTHNESS, 
      DEC_ROT_SMOOTHNESS, 
      INC_TRA_SMOOTHNESS, 
      DEC_TRA_SMOOTHNESS, 

      EXIT
    };

    GLFWwindow* m_window;

    const Graphics_scene* m_scene;
    const char* m_title;
    bool m_drawVertices;
    bool m_drawEdges;
    bool m_drawFaces;
    bool m_drawRays;
    bool m_drawLines;
    bool m_useMonoColor;
    bool m_inverseNormal;
    bool m_flatShading;

    bool m_areBuffersInitialized = false;

    bool m_drawWorldAxis = true;
    bool m_drawXYGrid = false;

    bool m_isOpengl4_3 = false;

    bool m_isFullscreen = false;
    
    Line_renderer m_worldAxisRenderer; 
    Line_renderer m_XYGridRenderer; 
    Line_renderer m_XYAxisRenderer; 
    Line_renderer m_clippingPlaneRenderer;

    float m_sizePoints = CGAL_SIZE_POINTS;
    float m_sizeEdges = CGAL_SIZE_EDGES;
    float m_sizeRays = CGAL_SIZE_RAYS;
    float m_sizeLines = CGAL_SIZE_LINES;

    CGAL::IO::Color m_facesMonoColor = CGAL_FACES_MONO_COLOR;
    CGAL::IO::Color m_verticeMonoColor = CGAL_VERTICES_MONO_COLOR;
    CGAL::IO::Color m_edgesMonoColor = CGAL_EDGES_MONO_COLOR;
    CGAL::IO::Color m_raysMonoColor = CGAL_RAYS_MONO_COLOR;
    CGAL::IO::Color m_linesMonoColor = CGAL_LINES_MONO_COLOR;

    vec4f m_lightPosition = CGAL_LIGHT_POSITION;
    vec4f m_ambientColor = CGAL_AMBIENT_COLOR;
    vec4f m_diffuseColor = CGAL_DIFFUSE_COLOR;
    vec4f m_specularColor = CGAL_SPECULAR_COLOR;

    float m_shininess = CGAL_SHININESS;

    vec4f m_clipPlane{0, 0, 1, 0};
    vec4f m_pointPlane{0, 0, 0, 1};

    mat4f m_modelViewMatrix;
    mat4f m_modelViewProjectionMatrix;

    Shader m_plShader;
    Shader m_faceShader; 
    Shader m_planeShader;
    Shader m_lineShader;

    vec2i m_windowSize{CGAL_WINDOW_WIDTH_INIT, CGAL_WINDOW_HEIGHT_INIT};
    vec2i m_oldWindowSize;
    vec2i m_oldWindowPosition;

    Camera m_camera;

    /***************CLIPPING PLANE****************/

    Display_mode m_displayModeEnum = Display_mode::CLIPPING_PLANE_OFF;
    Constraint_axis m_constraintAxisEnum = Constraint_axis::NO_CONSTRAINT;
    
    std::vector<float> m_clippingPlaneVertices;

    bool m_drawClippingPlane = true;                                                // will be toggled when alt+c is pressed, which is used for indicating whether or not to render the clipping plane ;
    float m_clippingPlaneTransparency = CGAL_CLIPPING_PLANE_RENDERING_TRANSPARENCY; // to what extent the transparent part should be rendered;
    float m_clippingPlaneMoveSpeed = CGAL_CLIPPING_PLANE_MOVE_SPEED;
    float m_clippingPlaneRotationSpeed = CGAL_CLIPPING_PLANE_ROT_SPEED;


    vec3f m_constraintAxis{1., 0., 0.};

    mat4f m_clippingMatrix = mat4f::Identity();

    /*********************/

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