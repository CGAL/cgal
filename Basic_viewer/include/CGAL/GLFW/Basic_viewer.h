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

#pragma once

#include <CGAL/Graphics_scene.h>
#include <CGAL/Basic_shaders.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>

#include "Shader.h"
#include "Input.h"
#include "Bv_Settings.h"
#include "math.h"
#include "Camera.h"
#include "Line_renderer.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace CGAL {
namespace GLFW {
  enum RenderMode
  {                   // rendering mode
    DRAW_ALL = -1,    // draw all
    DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
    DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
  };

  enum CAM_MODE
  {
    PERSPECTIVE,
    ORTHOGRAPHIC
  };
  enum CAM_ROTATION_MODE
  {
    OBJECT,
    FREE
  };

  enum ClippingMode
  { // clipping mode
    CLIPPING_PLANE_OFF = 0,
    CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF,
    CLIPPING_PLANE_SOLID_HALF_WIRE_HALF,
    CLIPPING_PLANE_SOLID_HALF_ONLY,
    CLIPPING_PLANE_END_INDEX
  };

  const int windowSamples = WINDOW_SAMPLES;

  void glfwErrorCallback(int error, const char *description);

  class Basic_viewer : public Input
  {
  public:
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;

  public:
    Basic_viewer(const Graphics_scene *graphicScene,
                 const char *title = "",
                 bool drawVertices = false,
                 bool drawEdges = true,
                 bool drawFaces = true,
                 bool useMonoColor = false,
                 bool inverseNormal = false,
                 bool drawRays = true,
                 bool drawLines = true);

    void show();
    void make_screenshot(const std::string &pngpath);

    /***** Getter & Setter ****/

    // Setter Section
    inline void position(const vec3f &pos) { m_camPosition = pos; }
    inline void forward(const vec3f &dir) { m_camForward = dir; }
    inline void set_scene(const Graphics_scene *scene)
    {
      m_scene = scene;
      m_areBuffersInitialized = false;
    }
    inline void window_size(const vec2f &size)
    {
      window_size_callback(m_window, size.x(), size.y());
    }

    inline void exit_app();

    inline void vertices_mono_color(const CGAL::IO::Color &c) { m_verticeMonoColor = c; }
    inline void edges_mono_color(const CGAL::IO::Color &c) { m_edgesMonoColor = c; }
    inline void rays_mono_color(const CGAL::IO::Color &c) { m_raysMonoColor = c; }
    inline void lines_mono_color(const CGAL::IO::Color &c) { m_linesMonoColor = c; }
    inline void faces_mono_color(const CGAL::IO::Color &c) { m_facesMonoColor = c; }

    inline void size_points(const float size) { m_sizePoints = size; }
    inline void size_edges(const float size) { m_sizeEdges = size; }
    inline void size_rays(const float size) { m_sizeRays = size; }
    inline void size_lines(const float size) { m_sizeLines = size; }

    inline void light_position(const vec4f &pos) { m_lightPosition = pos; }
    inline void light_ambient(const vec4f &color) { m_ambient = color; }
    inline void light_diffuse(const vec4f &color) { m_diffuse = color; }
    inline void light_specular(const vec4f &color) { m_specular = color; }
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
    inline vec3f position() const { return m_camPosition; }
    inline vec3f forward() const { return m_camForward; }

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
    inline vec4f light_ambient() const { return m_ambient; }
    inline vec4f light_diffuse() const { return m_diffuse; }
    inline vec4f light_specular() const { return m_specular; }
    inline float light_shininess() const { return m_shininess; }

    inline bool draw_vertices() const { return m_drawVertices; }
    inline bool draw_edges() const { return m_drawEdges; }
    inline bool draw_rays() const { return m_drawRays; }
    inline bool draw_lines() const { return m_drawLines; }
    inline bool draw_faces() const { return m_drawFaces; }
    inline bool use_mono_color() const { return m_useMonoColor; }
    inline bool inverse_normal() const { return m_inverseNormal; }
    inline bool flat_shading() const { return m_flatShading; }

    inline bool clipping_plane_enable() const { return m_useClippingPlane != CLIPPING_PLANE_OFF; }
    inline bool is_orthograpic() const { return m_camera.is_orthographic(); }

    CGAL::Plane_3<Local_kernel> clipping_plane() const;

  private:
    static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);
    static void cursor_callback(GLFWwindow *window, double xpos, double ypo);
    static void mouse_btn_callback(GLFWwindow *window, int button, int action, int mods);
    static void window_size_callback(GLFWwindow *window, int width, int height);
    static void scroll_callback(GLFWwindow *window, double xoffset, double yoffset);

    static GLFWwindow *create_window(int width, int height, const char *title, bool hidden = false);
    static void error_callback(int error, const char *description);

    void compile_shaders();
    void load_buffer(int i, int location, int gsEnum, int dataCount);
    void load_buffer(int i, int location, const std::vector<float> &vector, int dataCount);
    void init_buffers();
    void load_scene();

    void update_frame();
    void update_uniforms();

    void set_face_uniforms();
    void set_pl_uniforms();
    void set_clipping_uniforms();

    void render_scene();
    void draw_faces();
    void draw_rays();
    void draw_lines();

    void draw_faces_(RenderMode mode);
    void draw_vertices(RenderMode mode);
    void draw_edges(RenderMode mode);

    void generate_clipping_plane();
    void render_clipping_plane();

    void init_keys_actions();

    void start_action(ActionEnum action) override;
    void action_event(ActionEnum action) override;
    void end_action(ActionEnum action) override;

    void double_click_event(int btn) override;
    void scroll_event() override;

    void translate(const vec3f dir);
    void mouse_rotate();
    void mouse_translate();
    void set_cam_mode(CAM_MODE mode);
    void switch_rotation_mode();

    vec2f to_ndc(double, double);
    vec3f mapping_cursor_toHemisphere(double x, double y);
    mat4f get_rotation(vec3f const &start, vec3f const &end);
    void rotate_clipping_plane();

    void translate_clipping_plane();
    void translate_clipping_plane_cam_dir();
    // void translate_clipping_plane_n_dir();

    void switch_axis(int axis);

    // void scroll_event(float z);
    void fullscreen();
    void screenshot(const std::string &pngpath);

    void print_help();

    vec4f color_to_vec4(const CGAL::IO::Color &c) const;

    void draw_world_axis();
    void draw_xy_grid();

  private:
    GLFWwindow *m_window;
    const Graphics_scene *m_scene;
    const char *m_title;
    bool m_drawVertices;
    bool m_drawEdges;
    bool m_drawRays;
    bool m_drawLines;
    bool m_drawFaces;
    bool m_areBuffersInitialized = false;
    bool m_flatShading = true;
    bool m_useMonoColor;
    bool m_inverseNormal;

    bool m_draw_world_axis = true;
    bool m_draw_xy_grid = false;
    Line_renderer m_corner_world_axis; 
    Line_renderer m_world_axis; 

    float m_sizePoints = SIZE_POINTS;
    float m_sizeEdges = SIZE_EDGES;
    float m_sizeRays = SIZE_RAYS;
    float m_sizeLines = SIZE_LINES;

    double m_delta_time = 0; 

    CGAL::IO::Color m_facesMonoColor = FACES_MONO_COLOR;
    CGAL::IO::Color m_verticeMonoColor = VERTICES_MONO_COLOR;
    CGAL::IO::Color m_edgesMonoColor = EDGES_MONO_COLOR;
    CGAL::IO::Color m_raysMonoColor = RAYS_MONO_COLOR;
    CGAL::IO::Color m_linesMonoColor = LINES_MONO_COLOR;

    vec4f m_lightPosition = LIGHT_POSITION;
    vec4f m_ambient = AMBIENT_COLOR;
    vec4f m_diffuse = DIFFUSE_COLOR;
    vec4f m_specular = SPECULAR_COLOR;
    float m_shininess = SHININESS;

    vec4f m_clip_plane{0, 0, 1, 0};
    vec4f m_point_plane{0, 0, 0, 1};

    mat4f m_mv;
    mat4f m_mvp;
    bool m_is_opengl_4_3 = false;

    Shader m_pl_shader, m_face_shader, m_plane_shader, m_line_shader;

    /******* CAMERA ******/

    float m_cam_speed = CAM_MOVE_SPEED;
    float m_cam_rotation_speed = CAM_ROT_SPEED;
    float m_scene_rotation_speed = SCENE_ROT_SPEED;

    mat4f m_cam_projection;
    vec3f m_camPosition{0, 0, -5};
    vec2f m_cam_view{0, 0};
    vec3f m_camForward{0, 0, 1};
    float m_cam_orth_zoom = 1.0f;

    vec2f m_scene_view{0.0f, 0.0f};
    mat4f m_scene_rotation = mat4f::Identity();

    vec2i m_window_size{WINDOW_WIDTH_INIT, WINDOW_HEIGHT_INIT};
    vec2i m_old_window_size;
    vec2i m_old_window_pos;

    bool m_is_fullscreen = false;
    CAM_MODE m_camMode = PERSPECTIVE;
    CAM_ROTATION_MODE m_cam_rotation_mode = OBJECT;

    Camera m_camera;

    /***************CLIPPING PLANE****************/

    ClippingMode m_useClippingPlane = CLIPPING_PLANE_OFF;
    std::vector<float> m_array_for_clipping_plane;

    bool m_clipping_plane_rendering = true;                                                // will be toggled when alt+c is pressed, which is used for indicating whether or not to render the clipping plane ;
    float m_clipping_plane_rendering_transparency = CLIPPING_PLANE_RENDERING_TRANSPARENCY; // to what extent the transparent part should be rendered;
    float m_clipping_plane_move_speed = CLIPPING_PLANE_MOVE_SPEED;
    float m_clipping_plane_rot_speed = CLIPPING_PLANE_ROT_SPEED;

    int m_cstr_axis_enum = NO_AXIS;
    vec3f m_cstr_axis{1., 0., 0.};

    mat4f m_clipping_matrix = mat4f::Identity();

    enum Axis
    {
      NO_AXIS = 0,
      X_AXIS,
      Y_AXIS,
      Z_AXIS,
      NB_AXIS_ENUM
    };

    /*********************/

    enum Actions
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

      EXIT
    };

    enum VAO_TYPES
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
      VAO_CLIPPING_PLANE,
      NB_VAO_BUFFERS
    };

    GLuint m_vao[NB_VAO_BUFFERS];

    static const unsigned int NB_GL_BUFFERS = (Graphics_scene::END_POS - Graphics_scene::BEGIN_POS) +
                                              (Graphics_scene::END_COLOR - Graphics_scene::BEGIN_COLOR) + 3; // +2 for normals (mono and color), +1 for clipping plane

    GLuint m_vbo[NB_GL_BUFFERS]; // +1 for the vbo buffer of clipping plane
  };
} // end namespace GLFW
  using GLFW::Basic_viewer;

  inline void draw_graphics_scene(const Graphics_scene &graphics_scene, const char *title="CGAL Basic Viewer (GLFW)")
  {
    Basic_viewer basic_viewer(&graphics_scene, title);
    basic_viewer.show();
  }
} // end namespace CGAL 

#include "Basic_viewer_impl.h"
