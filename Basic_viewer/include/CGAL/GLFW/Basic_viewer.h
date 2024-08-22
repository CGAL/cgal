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
#include <memory>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

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
#include "internal/buffer/VAO.h"

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
      const char* title = "CGAL Basic Viewer (GLFW)",
      bool drawVertices = false,
      bool drawEdges = true,
      bool drawFaces = true,
      bool drawRays = true,
      bool drawLines = true,
      bool useDefaultColor = false,
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
      m_Scene = scene;
      m_AreBuffersInitialized = false;
    }

    inline void window_size(const vec2f &size) { window_size_callback(m_Window, size.x(), size.y()); }

    inline void size_points(const float size) { m_SizeVertices = size; }
    inline void size_edges(const float size) { m_SizeEdges = size; }
    inline void size_rays(const float size) { m_SizeRays = size; }
    inline void size_lines(const float size) { m_SizeLines = size; }

    inline void light_position(const vec4f& pos) { m_LightPosition = pos; }
    inline void light_ambient(const vec4f& color) { m_AmbientColor = color; }
    inline void light_diffuse(const vec4f& color) { m_DiffuseColor = color; }
    inline void light_specular(const vec4f& color) { m_SpecularColor = color; }
    inline void light_shininess(const float shininess) { m_Shininess = shininess; }

    inline void draw_vertices(bool b) { m_DrawVertices = b; }
    inline void draw_edges(bool b) { m_DrawEdges = b; }
    inline void draw_rays(bool b) { m_DrawRays = b; }
    inline void draw_lines(bool b) { m_DrawLines = b; }
    inline void draw_faces(bool b) { m_DrawFaces = b; }
    inline void use_mono_color(bool b) { m_UseDefaultColor = b; }
    inline void inverse_normal(bool b) { m_InverseNormal = b; }
    inline void flat_shading(bool b) { m_FlatShading = b; }
    inline void draw_world_axis(bool b) { m_DrawWorldAxis = b; }
    inline void draw_xy_grid(bool b) { m_DrawXYGrid = b; }

    inline void display_mode(DisplayMode mode) { m_DisplayMode = mode; }
    inline void draw_clipping_plane(bool b) { m_DrawClippingPlane = b; }

    inline void two_dimensional() { m_Camera.set_orthographic(); } 

    inline void azerty_layout() { set_keyboard_layout(KeyboardLayout::AZERTY); }
    inline void animation_duration(std::chrono::milliseconds duration) { m_AnimationController.set_duration(duration); }

    inline void scene_radius(float radius) { m_Camera.set_radius(radius); }
    inline void scene_center(const vec3f& center) { m_Camera.set_center(center); }
    inline void scene_center(float x, float y, float z) { m_Camera.set_center({x, y, z}); }
    inline void camera_position(const vec3f& position) { m_Camera.set_position(position); }
    inline void camera_position(float x, float y, float z) { m_Camera.set_position({x, y, z}); }
    inline void camera_orientation(const vec3f& forward, float upAngle) { m_Camera.set_orientation(forward, upAngle); }
    inline void zoom(float zoom) { m_Camera.set_size(zoom); }
    
    inline void align_camera_to_clipping_plane() { m_Camera.align_to_plane(m_ClippingPlane.get_normal()); }
    inline void clipping_plane_orientation(const vec3f& normal) { m_ClippingPlane.set_orientation(normal); } 
    inline void clipping_plane_translate_along_normal(float t) { m_ClippingPlane.translation(t*.1); } 
    inline void clipping_plane_translate_along_camera_forward(float t) { m_ClippingPlane.translation(m_Camera.get_forward(), t*.1); } 

    // Getter section
    inline vec3f position() const { return m_Camera.get_position(); }
    inline vec3f forward() const { return m_Camera.get_forward(); }
    inline vec3f right() const { return m_Camera.get_right(); }
    inline vec3f up() const { return m_Camera.get_up(); }

    inline const CGAL::IO::Color& vertices_mono_color() const { return m_Scene->get_default_color_point(); }
    inline const CGAL::IO::Color& edges_mono_color() const { return m_Scene->get_default_color_segment(); }
    inline const CGAL::IO::Color& rays_mono_color() const { return m_Scene->get_default_color_ray(); }
    inline const CGAL::IO::Color& lines_mono_color() const { return m_Scene->get_default_color_line(); }
    inline const CGAL::IO::Color& faces_mono_color() const { return m_Scene->get_default_color_face(); }

    inline float size_points() const { return m_SizeVertices; }
    inline float size_edges() const { return m_SizeEdges; }
    inline float size_rays() const { return m_SizeRays; }
    inline float size_lines() const { return m_SizeLines; }

    inline vec4f light_position() const { return m_LightPosition; }
    inline vec4f light_ambient() const { return m_AmbientColor; }
    inline vec4f light_diffuse() const { return m_DiffuseColor; }
    inline vec4f light_specular() const { return m_SpecularColor; }
    inline float light_shininess() const { return m_Shininess; }

    inline bool draw_vertices() const { return m_DrawVertices; }
    inline bool draw_edges() const { return m_DrawEdges; }
    inline bool draw_rays() const { return m_DrawRays; }
    inline bool draw_lines() const { return m_DrawLines; }
    inline bool draw_faces() const { return m_DrawFaces; }
    inline bool use_mono_color() const { return m_UseDefaultColor; }
    inline bool inverse_normal() const { return m_InverseNormal; }
    inline bool flat_shading() const { return m_FlatShading; }

    inline bool clipping_plane_enable() const { return m_DisplayMode != DisplayMode::CLIPPING_PLANE_OFF; }
    inline bool is_orthograpic() const { return m_Camera.is_orthographic(); }
    inline bool is_two_dimensional() const { return !is_orthograpic() && m_Scene->is_two_dimensional(); }

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

    void update_face_uniforms();
    void update_sphere_uniforms();
    void update_cylinder_uniforms();
    void update_pl_uniforms(const vec3f& defaultColor={0,0,0});
    void update_line_uniforms(float size, const vec3f& defaultColor={0,0,0});
    void update_clipping_uniforms();
    void update_world_axis_uniforms();
    void update_XY_axis_uniforms();
    void update_XY_grid_uniforms();
    void update_normals_uniforms();
    void update_triangles_uniforms();

    void generate_grid(Line_renderer& renderer, float size, int nbSubdivisions=10) const;
    void generate_grid(Line_renderer& renderer, const vec3f& color, float size, int nbSubdivisions=10) const;

    void render_scene(const float deltaTime);

    void draw(const float deltaTime=0);
    
    void draw_rays();
    void draw_edges();
    void draw_lines();
    void draw_vertices();

    void draw_faces();
    void draw_faces_bis(RenderingMode mode);

    void draw_world_axis();
    void draw_xy_grid();
    void draw_normals();
    void draw_triangles();

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

    void increase_size_edge(const float deltaTime);
    void decrease_size_edge(const float deltaTime);
    void increase_size_vertex(const float deltaTime);
    void decrease_size_vertex(const float deltaTime);

    void switch_display_mode();

    void fullscreen();
    void capture_screenshot(const std::string& filePath);

    vec4f color_to_normalized_vec4(const CGAL::IO::Color& c) const;
    vec3f color_to_normalized_vec3(const CGAL::IO::Color& c) const;

    void reset_camera_and_clipping_plane();
    void reset_clipping_plane();

    void change_pivot_point(); 
    void check_geometry_feature_availability(); 

    bool need_update() const; 
    
    std::vector<float> aggregating_data(int monoEnum, int coloredEnum) const;
    std::vector<float> aggregating_color_data(int monoEnum, int coloredEnum, const CGAL::IO::Color& monoColor) const;
 
  private:
    GLFWwindow* m_Window;

    const Graphics_scene* m_Scene;
    const char* m_Title;
    bool m_DrawVertices { false };
    bool m_DrawEdges { true };
    bool m_DrawFaces { true };
    bool m_DrawRays { true };
    bool m_DrawLines { true };
    bool m_DrawCylinderEdge { false };
    bool m_DrawSphereVertex { false };
    bool m_DrawTriangles { false };

    bool m_DisplayFaceNormal { false };
    bool m_UseNormalMonoColor { false };
    bool m_UseDefaultColor { false };

    bool m_InverseNormal { false };
    bool m_FlatShading { true };

    bool m_PrintApplicationState { true };
    bool m_AreBuffersInitialized { false };

    bool m_DrawWorldAxis { true };
    bool m_DrawXYGrid { false };
    bool m_DrawNormals { false };

    bool m_IsOpengl4_3 { false };
    bool m_IsFullscreen { false };

    bool m_GeometryFeatureEnabled { true };    
    
    Line_renderer m_WorldAxisRenderer; 
    Line_renderer m_XYGridRenderer; 
    Line_renderer m_XYAxisRenderer; 

    float m_SizeVertices { CGAL_SIZE_POINTS };
    float m_SizeEdges    { CGAL_SIZE_EDGES  };
    float m_SizeRays     { CGAL_SIZE_RAYS   };
    float m_SizeLines    { CGAL_SIZE_LINES  };
    float m_SizeNormals  { CGAL_SIZE_NORMALS};

    float m_NormalHeightFactor { CGAL_NORMAL_HEIGHT_FACTOR };

    CGAL::IO::Color m_NormalsMonoColor  = CGAL_NORMALS_MONO_COLOR;

    vec4f m_LightPosition { CGAL_LIGHT_POSITION };
    vec4f m_AmbientColor  { CGAL_AMBIENT_COLOR  };
    vec4f m_DiffuseColor  { CGAL_DIFFUSE_COLOR  };
    vec4f m_SpecularColor { CGAL_SPECULAR_COLOR };

    float m_Shininess { CGAL_SHININESS };

    vec3f m_DefaultColorRay;
    vec3f m_DefaultColorFace;
    vec3f m_DefaultColorLine;
    vec3f m_DefaultColorPoint;
    vec3f m_DefaultColorSegment;

    vec4f m_ClipPlane { 0, 0, 1, 0 };
    vec4f m_PointPlane { 0, 0, 0, 1 };

    mat4f m_ModelMatrix          { mat4f::Identity() };
    mat4f m_ViewMatrix           { mat4f::Identity() };
    mat4f m_ProjectionMatrix     { mat4f::Identity() };
    mat4f m_ViewProjectionMatrix { mat4f::Identity() };

    std::shared_ptr<Shader> m_ShaderFace; 
    std::shared_ptr<Shader> m_ShaderSphere;
    std::shared_ptr<Shader> m_ShaderPoint;
    std::shared_ptr<Shader> m_ShaderLine;
    std::shared_ptr<Shader> m_ShaderPl;
    std::shared_ptr<Shader> m_ShaderCylinder;
    std::shared_ptr<Shader> m_ShaderPlane;
    std::shared_ptr<Shader> m_ShaderGrid;
    std::shared_ptr<Shader> m_ShaderNormal;
    std::shared_ptr<Shader> m_ShaderArrow;
    std::shared_ptr<Shader> m_ShaderTriangles;

    vec2i m_WindowSize { CGAL_WINDOW_WIDTH_INIT, CGAL_WINDOW_HEIGHT_INIT };
    vec2i m_OldWindowSize { CGAL_WINDOW_WIDTH_INIT, CGAL_WINDOW_HEIGHT_INIT };

    vec2i m_OldWindowPosition { 0, 0 }; 

    float m_AspectRatio { 1.f };

    float m_DeltaTime { 0 };

    std::pair<vec3f, vec3f> m_BoundingBox; 

    Camera m_Camera;

    /***************CLIPPING PLANE****************/

    Clipping_plane m_ClippingPlane;

    DisplayMode m_DisplayMode { DisplayMode::CLIPPING_PLANE_OFF };
    
    bool m_DrawClippingPlane { true };  // will be toggled when alt+c is pressed, which is used for indicating whether or not to render the clipping plane ;

    /*********************/

    Animation_controller m_AnimationController;

    enum VAOEnum 
    {
      VAO_POINTS = 0,
      VAO_SEGMENTS,
      VAO_RAYS,
      VAO_LINES,
      VAO_FACES,
      NB_VAO_BUFFERS
    };

    unsigned int m_VAO[NB_VAO_BUFFERS];

    static const unsigned int NB_GL_BUFFERS = (Graphics_scene::END_POS - Graphics_scene::BEGIN_POS) +
                                              (Graphics_scene::END_COLOR - Graphics_scene::BEGIN_COLOR) + 2; // +2 for normals (mono and color)

    unsigned int m_VBO[NB_GL_BUFFERS]; // +1 for the vbo buffer of clipping plane

    enum ActionEnum
    {
      /*APPLICATION*/
      EXIT,
      PRINT_APPLICATION_STATE,

      /*WINDOW*/
      FULLSCREEN,
      SCREENSHOT,

      /*SCENE*/
      NORMALS_DISPLAY,
      TRIANGLES_DISPLAY,
      VERTICES_DISPLAY,
      FACES_DISPLAY,
      EDGES_DISPLAY,
      INVERSE_NORMAL,
      SHADING_MODE,

      FACE_NORMALS_DISPLAY,
      CYLINDER_EDGE_DISPLAY,
      SPHERE_VERTEX_DISPLAY,

      NORMALS_MONO_COLOR,
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

      WORLD_AXIS_DISPLAY, 
      XY_GRID_DISPLAY,

      SHOW_ENTIRE_SCENE,

      /*CAMERA*/
      ROTATE_CAMERA,
      TRANSLATE_CAMERA,
      
      UP,
      LEFT,
      RIGHT,
      DOWN,
      FORWARD,
      BACKWARDS,

      SWITCH_CAMERA_MODE,
      SWITCH_CAMERA_TYPE,
      SWITCH_CAMERA_CONSTRAINT_AXIS,

      RESET_CAMERA_AND_CP,

      INC_CAMERA_TRANSLATION_SPEED,
      DEC_CAMERA_TRANSLATION_SPEED,
      INC_CAMERA_ROTATION_SPEED,
      DEC_CAMERA_ROTATION_SPEED,

      INC_CAMERA_ROTATION_SMOOTHNESS, 
      DEC_CAMERA_ROTATION_SMOOTHNESS, 
      INC_CAMERA_TRANSLATION_SMOOTHNESS, 
      DEC_CAMERA_TRANSLATION_SMOOTHNESS, 

      CHANGE_PIVOT_POINT, 

      /*CLIPPING PLANE*/
      ROTATE_CLIPPING_PLANE,
      TRANSLATE_CLIPPING_PLANE,
      TRANSLATE_CP_ALONG_CAMERA_DIRECTION,
      TRANSLATE_CP_ALONG_NORMAL_DIRECTION,

      SWITCH_CLIPPING_PLANE_MODE,
      SWITCH_CLIPPING_PLANE_DISPLAY,
      SWITCH_CP_CONSTRAINT_AXIS,

      INC_CP_TRANSLATION_SPEED,
      INC_CP_ROTATION_SPEED,
      DEC_CP_TRANSLATION_SPEED,
      DEC_CP_ROTATION_SPEED,

      RESET_CLIPPING_PLANE,

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