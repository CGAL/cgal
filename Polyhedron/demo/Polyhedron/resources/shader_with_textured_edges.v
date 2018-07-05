#version 150
 in vec4 vertex;
 in vec2 v_texCoord; 
 uniform vec3 color_lines; 
 uniform mat4 mvp_matrix; 
 out vec3 fColors; 
 out vec2 f_texCoord; 
  
  uniform float point_size;
  void main(void)
  {
    gl_PointSize = point_size;
    f_texCoord = v_texCoord; 
    fColors = color_lines; 
    gl_Position = mvp_matrix * vertex; 
 }  
