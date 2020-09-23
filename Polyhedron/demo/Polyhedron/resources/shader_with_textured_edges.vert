#version 150
 in vec4 vertex;
 in vec2 v_texCoord; 
 uniform mat4 mvp_matrix; 
 out vec2 f_texCoord; 
  
  uniform float point_size;
  void main(void)
  {
    gl_PointSize = point_size;
    f_texCoord = v_texCoord; 
    gl_Position = mvp_matrix * vertex; 
 }  
