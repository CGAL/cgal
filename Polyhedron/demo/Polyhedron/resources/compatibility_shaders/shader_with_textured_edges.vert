//#version 100 
 attribute highp vec4 vertex;
 attribute highp vec2 v_texCoord; 
 uniform highp vec3 color_lines; 
 uniform highp mat4 mvp_matrix; 
 varying highp vec3 fColors; 
 varying highp vec2 f_texCoord; 
  
  uniform highp float point_size;
  void main(void)
  {
    gl_PointSize = point_size;
    f_texCoord = v_texCoord; 
    fColors = color_lines; 
    gl_Position = mvp_matrix * vertex; 
 }  
