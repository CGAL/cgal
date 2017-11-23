#version 120
 attribute highp vec4 vertex;
 attribute highp vec2 v_texCoord; 
 uniform highp vec3 colors;
 uniform highp mat4 mvp_matrix; 
 varying highp vec3 fColors; 
 varying highp vec2 f_texCoord; 
  
 void main(void) 
 { 
    f_texCoord = v_texCoord; 
    fColors = colors;
    gl_Position = mvp_matrix * vertex; 
 }  
