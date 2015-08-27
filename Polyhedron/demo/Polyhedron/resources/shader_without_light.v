#version 120
attribute highp vec4 vertex;
attribute highp vec3 colors;
uniform highp mat4 mvp_matrix;
uniform highp mat4 f_matrix;
varying highp vec4 color; 
void main(void)
{
   color = vec4(colors, 1.0); 
   gl_Position = mvp_matrix * f_matrix * vertex;
}
 
