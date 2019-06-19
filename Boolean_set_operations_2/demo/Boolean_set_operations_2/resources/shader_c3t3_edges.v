#version 120
attribute highp vec4 vertex;
attribute highp vec3 colors;
uniform highp mat4 mvp_matrix;
uniform highp vec4 cutplane;
varying highp vec4 color; 
void main(void)
{
  color = vec4(colors, vertex.x * cutplane.x  + vertex.y * cutplane.y  + vertex.z * cutplane.z  +  cutplane.w);
  gl_Position = mvp_matrix * vertex;
}
