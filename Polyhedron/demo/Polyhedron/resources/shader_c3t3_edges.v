#version 120
attribute highp vec4 vertex;
attribute highp vec3 colors;
uniform highp mat4 mvp_matrix;
uniform highp vec4 cutplane;
varying highp vec4 color; 
void main(void)
{
  if(vertex.x * cutplane.x  + vertex.y * cutplane.y  + vertex.z * cutplane.z  +  cutplane.w > 0)
  {
    color = vec4(0.0, 1.0, 1.0, 0.0);
  } 
  else 
  {
    color = vec4(colors, 1.0);
  }
  gl_Position = mvp_matrix * vertex;
}
 
