#version 150
in vec4 vertex;
in vec3 colors;
uniform  mat4 mvp_matrix;
uniform  vec4 cutplane;
out  vec4 color; 
uniform  float point_size;
void main(void)
{
  gl_PointSize = point_size;
  color = vec4(colors, vertex.x * cutplane.x  + vertex.y * cutplane.y  + vertex.z * cutplane.z  +  cutplane.w);
  gl_Position = mvp_matrix * vertex;
}
