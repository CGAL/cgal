#version 150
in vec4 vertex;
in vec3 colors;
in vec2 subdomain_in;
uniform  mat4 mvp_matrix;
uniform  vec4 cutplane;
out  vec4 color; 
uniform  float point_size;
flat out vec2 subdomain_out;
void main(void)
{
  subdomain_out = subdomain_in;
  gl_PointSize = point_size;
  color = vec4(colors, vertex.x * cutplane.x  + vertex.y * cutplane.y  + vertex.z * cutplane.z  +  cutplane.w);
  gl_Position = mvp_matrix * vertex;
}
