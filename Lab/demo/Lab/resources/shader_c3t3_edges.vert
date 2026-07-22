#version 150
in vec4 vertex;
in vec3 colors;
in vec2 subdomain_in;
uniform  mat4 mvp_matrix;
uniform  vec4 cutplane;
out  vec4 color; 
out float dist[6];
uniform  float point_size;
flat out vec2 subdomain_out;
uniform bool is_clipbox_on;
uniform mat4 clipbox1;
uniform mat4 clipbox2;

void compute_distances(void)
{
  for(int i=0; i<3; ++i)
  {
    dist[i]=
    clipbox1[i][0]*vertex.x+
    clipbox1[i][1]*vertex.y+
    clipbox1[i][2]*vertex.z +
    clipbox1[i][3];
    dist[i+3]=
    clipbox2[i][0]*vertex.x+
    clipbox2[i][1]*vertex.y+
    clipbox2[i][2]*vertex.z +
    clipbox2[i][3];
  }
}

void main(void)
{
  subdomain_out = subdomain_in;
  gl_PointSize = point_size;
  color = vec4(colors, vertex.x * cutplane.x  + vertex.y * cutplane.y  + vertex.z * cutplane.z  +  cutplane.w);
  if(is_clipbox_on)
   compute_distances();
  gl_Position = mvp_matrix * vertex;
}
