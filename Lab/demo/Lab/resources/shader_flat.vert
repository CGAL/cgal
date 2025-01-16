#version 150

in vec4 vertex;
in vec3 normals;
in vec4 colors;
out VS_OUT
{
  vec4 fP;
  vec4 out_color;
  float dist[6];
  vec4 vertex;
}vs_out;

uniform mat4 mvp_matrix;
uniform mat4 mv_matrix;
uniform bool is_clipbox_on;
uniform highp mat4 clipbox1;
uniform highp mat4 clipbox2;
uniform highp float point_size;

void compute_distances(void)
{
  for(int i=0; i<3; ++i)
  {
    vs_out.dist[i]=
    clipbox1[i][0]*vertex.x+
    clipbox1[i][1]*vertex.y+
    clipbox1[i][2]*vertex.z +
    clipbox1[i][3];
    vs_out.dist[i+3]=
    clipbox2[i][0]*vertex.x+
    clipbox2[i][1]*vertex.y+
    clipbox2[i][2]*vertex.z +
    clipbox2[i][3];
  }
}

void main(void)
{
  gl_PointSize = point_size;
  if(is_clipbox_on)
    compute_distances();
   vs_out.out_color=colors;
   vs_out.fP = mv_matrix * vertex;
   vs_out.vertex = vertex;
   gl_Position = mvp_matrix * vertex;
}
