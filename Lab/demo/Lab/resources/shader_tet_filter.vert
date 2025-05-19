#version 150

in vec4 vertex;
in vec3 normals;
in vec4 colors;
in float distance;
out VS_OUT
{
  vec4 fP;
  vec4 out_color;
  float dist[6];
  vec4 vertex;
  vec3 normal;
  flat float stat_value;
}vs_out;

uniform mat4 mvp_matrix;
uniform mat4 mv_matrix;
uniform mat4 norm_matrix;
uniform bool is_clipbox_on;
uniform highp mat4 clipbox1;
uniform highp mat4 clipbox2;

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
  if(is_clipbox_on)
    compute_distances();
   vs_out.out_color=colors;
   vs_out.fP = mv_matrix * vertex;
   vs_out.vertex = vertex;
   vs_out.stat_value= distance;

   mat3 norm_matrix_3;
   norm_matrix_3[0] = norm_matrix[0].xyz;
   norm_matrix_3[1] = norm_matrix[1].xyz;
   norm_matrix_3[2] = norm_matrix[2].xyz;
   vs_out.normal = norm_matrix_3* normals;

   gl_Position = mvp_matrix * vertex;
}
