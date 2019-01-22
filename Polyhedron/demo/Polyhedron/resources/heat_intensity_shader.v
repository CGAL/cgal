#version 150
in vec4 vertex;
in vec3 normals;
in vec4 colors;
in float distance;
uniform mat4 mvp_matrix;
uniform mat4 mv_matrix; 
out vec4 fP; 
out vec3 fN; 
out vec4 color;
out float dist[6];
out float out_distance;

uniform bool is_clipbox_on;
uniform mat4 clipbox1;
uniform mat4 clipbox2;
uniform float point_size;

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
   gl_PointSize = point_size;
   color = colors;
   //
   if(is_clipbox_on)
    compute_distances();
   fP = mv_matrix * vertex;
   mat3 mv_matrix_3;
   mv_matrix_3[0] = mv_matrix[0].xyz;
   mv_matrix_3[1] = mv_matrix[1].xyz;
   mv_matrix_3[2] = mv_matrix[2].xyz;
   fN = mv_matrix_3* normals; 
   out_distance = distance;
   gl_Position = mvp_matrix * vertex;
}
