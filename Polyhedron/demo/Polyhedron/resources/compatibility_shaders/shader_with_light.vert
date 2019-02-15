//#version 100 
attribute highp vec4 vertex;
attribute highp vec3 normals;
attribute highp vec4 colors;
uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix; 
uniform highp mat4 f_matrix; 
varying highp vec4 fP; 
varying highp vec3 fN; 
varying highp vec4 color;
varying highp float dist[6];
uniform bool is_clipbox_on;
uniform highp mat4 clipbox1;
uniform highp mat4 clipbox2;
uniform highp float point_size;

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
   gl_Position = mvp_matrix * f_matrix * vertex;
}
