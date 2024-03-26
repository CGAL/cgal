attribute highp vec4 vertex;
attribute highp vec3 colors;
attribute highp vec3 center;
attribute highp float radius;
uniform highp mat4 mvp_matrix;
uniform highp mat4 f_matrix;
varying highp vec4 color;
varying highp float dist[6];

void main(void)
{
  for(int i=0; i<6; ++i)
   dist[i] = 1.0;
  color = vec4(colors, 1.0);
  gl_Position =  mvp_matrix * f_matrix *
  vec4(radius*vertex.x + center.x, radius* vertex.y + center.y, radius*vertex.z + center.z, 1.0) ;
}
