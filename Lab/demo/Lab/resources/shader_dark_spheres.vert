#version 150
in vec4 vertex;
in vec3 colors;
in vec3 center;
in float radius;
uniform mat4 mvp_matrix;
uniform mat4 f_matrix;
out vec4 color;
out float dist[6];

void main(void)
{
  for(int i=0; i<6; ++i)
   dist[i] = 1.0;
  color = vec4(colors, 1.0);
  gl_Position =  mvp_matrix * f_matrix *
   vec4(radius*vertex.x + center.x, radius* vertex.y + center.y, 
   radius*vertex.z + center.z, 1.0) ;
}
