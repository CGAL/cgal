//#version 100 
attribute highp vec4 vertex;
attribute highp vec3 normals;
attribute highp vec2 v_texCoord; 

uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix; 
uniform highp mat4 f_matrix; 
uniform highp int is_two_side; 
uniform highp vec4 light_pos;  
uniform highp vec4 light_diff; 
uniform highp vec3 light_spec; 
uniform highp vec4 light_amb;  
uniform highp float spec_power; 
varying highp vec3 fColors; 
varying highp vec2 f_texCoord; 
       
uniform highp float point_size;
void main(void)
{
  gl_PointSize = point_size;
  highp vec4 P = mv_matrix * vertex;
  mat3 mv_matrix_3;
  mv_matrix_3[0] = mv_matrix[0].xyz;
  mv_matrix_3[1] = mv_matrix[1].xyz;
  mv_matrix_3[2] = mv_matrix[2].xyz;
  highp vec3 N = mv_matrix_3* normals;
  highp vec3 L = light_pos.xyz - P.xyz;
  N = normalize(N);
  L = normalize(L);
  highp vec3 diffuse;
  if(is_two_side == 1)
      diffuse = abs(dot(N,L)) * light_diff.xyz;
  else
      diffuse = max(dot(N,L), 0.0) * light_diff.xyz;
  f_texCoord = v_texCoord;
  fColors = vec3(1.0, 1.0, 1.0) * (light_amb.xyz + diffuse);
  gl_Position =  mvp_matrix * f_matrix * vertex;
}  
