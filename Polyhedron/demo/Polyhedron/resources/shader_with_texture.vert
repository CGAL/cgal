#version 150
in vec4 vertex;
in vec3 normals;
in vec2 v_texCoord; 

uniform mat4 mvp_matrix;
uniform mat4 mv_matrix; 
uniform mat4 f_matrix; 
uniform int is_two_side; 
uniform vec4 light_pos;  
uniform vec4 light_diff; 
uniform vec3 light_spec; 
uniform vec4 light_amb;  
uniform float spec_power; 
out vec3 fColors; 
out vec2 f_texCoord; 
       
uniform float point_size;
void main(void)
{
  gl_PointSize = point_size;
   vec4 P = mv_matrix * vertex; 
   mat3 mv_matrix_3;                    
   mv_matrix_3[0] = mv_matrix[0].xyz;   
   mv_matrix_3[1] = mv_matrix[1].xyz;   
   mv_matrix_3[2] = mv_matrix[2].xyz;   
   vec3 N = mv_matrix_3* normals; 
   vec3 L = light_pos.xyz - P.xyz; 
   N = normalize(N); 
   L = normalize(L); 
   vec3 diffuse;
   if(is_two_side == 1) 
       diffuse = abs(dot(N,L)) * light_diff.xyz; 
   else 
       diffuse = max(dot(N,L), 0.0) * light_diff.xyz; 
   f_texCoord = v_texCoord; 
   fColors = vec3(1.0, 1.0, 1.0) * (light_amb.xyz + diffuse);
   gl_Position =  mvp_matrix * f_matrix * vertex; 
}  
