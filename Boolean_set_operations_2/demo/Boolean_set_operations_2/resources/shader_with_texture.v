#version 120
attribute highp vec4 vertex;
attribute highp vec3 normal;
attribute highp vec3 color_facets;
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
       
void main(void) 
{ 
   vec4 P = mv_matrix * vertex; 
   vec3 N = mat3(mv_matrix)* normal; 
   vec3 L = light_pos.xyz - P.xyz; 
   N = normalize(N); 
   L = normalize(L); 
   vec3 diffuse;
   if(is_two_side == 1) 
       diffuse = abs(dot(N,L)) * light_diff.xyz; 
   else 
       diffuse = max(dot(N,L), 0.0) * light_diff.xyz; 
   f_texCoord = v_texCoord; 
   fColors = color_facets * (light_amb.xyz + diffuse);
   gl_Position =  mvp_matrix * f_matrix * vertex; 
}  
