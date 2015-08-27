#version 120
varying highp vec4 color;
varying highp vec4 fP; 
varying highp vec3 fN; 
uniform highp vec4 light_pos;  
uniform highp vec4 light_diff; 
uniform highp vec4 light_spec; 
uniform highp vec4 light_amb;  
uniform highp float spec_power ; 
uniform int is_two_side; 
void main(void) { 
   highp vec3 L = light_pos.xyz - fP.xyz; 
   highp vec3 V = -fP.xyz; 
   highp vec3 N; 
   if(fN == highp vec3(0.0,0.0,0.0))
       N = highp vec3(0.0,0.0,0.0); 
   else 
       N = normalize(fN); 
   L = normalize(L); 
   V = normalize(V); 
   highp vec3 R = reflect(-L, N); 
  vec4 diffuse; 
   if(is_two_side == 1) 
       diffuse = abs(dot(N,L)) * light_diff * color; 
   else 
       diffuse = max(dot(N,L), 0.0) * light_diff * color; 
   highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; 
gl_FragColor = vec4((color*light_amb).xyz + diffuse.xyz + specular.xyz,1); 
}  
