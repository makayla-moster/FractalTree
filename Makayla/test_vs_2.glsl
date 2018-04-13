#version 410

layout (location = 0) in vec3 vertex_position;
layout (location = 1) in vec3 vertex_normal;
//uniform mat4 projection_mat, view_mat, model_mat;
uniform mat4 rotateX, scale;
out vec3 position_eye, normal_eye;

void main () {
    position_eye = vertex_position; //vec3 (view_mat * model_mat * vec4 (vertex_position, 1.0));
    normal_eye = vertex_normal; // vec3 (view_mat * model_mat * vec4 (vertex_normal, 0.0));
    gl_Position = rotateX * scale * vec4 (position_eye, 1.0);
}