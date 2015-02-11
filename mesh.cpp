#include "mesh.h"
#include <vector>

const int Max_BSP_Depth = 15;
const int Min_Num_Triangles = 1000;//30;
const int Dimension = 3;

void Mesh::add(Vertex* v)
{
    v->index = vertices.size();
    vertices.push_back(v);
}

void Mesh::add(Edge* e)
{
    e->index = edges.size();
    edges.push_back(e);
}

void Mesh::add(Face* f)
{
    f->index = faces.size();
    faces.push_back(f);
}

Vertex* Mesh::get_vertex(int i)
{
    return vertices[i];
}

Edge* Mesh::get_edge(int i)
{
    return edges[i];
}

Face* Mesh::get_face(int i)
{
    return faces[i];
}

size_t Mesh::vertex_count() const
{
    return vertices.size();
}

size_t Mesh::edge_count() const
{
    return edges.size();
}

size_t Mesh::face_count() const
{
    return faces.size();
}

Edge* get_edge(Vertex* v0, Vertex* v1) 
{
    for(int i = 0; i < v0->f.size(); i++) {
        Face* fi = v0->f[i];
        for(int j = 0; j < fi->e.size(); j++) {
            Edge* ej = fi->e[j];
            if((ej->v[0] == v0 && ej->v[1] == v1) || (ej->v[0] == v1 && ej->v[1] == v0)) {
                return ej;
            }
        }
    }
    return NULL;
}

Edge* create_edge(Vertex* v0, Vertex* v1, Face* f, Mesh& mesh)
{                                              
    Edge* e = get_edge(v0,v1);
    if(e == NULL) {
        e = new Edge;
        e->v[0] = v0;
        e->v[1] = v1;
        e->f[0] = f;
        e->intersected = false;
        mesh.add(e);
    } else {
        e->f[1] = f;
    }
    return e;
} 

void normalize(Mesh& mesh)
{
    for(int i = 0; i < mesh.vertex_count(); i++) {
        Vertex* v = mesh.get_vertex(i);
        normalize(v->n);
    }
}

void compute_normals(Mesh& mesh)
{
    for(int i = 0; i < mesh.face_count(); i++) {
        Face* f = mesh.get_face(i);
        Vertex *v0 = f->v[0], *v1 = f->v[1], *v2 = f->v[f->v.size()-1];
        //Need at least 3 points to define a normal. Assume face is planar.
        double vec1[3] = {0,0,0};
        double vec2[3] = {0,0,0};
        for(int j = 0; j < 3; j++) {
            vec1[j] = v1->v[j] - v0->v[j];
            vec2[j] = v2->v[j] - v0->v[j];
        }
        cross_product(vec1, vec2, f->n);
        for(int j = 0; j < f->v.size(); j++) {
            Vertex* vj = f->v[j];
            vj->n[0] += f->n[0];
            vj->n[1] += f->n[1];
            vj->n[2] += f->n[2];
        }
        normalize(f->n);
    }
    normalize(mesh);
}

void center_on_screen(Mesh& mesh)
{
    double maxp[3], minp[3];
    if(mesh.vertex_count() == 0) {
        return;
    }
    Vertex* v = mesh.get_vertex(0);
    maxp[0] = minp[0] = v->v[0];
    maxp[1] = minp[1] = v->v[1];
    maxp[2] = minp[2] = v->v[2];
    for(int i = 1; i < mesh.vertex_count(); i++) {
        v = mesh.get_vertex(i);
        maxp[0] = maxp[0] < v->v[0] ? v->v[0] : maxp[0];
        maxp[1] = maxp[1] < v->v[1] ? v->v[1] : maxp[1];
        maxp[2] = maxp[2] < v->v[2] ? v->v[2] : maxp[2];
        minp[0] = minp[0] > v->v[0] ? v->v[0] : minp[0];
        minp[1] = minp[1] > v->v[1] ? v->v[1] : minp[1];
        minp[2] = minp[2] > v->v[2] ? v->v[2] : minp[2];
    }
    double center[3] = {(maxp[0] + minp[0])/2, (maxp[1] + minp[1])/2, (maxp[2] + minp[2])/2};
    double scale = 2.0/max(maxp[0] - minp[0], max(maxp[1] - minp[1], maxp[2] - minp[2]));
    for(int i = 0; i < mesh.vertex_count(); i++) {
        v = mesh.get_vertex(i);
        v->v[0] -= center[0];
        v->v[1] -= center[1];
        v->v[2] -= center[2];
        v->v[0] *= scale;
        v->v[1] *= scale;
        v->v[2] *= scale;
    }
}


void Box::extend(Face* f)
{
    for (int i = 0; i < f->v.size(); i++){
        for (int j = 0; j < 3; j++){
            if (!initialized || (f->v[i]->v[j] < bounds[j * 2]))
          bounds[j * 2] = f->v[i]->v[j];
            if (!initialized || (f->v[i]->v[j] > bounds[j * 2 + 1]))
          bounds[j * 2 + 1] = f->v[i]->v[j];
        }
        initialized = true;
    }
}

double get_midpoint(Face* f, int axis){
    double ans = 0.0;
    size_t size = f->v.size();
    for (size_t i = 0; i < size; i++){
        ans += (f->v[i]->v[axis]) / size;
    }
    return ans;
}

BSP_tree* build_BSP(vector<Face*> BSP_faces, int depth){
    BSP_tree* tree;
    tree->box = new Box;
    tree->l_child = NULL;
    tree->r_child = NULL;
    //tree->axis = depth % Dimension;

    ////initializing
   
    double midpoint;
    //size_t faces_size = BSP_faces.size();
    //if (size>Min_Num_Triangles){ //Terminating point
    //    for (int i = 0; i < size; i++){
    //        Face* current = faces[i];
    //        midpoint += (get_midpoint(current, tree->axis))/size;
    //        tree->box->extend(current);
    //    }// Spliting plane
    //    tree->split = midpoint;
    //    for (int i = 0; i < size; i++){
    //        Face* current = faces[i];
    //        if (get_midpoint(current,tree->axis) <= midpoint){
    //            tree->left.push_back(current);
    //        } else {
    //            tree->right.push_back(current);
    //        } 
    //    }
    //    tree->l_child = build_BSP(tree->left, depth+1);
    //    tree->r_child = build_BSP(tree->right, depth+1);
    //} else {
    //    printf("build the leaf at depth %d", depth);
    //    for (int i = 0; i < size; i++){
    //        tree->box->extend(faces[i]);
    //    }
    //    tree->left = faces;
    //    tree->right = faces;
    //}
    return tree;
}

