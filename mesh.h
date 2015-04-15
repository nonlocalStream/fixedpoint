#ifndef MESH_H
#define MESH_H

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include "mathutil.h"

//const int Dimension = 3;
using std::vector;
using std::max;

struct Face;

struct Vertex
{
    double v[3];
    double n[3];
    double t[2];
    bool has_normal, has_texture;
    int added; /* whether the vertex is added to the next mesh */
    size_t index;
    vector<Face*> f;
};

struct Edge
{
    Vertex* v[2];
    Face* f[2];
    Face* dual;
    double p[3];
    double dir[3];
    size_t index;
    bool intersected;
    bool is_ray;
};

struct Face
{
    vector<Vertex*> v;
    vector<Edge*> e;
    Edge* dual;
    double n[3];
    size_t index;
    bool visited;
};


class Mesh
{
public:
    void add(Vertex* v);
    void add(Edge* e);
    void add(Face* f);

    Vertex* get_vertex(int i);
    Edge* get_edge(int i);
    Face* get_face(int i);
    vector<Face*> get_faces();

    size_t vertex_count() const;
    size_t edge_count() const;
    size_t face_count() const;

private:
    vector<Vertex*> vertices;
    vector<Edge*> edges;
    vector<Face*> faces;
};

class Box
{

public:
    bool initialized;
    double bounds[6];//up, down, left, right, front, back; //not generalized enough consider: bounds[Dimension*2]
    void extend(Face* f);
    Box(){
      initialized = false;
    }
//private:
//    bool initialized;
};

class BSP_tree
{
public:
    Box* box;
    vector<Face*> left, right;
    double split;
    int axis;
    BSP_tree* l_child;
    BSP_tree* r_child;
};

Edge* get_edge(Vertex* v0, Vertex* v1);
Edge* create_edge(Vertex* v0, Vertex* v1, Face* f, Mesh& mesh);

void normalize(Mesh& mesh);
void compute_normals(Mesh& mesh);
void center_on_screen(Mesh& mesh);

//BSP_tree* build_BSP(vector<Face*> faces, size_t size, int depth);
double get_midpoint(Face* f, int axis);
void adjust_normals(Mesh& mesh);
void revert_normals(Mesh& mesh);
#endif
