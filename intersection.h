#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "mesh.h"

bool edge_triangle_intersection(Edge* e, Vertex* v0, Vertex* v1, Vertex* v2, double* n);
bool edge_face_intersection(Edge* e, Face* f);
//bool inBox(Edge *e, Box* box); 
#endif
