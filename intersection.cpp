#include "intersection.h"
#include "mathutil.h" 

bool edge_triangle_intersection(Edge* e, Vertex* v0, Vertex* v1, Vertex* v2, double* n)
{
    //Vertex* v0 = f->v[0];
    //Vertex* v1 = f->v[1];
    //Vertex* v2 = f->v[2];
    //double* n;
    //n = f->n;
    double p[3];
    double d = -dot(n, v0->v);
    
    if(norm(e->dir) < EPSILON) {
        return false;
    }    
 
    if (dot(n, e->dir) < EPSILON) {
        return false;
    }

    double t = -(dot(n, e->v[0]->v) + d) / dot(n, e->dir);
    if (t < 0) {
        return false;
    } else if(!e->is_ray && t > dist(e->v[0]->v, e->v[1]->v)) {
        return false;
    }

    p[0] = e->v[0]->v[0] + t * e->dir[0];
    p[1] = e->v[0]->v[1] + t * e->dir[1];
    p[2] = e->v[0]->v[2] + t * e->dir[2];
   
    double e0[3] = {v1->v[0] - v0->v[0], v1->v[1] - v0->v[1], v1->v[2] - v0->v[2]}; //v1 - v0
    double e1[3] = {v2->v[0] - v1->v[0], v2->v[1] - v1->v[1], v2->v[2] - v1->v[2]}; //v2 - v1
    double e2[3] = {v0->v[0] - v2->v[0], v0->v[1] - v2->v[1], v0->v[2] - v2->v[2]}; //v0 - v2
    
    double vp0[3] = {p[0] - v0->v[0], p[1] - v0->v[1], p[2] - v0->v[2]}; //p - v0
    double vp1[3] = {p[0] - v1->v[0], p[1] - v1->v[1], p[2] - v1->v[2]}; //p - v1
    double vp2[3] = {p[0] - v2->v[0], p[1] - v2->v[1], p[2] - v2->v[2]}; //p - v2

    double u0[3] = {v2->v[0] - v0->v[0], v2->v[1] - v0->v[1], v2->v[2] - v0->v[2]}; //v2 - v0
        
    double c[3];
    
    cross_product(e0, u0, c);
    double denom = dot(n, c);
    
    cross_product(e1, vp1, c);
    double alpha = dot(n, c)/denom;
    
    cross_product(e2, vp2, c);
    double beta = dot(n, c)/denom;
    
    cross_product(e0, vp0, c);
    double gamma = dot(n, c)/denom;

    bool in_triangle = alpha + beta + gamma < 1.05 && alpha + beta + gamma > 0.95 && alpha >= 0 && beta >= 0 && gamma >= 0;
       
    if(in_triangle) {
        e->p[0] = p[0];
        e->p[1] = p[1];
        e->p[2] = p[2];
        e->intersected = true;
    }
    return in_triangle;
}

bool edge_face_intersection(Edge* e, Face* f) {
    for (int i = 1; i < f->v.size()-1; i++) {
        if (edge_triangle_intersection(e, f->v[0], f->v[i], f->v[i+1], f->n)) {
            e->intersected_face = f;
            return true;
    }
    }
    return false;
}

//bool inBox(Edge *e, Box* box) {
//    double d[3];
//    double x1, y1, z1, x2, y2, z2;
//    double left = box->bounds[0];
//    double right = box->bounds[1];
//    double front = box->bounds[2];
//    double back = box->bounds[3];
//    double bottom = box->bounds[4];
//    double top = box->bounds[5];
//
//    if (!e->is_ray) {
//            for (int i = 0; i++; i < 3) {
//                d[i] = e->v[0]->v[i] - e->v[1]->v[i];
//            }
//            x1 = d[0];
//            y1 = d[1];
//            z1 = d[2];
//        if (!(((e->v[0]->v[2] > top)&&(e->v[1]->v[2] > top))||((e->v[0]->v[2] < bottom)&&(e->v[1]->v[2] < bottom)))) {
//            // Check intersection with top side
//            z2 = top;
//            x2 = x1 / z1 * z2;
//            y2 = y1 / z1 * z2;
//            if ((left <= x2 <= right) && (front <= y2 <= back)) {
//                return true;
//            } 
//            // Check intersection with bottom side
//            z2 = bottom;
//            x2 = x1 / z1 * z2;
//            y2 = y1 / z1 * z2;
//            if ((left <= x2 <= right) && (front <= y2 <= back)) {
//                return true;
//            } 
//        }
//        if (!(((e->v[0]->v[1] > back)&&(e->v[1]->v[1] > back))||((e->v[0]->v[1] < front)&&(e->v[1]->v[1] < front)))) {
//            // Check intersection with left side
//            x2 = left;
//            y2 = y1 / x1 * x2;
//            z2 = z1 / x1 * x2;
//            if ((top <= z2 <= bottom) && (front <= y2 <= back)) {
//                return true;
//            } 
//            // Check intersection with right side
//            x2 = right;
//            y2 = y1 / x1 * x2;
//            z2 = z1 / x1 * x2;
//            if ((top <= z2 <= bottom) && (front <= y2 <= back)) {
//                return true;
//            } 
//        }
//        if (!(((e->v[0]->v[0] > right)&&(e->v[1]->v[0] > right))||((e->v[0]->v[0] < left)&&(e->v[1]->v[0] < left)))) {
//            // Check intersection with front side
//            y2 = front;
//            x2 = x1 / y1 * y2;
//            z2 = z1 / y1 * y2;
//            if ((left <= x2 <= right) && (top <= z2 <= bottom)) {
//                return true;
//            } 
//            // Check intersection with front side
//            y2 = back;
//            x2 = x1 / y1 * y2;
//            z2 = z1 / y1 * y2;
//            if ((left <= x2 <= right) && (top <= z2 <= bottom)) {
//                return true;
//            } 
//        }
//    } else {
//        x1 = e->dir[0];
//        y1 = e->dir[1];
//        z1 = e->dir[2];
//        if (!(((e->v[0]->v[2] > top)&&(e->dir[2] > 0))||((e->v[0]->v[2] < bottom)&&(e->dir[2] < 0)))) {
//            // Check intersection with top side
//            z2 = top;
//            x2 = x1 / z1 * z2;
//            y2 = y1 / z1 * z2;
//            if ((left <= x2 <= right) && (front <= y2 <= back)) {
//                return true;
//            } 
//            // Check intersection with bottom side
//            z2 = bottom;
//            x2 = x1 / z1 * z2;
//            y2 = y1 / z1 * z2;
//            if ((left <= x2 <= right) && (front <= y2 <= back)) {
//                return true;
//            } 
//        }
//        if (!(((e->v[0]->v[1] > back)&&(e->dir[1] > 0))||((e->v[0]->v[1] < front)&&(e->dir[1] < front)))) {
//            // Check intersection with left side
//            x2 = left;
//            y2 = y1 / x1 * x2;
//            z2 = z1 / x1 * x2;
//            if ((top <= z2 <= bottom) && (front <= y2 <= back)) {
//                return true;
//            } 
//            // Check intersection with right side
//            x2 = right;
//            y2 = y1 / x1 * x2;
//            z2 = z1 / x1 * x2;
//            if ((top <= z2 <= bottom) && (front <= y2 <= back)) {
//                return true;
//            } 
//        }
//        if (!(((e->v[0]->v[0] > right)&&(e->dir[0] > 0))||((e->v[0]->v[0] < left)&&(e->dir < 0)))) {
//            // Check intersection with front side
//            y2 = front;
//            x2 = x1 / y1 * y2;
//            z2 = z1 / y1 * y2;
//            if ((left <= x2 <= right) && (top <= z2 <= bottom)) {
//                return true;
//            } 
//            // Check intersection with front side
//            y2 = back;
//            x2 = x1 / y1 * y2;
//            z2 = z1 / y1 * y2;
//            if ((left <= x2 <= right) && (top <= z2 <= bottom)) {
//                return true;
//            } 
//        }
//        
//          return true;        
//    }
//    return false;
//}
