#include <iostream>
#include <cstdlib>
#include <algorithm>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "io.h"
#include "mesh.h"
#include "mathutil.h"
#include "glcontext.h"
#include "intersection.h"
#ifndef __PYRAMID__
#define __PYRAMID__
#include "pyramid.c"
#include <time.h>
#include <stdlib.h>
#endif

const int Max_BSP_Depth = 15;
const int Min_Num_Triangles = 3;
const int Dimension = 3;

using namespace std;


Context ctx;
Mesh mesh;
Mesh vor;
Mesh del;
vector<Mesh> meshes;
int curr_mesh = 0;
int countNum = 0;

GLUquadricObj *quadric = gluNewQuadric();
bool draw_intersection = false;
bool draw_disappear_and_appear = false;
bool draw_voronoi = false;

void clear_intersection()
{
    for(size_t i = 0; i < vor.edge_count(); i++) {
        Edge* e = vor.get_edge(i);
        e->intersected = false;
    }
}

void print_intersection(FILE* f)
{
    for(size_t i = 0; i < vor.edge_count(); i++) {
        Edge* e = vor.get_edge(i);
        if (e->intersected){
            fprintf(f, "%d th is TRUE: %f, %f, %f\n", i, e->p[0], e->p[1], e->p[2]);
        } else {
            fprintf(f, "%d th is FALSE\n", i);
        }
    }
}

int count_triangles(BSP_tree* bsp)
{
    if ((bsp->l_child == NULL)&&(bsp->r_child == NULL)) {
        return bsp->left.size();
    } else {
        return count_triangles(bsp->l_child)+count_triangles(bsp->r_child);
    }
}

void create_mesh()
{
    Mesh next;
   for(size_t i = 0; i < vor.edge_count(); i++) {
        Edge* e = vor.get_edge(i);
        next.add(e);
        if(e->intersected) {
            Face* f = e->dual;
            next.add(f);
            for(size_t j = 0; j < f->v.size(); j++) {
                Vertex* v = f->v[j];
                //if (!v->added){
                //    v->added = true;
                //    next.add(v);
                //}
            }
        }
    }
   for(size_t j = 0; j < del.vertex_count();j++) {
        Vertex *v = del.get_vertex(j);
        v->added = 0;
    }
     for(size_t j = 0; j < del.vertex_count();j++) {
        Vertex *v = del.get_vertex(j);
        v->added++;
        next.add(v);
    }
    compute_normals(next);
    meshes.push_back(next);
    curr_mesh = meshes.size()-1; 
    mesh = next;
}

void intersection()
{
    for(size_t i = 0; i < vor.edge_count(); i++) {
        Edge* e = vor.get_edge(i);
        for(size_t j = 0; j < meshes[curr_mesh].face_count(); j++) {
            Face* f = meshes[curr_mesh].get_face(j);
            if(!e->intersected) {
                bool hit = edge_face_intersection(e,f);
                if(!hit) {
                    f->n[0] *= -1, f->n[1] *= -1, f->n[2] *= -1;
                    edge_face_intersection(e,f);
                    f->n[0] *= -1, f->n[1] *= -1, f->n[2] *= -1;
                }
            }
        }
    }
}

bool edge_side_intersection(vertex x,
                            vertex w,
                            vertex y,
                            vertex z,
                            vertex a,
                            vertex b)
{
    // this ensures a and b are on opposite sides of the rectangle x,w,y,z (counterclockwise) or one of themis on the rectangle
    if (orient3d(x,w,y,a)*orient3d(x,w,y,b) > 0) {
      return false;
    }
    /* Orient3D(x,w,a,b), Orient3D(w,y,a,b)
     * Orient3D(y,z,a,b), Orient3D(z,x,a,b)
     * all have the same signs
     * which guarantee that line a-b intersect with the rectangle
     */
    if (orient3d(x,w,a,b)*orient3d(w,y,a,b) < 0) {
      return false;
    }
    if (orient3d(w,y,a,b)*orient3d(y,z,a,b) < 0) {
      return false;
    }
    if (orient3d(y,z,a,b)*orient3d(z,x,a,b) < 0) {
      return false;
    }
    return true;
}
bool ray_side_intersection(vertex x,
                            vertex w,
                            vertex y,
                            vertex z,
                            vertex a,
                            vertex b)
{
    // this ensures a and a+(infinity) * b are on opposite sides of the rectangle x,w,y,z (counterclockwise) or one of triangle on the rectangle
    if (edge_side_intersection(x,w,y,z,a,b))
        return true;


    //if (abs(orient3d(a,x,w,y)) > 0 &&// >EPSILON  &&
    //   (abs(orient3d(a,w,y,z))<=abs(orient3d(b,w,y,z)))){//+EPSILON)) {
    //  return false;
    //}
    /* Orient3D(x,w,a,b), Orient3D(w,y,a,b)
     * Orient3D(y,z,a,b), Orient3D(z,x,a,b)
     * all have the same signs
     * which guarantee that line a-b intersect with the rectangle
     */
    if (orient3d(x,w,a,b)*orient3d(w,y,a,b) < 0) {
      return false;
    }
    if (orient3d(w,y,a,b)*orient3d(y,z,a,b) < 0) {
      return false;
    }
    if (orient3d(y,z,a,b)*orient3d(z,x,a,b) < 0) {
      return false;
    }
    return true;
}

bool inBox(Edge *e, Box* box) {
  double* bounds;
  bounds = box->bounds;
  REAL v000[3]= {bounds[0], bounds[2], bounds[4]};
  REAL v001[3]= {bounds[0], bounds[2], bounds[5]};
  REAL v010[3]= {bounds[0], bounds[3], bounds[4]};
  REAL v011[3]= {bounds[0], bounds[3], bounds[5]};
  REAL v100[3]= {bounds[1], bounds[2], bounds[4]};
  REAL v101[3]= {bounds[1], bounds[2], bounds[5]};
  REAL v110[3]= {bounds[1], bounds[3], bounds[4]};
  REAL v111[3]= {bounds[1], bounds[3], bounds[5]};
  if (e->is_ray) {
      REAL e_v0[3];
      REAL e_v1[3];
      for (int i = 0; i < 3; i++) {
        e_v0[i] = e->v[0]->v[i];
        e_v1[i] = e->v[0]->v[i] + e->dir[i];
      }
      int inside = 1;
      for (int i=0; i<3; i++){
          if ((e_v0[i] < bounds[i*2]) || (e_v0[i] > bounds[i*2+1])) inside = 0;
          if ((e_v1[i] < bounds[i*2]) || (e_v1[i] > bounds[i*2+1])) inside = 0;
      }
      
      if (inside) return true;
 
      //check if intersect with 1 side of the box
      bool result = ray_side_intersection(v001, v101, v111, v011, e_v0, e_v1) || //up
             ray_side_intersection(v000, v010, v011, v001, e_v0, e_v1) || //left
             ray_side_intersection(v100, v110, v111, v101, e_v0, e_v1) || //right
             ray_side_intersection(v000, v100, v101, v001, e_v0, e_v1) || //front
             ray_side_intersection(v010, v110, v111, v011, e_v0, e_v1); //back
      return result;
  } else {
      REAL e_v0[3];
      REAL e_v1[3];
      for (int i = 0; i < 3; i++) {
        e_v0[i] = e->v[0]->v[i];
        e_v1[i] = e->v[1]->v[i];
      }

      //check if completely inside the box
      int inside = 1;
      for (int i=0; i<3; i++){
          if ((e_v0[i] < bounds[i*2]) || (e_v0[i] > bounds[i*2+1])) inside = 0;
          if ((e_v1[i] < bounds[i*2]) || (e_v1[i] > bounds[i*2+1])) inside = 0;
      }
      
      if (inside) return true;
      
      //check if intersect with 1 side of the box
      return edge_side_intersection(v001, v101, v111, v011, e_v0, e_v1) || //up
             edge_side_intersection(v000, v100, v110, v010, e_v0, e_v1) || //down
             edge_side_intersection(v000, v010, v011, v001, e_v0, e_v1) || //left
             edge_side_intersection(v100, v110, v111, v101, e_v0, e_v1) || //right
             edge_side_intersection(v000, v100, v101, v001, e_v0, e_v1) || //front
             edge_side_intersection(v010, v110, v111, v011, e_v0, e_v1); //back
  }
  return true;
}

double quickselect(double * a, 
                   int l,
                   int r,
                   int k) 
{
  //For (int i = l; i<=r; i++) {
  //    printf("%f, ", a[i]);
  //}
  //printf("\n");
  if (l < r) {
    int pivot_index = rand() % (r-l+1) + l;
    double pivot = a[pivot_index];
    double temp;
    a[pivot_index] = a[r];
    a[r] = pivot;
    int i = l - 1;
    int j = r;
    do {
        do { i++; } while (a[i] < pivot);
        do { j--; } while ((pivot < a[j]) && (j > l));
        if (i < j) {
            temp = a[i];
            a[i] = a[j];
            a[j] = temp;
        }
    } while (i < j);
    a[r] = a[i];
    a[i] = pivot;
    if (i == k) {
        return a[i];
    } else if (k < i) {
        return quickselect(a, l, i - 1, k);
    } else {
        return quickselect(a, i + 1, r, k);
    }
  }
  else {
    if (l == k) return a[l]; else return a[r];
  }
}

BSP_tree* build_BSP(const vector<Face*> & faces, size_t size, int depth) {
    srand(time(NULL));
    BSP_tree* tree = new BSP_tree;
    tree->box = new Box;
    tree->l_child = NULL;
    tree->r_child = NULL;
    tree->axis = depth % Dimension;

    //////initializing
   
    double midpoint = 0;
    double a[size * 3];

    if ((size>Min_Num_Triangles)&&(depth < Max_BSP_Depth)){ //Terminating point
        for (int i = 0; i < size; i++){
            Face* current = faces[i];
            //midpoint += (get_midpoint(current, tree->axis))/(double)size;
           
            a[i*3] = current->v[0]->v[tree->axis];
            a[i*3+1] = current->v[1]->v[tree->axis];
            a[i*3+2] = current->v[2]->v[tree->axis];
            tree->box->extend(current);
        }// Spliting plane
        //for (int i = 0; i < 60; i++) {
        //    printf("a[%d] = %f, ", i, a[i]);
        //}
        midpoint = quickselect(a, 0, 3*size-1, 3*size/2);
        //printf("\n midpoint is %f\n", midpoint);
        tree->split = midpoint;
        size_t left_size = 0;
        size_t right_size = 0;

        for (int i = 0; i < size; i++){
            Face* current = faces[i];
            if (get_midpoint(current,tree->axis) <= midpoint){
                tree->left.push_back(current);
                left_size++;
            } else {
                tree->right.push_back(current);
                right_size++;
            } 
        }
        //printf("spliting at depth %d with %d left tris and %d right tris with threshold %f\n", depth, left_size, right_size, tree->split);
        tree->l_child = build_BSP(tree->left, left_size, depth+1);
        tree->r_child = build_BSP(tree->right, right_size, depth+1);
    } else {
        //printf("build the leaf at depth %d\n", depth);
        for (int i = 0; i < size; i++){
            tree->box->extend(faces[i]);
        }
        tree->left = faces;
        tree->right = faces;
    }
    return tree;
}

void intersection_with_BSP_helper(Edge* e, BSP_tree* bsp, int depth){
    double * bounds;
    if ((bsp->l_child == NULL)&&(bsp->r_child == NULL)) {
        //printf("reach the leaf!\n");
        for (int i = 0; i < bsp->left.size(); i++) {
            Face* f = bsp->left[i];
            if(!e->intersected) {
                bool hit = edge_face_intersection(e,f);
                if(!hit) {
                    f->n[0] *= -1, f->n[1] *= -1, f->n[2] *= -1;
                    edge_face_intersection(e,f);
                    f->n[0] *= -1, f->n[1] *= -1, f->n[2] *= -1;
                }
            }
        }
    } else {
    bounds = bsp->l_child->box->bounds;
    //printf("left reaching depth %d---1st: %f~%f; 2nd: %f~%f; 3rd: %f~%f\n", depth, bounds[0], bounds[1],bounds[2], bounds[3],bounds[4], bounds[5]);
    bounds = bsp->r_child->box->bounds;
    //printf("reaching depth %d---1st: %f~%f; 2nd: %f~%f; 3rd: %f~%f\n", depth, bounds[0], bounds[1],bounds[2], bounds[3],bounds[4], bounds[5]);
        if (inBox(e, bsp->l_child->box)){
          //printf("intersect with left!\n");
          intersection_with_BSP_helper(e, bsp->l_child, depth+1);
        }
        if (inBox(e, bsp->r_child->box)){
          //printf("intersect with right!\n");
          intersection_with_BSP_helper(e, bsp->r_child, depth+1);
        }
    }
}


void intersection_with_BSP()
{
    countNum++;
    std::vector<Face*> triangles;
    size_t size = 0;
    //cout << "face_count: " << meshes[curr_mesh].face_count() << endl;
    for(size_t j = 0; j < meshes[curr_mesh].face_count(); j++) {
        Face* f = meshes[curr_mesh].get_face(j);
        triangles.push_back(f);
        size++;
    }
    //cout << "sizes: " << size << endl;
    BSP_tree* bsp = build_BSP(triangles, size, 0);
    cout << "finish building BSP!" << endl; 
    for(size_t i = 0; i < vor.edge_count(); i++) {
        Edge* e = vor.get_edge(i);
    //printf("edge (%f,%f,%f)-(%f,%f,%f)\n", e->v[0]->v[0], e->v[0]->v[1], e->v[0]->v[2],e->v[1]->v[0], e->v[1]->v[1], e->v[1]->v[2]);
        intersection_with_BSP_helper(e, bsp, 0);
    }
}


void draw_point(double* v, double r, double g, double b)
{
    glColor3f(r,g,b);
    glPushMatrix();
    glTranslatef(v[0], v[1], v[2]);
    gluSphere(quadric, 0.02, 10, 10);
    glPopMatrix();
    glColor3f(1,1,1);
}

void display()
{
    float mat_diffuse[4] = {.5,.5,.5,0.5};
    float mat_ambient[4] = {0,0,0,1};
    float mat_specular[4] = {0,0,0,1};
    float mat_shininess = 1;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, ctx.wireframe ? GL_LINE : GL_FILL);
    
    //glEnable(GL_LIGHT0);
    //glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ctx.global_ambient);
    //glLightfv1GL_LIGHT0, GL_POSITION, ctx.light_pos);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &mat_shininess);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);
    //glColorMaterial(GL_FRONT_AND_BACK, GL_EMISSION);
    glColorMaterial(GL_FRONT_AND_BACK, GL_SPECULAR);
    glColor3f(1, 1, 1);

    glPushMatrix();
    glLoadIdentity();
    glTranslatef(0, 0, -4);
    glMultMatrixf(ctx.mat);
    
    float delta, r,g, b;
    for(int i = 0; i < meshes[curr_mesh].face_count(); i++) {
        Face* f = meshes[curr_mesh].get_face(i);
        glPolygonMode(GL_FRONT_AND_BACK, ctx.wireframe ? GL_LINE : GL_FILL);
        //glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        glBegin(GL_POLYGON);
        if (curr_mesh > 0) {
            if (f->fix_level > 0) {
                delta = (float)(f->fix_level)/5;
                b = (1.0+delta)/3;
                g = (1.0-delta/2)/3;
                r = (1.0-delta/2)/3;
                glColor3f(r,g,b);
            } else {
                glColor3f(1, 1, 1);
            }
        }

        for(int j = 0; j < f->v.size(); j++) {
            Vertex* v = f->v[j];
            glNormal3dv(v->n);
            glTexCoord2dv(v->t);
            glVertex3dv(v->v);  
        } 
        glEnd();
    } 

    //Draw Voronoi edges
    if(draw_voronoi) {
        glDisable(GL_LIGHTING);
        for(size_t i = 0; i < vor.edge_count(); i++) {
            Edge* e = vor.get_edge(i);
            //if (e->intersected) {
            //    glColor3f(1,0,0);
            //} else {
            //    glColor3f(1,1,1);
            //}
            glBegin(GL_LINES);
            glVertex3dv(e->v[0]->v);
            glVertex3dv(e->v[1]->v);
            glEnd();    
        }
        glEnable(GL_LIGHTING);
    }

    if(draw_intersection) {
        //Draw intersection points
        //mat_diffuse[0] = 1, mat_diffuse[1] = 0;
        //glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
        for(size_t i = 0; i < vor.edge_count(); i++) {
            Edge* e = vor.get_edge(i);
            if(e->intersected) {
                draw_point(e->p,1,0,0);
            }
        //mat_diffuse[0] = 1, mat_diffuse[1] = 0;
        //glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
        }
    }

    if((curr_mesh > 0) && draw_disappear_and_appear){
        glColor3f(1,0,0);
        for(size_t i = 0; i < meshes[curr_mesh - 1].face_count(); i++) {
            Face* f = meshes[curr_mesh - 1].get_face(i);
            if (f->to_disappear){
                glPolygonMode(GL_FRONT_AND_BACK, ctx.wireframe ? GL_LINE : GL_FILL);
                //glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
                mat_diffuse[0] = 1, mat_diffuse[1] = 0;
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
                glBegin(GL_POLYGON);
                for(int j = 0; j < f->v.size(); j++) {
                    Vertex* v = f->v[j];
                    glNormal3dv(v->n);
                    glTexCoord2dv(v->t);
                    glVertex3dv(v->v);  
                } 
                glEnd();
                mat_diffuse[0] = 1, mat_diffuse[1] = 0;
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
            }
        }
        glColor3f(0,1,0);
        for(size_t i = 0; i < meshes[curr_mesh].face_count(); i++) {
            Face* f = meshes[curr_mesh].get_face(i);
            if (f->newly_appear){
                glDisable(GL_LIGHTING);
                glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
                glBegin(GL_POLYGON);
                for(int j = 0; j < f->v.size(); j++) {
                    Vertex* v = f->v[j];
                    glNormal3dv(v->n);
                    glTexCoord2dv(v->t);
                    glVertex3dv(v->v);  
                } 
                glEnd();
                glEnable(GL_LIGHTING);
            }
        }
        glColor3f(1,1,1);
    }

//    if(draw_BSP) {
//        Box* b;
//        double left,up,front = 100,100,100;
//        double right,down,back = 300,300,300;
//        b 
//        for(int i = 0; i < meshes[curr_mesh].face_count(); i++) {
//            Face* f = meshes[curr_mesh].get_face(i);
//            glBegin(GL_POLYGON);
//            for(int j = 0; j < f->v.size(); j++) {
//                Vertex* v = f->v[j];
//                glNormal3dv(v->n);
//                glTexCoord2dv(v->t);
//                glVertex3dv(v->v);  
//            } 
//            glEnd();
//    } 


    glPopMatrix();
    glutSwapBuffers();
}

void reshape(int w, int h)
{
    ctx.window[0] = -static_cast<double>(w)/h;
    ctx.window[1] = -ctx.window[0];
    ctx.window[2] = 1.0;
    ctx.window[3] = -1.0;

    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(ctx.fovy, static_cast<double>(w)/h, ctx.znear, ctx.zfar);
    glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y)
{
    FILE *f;
    double a[10] = {10,60,20,40,80,30,90,100,50,70};
    switch(key) {
        /*case 'i':
            cout << "normal intersection" << endl;
            intersection();
        break;
        */
        case 'i':
            cout << "BSP intersection" << endl;
            clear_intersection();
            intersection_with_BSP();
            f = fopen("BSPintsec","w+");
            print_intersection(f);
            fclose(f);
        break;
        case 'p':
            draw_intersection = !draw_intersection;
        break;
        case 'd':
            draw_disappear_and_appear = !draw_disappear_and_appear;
        break;
        case 'v':
            draw_voronoi = !draw_voronoi;
        break;
        case 'r':
            revert_normals(meshes[curr_mesh]);
        break;
        case 'a':
            adjust_normals(meshes[curr_mesh]);
        break;
        case 'b':
            if (curr_mesh > 0){
                curr_mesh--;
                mesh = meshes[curr_mesh];
                update_color_var(meshes[curr_mesh], meshes[curr_mesh-1], curr_mesh == 0);
                clear_intersection();
                intersection_with_BSP();
                cout << "Face count: " <<  meshes[curr_mesh].face_count() << endl;
            }
            cout << "Current Frame:" << curr_mesh <<endl;
        break;  
        case 'n':
            clear_intersection();
            intersection_with_BSP();
            curr_mesh++;
            if (curr_mesh >= meshes.size()){
                create_mesh();
                adjust_normals(meshes[curr_mesh]);
            }
               
            update_color_var(meshes[curr_mesh], meshes[curr_mesh-1], curr_mesh == 0);
            cout << "Face count: " <<  meshes[curr_mesh].face_count() << endl;
            cout << "Current Frame:" << curr_mesh <<endl;
        break;  
        case 'w':
            ctx.wireframe = !ctx.wireframe;
        break;
        case 27: //ESC key
            exit(0);
        break;
        default:
            cout << "Unassigned character: " << key << endl;
        break;
    }
}

void menu()
{
    cout << "\033[1;31m   Fixed Point Menu\033[0m" << endl;
    cout << "\033[1;31mOption                                    Key\033[0m" << endl;
    cout << "Compute Intersections                      i" << endl;
    cout << "Compute Intersections with BSP             s" << endl;
    cout << "Draw Intersections                         p" << endl;
    cout << "Draw Voronoi Edges                         v" << endl;
    cout << "Iterate RDT Operator: next                 n" << endl;
    cout << "Iterate RDT Operator: back                 b" << endl;
    cout << "Toggle Wireframe                           w" << endl;
    cout << "Revert all normals                         r" << endl;
    cout << "Adjust normals                             a" << endl;
    cout << "Quit                                       ESC" << endl;
    cout << "Quit                                       ESC" << endl;
}

void world_coords(int x, int y, double* p)
{
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    p[0] = static_cast<double>(x-viewport[0])/viewport[2];
    p[1] = static_cast<double>(y-viewport[1])/viewport[3];

    p[0] = ctx.window[0] + p[0]*(ctx.window[1]-ctx.window[0]);
    p[1] = ctx.window[2] + p[1]*(ctx.window[3]-ctx.window[2]);
    p[2] = ctx.znear;
}

void mouse(int button, int state, int x, int y)
{
    int cursor = GLUT_CURSOR_RIGHT_ARROW;
    if(state == GLUT_DOWN) {
        if(button == GLUT_LEFT_BUTTON) {
            cursor = GLUT_CURSOR_CYCLE;
            ctx.mouse_button[0] = true;
        } else if(button == GLUT_MIDDLE_BUTTON) {
            cursor = GLUT_CURSOR_CROSSHAIR;
            ctx.mouse_button[1] = true;
        } else if(button == GLUT_RIGHT_BUTTON) {
            cursor = GLUT_CURSOR_UP_DOWN;
            ctx.mouse_button[2] = true;
        }
    } else {
        ctx.mouse_button[0] = ctx.mouse_button[1] = ctx.mouse_button[2] = false;
    }
    glutSetCursor(cursor);
    ctx.mouse_pos[0] = x;
    ctx.mouse_pos[1] = y;
    world_coords(x, y, ctx.drag);
}

void motion(int x, int y)
{
    int dx = x - ctx.mouse_pos[0];
    int dy = y - ctx.mouse_pos[1];
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    if(dx == 0 && dy == 0) {
        return;
    } else if(ctx.mouse_button[0]) {
        double angle = sqrt(dy*dy + dx*dx)/static_cast<double>(viewport[2]+1)*180.0;
        double rx = ctx.matinv[0]*dy + ctx.matinv[4]*dx;
        double ry = ctx.matinv[1]*dy + ctx.matinv[5]*dx;
        double rz = ctx.matinv[2]*dy + ctx.matinv[6]*dx;
        glRotatef(angle, rx, ry, rz);
    } else if(ctx.mouse_button[1]) {
        double p[3];
        world_coords(x, y, p);
        glLoadIdentity();
        glTranslatef(p[0] - ctx.drag[0], p[1] - ctx.drag[1], p[2] - ctx.drag[2]);
        glMultMatrixf(ctx.mat);
        ctx.drag[0] = p[0], ctx.drag[1] = p[1], ctx.drag[2] = p[2];
    } else if(ctx.mouse_button[2]) {
        glLoadIdentity();
        glTranslatef(0, 0, dy*0.01);
        glMultMatrixf(ctx.mat);
    }

    ctx.mouse_pos[0] = x;
    ctx.mouse_pos[1] = y;
    
    glGetFloatv(GL_MODELVIEW_MATRIX, ctx.mat);
    inverse(ctx.mat, ctx.matinv);
    glutPostRedisplay();
}

void idle()
{
    glutPostRedisplay();
}

void convert(tetrahedralizeio* out, tetrahedralizeio* vorout)
{
    for(size_t i = 0; i < vorout->numberofpoints; i++) {
        Vertex* v = new Vertex;
        v->v[0] = vorout->pointlist[3*i];
        v->v[1] = vorout->pointlist[3*i+1];
        v->v[2] = vorout->pointlist[3*i+2];
        v->has_normal = v->has_texture = false;
        vor.add(v);
    }
    for(size_t i = 0; i < vorout->numberofedges; i++) {
        Edge* e = new Edge;
        long vi = vorout->edgelist[2*i];
        long ui = vorout->edgelist[2*i+1];
        Vertex *v, *u;    
        
        if(vi == -1) {
            long tmp = vi;
            vi = ui;
            ui = tmp;
        }
         
        if (ui == -1) {
            double t = 10;
            v = vor.get_vertex(vi);
            e->v[0] = v;
            e->dir[0] = vorout->normlist[3*i];
            e->dir[1] = vorout->normlist[3*i+1];
            e->dir[2] = vorout->normlist[3*i+2];
            normalize(e->dir);
            
            u = new Vertex;
            u->v[0] = v->v[0] + t * e->dir[0];
            u->v[1] = v->v[1] + t * e->dir[1];
            u->v[2] = v->v[2] + t * e->dir[2];
            u->index = mesh.vertex_count();
            u->has_normal = u->has_texture = false;
            vor.add(u);

            e->v[1] = u;
            e->is_ray = true;
        } else {
            v = vor.get_vertex(vi);
            u = vor.get_vertex(ui);
            e->v[0] = v;
            e->v[1] = u;
            e->dir[0] = u->v[0] - v->v[0];
            e->dir[1] = u->v[1] - v->v[1];
            e->dir[2] = u->v[2] - v->v[2];

            if(norm(e->dir) > EPSILON) {
                normalize(e->dir);
            }
            e->is_ray = false;
        }
        e->intersected = false;
        vor.add(e);
    }
    
    for(size_t i = 0; i < out->numberofpoints; i++) {
        Vertex* v = new Vertex;
        v->v[0] = out->pointlist[3*i];
        v->v[1] = out->pointlist[3*i+1];
        v->v[2] = out->pointlist[3*i+2];
        v->has_normal = v->has_texture = false;
        del.add(v);
    }
    for(size_t i = 0; i < out->numberoffaces; i++) {
        Face* f = new Face;
        Vertex* v = del.get_vertex(out->facelist[3*i]);
        Vertex* u = del.get_vertex(out->facelist[3*i+1]);
        Vertex* w = del.get_vertex(out->facelist[3*i+2]);
        f->v.push_back(v);
        f->v.push_back(u);
        f->v.push_back(w);

        Edge* e = vor.get_edge(i);
        f->dual = e;
        e->dual = f;

        del.add(f);
    }
    compute_normals(del);
}

void delaunay()
{
    tetrahedralizeio* in = new tetrahedralizeio;
    tetrahedralizeio* out = new tetrahedralizeio;
    tetrahedralizeio* vorout = new tetrahedralizeio;
    
    in->pointlist = new REAL[mesh.vertex_count()*3];
    in->pointattributelist = NULL;
    in->pointmarkerlist = NULL;
    in->numberofpoints = mesh.vertex_count();
    in->numberofpointattributes = 0;

    out->pointlist = NULL;
    out->edgelist = NULL;
    out->facelist = NULL;

    vorout->pointlist = NULL;
    vorout->edgelist = NULL;
    vorout->normlist = NULL;
    vorout->facelist = NULL;
    
    vorout->numberofpointattributes = 0;

    for(size_t i = 0; i < mesh.vertex_count(); i++) {
        Vertex* v = mesh.get_vertex(i);
        in->pointlist[3*i] = v->v[0];
        in->pointlist[3*i+1] = v->v[1];
        in->pointlist[3*i+2] = v->v[2];
    }
    tetrahedralize("vzefQ", in, out, vorout);
    convert(out, vorout);
    
    delete in;
    delete out;
    delete vorout;
}

void init()
{
    glClearColor(0.0, 0.71, 1.0, 1.0);
    //glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel(GL_SMOOTH);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ctx.ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, ctx.diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, ctx.specular);
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(ctx.fovy, 1, ctx.znear, ctx.zfar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(ctx.width, ctx.height);
    glutCreateWindow("Del");
    init();
    menu();

    mesh = read_obj(argv[1]);
    center_on_screen(mesh);
    meshes.push_back(mesh);
    
    delaunay();
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(idle);

    glutMainLoop();

    return 0;
}
