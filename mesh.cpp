#include "mesh.h"
#include "assert.h"
#include <map>

const double THRE = 0.5;
const long P = 100663319;
const int ADJ_DEGREE = 20;
using namespace std;

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
            vj->n[0] = 0;
            vj->n[1] = 0;
            vj->n[2] = 0;
        }

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


void revert_face(Face* f);
void compute_adjacent(Mesh& mesh, int adj_face[][ADJ_DEGREE+1]);

//void empty(){
//    int i;
//    for (i = 0; i < P; i++){
//    //}
//}
//long hash(int i, int j, int v_size){
//    return (i*v_size+j) % P;
//}
//
//bool find(int i, int j, )

void compute_adjacent(Mesh& mesh, int adj_face[][ADJ_DEGREE+1]) {
  int i,j,size;
  int v_size = mesh.vertex_count();
  int f_size = mesh.face_count();
  long index;
  map<long, int> face_for_e;

  //int face_for_e[v_size][v_size];
  int me, other;
  //for (i = 0; i<v_size; i++){
  //    for (j=0; j<v_size; j++){
  //        face_for_e[i][j]= -1;
  //    }
  //}
  for (i = 0; i< f_size; i++){
      adj_face[i][0] = 0;
  }
  for (i = 0; i< mesh.face_count(); i++){
    Face * f = mesh.get_face(i);
    int v0= f->v[2]->index; 
    int v1;
    for (j = 0; j< 3; j++){
      v1 = f->v[j]->index;
      int smaller, greater;
      if (v0 < v1) {
        smaller = v0;
        greater = v1;
      } else {
        smaller = v1;
        greater = v0;
      }
      index = smaller*v_size+greater;
      if (face_for_e.find(index)==face_for_e.end()) {
          face_for_e[index] = f->index;
      } else {

        other = face_for_e[index];
         
 
        //other = face_for_e[smaller][greater];
        me = f->index;
        //printf("other: %d\n", other);
        //printf("me: %d\n", me);
        //if (adj_face[me][0] <3){
          adj_face[me][0]++;
          size = adj_face[me][0];
          adj_face[me][size] = other;
        //}
        //if (adj_face[other][0] <3){
          adj_face[other][0]++;
          size = adj_face[other][0];
          adj_face[other][size] = me;
        //}
      }
      v0 = v1;
    }
  }
  for (i=0; i < v_size; i++){
      Vertex * v = mesh.get_vertex(i);
      //printf("vertex %d:(%f,%f,%f)\n", v->index, v->v[0], v->v[1], v->v[2]);
  }
  /*
  for (i=0; i < v_size; i++){
      for (j=0; j < v_size; j++)
        if (face_for_e.find(i*v_size+j)!=face_for_e.end())
            printf("face_for_e(%d,%d)=%d; ",i,j,face_for_e[i*v_size+j]);
      printf("\n");
  }
      */
}


/* reverse the direction of f2 if it 
 * doesn't agree with f's
 *       v1
 *     / | \
 *    /  |  \
 * v2 \n1|n2 / v3
 *     \ |  /
 *      \| /
 *      v0
 */

void fix_adjacent_normal_new(Face* f1, Face* f2){
  /* v0, v1 are the shared vertex between f1,f2 */
  /* return true if it did revert the face
   * false if no change
   */
  int i,j;
  int v0_1 = -1;
  int v0_2 = -1;
  int v1_1, v1_2;
  for (i = 0; i < 3; i++){
      int me = f1->v[i]->index;
      for (j = 0; j < 3; j++){
          if (f2->v[j]->index == me) {
              if (v0_1 == -1) {
                v0_1 = i;
                v0_2 = j;
              } else {
                v1_1 = i;
                v1_2 = j;
              }
          }
      }
  }
  bool clockwise1, clockwise2; 
  /* clockwise v0->v1->v2->v0
   * countclockwise v0->v2->v1->v0 */
  if ((v0_1 == v1_1-1) || ((v0_1 == 2) && (v1_1 == 0))){
      clockwise1 = true;
  } else {
      clockwise1 = false;
  }
  if ((v0_2 == v1_2-1) || ((v0_2 == 2) && (v1_2 == 0))){
      clockwise2 = true;
  } else {
      clockwise2 = false;
  }
  if (clockwise1 == clockwise2) {
    //printf("face%d is reverted\n", f2->index);
    revert_face(f2);
  }
}

void dfs_from_face(Mesh& mesh, Face* f, int adj_face[][1+ADJ_DEGREE]){
    f->visited = true;
    int me = f->index;
    //printf("face %d's adjacent to: ", me);
    for (int i = 1; i <= adj_face[me][0]; i++) {
        int other = adj_face[me][i];
        Face * f2 = mesh.get_face(other);
        //printf("face %d, ", other);
        if (!f2->visited) {
            fix_adjacent_normal_new(f, f2);
            dfs_from_face(mesh, f2, adj_face);
        }
    }
    //printf("\n");
}
bool need_fix_adjacent(Face* f1, Face* f2){
  /* return true if f2 have wrong normal direction as regard to f1
   * false if no change
   */
  int i,j;
  int v0_1 = -1;
  int v0_2 = -1;
  int v1_1, v1_2;
  for (i = 0; i < 3; i++){
      int me = f1->v[i]->index;
      for (j = 0; j < 3; j++){
          if (f2->v[j]->index == me) {
              if (v0_1 == -1) {
                v0_1 = i;
                v0_2 = j;
              } else {
                v1_1 = i;
                v1_2 = j;
              }
          }
      }
  }
  bool clockwise1, clockwise2; 
  /* clockwise v0->v1->v2->v0
   * countclockwise v0->v2->v1->v0 */
  if ((v0_1 == v1_1-1) || ((v0_1 == 2) && (v1_1 == 0))){
      clockwise1 = true;
  } else {
      clockwise1 = false;
  }
  if ((v0_2 == v1_2-1) || ((v0_2 == 2) && (v1_2 == 0))){
      clockwise2 = true;
  } else {
      clockwise2 = false;
  }
  return (clockwise1 == clockwise2);
}

void check_normals(Mesh& mesh, int adj_face[][1+ADJ_DEGREE]){
    for (int i = 0; i < mesh.face_count(); i++){
        int me, other;
        Face * f = mesh.get_face(i);
        me = f->index;
        bool error = 0;
        for (int j = 1; j <= adj_face[me][0]; j++) {
          other = adj_face[me][j];
          Face * f2 = mesh.get_face(other);
          //printf("face %d, ", other);
          if (need_fix_adjacent(f, f2)){
              error = 1;
              printf("face%d is opposite to face%d\n", f2->index, f->index);
              printf("face%d: (%d,%d,%d)\n", f2->index, f2->v[0]->index, f2->v[1]->index, f2->v[2]->index);
              printf("face%d: (%d,%d,%d)\n", f->index, f->v[0]->index, f->v[1]->index, f->v[2]->index);
          }
        }
        if (error){
          printf("face%d is adjacent to:", me);
          for (int j = 1; j <= adj_face[me][0]; j++) {
            other = adj_face[me][j];
            printf("face%d, ", other);
          }
          printf("\n\n");
        }
   }
}

void adjust_normals(Mesh& mesh)
{
    int f_size = mesh.face_count();
    int adj_face[f_size][ADJ_DEGREE+1];
    compute_adjacent(mesh, adj_face);
    for (int i = 0; i < f_size; i++) {
        if (adj_face[i][0] > 3) {
            //printf("face %d's adjacent: ", i);
                Face * f = mesh.get_face(i);
                int k;
                for (k=0; k < 3; k++){
                      Vertex * v = f->v[k];
                      //printf("vertex %d:(%f,%f,%f)\n", v->index, v->v[0], v->v[1], v->v[2]);
                  }
            for (int j = 0; j <= adj_face[i][0]; j++){
                //printf("%d, \n", adj_face[i][j]);
                Face * f = mesh.get_face(adj_face[i][j]);
                int k;
                for (k=0; k < 3; k++){
                      Vertex * v = f->v[k];
                      //printf("vertex %d:(%f,%f,%f)\n", v->index, v->v[0], v->v[1], v->v[2]);
                  }

            }
            //printf("\n");
        }
    }
    
     
    for (int i = 0; i < f_size; i++){
        mesh.get_face(i)->visited = false; 
    }
    for (int i = 0; i < f_size; i++){
      if (!mesh.get_face(i)->visited) {
        dfs_from_face(mesh, mesh.get_face(i), adj_face);
      }
    }
    compute_normals(mesh);
    //check_normals(mesh, adj_face);
}

void revert_face(Face* f){
    revert(f->n);
    Vertex * tmp = f->v[0];
    f->v[0] = f->v[1];
    f->v[1] = tmp;
}

void revert_normals(Mesh& mesh) {
    //printf("hello! %d", mesh.face_count());
    for (int i = 0; i < mesh.face_count(); i++){
        Face * f = mesh.get_face(i);
        //printf("%f--->", mesh.get_face(i)->n[0]);
        revert_face(f);
        //revert((mesh.get_face(i))->n);
        //printf("%f\n", mesh.get_face(i)->n[0]);
    }
    //for (int i = 0; mesh.vertex_count(); i++){
    //    //Vertex * v = mesh.get_vertex(i);
    //    //printf("%f--->", mesh.get_face(i)->n[0]);
    //    revert((mesh.get_vertex(i))->n);
    //    //printf("%f\n", mesh.get_face(i)->n[0]);
    //}
}

void update_color_var(Mesh &mesh, Mesh &prev_mesh, bool first){
    if (first) {
        for (int i = 0; i < mesh.face_count(); i++){
            mesh.get_face(i)->to_disappear = false;
            mesh.get_face(i)->newly_appear = false;
            mesh.get_face(i)->fix_level = 0;
        }
        return;
    }
    for (int i = 0; i < prev_mesh.face_count(); i++){
        prev_mesh.get_face(i)->to_disappear = true;
    }
    for (int i = 0; i < mesh.face_count(); i++){
        mesh.get_face(i)->to_disappear = false;
        mesh.get_face(i)->newly_appear = true;
    }
    for (int i = 0; i < prev_mesh.face_count(); i++){
        prev_mesh.get_face(i)->newly_appear = false;
    }
    for (int i = 0; i < mesh.face_count(); i++){
        Face * f = mesh.get_face(i);
        if (f->dual){
            int count = 5;
            Face * walk = f;
            while ((walk) && (walk->dual) && (walk->dual->intersected_face != walk) && (count > 0)){
                walk = walk->dual->intersected_face;
                count--;
            }
            if (count <= 0){
                printf("error: 5 is not the max depth of dependencies");
            } else if (walk && (walk->dual) && (walk->dual->intersected_face == walk)) {
                f->fix_level = count;
            //printf("count:%d\n",count);
            } else {
                f->fix_level = 0;
            }

        } else {
            f->fix_level = 0;
        }
    }

}
