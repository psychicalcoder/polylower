#include "mesh.hpp"
#include "triangulation2d.hpp"
#include <cmath>
#include <cstdio>
namespace polylower {

  Mesh::Mesh(int num_vertices, int num_surfaces) :
    vertices(num_vertices), faces(num_surfaces) {
  }
  
  Mesh::Mesh(vector<Vertex> &&vs, vector<vector<int> > &&fs) :
    vertices(vs), faces(fs) {}
  
  int Mesh::num_faces() const { return faces.size(); }  
  int Mesh::num_vertex() const { return vertices.size() - 1; }

  std::vector<vector<int>>& Mesh::get_faces() { return this->faces; }
  std::vector<Vertex>& Mesh::get_vertices() { return this->vertices; }
  
  std::vector<vector<int>> Mesh::get_faces() const { return this->faces; }
  std::vector<Vertex> Mesh::get_vertices() const { return this->vertices; }

  
  void Mesh::triangulation() {
    int numfaces = this->num_faces();
    for (int i = 0; i < numfaces; i++) {
      // printf("face %d\n", i+1);
      size_t v_num = faces[i].size();
      if (v_num > 3) {
        vector<int> vns = faces[i];
        Vertex v0 = vertices[vns[0]];
        Vertex v1 = unit(vertices[vns[1]] - v0);
        Vertex v2 = unit(
                         cross(cross(v1,
                                     vertices[vns[2]] - v0),
                               v1
                               )
                         );
        
        Point *p = new Point[v_num];
        for (int j = 0; j < v_num; j++) {
          Vertex t = vertices[vns[j]] - v0;
          p[j].x = dot(t, v1);
          p[j].y = dot(t, v2);
          // printf("  vertex %d -> (%f, %f)\n", j, p[j].x, p[j].y);
        }
        vector<int> tris = triangulation2d(p, v_num);
        for (int i = 3; i < tris.size(); i += 3) {
          this->faces.push_back(vector<int>(3));
          this->faces.back()[0] = vns[tris[i]];
          this->faces.back()[1] = vns[tris[i+1]];
          this->faces.back()[2] = vns[tris[i+2]];
        }
        this->faces[i].resize(3);
        this->faces[i][0] = vns[tris[0]];
        this->faces[i][1] = vns[tris[1]];
        this->faces[i][2] = vns[tris[2]];
      }
    }
  }

  Vertex Vertex::operator-(const Vertex &subtrahend) const {
    double xx = this->pos_x - subtrahend.pos_x;
    double yy = this->pos_y - subtrahend.pos_y;
    double zz = this->pos_z - subtrahend.pos_z;
    return Vertex(xx, yy, zz);
  }

  bool Vertex::operator==(const Vertex &rhs) const {
    return this->pos_x == rhs.pos_x && this->pos_y == rhs.pos_y && this->pos_z == rhs.pos_z;
  }

  Vertex cross(const Vertex &u, const Vertex &v) {
    double x = u.y() * v.z() - u.z() * v.y();
    double y = u.z() * v.x() - u.x() * v.z();
    double z = u.x() * v.y() - u.y() * v.x();
    return Vertex(x, y, z);
  }

  Vertex unit(const Vertex &v) {
    double len = sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
    return Vertex(v.x() / len, v.y() / len, v.z() / len);
  }

  double dot(const Vertex &u, const Vertex &v) {
    return u.x() * v.x() + u.y() * v.y() + u.z() * v.z();
  }

}; // namespace polylower
