#include "mesh.hpp"

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

  Vertex Vertex::operator-(const Vertex &subtrahend) const {
    double xx = this->pos_x - subtrahend.pos_x;
    double yy = this->pos_y - subtrahend.pos_y;
    double zz = this->pos_z - subtrahend.pos_z;
    return Vertex(xx, yy, zz);
  }
  
}; // namespace polylower
