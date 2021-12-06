#ifndef _POLYLOWER_MESH_H
#define _POLYLOWER_MESH_H

#include <vector>

namespace polylower {

  using namespace std;
  
  class Vertex {
  public:
    Vertex(): pos_x(0.0f), pos_y(0.0f), pos_z(0.0f) {}
    Vertex(float pos_x, float pos_y, float pos_z):
      pos_x(pos_x), pos_y(pos_y), pos_z(pos_z) {}
    float x() const { return this->pos_x; }
    float y() const { return this->pos_y; }
    float z() const { return this->pos_z; }

    Vertex operator-(const Vertex &subtrahend) const;

  private:
    float pos_x, pos_y, pos_z;
  };
  
  class Mesh {
  private:
    vector<Vertex> vertices;
    vector<vector<int> > faces;
  public:
    Mesh(int, int);
    Mesh(vector<Vertex> &&, vector<vector<int> > &&);
    Mesh() = default;
    int num_vertex() const ;
    int num_faces() const ;
    vector<Vertex>& get_vertices();
    vector<vector<int>>& get_faces();
    
    vector<Vertex> get_vertices() const ;
    vector<vector<int>> get_faces() const ;
  };

}

#endif
