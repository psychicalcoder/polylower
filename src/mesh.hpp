#ifndef _POLYLOWER_MESH_H
#define _POLYLOWER_MESH_H

#include <vector>

namespace polylower {

  using namespace std;
  
  class Vertex {
  public:
    Vertex(): pos_x(0.0), pos_y(0.0), pos_z(0.0) {}
    Vertex(double pos_x, double pos_y, double pos_z):
      pos_x(pos_x), pos_y(pos_y), pos_z(pos_z) {}
    double x() const { return this->pos_x; }
    double y() const { return this->pos_y; }
    double z() const { return this->pos_z; }

    
    void set_x(double x) { this->pos_x = x; }
    void set_y(double y) { this->pos_y = y; }
    void set_z(double z) { this->pos_z = z; }

    Vertex operator-(const Vertex &subtrahend) const;
    bool operator==(const Vertex &rhs) const;

  private:
    double pos_x, pos_y, pos_z;
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

    void triangulation();
  };

  Vertex unit(const Vertex& v);
  Vertex cross(const Vertex &u, const Vertex &v);
  double dot(const Vertex &u, const Vertex &v);

}

#endif
