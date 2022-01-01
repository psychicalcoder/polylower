#ifndef POLYLOWER_MESH_H
#define POLYLOWER_MESH_H

#include <vector>
#include <list>
#include <queue>
#include <array>
#include "vec.hpp"
#include "util.hpp"
using namespace std;

struct edge {
  int u, v;
  double cost;
  vec target;
  int time_stamp;
  edge(int u, int v, double cost, const vec &target, int time_stamp)
    : u(u), v(v), cost(cost), target(target), time_stamp(time_stamp) {}
  bool operator<(const edge &rhs) const {
    return cost < rhs.cost;
  }
};


class Mesh {
private:
  vector<vec> V;
  vector<array<int, 3>> F;
  vector<Mat4> Q;
  vector<int> validV;
  vector<int> validF;
  vector<int> timestamps;
  vector<list<int>> facelist;
  priority_queue<edge> edges;
  void initQ();
  double calcCost(const vec& target, const Mat4 &q);
  bool faceflip(const vec& target, int u, int v);
  void selectEdges();
  void addedge(int u, int v, int timestamp);
  bool is_edge_valid(const edge &e) const;
  void remove_vertex(int v);
  void refresh();
  
public:
  Mesh() {}
  ~Mesh() {}
  void simplify(double rate);
  void load(const string &filename);
  void dump(const string &filename);
};

#endif
