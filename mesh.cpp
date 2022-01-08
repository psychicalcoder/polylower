#include "mesh.hpp"
#include "vec.hpp"
#include "util.hpp"
#include <stdexcept>
#include <vector>
#include <array>
#include <set>
#include <fstream>
#include <sstream>
#include <cstdio>

#include <pybind11/pybind11.h>

using namespace std;
namespace py = pybind11;

void Mesh::initQ() {
  vector<array<int, 3>>::iterator it = F.begin();
  Q.clear();
  Q.resize(V.size());
  while (it != F.end()) {
    array<int, 3> &f = (*it);
    vec &v0 = V[f[0]];
    vec v1 = V[f[1]] - v0;
    vec v2 = V[f[2]] - v0;
    vec n = (v1*v2).unit();
    double abcd[4] = {n.x, n.y, n.z,
      -(n.x * v0.x + n.y * v0.y + n.z * v0.z)};

    // printf("abcd %lf %lf %lf %lf\n", abcd[0], abcd[1], abcd[2], abcd[3]);
    
    for (int k = 0; k < 3; k++) {
      Mat4 &q = Q[f[k]];
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          q.m[i*4+j] += abcd[i] * abcd[j];
        }
      }
      // printf("vertex %d\n", f[k]);
      // for (int i = 0; i < 4; i++) {
      //   for (int j = 0; j < 4; j++) {
      //     printf(" %lf",  Q[f[k]].m[i*4+j]);
      //   }
      //   putchar('\n');
      // }
    }
    
    it++;
  }
}

double Mesh::calcCost(const vec& target, const Mat4 &q) {
  const double x[4] = { target.x, target.y, target.z, 1.0 };
  // printf("x %lf %lf %lf %lf\n", x[0], x[1], x[2], x[3]);
  // printf("q");
  // for (int i = 0; i < 16; i++) printf(" %lf", q.m[i]);
  // printf("\n");
  double cost = 0.0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      cost += x[i] * x[j] * q.m[i*4+j];
    }
  }
  // printf("cost %lf\n", cost);
  return cost;
}

bool Mesh::faceflip(const vec &target, int u, int v) const {
  for (int f : facelist[u]) {
    if (F[f][0] == v || F[f][1] == v || F[f][2] == v) continue;
    if (F[f][0] == u) {
      vec oldn = ((V[F[f][1]] - V[F[f][0]]) * (V[F[f][2]] - V[F[f][0]])).unit();
      vec newn = ((V[F[f][1]] - target) * (V[F[f][2]] - target)).unit();
      if ((oldn ^ newn) < 0.0) return true;
    } else if (F[f][1] == u) {
      vec oldn = ((V[F[f][1]] - V[F[f][0]]) * (V[F[f][2]] - V[F[f][0]])).unit();
      vec newn = ((target - V[F[f][0]]) * (V[F[f][2]] - V[F[f][0]])).unit();
      if ((oldn ^ newn) < 0.0) return true;
    } else {
      vec oldn = ((V[F[f][1]] - V[F[f][0]]) * (V[F[f][2]] - V[F[f][0]])).unit();
      vec newn = ((V[F[f][1]] - V[F[f][0]]) * (target - V[F[f][0]])).unit();
      if ((oldn ^ newn) < 0.9) return true;
    }
  }
  for (int f : facelist[v]) {
    if (F[f][0] == u || F[f][1] == u || F[f][2] == u) continue;
    if (F[f][0] == v) {
      vec oldn = ((V[F[f][1]] - V[F[f][0]]) * (V[F[f][2]] - V[F[f][0]])).unit();
      vec newn = ((V[F[f][1]] - target) * (V[F[f][2]] - target)).unit();
      if ((oldn ^ newn) < 0.0) return true;
    } else if (F[f][1] == v) {
      vec oldn = ((V[F[f][1]] - V[F[f][0]]) * (V[F[f][2]] - V[F[f][0]])).unit();
      vec newn = ((target - V[F[f][0]]) * (V[F[f][2]] - V[F[f][0]])).unit();
      if ((oldn ^ newn) < 0.0) return true;
    } else {
      vec oldn = ((V[F[f][1]] - V[F[f][0]]) * (V[F[f][2]] - V[F[f][0]])).unit();
      vec newn = ((V[F[f][1]] - V[F[f][0]]) * (target - V[F[f][0]])).unit();
      if ((oldn ^ newn) < 0.0) return true;
    }
  }
  return false;
}

void Mesh::addedge(int u, int v, int timestamp) {
  Mat4 Qbar = Q[u] + Q[v];
  vec target;
  int ret = Qbar.solve(target);
  if (ret < 0) {
    target.x = (V[u].x + V[v].x) / 2.0;    
    target.y = (V[u].y + V[v].y) / 2.0;
    target.z = (V[u].z + V[v].z) / 2.0;
  }
  double cost = calcCost(target, Qbar);
  if (!faceflip(target, u, v)) {
    // printf("cost %d %d = %lf\n", u, v, cost);
    edges.push(edge(u, v, cost, target, timestamp));
    // printf("edge %d %d cost %lf target %.6f %.6f %.6f\n", u, v, cost, target.x, target.y, target.z);
  }
}

void Mesh::selectEdges() {
  timestamps.clear();
  timestamps.resize(V.size(), 0);
  set<pair<int, int>> es;
  for (array<int, 3> &f : F) {
    for (int i = 0; i < 3; i++) { 
      pair<int, int> ep = f[i] > f[(i+1)%3] ? make_pair(f[(i+1)%3], f[i]) : make_pair(f[i], f[(i+1)%3]);
      if (!es.count(ep)) {
        es.insert(ep);
        addedge(ep.first, ep.second, 0);
      }
    }
  }
}

inline bool Mesh::is_edge_valid(const edge &e) const {
  return
    validV[e.u] &&
    validV[e.v] &&
    e.time_stamp >= timestamps[e.u] &&
    e.time_stamp >= timestamps[e.v]; //&& !faceflip(e.target, e.u, e.v);
}

void Mesh::simplify(double rate) {
  if (rate > 1.0) return;
  if (rate < 0.0) return;
  int num_iter = V.size() * (1 - rate);

  initQ();
  selectEdges();

  int current_time = 1;
  
  while (num_iter--) {
    while (!edges.empty() && !is_edge_valid(edges.top()))
      edges.pop();
    if (edges.empty()) break;
    
    edge e = edges.top(); edges.pop();

    // printf("contract %d %d -> (%lf, %lf, %lf)\n",
    //       e.u, e.v, e.target.x, e.target.y, e.target.z);

    Q[e.u] += Q[e.v];
    validV[e.v] = 0;
    timestamps[e.u] = current_time;
    V[e.u] = e.target;
    // printf("new %lf %lf %lf\n", V[e.u].x, V[e.u].y, V[e.u].z); 

    set<pair<int, int>> es;
    list<int> involvingface;

    for (int f_num : facelist[e.u]) {
      array<int, 3> &f = F[f_num];
      if (f[0] == e.v || f[1] == e.v || f[2] == e.v) {
        involvingface.push_back(f_num);
        continue;
      }
      pair<int, int> ep1, ep2;
      if (f[0] == e.u) {
        ep1 = f[1] > f[0] ? make_pair(f[0], f[1]) : make_pair(f[1], f[0]);
        ep2 = f[2] > f[0] ? make_pair(f[0], f[2]) : make_pair(f[2], f[0]);
      } else if (f[1] == e.u) {
        ep1 = f[0] > f[1] ? make_pair(f[1], f[0]) : make_pair(f[0], f[1]);
        ep2 = f[2] > f[1] ? make_pair(f[1], f[2]) : make_pair(f[2], f[1]);
      } else {
        ep1 = f[0] > f[2] ? make_pair(f[2], f[0]) : make_pair(f[0], f[2]);
        ep2 = f[1] > f[2] ? make_pair(f[2], f[1]) : make_pair(f[1], f[2]);   
      }
      if (!es.count(ep1)) {
        addedge(ep1.first, ep1.second, current_time);
        es.insert(ep1);
      }
      if (!es.count(ep2)) {
        addedge(ep2.first, ep2.second, current_time);
        es.insert(ep2);
      }
    }

    for (int f_num : facelist[e.v]) {
      array<int, 3> &f = F[f_num];
      if (f[0] == e.u || f[1] == e.u || f[2] == e.u) {
        // involvingface.push_back(f_num);
        // printf("hit2\n");
        continue;
      }
      pair<int, int> ep1, ep2;
      if (f[0] == e.v) {
        f[0] = e.u;
        ep1 = f[1] > f[0] ? make_pair(f[0], f[1]) : make_pair(f[1], f[0]);
        ep2 = f[2] > f[0] ? make_pair(f[0], f[2]) : make_pair(f[2], f[0]);
      } else if (f[1] == e.v) {
        f[1] = e.u;
        ep1 = f[0] > f[1] ? make_pair(f[1], f[0]) : make_pair(f[0], f[1]);
        ep2 = f[2] > f[1] ? make_pair(f[1], f[2]) : make_pair(f[2], f[1]);
      } else {
        f[2] = e.u;
        ep1 = f[0] > f[2] ? make_pair(f[2], f[0]) : make_pair(f[0], f[2]);
        ep2 = f[1] > f[2] ? make_pair(f[2], f[1]) : make_pair(f[1], f[2]);   
      }
      if (!es.count(ep1)) {
        addedge(ep1.first, ep1.second, current_time);
        es.insert(ep1);
      }
      if (!es.count(ep2)) {
        addedge(ep2.first, ep2.second, current_time);
        es.insert(ep2);
      }
      facelist[e.u].push_back(f_num);
    }

    for (int f_num : involvingface) {
      facelist[F[f_num][0]].remove(f_num);
      facelist[F[f_num][1]].remove(f_num);
      facelist[F[f_num][2]].remove(f_num);
      validF[f_num] = 0;
    }

    current_time++;
  }
  refresh();
}

void Mesh::aggregation(int epoch) {
  initQ();
  selectEdges();

  int current_time = 1;
  
  while (epoch--) {
    while (!edges.empty() && !is_edge_valid(edges.top()))
      edges.pop();
    if (edges.empty()) break;
    
    edge e = edges.top(); edges.pop();

    printf("contract %d (%lf, %lf, %lf) + %d (%lf, %lf, %lf) -> (%lf, %lf, %lf)\n",
           e.u, V[e.u].x, V[e.u].y, V[e.u].z,
           e.v, V[e.v].x, V[e.v].y, V[e.v].z,
           e.target.x, e.target.y, e.target.z);

    Q[e.u] += Q[e.v];
    validV[e.v] = 0;
    timestamps[e.u] = current_time;
    V[e.u] = e.target;
    // printf("new %lf %lf %lf\n", V[e.u].x, V[e.u].y, V[e.u].z); 

    set<pair<int, int>> es;
    list<int> involvingface;

    for (int f_num : facelist[e.u]) {
      array<int, 3> &f = F[f_num];
      if (f[0] == e.v || f[1] == e.v || f[2] == e.v) {
        involvingface.push_back(f_num);
        continue;
      }
      pair<int, int> ep1, ep2;
      if (f[0] == e.u) {
        ep1 = f[1] > f[0] ? make_pair(f[0], f[1]) : make_pair(f[1], f[0]);
        ep2 = f[2] > f[0] ? make_pair(f[0], f[2]) : make_pair(f[2], f[0]);
      } else if (f[1] == e.u) {
        ep1 = f[0] > f[1] ? make_pair(f[1], f[0]) : make_pair(f[0], f[1]);
        ep2 = f[2] > f[1] ? make_pair(f[1], f[2]) : make_pair(f[2], f[1]);
      } else {
        ep1 = f[0] > f[2] ? make_pair(f[2], f[0]) : make_pair(f[0], f[2]);
        ep2 = f[1] > f[2] ? make_pair(f[2], f[1]) : make_pair(f[1], f[2]);   
      }
      if (!es.count(ep1)) {
        addedge(ep1.first, ep1.second, current_time);
        es.insert(ep1);
      }
      if (!es.count(ep2)) {
        addedge(ep2.first, ep2.second, current_time);
        es.insert(ep2);
      }
    }

    for (int f_num : facelist[e.v]) {
      array<int, 3> &f = F[f_num];
      if (f[0] == e.u || f[1] == e.u || f[2] == e.u) {
        // involvingface.push_back(f_num);
        // printf("hit2\n");
        continue;
      }
      pair<int, int> ep1, ep2;
      if (f[0] == e.v) {
        f[0] = e.u;
        ep1 = f[1] > f[0] ? make_pair(f[0], f[1]) : make_pair(f[1], f[0]);
        ep2 = f[2] > f[0] ? make_pair(f[0], f[2]) : make_pair(f[2], f[0]);
      } else if (f[1] == e.v) {
        f[1] = e.u;
        ep1 = f[0] > f[1] ? make_pair(f[1], f[0]) : make_pair(f[0], f[1]);
        ep2 = f[2] > f[1] ? make_pair(f[1], f[2]) : make_pair(f[2], f[1]);
      } else {
        f[2] = e.u;
        ep1 = f[0] > f[2] ? make_pair(f[2], f[0]) : make_pair(f[0], f[2]);
        ep2 = f[1] > f[2] ? make_pair(f[2], f[1]) : make_pair(f[1], f[2]);   
      }
      if (!es.count(ep1)) {
        addedge(ep1.first, ep1.second, current_time);
        es.insert(ep1);
      }
      if (!es.count(ep2)) {
        addedge(ep2.first, ep2.second, current_time);
        es.insert(ep2);
      }
      facelist[e.u].push_back(f_num);
    }

    for (int f_num : involvingface) {
      facelist[F[f_num][0]].remove(f_num);
      facelist[F[f_num][1]].remove(f_num);
      facelist[F[f_num][2]].remove(f_num);
      validF[f_num] = 0;
    }
    current_time++;
  }
  refresh();
}

void Mesh::load(const string &filename) {
  ifstream obj(filename);
  if (!obj.is_open()) {
    throw invalid_argument("Failed to open file");
  }
  V.clear();
  F.clear();
  vector<vector<int>> faces;
  char *buf = new char[256];
  while (!obj.eof()) {
    obj.getline(buf, 256);
    string attr;
    stringstream ss(buf);
    ss >> attr;
    if (attr == "v") {
      double x, y, z;
      ss >> x >> y >> z;
      V.push_back(vec(x, y, z));
    } else if (attr == "f") {
      string vstr;
      vector<int> vlist;
      while (ss >> vstr) {
        int v_idx;
        stringstream vstm;
        vstm << vstr;
        vstm >> v_idx;
        vlist.push_back(v_idx-1);
      }
      faces.push_back(vlist);
    }
  }
  delete [] buf;
  obj.close();
  
  triangularize(V, faces, F);
  validV.clear();
  validV.resize(V.size(), 1);
  validF.clear();
  validF.resize(F.size(), 1);
  refresh();
  
}

void Mesh::dump(const string &filename) {
  ofstream obj(filename);
  obj << "# n_vertices " << V.size() << '\n';
  obj << "# n_faces " << F.size() << '\n';
  vector<vec>::const_iterator it = V.cbegin();
  while (it != V.cend()) {
    obj << "v " << it->x << ' ' << it->y << ' ' << it->z << '\n';
    it++;
  }

  for (const array<int, 3> &f : F) {
    obj << "f " << f[0]+1 << ' ' << f[1]+1 << ' ' << f[2]+1 << '\n';
  }
  obj.flush();
  obj.close();
}

void Mesh::refresh() {
  int newid = 0;
  vector<int> vlink(V.size(), 0);
  for (int oldid = 0; oldid < V.size(); oldid++) {
    if (validV[oldid]) {
      vlink[oldid] = newid;
      V[newid] = V[oldid];
      newid++;
    }
  }
  V.resize(newid);
  newid = 0;
  for (int oldid = 0; oldid < F.size(); oldid++) {
    if (validF[oldid]) {
      F[newid][0] = vlink[F[oldid][0]];
      F[newid][1] = vlink[F[oldid][1]];
      F[newid][2] = vlink[F[oldid][2]];
      newid++;
    }
  }
  F.resize(newid);

  facelist.clear();
  facelist.resize(V.size());
  int fc = 0;
  for (array<int, 3> &f : F) {
    facelist[f[0]].push_back(fc);
    facelist[f[1]].push_back(fc);
    facelist[f[2]].push_back(fc);
    fc++;
  }
  validV.resize(V.size());
  fill(validV.begin(), validV.end(), 1);
  validF.resize(F.size());
  fill(validF.begin(), validF.end(), 1);
  while (!edges.empty()) edges.pop();
  timestamps.resize(V.size());
  fill(timestamps.begin(), timestamps.end(), 0);
}

PYBIND11_MODULE(polylower, m) {
  m.doc() = "Polylower";
  py::class_<Mesh>(m, "Mesh")
    .def(py::init<>())
    .def("load", &Mesh::load)
    .def("dump", &Mesh::dump)
    .def("aggregation", &Mesh::aggregation)
    .def("simplify", &Mesh::simplify);
}
