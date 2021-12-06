#include "obj.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <exception>
#include <stdexcept>

using namespace std;

namespace polylower {

  Mesh load_mesh(const char *path) {
    ifstream obj(path);
    if (!obj.is_open()) {
      throw invalid_argument("Failed to open file!");
    }
    vector<Vertex> vertices;
    vector<vector<int> > faces;
    const size_t buff_len = 256;
    char *buff = new char[buff_len];
    vertices.push_back(Vertex());
    while (!obj.eof()) {
      obj.getline(buff, buff_len);
      string attr;
      stringstream ss(buff);
      ss >> attr;
      if (attr == "v") {
	double x, y, z;
	ss >> x >> y >> z;
	vertices.push_back(Vertex(x, y, z));
      } else if (attr == "f") {
	string vstr;
	vector<int> vlist;
	while (ss >> vstr) {
	  int v_idx;
	  stringstream vstm;
	  vstm << vstr;
	  vstm >> v_idx;
	  vlist.push_back(v_idx);
	}
	faces.push_back(vlist);
      }
    }
    delete [] buff;
    obj.close();
  
    return Mesh(std::move(vertices), std::move(faces));
  }

  void save_mesh(const char *path, const Mesh &mesh) {
    ofstream obj(path);
    obj << "# n_vertices " << mesh.num_vertex() << '\n';
    obj << "# n_faces " << mesh.num_faces() << '\n';
    const vector<Vertex> &vs = mesh.get_vertices();
    vector<Vertex>::const_iterator itr = vs.cbegin();
    itr++;
    while (itr != vs.cend()) {
      obj << "v " << itr->x() << " " << itr->y() << " " << itr->z() << '\n';
      itr++;
    }

    for (const vector<int> &f : mesh.get_faces()) {
      obj << "f";
      for (const int v : f) {
	obj << ' ' << v;
      }
      obj << '\n';
    }
    obj.flush();
    obj.close();
  }

}
