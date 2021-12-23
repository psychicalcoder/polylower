#include "polylower.hpp"
#include "mesh.hpp"
#include <vector>
#include <queue>
#include <list>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <sstream>
#define EPS 1e-6
using namespace std;
namespace py = pybind11;

namespace polylower {

  struct Mat4 {
    double m[4][4];
  };

  struct NotSolvableException : public exception {
    const char *info;
    NotSolvableException() : info("The linear system is not solvable") {}
    NotSolvableException(const char *info) : info(info) {}
  };
 
  vector<Mat4> calcQ(const Mesh &mesh) {
    vector<Mat4> Q(mesh.num_vertex()+1);
    vector<Mat4> faceQ(mesh.num_faces());
    int face_i = 0;
    const vector<Vertex> &V = mesh.get_vertices();
    const vector<vector<int>> &F = mesh.get_faces();
    list<int> *v2p = new list<int>[mesh.num_vertex()+1];
    for (int i = 0; i < F.size(); i++) {
      const vector<int> &vs = F[i];
      Vertex v1 = V[vs[1]] - V[vs[0]];
      Vertex v2 = V[vs[2]] - V[vs[0]];
      Vertex norm = unit(cross(v1, v2));
      double abcd[4];
      abcd[0] = norm.x();
      abcd[1] = norm.y();
      abcd[2] = norm.z();
      abcd[3] = - ( norm.x() * V[vs[0]].x() +
		    norm.y() * V[vs[0]].y() +
		    norm.z() * V[vs[0]].z() );
      for (size_t row = 0; row < 4; row++) {
	for (size_t col = 0; col < 4; col++) {
	  faceQ[i].m[row][col] = abcd[row] * abcd[col];
	}
      }
      for (int vi : vs) {
	v2p[vi].push_back(i);
      }
    }
    for (int i = 1; i <= mesh.num_vertex(); i++) {
      for (int f : v2p[i]) {
	for (size_t row = 0; row < 4; row++) {
	  for (size_t col = 0; col < 4; col++) {
	    Q[i].m[row][col] += faceQ[f].m[row][col];
	  }
	}
      }
    }

    return Q;
  }

  Vertex gaussian4(const Mat4 &stm) {
    double **T;
    T = new double*[4];
    for (int i = 0; i < 4; i++) {
      T[i] = new double[5];
      for (int j = 0; j < 4; j++) T[i][j] = stm.m[i][j];
      T[i][4] = 0.0;
    }
    T[3][4] = 1.0;

    // gaussian elimination
    int row = 0, col = 0;
    while (row < 4 && col < 5) {
      int maxi = row;
      for (int k = row+1; k < 4; k++) {
	if (abs(T[k][col]) > abs(T[maxi][col])) {
	  maxi = k;
	}
      }
      if (abs(T[maxi][col]) > EPS) {
	swap(T[maxi], T[row]);
	for (int k = row+1; k < 4; k++) {
	  double scale = T[k][col] / T[row][col];
	  for (int j = row; j < 5; j++) {
	    T[k][j] -=  scale * T[row][j];
	  }
	}
	row++;
      }
      col++;
    }

    if (row != 4) throw NotSolvableException();
    // jordan 
    for (int i = 3; i >= 0; i--) {
      T[i][4] /= T[i][i];
      T[i][i] = 1.0;
      if (abs(T[i][i] - T[i][4]) > EPS)
	throw NotSolvableException();
      for (int j = i-1; j >= 0; j--) {
        T[j][4] -= T[j][i] * T[i][4];
	T[j][i] = 0.0;
      }
    }

    for (int i = 0; i < 4; i++) delete [] T[i];
    delete [] T;
    
    return Vertex(T[0][4], T[1][4], T[2][4]);
  }
  
  Mesh simplify(const Mesh &mesh, double ratio) {
    int num_iter = mesh.num_vertex() - int(mesh.num_vertex() * ratio);

    vector<Mat4> Q = calcQ(mesh);
    
    return Mesh();
  }

} // namespace polylower

PYBIND11_MODULE(polylower, m) {
  m.doc() = "Polylower";
  m.def("unit", &polylower::unit, "return an unit vector");
  m.def("dot", &polylower::dot, "dot product");
  m.def("cross", &polylower::cross, "cross product");
  py::class_<polylower::Vertex>(m, "Vertex")
    .def(py::init<>())
    .def(py::init<double, double, double>())
    .def_property("x", &polylower::Vertex::x, &polylower::Vertex::set_x)
    .def_property("y", &polylower::Vertex::y, &polylower::Vertex::set_y)
    .def_property("z", &polylower::Vertex::z, &polylower::Vertex::set_z)
    .def("__repr__", [](const polylower::Vertex &v) {
      std::stringstream ss;
      ss << "<polylower.Vertex (" << v.x() << ", " << v.y() << ", " << v.z() << ")>";
      return ss.str();
    })
    .def(py::self - py::self)
    .def(py::self == py::self);
}
