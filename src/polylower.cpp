#include "polylower.hpp"
#include <vector>
#include <queue>
#include <cmath>
#include <list>
#include <stdexcept>
#define EPS 1e-6
using namespace std;

namespace polylower {

  struct Mat4 {
    double m[4][4];
  };

  struct NotSolvableException : public exception {
    const char *info;
    NotSolvableException() : info("The linear system is not solvable") {}
    NotSolvableException(const char *info) : info(info) {}
  };

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
  
}
