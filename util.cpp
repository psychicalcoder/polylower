#include "util.hpp"
#include "triangulation2d.hpp"
#include <cmath>
#include <utility>

using namespace std;
using namespace polylower;

int Mat4::solve(vec &ans) {
  double **T = new double*[4];
  for (int i = 0; i < 3; i++) {
    T[i] = new double[5];
    for (int j = 0; j < 4; j++)
      T[i][j] = m[i*4 + j];
  }
  T[3] = new double[5];
  T[3][0] = T[3][1] = T[3][2] = 0.0;
  T[3][4] = 1.0;

  int row = 0, col = 0;
  while (row < 4 && col < 4) {
    int maxi = row;
    for (int k = row+1; k < 4; k++) {
      if (fabs(T[k][col]) > fabs(T[maxi][col])) {
        maxi = k;
      }
    }
    if (fabs(T[maxi][col]) > PL_EPS) {
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

  if (row != 4) {
    for (int i = 0; i < 4; i++) delete T[i];
    delete [] T;
    return -1;
  }

  ans.x = T[0][4];
  ans.y = T[1][4];
  ans.z = T[2][4];

  for (int i = 0; i < 4; i++) delete T[i];
  delete [] T;
  
  return 0;
}

void triangularize(const vector<vec> &vertices, const vector<vector<int>> &faces, vector<array<int, 3>> &tfaces) {
  int facenum = faces.size();
  for (int i = 0; i < facenum; i++) {
    int v_num = faces[i].size();
    if (v_num > 3) {
      vector<int> vns = faces[i];
      vec v0 = vertices[vns[0]];
      vec v1 = (vertices[vns[1]] - v0).unit();
      vec v2 = ((v1*(vertices[vns[2]]-v0))*v1).unit();
      Point *p = new Point[v_num];
      for (int j = 0; j < v_num; j++) {
        vec t = vertices[vns[j]] - v0;
        p[j].x = t ^ v1;
        p[j].y = t ^ v2;
      }
      vector<int> tris = triangulation2d(p, v_num);
      for (int i = 0; i < tris.size(); i += 3) {
        tfaces.push_back(array<int, 3>{vns[tris[i]], vns[tris[i+1]], vns[tris[i+2]]});
      }
    } else {
      tfaces.push_back(array<int, 3>{faces[i][0], faces[i][1], faces[i][2]});
    }
  }
}
