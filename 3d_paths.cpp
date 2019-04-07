#include <cstring>
#include <iostream>

int grid[2][2][2];

void solve(int c, int i, int j, int k) {
  if (c == 9) {
    std::cout << "solution" << std::endl;
    for (int pi = 0; pi < 2; pi++) {
      for (int pj = 0; pj < 2; pj++) {
	for (int pk = 0; pk < 2; pk++) {
	  std::cout << grid[pi][pj][pk] << '\t';
	}
	std::cout << std::endl;
      }
    }
    return;
  }
  int dirs[6][3] = {
      {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}
  };
  for (int dir = 0; dir < 6; dir++) {
    int ni = i + dirs[dir][0];
    int nj = j + dirs[dir][1];
    int nk = k + dirs[dir][2];
    if (ni < 0 || ni >= 2 || nj < 0 || nj >= 2 || nk < 0 || nk >= 2 || grid[ni][nj][nk] != 0) {
      continue;
    }
    grid[ni][nj][nk] = c;
    solve(c + 1, ni, nj, nk);
    grid[ni][nj][nk] = 0;
  }
}

int main() {
  std::memset(grid, 0, sizeof(grid));
  grid[0][0][0] = 1;
  grid[0][0][1] = 2;
  solve(3, 0, 0, 1);
  return 0;
}
