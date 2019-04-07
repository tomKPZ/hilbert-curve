#include <cmath>
#include <cstdint>
#include <iostream>

constexpr int N = 4;
float sqrt2 = std::sqrt(2);

float cache[N][N][N][N];

float solution_error = 999999999;
int current[N][N];

void print_solution() {
  std::cout << "Error: " << solution_error << '\n';
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      std::cout << current[i][j] << '\t';
    }
    std::cout << '\n';
  }
}

void copy_current_to_solution(float current_error) {
  solution_error = current_error;
  print_solution();
}

float error(int depth, int i1, int j1) {
  float err = 0;
  for (int i2 = 0; i2 < N; i2++) {
    for (int j2 = 0; j2 < N; j2++) {
      if (current[i2][j2] == -1)
	continue;

      float expected = cache[i1][j1][i2][j2];
      float actual = std::abs(depth - current[i2][j2]);

      // Expected: [0, sqrt(2)*(N-1)]
      // Actual:   [0, N^2 - 1]
      err += std::abs((N*N - 1) * expected - sqrt2 * (N - 1) * actual);
    }
  }
  return err;
}

template <int depth> void solve(float current_error) {
  if (current_error >= solution_error)
    return;

  if constexpr (depth == N * N) {
    int err = current_error;
    if (err < solution_error) {
      solution_error = err;
      copy_current_to_solution(current_error);
    }
    return;
  } else {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (current[i][j] != -1)
          continue;

	float new_error = error(depth, i, j);
        current[i][j] = depth;
        solve<depth + 1>(current_error + new_error);
        current[i][j] = -1;
      }
    }
  }
}

void solve() {
  current[0][0] = 0;
  for (int i = 0; i < N; i++) {
    for (int j = i; j < N; j++) {
      if (current[i][j] != -1)
	continue;
	      
      float new_error = error(1, i, j);
      current[i][j] = 1;
      solve<2>(new_error);
      current[i][j] = -1;
    }
  }
  current[0][0] = -1;

  current[0][1] = 0;
  solve<1>(0);
  current[0][1] = -1;

  current[1][1] = 0;
  for (int i = 0; i < N; i++) {
    for (int j = i; j < N; j++) {
      if (current[i][j] != -1)
	continue;
	      
      float new_error = error(1, i, j);
      current[i][j] = 1;
      solve<2>(new_error);
      current[i][j] = -1;
    }
  }
  current[1][1] = -1;
}

void init_current() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      current[i][j] = -1;
    }
  }
}

void init_cache() {
  for (int i1 = 0; i1 < N; i1++) {
    for (int j1 = 0; j1 < N; j1++) {
      for (int i2 = 0; i2 < N; i2++) {
	for (int j2 = 0; j2 < N; j2++) {
	  int di = std::abs(i2 - i1);
	  int dj = std::abs(j2 - j1);
	  cache[i1][j1][i2][j2] = std::sqrt(di*di + dj*dj);
	}
      }
    }
  }
}

int main() {
  init_current();
  init_cache();
  // solve<0>(0);
  solve();
  return 0;
}
