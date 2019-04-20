#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>

void rotate(int N, int K) {
  for (int write = 0; write < N;) {
    for (int read = K; read < N; read++) {
      if (K == write) {
        K = read;
      }
      std::cout << write << ' ' << read << std::endl;
      write++;
    }
  }
}

void rotate_arr(int* arr, int N, int K) {
  for (int write = 0; write < N;) {
    for (int read = K; read < N; read++) {
      if (K == write) {
        K = read;
      }
      std::swap(arr[write], arr[read]);
      write++;
    }
  }
}

void verify() {
  for (int i = 0; i < 100; i++) {
    int arr1[i];
    int arr2[i];
    for (int j = 0; j < i; j++) {
      std::iota(arr1, arr1 + i, 0);
      std::iota(arr2, arr2 + i, 0);
      rotate_arr(arr1, i, j);
      std::rotate(arr2, arr2 + j, arr2 + i);
      for (int k = 0; k < i; k++) {
        assert(arr1[k] == arr2[k]);
      }
    }
  }
}

int main() {
  verify();
  return 0;
}
