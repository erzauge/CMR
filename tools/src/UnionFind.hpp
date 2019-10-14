#pragma once
#include <vector>

class UnionFind{
private:
  std::vector<long> node;
  long N;
public:
  UnionFind(long N_n);
  void add(long x);
  long Find(long x);
  void Union(long a, long b);
  std::vector<long> ClusterSize();
  long Label(long *label);
  long size();
  bool Singleton(long x);
};
