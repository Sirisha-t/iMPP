#ifndef UNIONFIND_H_
#define UNIONFIND_H_
# include <unordered_map>
# include <stack>
# include <iostream>
struct node
{
  node* parent;
  node() {parent = nullptr;};
};

class DisjointSet
{
private:
  std::unordered_map<int, node*> entry;
  void deleteNode();
public:
  DisjointSet() {};
  ~DisjointSet() {deleteNode(); };
  void insert(int x);
  node* find(int x);
  void unionNode(int x, int y);
  void visual();
};
#endif
