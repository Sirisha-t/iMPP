# include "unionFind.hpp"
void DisjointSet::deleteNode()
{
  std::unordered_map<node*, int> book;
  for(auto x : entry)
  {
    node* this_node = x.second;
    book[this_node] = 0;
    while(this_node->parent != this_node)
    {
      this_node = this_node->parent;
      book[this_node] = 0;
    }
  }
  for(auto x : book)
  {
    delete x.first;
  }
}

void DisjointSet::insert(int x)
{
  node* new_node = new node();
  new_node->parent = new_node;
  entry[x] = new_node;
}

void DisjointSet::unionNode(int x, int y)
{
  if(entry.count(x) == 0) {insert(x); unionNode(x, y);}
  if(entry.count(y) == 0) {insert(y); unionNode(x, y);}
  node* node1 = find(x);
  node* node2 = find(y);
  if(node1 != node2)  node2->parent = find(x);
}

node* DisjointSet::find(int x)
{
  if(entry.count(x) == 0)
  {
    insert(x);
    return entry[x];
  }
  node* this_node = entry[x];
  std::stack<node*> to_compress;
  while(this_node->parent != this_node)
  {
    to_compress.push(this_node);
    this_node = this_node->parent;
  }
  while(!to_compress.empty())
  {
    node* top = to_compress.top();
    to_compress.pop();
    top->parent = this_node;
  }
  return this_node;
}

void DisjointSet::visual()
{
  for(auto x : entry)
  {
    std::cout << x.first << ":" << x.second << "->" << x.second->parent << "\n";
  }
}
