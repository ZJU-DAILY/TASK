
//#pragma GCC optimize(2)

#include "mtrie.h"

MTrieNode::~MTrieNode() {
  if (child) { delete child; child = NULL; }
  if (next) { delete next; next = NULL; }
}

MTrieNode* MTrieNode::insertChild(TrieNode* node, pair<int, int> ed) {
  MTrieNode *cnode = child;
  while (cnode) {
    if (node == cnode->tnode) break;
    cnode = cnode->next;
  }
  if (cnode == NULL) {
    cnode = new MTrieNode();
    cnode->parent = this;
    cnode->tnode = node;
    cnode->next = child;
    child = cnode;
  }
  cnode->eds.push_back(ed);
  return cnode;
}
