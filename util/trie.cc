//#pragma GCC optimize(2)

#include "trie.h"

TrieNode::~TrieNode() {
    if (child) {
        delete child;
        child = NULL;
    }
    if (next) {
        delete next;
        next = NULL;
    }
}

void TrieNode::getRecords(vector<int> &recs) {
    queue<TrieNode *> q;
    q.push(this);
    while (!q.empty()) {
        TrieNode *node = q.front();
        q.pop();
        if (!node->rids.empty())
            recs.insert(recs.end(), node->rids.begin(), node->rids.end());
        TrieNode *cnode = node->child;
        while (cnode) {
            q.push(cnode);
            cnode = cnode->next;
        }
    }
}

int TrieNode::preorder(int id, Trie *trie) {
    this->id = id;
    if (!rids.empty()) {  // if the node corresponds to a string in the dataset, add the rid and node id information into the trie
        for (auto it = rids.begin(); it != rids.end(); it++)
            trie->ids.push_back(make_pair(id, *it));
    }
    TrieNode *node = this->child;
    while (node) {
        id = id + 1;
        id = node->preorder(id, trie);
        node = node->next;
    }
    this->last = id;
    return id;
}

int Trie::buildIdx() {
    return root->preorder(1, this) + 1;
}

void Trie::bfsIndex() {  // level traversal
    queue<TrieNode *> q;
    q.push(root);
    while (!q.empty()) {
        TrieNode *node = q.front();
        q.pop();
        index[node->key][node->depth].push_back(node);
        TrieNode *cnode = node->child;
        while (cnode) {
            q.push(cnode);
            cnode = cnode->next;
        }
    }
}

/* build length index*/
/*
for (auto it = index.begin(); it != index.end(); it++) {
  auto first = it->second.begin();
  int current = (*first)->depth;
  auto low = lower_bound(first, it->second.end(), current + 1, comp);
  int prev = 0;
  while (true) {
    index_depth[it->first].push_back(Triple(current, prev, distance(first, low)));
    if (low == it->second.end()) break;
    prev = distance(first, low);
    current = (*low)->depth;
    low = lower_bound(low, it->second.end(), current + 1, comp);
  }
}
*/

TrieNode *TrieNode::insertChild(char ch) {
    TrieNode *node = child;
    while (node) {
        if (ch == node->key) break;
        node = node->next;
    }
    if (node == NULL) {
        node = new TrieNode();
        node->parent = this;
        node->key = ch;
        node->depth = depth + 1;
        node->next = child;
        child = node;
    }
    return node;
}

TrieNode *Trie::append(const char *str, const int rid) {
    TrieNode *node = root;
    while (*str) {
        node = node->insertChild(*str);
        ++str;
    }
    node->rids.push_back(rid);
    return node;
}

int Trie::appendNewKey(string keyword, vector<TrieNode*> &path) {
    int insertPos = 1;
    TrieNode *node = root;
    path.push_back(node);
    for (int i = 0; i < keyword.size(); i++) {
        TrieNode *node1 = node->child;
        TrieNode *node2 = node->child;
        while (node2) {
            if (keyword[i] == node2->key)
                break;
            node1 = node2;
            node2 = node2->next;
        }
        if (node2 == NULL) {
            node2 = new TrieNode();
            node2->parent = node;
            node2->key = keyword[i];
            node2->depth = node->depth + 1;
            if (node1 != NULL)
                node1->next = node2;
            else
                node->child = node2;
        } else {
            insertPos++;
        }
        node = node2;
        path.push_back(node);
    }
    return insertPos;
}

TrieNode *Trie::match(string str) {
    TrieNode *activeNode = root;
    for (int i = 0; i < str.length(); i++) {
        if (activeNode == NULL)
            break;
        vector<TrieNode *> canNodes = index[str[i]][i + 1];
        bool isMatch = false;  // whether the corresponding node is matched
        for (int j = 0; j < canNodes.size(); j++) {
            if (canNodes[j]->id >= activeNode->id && canNodes[j]->id <= activeNode->last) {
                activeNode = canNodes[j];
                isMatch = true;
                break;
            }
        }
        if (!isMatch)
            activeNode = NULL;
    }
    return activeNode;
}

TrieNode *Trie::match(TrieNode *activeNode, char ch) {
    if (activeNode == NULL)
        return NULL;
    vector<TrieNode*> canNodes = index[ch][activeNode->depth + 1];
    bool isMatch = false;
    for (int j = 0; j < canNodes.size(); j++) {
        if (canNodes[j]->id >= activeNode->id && canNodes[j]->id <= activeNode->last) {
            activeNode = canNodes[j];
            isMatch = true;
            break;
        }
    }
    if (!isMatch)
        activeNode = NULL;
    return activeNode;
}
