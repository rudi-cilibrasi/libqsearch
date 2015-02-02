#ifndef __LIBQS_H
#define __LIBQS_H

#include <stddef.h>
#include <string.h>

#define QST_NODE_COUNT(leaf_count) (4*leaf_count - 6)
#define QST_NODELIST_COUNT(leaf_count) (2*leaf_count - 2)
#define QST_PATH_LENGTH_COUNT(leaf_count) (QST_NODELIST_COUNT(leaf_count)*QST_NODELIST_COUNT(leaf_count))
#define QST_BYTE_SIZE(type_param, leaf_count) (sizeof(type_param) * (1 + QST_NODE_COUNT(leaf_count)))
#define QST_LEAF_COUNT(tree) (tree[-1])
#define QST_EMPTY_FLAG(type_param) (~((type_param)0))
#define QST_NLIST_SIZE(tree, listof) (listof < tree[-1] ? 1 : 3)
#define QST_NLIST_BASE(tree, listof) (listof < tree[-1] ? listof : tree[-1] + 3*(listof-tree[-1]))
#define QST_REMOVE_FROM(type_param, tree, listof, target)                      \
do {                                                                           \
           uint32_t qnb = QST_NLIST_BASE(tree, listof);                        \
           if (listof < tree[-1]) {                                            \
               tree[qnb] = QST_EMPTY_FLAG(type_param);                         \
               break;                                                          \
           }                                                                   \
           if (tree[qnb+0] == target) {                                        \
               tree[qnb+0]  = QST_EMPTY_FLAG(type_param);                      \
               break;                                                          \
           }                                                                   \
           if (tree[qnb+1] == target) {                                        \
               tree[qnb+1]  = QST_EMPTY_FLAG(type_param);                      \
           } else {                                                            \
               tree[qnb+2]  = QST_EMPTY_FLAG(type_param);                      \
           }                                                                   \
} while (0)

#define QST_REMOVE_FROM_BOTH(type_param, tree, listof, target)              \
do {                                                                        \
   QST_REMOVE_FROM(type_param, tree, listof, target);                       \
   QST_REMOVE_FROM(type_param, tree, target, listof);                       \
} while (0)

#define QST_CONNECT_ONE(type_param, tree, listof, target)                       \
do {                                                                            \
           uint32_t qnb = QST_NLIST_BASE(tree, listof);                         \
           if (listof < tree[-1]) {                                             \
             tree[qnb] = target;                                                \
             break;                                                             \
           }                                                                    \
           if (tree[qnb+0]==QST_EMPTY_FLAG(type_param)){                        \
             tree[qnb+0] = target;                                              \
             break;                                                             \
           }                                                                    \
           if (tree[qnb+1]==QST_EMPTY_FLAG(type_param)){                        \
             tree[qnb+1] = target;                                              \
             break;                                                             \
           }                                                                    \
           tree[qnb+2] = target;                                                \
} while(0)

#define QST_CONNECT_BOTH(type_param, tree, listof, target)                  \
   do {                                                                     \
     QST_CONNECT_ONE(type_param, tree, listof, target);                     \
     QST_CONNECT_ONE(type_param, tree, target, listof); }                   \
   while(0)

#define QST_DECLARE_PATH_LENGTH(type_param, path_name, leaf_count)          \
  type_param path_name ## _back[QST_PATH_LENGTH_COUNT(leaf_count) + 1] = { QST_NODELIST_COUNT(leaf_count) },           \
  *path_name = path_name ## _back + 1

#define QST_DECLARE_TREE(type_param, tree_name, leaf_count)                 \
  type_param tree_name ## _store[QST_NODE_COUNT(leaf_count) + 1] = { leaf_count },           \
  *tree_name = tree_name ## _store + 1

#define QST_DISINTEGRATE_TREE(tree) memset(tree, 0xff, QST_NODE_COUNT(tree[-1])*sizeof(tree[0]))

#define qstIsConnected(tree, i, j) (i == j ? 0 :                              \
          i > j ? qstIsConnectedB(tree, j, i) : qstIsConnectedB(tree, i, j))

#define qstIsConnectedB(tree, i, j)                                        \
   (        i < tree[-1] ? tree[i] == j :                                  \
    tree[QST_NLIST_BASE(tree, i)+0] == j ||                                \
    tree[QST_NLIST_BASE(tree, i)+1] == j ||                                \
    tree[QST_NLIST_BASE(tree, i)+2] == j)

#define qstWritePathMatrix(path, tree)          do {                       \
  path[-1] = QST_PATH_LENGTH_COUNT(tree[-1]);                              \
  uint32_t n = path[-1];                                                   \
  uint32_t i, j, k;                                                        \
  for (i = 0; i < n; i += 1) {                                             \
    for (j = 0; j < n; j += 1) {                                           \
      path[i*n + j] = n;                                                   \
    }                                                                      \
    path[i*n + i] = 0;                                                     \
  }                                                                        \
  for (i = 0; i < n; i += 1) {                                             \
    uint32_t neighbors_start  = QST_NLIST_BASE(tree, i);                   \
    uint32_t neighbors_length = QST_NLIST_SIZE(tree, i);                   \
    for (j = 0; j < neighbors_length; j += 1) {                            \
      path[i*n + tree[neighbors_start + j]] = 1;                           \
    }                                                                      \
  }                                                                        \
  for (k = 0; k < n; k += 1) {                                             \
    for (i = 0; i < n; i += 1) {                                           \
      for (j = 0; j < n; j += 1) {                                         \
        if (path[i*n+j] > path[i*n+k] + path[k*n+j]) {                     \
          path[i*n+j] = path[i*n+k] + path[k*n+j];                         \
        }                                                                  \
      }                                                                    \
    }                                                                      \
  }                                                                        \
}  while (0)

#define qstMakeFixedStartingTree(tree)          do {                       \
  int howManyLeaves = tree[-1];                                            \
  memset(tree, 0xff, 2*QST_NODE_COUNT(howManyLeaves));                     \
  uint32_t i, jTop = howManyLeaves - 2;                                    \
  for (i = 0; i < jTop; i += 1) {                                          \
    QST_CONNECT_BOTH(uint16_t, tree, i, howManyLeaves + i);                \
    if (i > 0) {                                                           \
      QST_CONNECT_BOTH(uint16_t, tree, i+howManyLeaves-1,howManyLeaves+i); \
    }                                                                      \
  }                                                                        \
  QST_CONNECT_BOTH(uint16_t, tree, howManyLeaves-2, howManyLeaves);        \
  QST_CONNECT_BOTH(uint16_t, tree, howManyLeaves-1, 2*howManyLeaves-3);    \
}  while (0)

#endif
