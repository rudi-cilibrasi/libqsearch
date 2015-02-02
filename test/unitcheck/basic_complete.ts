/* A complete test example */

#include <stdio.h>
#include <qsearch.h>

#test qsearch_test
  QST_DECLARE_TREE(uint16_t, tree, 4);
  ck_assert(tree != NULL);
  ck_assert(tree[-1] == 4);
  QST_DISINTEGRATE_TREE(tree);
  int i, j;
  for (i = 0; i < QST_LEAF_COUNT(tree); ++i) {
    for (j = 0; j < QST_LEAF_COUNT(tree); ++j) {
      ck_assert( ! qstIsConnected(tree, i, j));
    }
  }
  ck_assert( ! qstIsConnected(tree, 0, 4));
  ck_assert( ! qstIsConnected(tree, 4, 0));
  QST_CONNECT_BOTH(uint16_t, tree, 0, 4);
  ck_assert(qstIsConnected(tree, 0, 4));
  ck_assert(qstIsConnected(tree, 4, 0));
  QST_REMOVE_FROM_BOTH(uint16_t, tree, 0, 4);
  ck_assert( ! qstIsConnected(tree, 0, 4));
  ck_assert( ! qstIsConnected(tree, 4, 0));

#test qsearch_pathtest
  QST_DECLARE_TREE(uint16_t, tree, 4);
  QST_DECLARE_PATH_LENGTH(uint16_t, pathlen, 4);
  qstMakeFixedStartingTree(tree);
  ck_assert(pathlen != NULL);
  ck_assert(pathlen[-1] == 6);
  ck_assert(tree[-1] == 4);
  ck_assert(sizeof(pathlen_back) == sizeof(pathlen[0])*37);
  qstWritePathMatrix(pathlen, tree);

