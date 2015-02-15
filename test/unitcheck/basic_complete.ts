/* A complete test example */

#include <stdio.h>
#include <qsearch.h>
#include <stdlib.h>
#include <math.h>

#define MAX_LEAVES_TEST 12

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
  ck_assert(tree[4] == 0xffff);
  QST_CONNECT_BOTH(uint16_t, tree, 0, 4);
  ck_assert(qstIsConnected(tree, 0, 4));
  ck_assert(qstIsConnected(tree, 4, 0));
  ck_assert(tree[4] == 0);
  QST_REMOVE_FROM_BOTH(uint16_t, tree, 0, 4);
  ck_assert( ! qstIsConnected(tree, 0, 4));
  ck_assert( ! qstIsConnected(tree, 4, 0));

#test qsearch_pathtest
  QST_DECLARE_TREE(uint16_t, tree, 4);
  QST_DECLARE_PATH_LENGTH(uint16_t, pathlen, 4);
  QST_DECLARE_TRUNCATED_PATH_LENGTH(uint16_t, smallpathlen, 4);
  ck_assert(sizeof(tree_store) == sizeof(tree_store[0])*11);
  ck_assert(tree_store == tree - 1);
  qsInitializeTree((struct QSTree *) tree, 4);
  ck_assert(tree[-1] == 4);
  qstMakeFixedStartingTree((struct QSTree *) tree);
  ck_assert(tree[-1] == 4);
  ck_assert(pathlen != NULL);
  ck_assert(pathlen[-1] == 6);
  ck_assert(tree[-1] == 4);
  ck_assert(sizeof(pathlen_back) == sizeof(pathlen[0])*37);
  ck_assert(smallpathlen[-1] == 4);
  ck_assert(pathlen[-1] == 6);
  qstWritePathMatrix(pathlen, tree);
  ck_assert(pathlen[-1] == 6);
  qstWriteTruncatedPathMatrix(smallpathlen, pathlen);
  ck_assert(pathlen[-1] == 6);
  ck_assert(smallpathlen[-1] == 4);
  ck_assert(smallpathlen[1*4+1] == 0);
  ck_assert(pathlen[0] == smallpathlen[0]);

#test qsearch_iteratequartets
static int counter = 0;
void qcfunc(uint16_t *tree, int i, int a, int b, int c, int d) {
  counter += 1;
}
  QST_DECLARE_TREE(uint16_t, tree, 5);
  qstIterateQuartetsForTree(tree, qcfunc, 0);
  ck_assert(counter == 5);

#test qsearch_newtree_test
  int leaf_count;
  for (leaf_count = 4; leaf_count < 10; ++leaf_count) {
    struct QSTree *tree = qsNewTree(leaf_count);
    ck_assert(tree != NULL);
    int zeroIsGood;
    zeroIsGood = qsVerifyTree(tree);
    if (zeroIsGood != 0) {
      fprintf(stderr, "Tree verification error N in newtree test.\n");
      qsPrintTree(tree);
      exit(1);
    }
    qsFreeTree(tree);
  }

#test qsearch_inittree_test
  int leaf_count;
  for (leaf_count = 4; leaf_count < 10; ++leaf_count) {
    struct QSTree *tree = (struct QSTree *) (((uint16_t *)malloc(qsTreeAllocationSize(leaf_count)))+1);
    qsInitializeTree(tree, leaf_count);
    ck_assert(tree != NULL);
    qsFreeTree(tree);
  }

#test qsearch_scoretree_test
  int leaf_count;
  QST_DECLARE_PATH_LENGTH(uint16_t, pathlen, 10);
  QST_DECLARE_TRUNCATED_PATH_LENGTH(uint16_t, smallpathlen, 10);
  for (leaf_count = 4; leaf_count < 10; ++leaf_count) {
    struct QSTree *tree = qsNewTree(leaf_count);
    ck_assert(tree != NULL);
    qstWritePathMatrix(pathlen, tree);
    qstWriteTruncatedPathMatrix(smallpathlen, pathlen);
    double *distmatrix = calloc(leaf_count * leaf_count , sizeof(double));
    int i, j;
    for (i = 0; i < leaf_count; ++i) {
      for (j = 0; j < leaf_count; ++j) {
        double min = (i < j ? i : j);
        double max = (i > j ? i : j);
        double sum = (i + j) * 0.17 + min * min * 0.3 + max * max * max * 0.01;
        double result = fabs(sin(sum));
        distmatrix[i*leaf_count + j] = result;
      }
      distmatrix[i*leaf_count + i] = 0;
    }
    for (i = 0; i < 10; i += 1) {
      double score = qsScoreTree(tree, smallpathlen, distmatrix);
      ck_assert(score >= 0);
      ck_assert(score <= 1);
    }
    qsFreeTree(tree);
    free(distmatrix);
  }

#test qsearch_normalize_test
  int leaf_count;
  for (leaf_count = 4; leaf_count < 10; ++leaf_count) {
    struct QSTree *tree = qsNewTree(leaf_count);
    ck_assert(tree != NULL);
    int retval = qsNormalizeTree(tree);
    ck_assert(retval == 0);
    qsFreeTree(tree);
  }

#test qsearch_pathfromto_test
  int leaf_count;
  QST_DECLARE_PATH_LENGTH(uint16_t, pathlen, 10);
  uint16_t path[10];
  for (leaf_count = 4; leaf_count < 10; ++leaf_count) {
    struct QSTree *tree = qsNewTree(leaf_count);
    ck_assert(tree != NULL);
    qstWritePathMatrix(pathlen, tree);
    int i, j, k;
    for (i = 0; i < 2*leaf_count-2; i += 1) {
      for (k = 0; k < 2*leaf_count-2; k += 1) {
        int src = i;
        int dest = k;
        int path_length = qsPathFromTo(tree, pathlen, src , dest, path);
        ck_assert(path[0] == src);
        ck_assert(path[path_length-1] == dest);
        for (j = 1; j < path_length; ++j) {
          ck_assert(qsIsConnected(tree, path[j-1], path[j]));
        }
      }
    }
    qsFreeTree(tree);
  }

#test qsearch_iteratemutations_test
int mutationHandler(const struct QSTree *tree, const struct QSTree *nexttree, int sequence_number,
                       uint64_t mutation_code, void *obj) {

  int zeroIsGood;
  zeroIsGood = qsVerifyTree(nexttree);
  if (zeroIsGood != 0) {
    fprintf(stderr, "Tree verification error A in mutation test.\n");
    qsPrintTree(nexttree);
    exit(1);
  }
  return 0;
}

  int leaf_count;
  QST_DECLARE_PATH_LENGTH(uint16_t, pathlen, MAX_LEAVES_TEST);
  for (leaf_count = 4; leaf_count < MAX_LEAVES_TEST; ++leaf_count) {
    struct QSTree *tree = qsNewTree(leaf_count);
    ck_assert(tree != NULL);
    qstWritePathMatrix(pathlen, tree);
    qsIterateMutations(tree, pathlen, pathlen, mutationHandler);
    qsFreeTree(tree);
  }

#test qsearch_randommutations_test
  int leaf_count, i;
  for (leaf_count = 4; leaf_count < MAX_LEAVES_TEST; ++leaf_count) {
    struct QSTree *tree = qsNewTree(leaf_count);
    for (i = 0; i < leaf_count * 10; ++i) {
      qsApplyRandomMutation(tree);
      int zeroIsGood = qsVerifyTree(tree);
      if (zeroIsGood != 0) {
        fprintf(stderr, "Tree verification error M in random mutation test.\n");
        qsPrintTree(tree);
        exit(1);
      }
    }
    qsFreeTree(tree);
  }

#test qsearch_newrandomtree_test
  int leaf_count;
  for (leaf_count = 4; leaf_count < MAX_LEAVES_TEST; ++leaf_count) {
    struct QSTree *tree = qsNewRandomTree(leaf_count);
    int zeroIsGood = qsVerifyTree(tree);
    if (zeroIsGood != 0) {
      fprintf(stderr, "Tree verification error in random new tree test.\n");
      qsPrintTree(tree);
      exit(1);
    }
    qsFreeTree(tree);
  }

#test qsearch_treehash_test
  int leaf_count;
  for (leaf_count = 4; leaf_count < MAX_LEAVES_TEST; ++leaf_count) {
    struct QSTree *tree = qsNewRandomTree(leaf_count);
    int zeroIsGood = qsVerifyTree(tree);
    if (zeroIsGood != 0) {
      fprintf(stderr, "Tree verification error in random new tree test.\n");
      qsPrintTree(tree);
      exit(1);
    }
    char hval[17];
    hval[0] = 5;
    hval[1] = 0;
    qsTreeHashHex(tree, hval);
    qsFreeTree(tree);
    ck_assert(strlen(hval) == 16);
  }

#test qsearch_mcmc_test
  int leaf_count;
  QST_DECLARE_PATH_LENGTH(uint16_t, pathlen, 10);
  QST_DECLARE_TRUNCATED_PATH_LENGTH(uint16_t, smallpathlen, 10);
  for (leaf_count = 4; leaf_count < 10; ++leaf_count) {
    struct QSTree *tree = qsNewTree(leaf_count);
    ck_assert(tree != NULL);
    qstWritePathMatrix(pathlen, tree);
    qstWriteTruncatedPathMatrix(smallpathlen, pathlen);
    double *distmatrix = calloc(leaf_count * leaf_count , sizeof(double));
    int i, j;
    for (i = 0; i < leaf_count; ++i) {
      for (j = 0; j < leaf_count; ++j) {
        double min = (i < j ? i : j);
        double max = (i > j ? i : j);
        double sum = (i + j) * 0.17 + min * min * 0.3 + max * max * max * 0.01;
        double result = fabs(sin(sum));
        distmatrix[i*leaf_count + j] = result;
      }
      distmatrix[i*leaf_count + i] = 0;
    }
    for (i = 0; i < 10; i += 1) {
      const double beta = 1.0;
      const double score = qsStepMCMC(tree, distmatrix, beta);
      ck_assert(score >= 0);
      ck_assert(score <= 1);
    }
    qsFreeTree(tree);
    free(distmatrix);
  }

#test qsearch_mcmcsolver_test
  int leaf_count;
  for (leaf_count = 4; leaf_count < 6; ++leaf_count) {
    struct QSTree *tree;
    double *distmatrix = calloc(leaf_count * leaf_count , sizeof(double));
    int i, j;
    for (i = 0; i < leaf_count; ++i) {
      for (j = 0; j < leaf_count; ++j) {
        double min = (i < j ? i : j);
        double max = (i > j ? i : j);
        double sum = (i + j) * 0.17 + min * min * 0.3 + max * max * max * 0.01;
        double result = fabs(sin(sum));
        distmatrix[i*leaf_count + j] = result;
      }
      distmatrix[i*leaf_count + i] = 0;
    }
    const double score = qsSolveMCMC(&tree, leaf_count, distmatrix);
    ck_assert(score == 1);
    qsFreeTree(tree);
    free(distmatrix);
  }

