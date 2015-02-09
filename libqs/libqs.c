#include <stdlib.h>
#include <stdio.h>
#include "include/qsearch/libqs.h"

uint32_t qsTreeAllocationSize(uint32_t leaf_count) {
  return QST_BYTE_SIZE(uint16_t, leaf_count);
}

static void verifyLeafCount(uint32_t leaf_count) {
  if (leaf_count < 4 || leaf_count > 16000) {
    fprintf(stderr, "error:  must have 4 <= leaf_count <= 16000\n");
    exit(1);
  }
}

uint32_t qsInitializeTree(struct QSTree *tree, uint32_t leaf_count) {
  uint16_t *tr = (uint16_t *) tree;
  verifyLeafCount(leaf_count);
  tr[-1] = leaf_count;
  qstMakeFixedStartingTree(tree);
  return 0;
}

struct QSTree *qsNewTree(uint32_t leaf_count) {
  uint16_t *tr = calloc(QST_BYTE_SIZE(uint16_t, leaf_count), 1);
  struct QSTree *result;
  verifyLeafCount(leaf_count);
  tr += 1;
  result = (struct QSTree *) tr;
  tr[-1] = leaf_count;
  qstMakeFixedStartingTree(result);
  return result;
}

void qsFreeTree(struct QSTree *tree) {
  uint16_t *tr = (uint16_t *) tree;
  free(tr-1);
}

uint32_t qsCopyTreeOver(struct QSTree *destination, const struct QSTree *source) {
  memcpy(destination, source, QST_BYTE_SIZE(uint16_t, ((uint16_t *)source)[-1]));
  return 0;
}

uint32_t qsLeafCount(const struct QSTree *tree) {
  uint16_t *tr = (uint16_t *) tree;
  return tr[-1];
}

uint32_t qsNodeCount(const struct QSTree *tree) {
  uint16_t *tr = (uint16_t *) tree;
  return QST_NODELIST_COUNT(tr[-1]);
}

uint32_t qsIsConnected(const struct QSTree *tree, uint32_t a, uint32_t b) {
  uint16_t *tr = (uint16_t *) tree;
  return qstIsConnected(tr, a, b);
}

uint32_t qsVerifyTree(const struct QSTree *tree) {
  return 1;
}

uint32_t qsNormalizeTree(struct QSTree *tree) {
  int i, j, k;
  uint16_t *tr = (uint16_t *) tree;
  uint16_t tmp;
  int leaf_count = tr[-1];
  for (i = 0; i < leaf_count - 2; i += 1) {
    uint16_t *ind = &tr[leaf_count + i * 3];
    for (j = 0; j < 2; j += 1) {
      for (k = 0; k < 2; k += 1) {
        if (ind[k] > ind[k+1]) { tmp=ind[k]; ind[k]=ind[k+1]; ind[k+1]=tmp; }
      }
    }
  }
  return 0;
}


struct QSTScoreContext {
  const struct QSTree *tree;
  const double *distmatrix;   // same size as tree leaf count squared
  const uint16_t *pathmatrix; // truncated with -1 holding length
  double totmin, totmax, totcur;
};

static void qsScoreFuncPrivate(const struct QSTree *tree, struct QSTScoreContext *param, int a, int b, int c, int d) {
  double scores[3];
  int topos[3][4] = { { a, b, c, d }, { a, c, b, d }, { a, d, b, c } };
  int i;
  int leaf_count = ((uint16_t *) tree)[-1];
  for (i = 0; i < 3; ++i) {
    scores[i] = param->distmatrix[topos[i][0]*leaf_count+topos[i][1]] +
                param->distmatrix[topos[i][2]*leaf_count+topos[i][3]];
  }
  double minScore = scores[0], maxScore = scores[0];
  double s = scores[1];
  const uint16_t *pathmatrix = param->pathmatrix;
  if (s < minScore) { minScore = s; }
  if (s > maxScore) { maxScore = s; }
  s = scores[2];
  if (s < minScore) { minScore = s; }
  if (s > maxScore) { maxScore = s; }
  for (i = 0; i < 3; ++i) {
    int *q = &topos[i][0];
    if (pathmatrix[q[0]*leaf_count + q[1]] + pathmatrix[q[2]*leaf_count + q[3]] <
        pathmatrix[q[0]*leaf_count + q[2]] + pathmatrix[q[1]*leaf_count + q[3]]) {
      param->totcur += scores[i];
      param->totmin += minScore;
      param->totmax += maxScore;
    }
  }
}

double qsScoreTree(const struct QSTree *tree, const uint16_t *pathmatrix,
 const double *distmatrix) {
  struct QSTScoreContext sc;
  sc.tree = tree; sc.pathmatrix = pathmatrix; sc.distmatrix = distmatrix;
  sc.totmin = 0.0; sc.totmax = 0.0; sc.totcur = 0.0;
  qstIterateQuartetsForTree(tree, qsScoreFuncPrivate, &sc);
  double score = 1.0 - ((sc.totcur - sc.totmin) / (sc.totmax - sc.totmin));
  return score;
}

int qsTreeCompare(const struct QSTree *tree_a, const struct QSTree *tree_b) {
  uint16_t *tr_a = (uint16_t *) tree_a;
  uint16_t *tr_b = (uint16_t *) tree_b;
  int retval = memcmp(tr_a-1, tr_b-1, 2);
  if (retval != 0) {
    return retval;
  }
  return memcmp(tr_a, tr_b, qsTreeAllocationSize(tr_a[-1]) - 2);
}
