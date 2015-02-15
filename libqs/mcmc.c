#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "include/qsearch/libqs.h"

struct MCMCContext {
  const double *distmatrix;
  double beta;
  double total_weight;
  double nonmove_weight;
  double cutoff_weight;
  double last_score;
  uint64_t mutation_code;
};

static double scoreToWeight(double score, double beta) {
  double invprob = (1 - score) * beta;
  if (invprob < 0) { invprob = 0; }
  return exp(-invprob);
}


static int mutationAccumulator(const struct QSTree *tree,  const struct QSTree *nexttree, int sequence_number,
                       uint64_t mutation_code, void *obj) {
  struct MCMCContext *mcc = (struct MCMCContext *) obj;
  uint16_t *fullpathmatrix = qsNewFullPathMatrix(qsLeafCount(nexttree));
  uint16_t *pathmatrix = qsNewPathMatrix(qsLeafCount(nexttree));
  qstWritePathMatrix(fullpathmatrix, nexttree);
  qstWriteTruncatedPathMatrix(pathmatrix, fullpathmatrix);
  double score = qsScoreTree(nexttree, pathmatrix, mcc->distmatrix);
  mcc->last_score = score;
  qsFreePathMatrix(pathmatrix);
  qsFreeFullPathMatrix(fullpathmatrix);
  double weight = scoreToWeight(score, mcc->beta);
  mcc->total_weight += weight;
  if (mcc->cutoff_weight >= 0 && mcc->total_weight >= mcc->cutoff_weight) {
    mcc->mutation_code = mutation_code;
    return 1;
  }
  return 0;
}


#if 0
void qsApplyRandomMutation(struct QSTree *tree) {
  uint16_t *utree = (uint16_t *) tree;
  int counter = 0;
  uint64_t muta[2] = { 0, 0 };
  qstWritePathMatrix(fullpathmatrix, tree);
  qsIterateMutations(tree, fullpathmatrix, &counter, mutationCounter);
  muta[0] = rand() % counter; // TODO
#endif

double qsStepMCMC(struct QSTree *tree, const double *distmatrix, double beta) {
  struct MCMCContext *mcc = calloc(sizeof(struct MCMCContext), 1);
  uint16_t *fullpathmatrix = qsNewFullPathMatrix(qsLeafCount(tree));
  uint16_t *pathmatrix = qsNewPathMatrix(qsLeafCount(tree));
//  int leaf_count = qsLeafCount(tree);
//  int node_count = qsNodeCount(tree);
  qstWritePathMatrix(fullpathmatrix, tree);
  qstWriteTruncatedPathMatrix(pathmatrix, fullpathmatrix);
  double score = qsScoreTree(tree, pathmatrix, distmatrix);
  mcc->distmatrix = distmatrix;
  mcc->beta = beta;
  mcc->nonmove_weight = scoreToWeight(score, beta);
  mcc->total_weight = mcc->nonmove_weight;
  mcc->cutoff_weight = -1;
  qsIterateMutations(tree, fullpathmatrix, mcc, mutationAccumulator);
  double normf = (rand() % 1000000000) / 1000000000.0;
  mcc->cutoff_weight = normf * mcc->total_weight;
  mcc->total_weight = mcc->nonmove_weight;
  if (mcc->total_weight >= mcc->cutoff_weight) {
    return score;
  } else {
    qsIterateMutations(tree, fullpathmatrix, mcc, mutationAccumulator);
  }
  qsFreePathMatrix(pathmatrix);
  qsFreeFullPathMatrix(fullpathmatrix);
  score = mcc->last_score;
  free(mcc);
  return score;
}

static int areTreesEqual(struct QSTree **arr, int tree_count) {
  int i;
  for (i = 1; i < tree_count; ++i) {
    if (qsTreeCompare(arr[i], arr[0]) != 0) {
      return 0;
    }
  }
  return 1;
}

double qsSolveMCMC(struct QSTree **result, int leaf_count, const double *distmatrix) {
  struct QSTree *trees[10];
  int i;
  int tree_sizes[] = {5, 4, 4, 3, 3, 3};
  int tree_count = 2;
  if (leaf_count < 4) {
    fprintf(stderr, "Error, leaf_count must be at least 4.\n");
    exit(1);
  }
  if (leaf_count < 10) {
    tree_count = tree_sizes[leaf_count - 4];
  }
  for (i = 0; i < tree_count; ++i) {
    trees[i] = qsNewRandomTree(leaf_count);
  }
  int tree_pointer = 0;
  double score = 0;
  uint64_t itercount = 5;
  while (!areTreesEqual(trees, tree_count)) {
    itercount += 1;
    double lg = log(itercount);
    double beta = lg*lg*lg;
    tree_pointer = (tree_pointer + 1) % tree_count;
    score = qsStepMCMC(trees[tree_pointer], distmatrix, beta);
//    printf("score for %d = %f\n", tree_pointer, score);
    if (score == 1.0) {
      break;
    }
  }
  *result = qsNewCloneOf(trees[tree_pointer]);
  return score;
}
