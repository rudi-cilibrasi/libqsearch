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
  double base_score;
  uint64_t mutation_code;
};

double qsStepHillClimb(struct QSTree *tree, const double *distmatrix);

static double scoreToWeight(double score_delta, double beta) {
  const double LOWEND = -2;
  double invprob = -score_delta * beta;
  if (invprob < LOWEND) { invprob = LOWEND; }
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
  double weight = scoreToWeight(score - mcc->base_score, mcc->beta);
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
  //mcc->nonmove_weight = scoreToWeight(0, beta);
  mcc->nonmove_weight = 0;
  mcc->base_score = score;
  mcc->total_weight = mcc->nonmove_weight;
  mcc->cutoff_weight = -1;
  qsIterateMutations(tree, fullpathmatrix, mcc, mutationAccumulator);
  double normf = (rand() % 1000000000) / 1000000000.0;
  mcc->cutoff_weight = normf * mcc->total_weight;
  mcc->total_weight = mcc->nonmove_weight;
  if (mcc->total_weight >= mcc->cutoff_weight) {
    mcc->last_score = score;
  } else {
    qsIterateMutations(tree, fullpathmatrix, mcc, mutationAccumulator);
    qsApplyMutation(tree, fullpathmatrix, mcc->mutation_code);
  }
  qsFreePathMatrix(pathmatrix);
  qsFreeFullPathMatrix(fullpathmatrix);
  score = mcc->last_score;
  free(mcc);
  return score;
}

static int areTreesEqual(struct QSTree **arr, int tree_count, double *scores) {
  int i;
  for (i = 1; i < tree_count; ++i) {
    if (scores[0] != scores[i]) {
//      printf("Found tree %d and %d not equal\n", 0, i);
//      qsPrintTree(arr[0]);
//      qsPrintTree(arr[i]);
      return 0;
    }
  }
//  printf("All trees equal.\n");
  return 1;
}

double qsSolveHillClimb(struct QSTree **result, int leaf_count, const double *distmatrix) {
  struct QSTree *trees[10];
  int i;
  int tree_sizes[] = { 5, 4, 3, 3, 3 };
  int tree_count = 2;
  if (leaf_count < 4) {
    fprintf(stderr, "Error, leaf_count must be at least 4.\n");
    exit(1);
  }
  if (leaf_count < 9) {
    tree_count = tree_sizes[leaf_count - 4];
  }
  int tree_pointer = 0;
  double score = 0, scores[10];
  uint64_t itercount = 0;
  int finished = 0;
  while (!finished) {
    for (i = 0; i < tree_count; ++i) {
      trees[i] = qsNewRandomTree(leaf_count);
    }
    for (i = 0; i < 10; ++i) {
      scores[i] = -1;
    }
    scores[0] = score;
    while (!areTreesEqual(trees, tree_count, scores)) {
      itercount += 1;
      tree_pointer = (tree_pointer + 1) % tree_count;
      score = qsStepHillClimb(trees[tree_pointer], distmatrix);
      scores[tree_pointer] = score;
  //    printf("iter=%d tp=%d score=%f beta=%f\n", (int) itercount, tree_pointer, score, beta);
  //    printf("score for %d = %f\n", tree_pointer, score);
      if (score == 1.0) {
        break;
      }
    }
    if (score == 1.0 || areTreesEqual(trees, tree_count, scores)) {
      *result = qsNewCloneOf(trees[tree_pointer]);
      finished = 1;
    }
    for (i = 0; i < tree_count; ++i) {
      qsFreeTree(trees[i]);
    }
  }
  return score;
}
double qsSolveMCMC(struct QSTree **result, int leaf_count, const double *distmatrix) {
  struct QSTree *trees[10];
  int i;
  int tree_sizes[] = { 5, 4, 3, 3, 3 };
  int tree_count = 2;
  if (leaf_count < 4) {
    fprintf(stderr, "Error, leaf_count must be at least 4.\n");
    exit(1);
  }
  if (leaf_count < 9) {
    tree_count = tree_sizes[leaf_count - 4];
  }
  int tree_pointer = 0;
  double score = 0, scores[10];
  uint64_t itercount = 5;
  uint64_t itertop = 1024;
  double betalev = 5, betastep = 1;
  int finished = 0;
  while (!finished) {
    for (i = 0; i < tree_count; ++i) {
      trees[i] = qsNewRandomTree(leaf_count);
    }
    for (i = 0; i < 10; ++i) {
      scores[i] = -1;
    }
    scores[0] = score;
    while (!areTreesEqual(trees, tree_count, scores)) {
      itercount += 1;
      if (itercount == itertop) {
        itertop = (itertop*3)/2; itercount = 5;
        //betalev += betastep;
        betalev = (betalev*5)/3;
//        betastep += 0.5;
        break;
      }
      double nrm = (double) itercount / (double) itertop;
      double beta = betalev * (1.0-exp(-8*nrm));
      tree_pointer = (tree_pointer + 1) % tree_count;
      score = qsStepMCMC(trees[tree_pointer], distmatrix, beta);
      scores[tree_pointer] = score;
  //    printf("iter=%d tp=%d score=%f beta=%f\n", (int) itercount, tree_pointer, score, beta);
  //    printf("score for %d = %f\n", tree_pointer, score);
      if (score == 1.0) {
        break;
      }
    }
    if (score == 1.0 || areTreesEqual(trees, tree_count, scores)) {
      *result = qsNewCloneOf(trees[tree_pointer]);
      finished = 1;
    }
    for (i = 0; i < tree_count; ++i) {
      qsFreeTree(trees[i]);
    }
  }
  return score;
}

double qsStepHillClimb(struct QSTree *tree, const double *distmatrix) {
  uint16_t *fullpathmatrix = qsNewFullPathMatrix(qsLeafCount(tree));
  uint16_t *pathmatrix = qsNewPathMatrix(qsLeafCount(tree));
  qstWritePathMatrix(fullpathmatrix, tree);
  qstWriteTruncatedPathMatrix(pathmatrix, fullpathmatrix);
  double score = qsScoreTree(tree, pathmatrix, distmatrix);
  double candscore;
  for (;;) {
    struct QSTree *candtree = qsNewCloneOf(tree);
    do {
      qsApplyRandomMutation(candtree);
    } while (rand() % 2 == 0);
    qstWritePathMatrix(fullpathmatrix, candtree);
    qstWriteTruncatedPathMatrix(pathmatrix, fullpathmatrix);
    candscore = qsScoreTree(candtree, pathmatrix, distmatrix);
    if (candscore > score) {
      qsCopyTreeOver(tree, candtree);
      qsFreeTree(candtree);
      qsFreePathMatrix(pathmatrix);
      qsFreeFullPathMatrix(fullpathmatrix);
      printf("candscore=%f\n", candscore);
      return candscore;
    }
    qsFreeTree(candtree);
  }
}
