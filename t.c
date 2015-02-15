#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <qsearch.h>

int main(int argc, char **argv)
{
  int leaf_count = 15;
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
  const double score = qsSolveHillClimb(&tree, leaf_count, distmatrix);
  printf("Got score = %f with leaf_count = %d\n", score, leaf_count);
  qsPrintTree(tree);
  qsFreeTree(tree);
  free(distmatrix);
  return 0;
}

