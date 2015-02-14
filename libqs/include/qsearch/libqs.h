#ifndef __LIBQS_H
#define __LIBQS_H

#include <stdint.h>
#include <string.h>

struct QSTDirectedEdge {
  uint32_t from, to;
};

struct QSTUInt64Table;

struct QSTUInt64Table *qsNewUInt64Table(void);
void qsAddUInt64ToTable(struct QSTUInt64Table *hashtab, uint64_t val);
int qsIsUInt64InTable(const struct QSTUInt64Table *hashtab, uint64_t val);
void qsFreeUInt64Table(struct QSTUInt64Table *hashtab);

struct QSTree;

uint32_t qsTreeAllocationSize(uint32_t leaf_count);
uint32_t qsInitializeTree(struct QSTree *tree, uint32_t leaf_count);
struct QSTree *qsNewTree(uint32_t leaf_count);
struct QSTree *qsNewCloneOf(const struct QSTree *orig);
struct QSTree *qsNewRandomTree(uint32_t leaf_count);
void qsFreeTree(struct QSTree *tree);
uint32_t qsCopyTreeOver(struct QSTree *destination, const struct QSTree *source);
uint32_t qsLeafCount(const struct QSTree *tree);
uint32_t qsNodeCount(const struct QSTree *tree);
uint32_t qsIsConnected(const struct QSTree *tree, uint32_t a, uint32_t b);
uint32_t qsVerifyTree(const struct QSTree *tree);
double qsScoreTree(const struct QSTree *tree, const uint16_t *pathmatrix,
 const double *distmatrix);   // same size as tree leaf count squared
uint32_t qsNormalizeTree(struct QSTree *tree);
int qsTreeCompare(const struct QSTree *tree_a, const struct QSTree *tree_b);
int qsPathFromTo(const struct QSTree *tree, const uint16_t *fullpathmatrix, int a, int b, uint16_t *path_buffer);
void qsPrintTree(const struct QSTree *tree);
void qsIterateMutations(const struct QSTree *tree,
                        const uint16_t *fullpathmatrix,
                        void *obj,
  int (*mutationHandler)(const struct QSTree *tree, const struct QSTree *nexttree, int sequence_number,
                         uint64_t mutation_code, void *obj));
void qsApplyMutation(struct QSTree *tree,
                        const uint16_t *fullpathmatrix,
                        uint64_t mutation_code);
void qsApplyRandomMutation(struct QSTree *tree);
void qsScrambleTree(struct QSTree *tree);
uint64_t qsTreeHash(const struct QSTree *tree);
uint64_t qsTreeHashHex(const struct QSTree *tree, char hval[17]);

#define QST_NODE_COUNT(leaf_count) (4*leaf_count - 6)
#define QST_NODELIST_COUNT(leaf_count) (2*leaf_count - 2)
#define QST_PATH_LENGTH_COUNT(leaf_count) (QST_NODELIST_COUNT(leaf_count)*QST_NODELIST_COUNT(leaf_count))
#define QST_BYTE_SIZE(type_param, leaf_count) (sizeof(type_param) * (1 + QST_NODE_COUNT(leaf_count)))
#define QST_LEAF_COUNT(tree) (tree[-1])
#define QST_EMPTY_FLAG(type_param) ((type_param) ~((type_param)0))
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

#define QST_DECLARE_TRUNCATED_PATH_LENGTH(type_param, path_name, leaf_count) \
  type_param path_name ## _back[leaf_count * leaf_count + 1] = { leaf_count },            \
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

#define qstIterateQuartetsForTree(tree, qfunc, obj)                          \
  do {                                                                       \
    uint32_t leaf_count = ((uint16_t *)tree)[-1];                            \
    uint32_t a, b, c, d;                                                     \
    for (a = 0; a < leaf_count; a += 1) {                                    \
      for (b = a + 1; b < leaf_count; b += 1) {                              \
        for (c = b + 1; c < leaf_count; c += 1) {                            \
          for (d = c + 1; d < leaf_count; d += 1) {                          \
            qfunc(tree, obj, a, b, c, d);                                    \
          }                                                                  \
        }                                                                    \
      }                                                                      \
    }                                                                        \
  } while(0)

#define qstWriteTruncatedPathMatrix(smallpath, path)  do {                 \
  uint32_t sn = smallpath[-1], n = path[-1];                               \
  uint32_t i;                                                              \
  for (i = 0; i < sn; i += 1) {                                            \
    memcpy(&smallpath[i * sn], &path[i * n], sn * sizeof(smallpath[0]));   \
  }                                                                        \
} while(0)

#define qstWritePathMatrix(path, utree)          do {                      \
  const uint16_t *__ytree = (const uint16_t *) utree;                      \
  if (__ytree[-1] < 4 || __ytree[-1] > 16000) {                            \
    fprintf(stderr, "Error, bad __ytree length: %d\n", __ytree[-1]);       \
    exit(1);                                                               \
  }                                                                        \
  path[-1] = QST_NODELIST_COUNT(__ytree[-1]);                              \
  uint32_t n = path[-1];                                                   \
  uint32_t i, j, k;                                                        \
  for (i = 0; i < n; i += 1) {                                             \
    for (j = 0; j < n; j += 1) {                                           \
      path[i*n + j] = n;                                                   \
    }                                                                      \
    path[i*n + i] = 0;                                                     \
  }                                                                        \
  for (i = 0; i < n; i += 1) {                                             \
    uint32_t neighbors_start  = QST_NLIST_BASE(__ytree, i);                \
    uint32_t neighbors_length = QST_NLIST_SIZE(__ytree, i);                \
    for (j = 0; j < neighbors_length; j += 1) {                            \
      path[i*n + __ytree[neighbors_start + j]] = 1;                        \
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

#define qstMakeFixedStartingTree(utree)          do {                      \
  uint16_t *__ztree = (uint16_t *) utree;                                  \
  int howManyLeaves = __ztree[-1];                                         \
  if (__ztree[-1] < 4 || __ztree[-1] > 16000) {                            \
    fprintf(stderr, "Error, bad tree length: %d\n", __ztree[-1]);          \
    exit(1);                                                               \
  }                                                                        \
  memset(__ztree, 0xff, 2*QST_NODE_COUNT(howManyLeaves));                  \
  uint32_t i, jTop = howManyLeaves - 2;                                    \
  for (i = 0; i < jTop; i += 1) {                                          \
    QST_CONNECT_BOTH(uint16_t, __ztree, i, howManyLeaves + i);             \
    if (i > 0) {                                                           \
      QST_CONNECT_BOTH(uint16_t, __ztree, i+howManyLeaves-1,howManyLeaves+i); \
    }                                                                      \
  }                                                                        \
  QST_CONNECT_BOTH(uint16_t, __ztree, howManyLeaves-2, howManyLeaves);     \
  QST_CONNECT_BOTH(uint16_t, __ztree, howManyLeaves-1, 2*howManyLeaves-3); \
  qsNormalizeTree(utree);                                                  \
}  while (0)


#define qsFullPathLengthByteSize(leafcount) (sizeof(uint16_t)*(QST_PATH_LENGTH_COUNT(leafcount)+1))
#define qsNewFullPathMatrix(leafcount) (((uint16_t *) calloc(1, qsFullPathLengthByteSize(leafcount))) + 1)
#define qsFreeFullPathMatrix(ptr) free(((uint16_t *) ptr)-1)

#endif
