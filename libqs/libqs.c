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
  qsNormalizeTree(tree);
  return 0;
}

struct QSTree *qsNewRandomTree(uint32_t leaf_count) {
  struct QSTree *tree = qsNewTree(leaf_count);
  int i;
  for (i = 0; i < 10 * leaf_count; ++i) {
    qsApplyRandomMutation(tree);
  }
  return tree;
}

struct QSTree *qsNewCloneOf(const struct QSTree *orig) {
  const uint16_t *tr = (uint16_t *) orig;
  int len;
  uint16_t *ntr = calloc(len=QST_BYTE_SIZE(uint16_t, tr[-1]), 1);
  memcpy(ntr, tr-1, len);
  return (struct QSTree *) (ntr+1);
}

struct QSTree *qsNewTree(uint32_t leaf_count) {
  uint16_t *tr = calloc(QST_BYTE_SIZE(uint16_t, leaf_count), 1);
  struct QSTree *result;
  verifyLeafCount(leaf_count);
  tr += 1;
  result = (struct QSTree *) tr;
  tr[-1] = leaf_count;
  qstMakeFixedStartingTree(result);
  qsNormalizeTree(result);
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
  uint16_t histo[32768];
  uint16_t *tr = (uint16_t *) tree;
  int leaf_count = tr[-1];
  memset(histo, 0, sizeof(histo[0]) * (2 * leaf_count - 2));
  int i;
  for (i = 0; i < QST_NODE_COUNT(leaf_count); ++i) {
    int cur = tr[i];
    histo[cur] += 1;
  }
  for (i = 0; i < leaf_count; ++i) {
    if (histo[i] != 1) {
      printf("bad leaf %d has count %d.\n", i, histo[i]);
      return 1;
    }
  }
  for (i = 0; i < leaf_count - 2; ++i) {
    if (histo[i + leaf_count] != 3) {
      printf("bad kernel %d has count %d.\n", i + leaf_count, histo[i+leaf_count]);
      return 1;
    }
  }
  return 0;
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
    int consistent =
     (pathmatrix[q[0]*leaf_count + q[1]] + pathmatrix[q[2]*leaf_count + q[3]] <
        pathmatrix[q[0]*leaf_count + q[2]] + pathmatrix[q[1]*leaf_count + q[3]]);
    if ( consistent ) {
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

int qsPathFromTo(const struct QSTree *tree, const uint16_t *fullpathmatrix, int a, int b, uint16_t *path_buffer) {
  int path_length = 0;
  uint16_t *utree = (uint16_t *) tree;
  int pwidth = qsNodeCount(tree);
  while (a != b) {
    int nlist = QST_NLIST_BASE(utree, a);
    int nsize = QST_NLIST_SIZE(utree, a);
    int dist[3];
    int i;
    for (i = 0; i < nsize; ++i) {
      int neighbor = utree[nlist+i];
      dist[i] = fullpathmatrix[pwidth*b+neighbor];
    }
    int mindist = dist[0];
    int mini = 0;
    if (nsize > 1) {
      if (dist[1] < mindist) { mindist = dist[1]; mini = 1; }
      if (dist[2] < mindist) { mindist = dist[2]; mini = 2; }
      if (dist[0] == mindist && dist[1] == mindist && dist[2] == mindist) {
        fprintf(stderr, "Error in path length matrix.\n");
        exit(1);
      }
    }
    path_buffer[0] = a;
    path_buffer += 1;
    a = utree[nlist+mini];
    path_length += 1;
  }
  path_buffer[0] = a;
  path_length += 1;
  return path_length;
}

void qsPrintTree(const struct QSTree *tree) {
  uint16_t *utree = (uint16_t *) tree;
  int i, j;
  for (i = 0; i < 2*utree[-1] - 2; i += 1) {
    printf("[%d, [", i);
    int nlist = QST_NLIST_BASE(utree, i);
    int nsize = QST_NLIST_SIZE(utree, i);
    for (j = 0; j < nsize; ++j) {
      int neighbor = utree[nlist+j];
      printf("%d",neighbor);
      if (j != nsize - 1) {
        printf(", ");
      }
    }
    printf("]]");
    if (i != 2*utree[-1]-1) {
      printf(",");
    }
    printf("\n");
  }
}

static int isNewMutation(const struct QSTree *tree,  const uint16_t *fullpathmatrix, uint64_t mut, struct QSTUInt64Table *old_trees, struct QSTree **holder) {
  int result;
  if (*holder)
    qsFreeTree(*holder);
  *holder = qsNewCloneOf(tree);
  qsApplyMutation(*holder, fullpathmatrix, mut);
  uint64_t hval = qsTreeHash(*holder);
  result = qsIsUInt64InTable(old_trees, hval) ? 0 : 1;
  if (result) {
    qsAddUInt64ToTable(old_trees, hval);
  }
  return result;
}

void qsIterateMutations(const struct QSTree *tree,
                        const uint16_t *fullpathmatrix,
                        void *obj,
  int (*mutationHandler)(const struct QSTree *tree, const struct QSTree *nexttree,  int sequence_number,
                         uint64_t mutation_code, void *obj)) {
  int leaf_count = qsLeafCount(tree);
  int node_count = qsNodeCount(tree);
  int kern_count = leaf_count - 2;
  int i, j, terminate;
  uint16_t mut[4];
  struct QSTUInt64Table *old_trees = qsNewUInt64Table();
  uint16_t *utree = (uint16_t *) tree;
  uint64_t *m64 = (uint64_t *) &mut[0];
  mut[0] = 0;
  mut[3] = 0;
  int seqno = 0;
  struct QSTree *holder = NULL;
  qsAddUInt64ToTable(old_trees, qsTreeHash(tree));
  for (i = 0; i < leaf_count; ++i) {
    mut[1] = i;
    for (j = 0; j < leaf_count; ++j) {
      if (i == j) { continue; }
      if (fullpathmatrix[i*node_count + j] <= 2) { continue; }
      mut[2] = j;
      if (isNewMutation(tree, fullpathmatrix, *m64, old_trees, &holder)) {
        terminate = mutationHandler(tree, holder, seqno, *m64, obj);
        if (terminate) { goto done; }
        seqno++;
      }
    }
  }
  mut[0] = 1;
  uint16_t path_buffer[16383];
  for (i = 0; i < node_count; ++i) {
    int jpre;
    mut[1] = i;
    for (jpre = 0; jpre < kern_count; ++jpre) {
      int j = jpre + leaf_count;
      if (i == j) { continue; }
      if (fullpathmatrix[i*node_count + j] <= 2) { continue; }
      int nlist = QST_NLIST_BASE(utree, j);
      int nsize = QST_NLIST_SIZE(utree, j);
      mut[2] = j;
      int pathlen = qsPathFromTo(tree, fullpathmatrix, i, j, path_buffer);
      int mi;
      for (mi = 0; mi < nsize; mi++) {
        int m3 = utree[nlist + mi];
        if (m3 == path_buffer[pathlen-2] || m3 == path_buffer[1])
          continue;
        mut[3] = m3;
        if (isNewMutation(tree, fullpathmatrix, *m64, old_trees, &holder)) {
          terminate = mutationHandler(tree, holder, seqno, *m64, obj);
          if (terminate) { goto done; }
          seqno++;
        }
      }
    }
  }
  mut[0] = 2;
  int ipre;
  for (ipre = 0; ipre < kern_count; ++ipre) {
    int i = ipre + leaf_count;
    int jpre;
    mut[1] = i;
    for (jpre = ipre; jpre < kern_count; ++jpre) {
      int j = jpre + leaf_count;
      if (i == j) { continue; }
      if (fullpathmatrix[i*node_count + j] <= 2) { continue; }
      mut[2] = j;
      mut[3] = 0;
      if (isNewMutation(tree, fullpathmatrix, *m64, old_trees, &holder)) {
        terminate = mutationHandler(tree, holder, seqno, *m64, obj);
        if (terminate) { goto done; }
        seqno++;
      }
    }
  }
  done:
    qsFreeUInt64Table(old_trees);
    qsFreeTree(holder);
    return;
}

static int mutationCounter(const struct QSTree *tree, const struct QSTree *nexttree, int sequence_number,
                       uint64_t mutation_code, void *obj) {
  int *iptr = (int *) obj;
  (*iptr) += 1;
  return 0;
}

static int mutationExtractor(const struct QSTree *tree,  const struct QSTree *nexttree, int sequence_number,
                       uint64_t mutation_code, void *obj) {
  uint64_t *mptr = (uint64_t *) obj;
  if (sequence_number == mptr[0]) {
    mptr[1] = mutation_code;
  }
  return 0;
}


void qsApplyRandomMutation(struct QSTree *tree) {
  uint16_t *utree = (uint16_t *) tree;
  int counter = 0;
  uint64_t muta[2] = { 0, 0 };
  uint16_t *fullpathmatrix = qsNewFullPathMatrix(utree[-1]);
  qstWritePathMatrix(fullpathmatrix, tree);
  qsIterateMutations(tree, fullpathmatrix, &counter, mutationCounter);
  muta[0] = rand() % counter; // TODO
  qsIterateMutations(tree, fullpathmatrix, muta, mutationExtractor);
  if (muta[1] != 0) {
    qsApplyMutation(tree, fullpathmatrix, muta[1]);
  } else {
    fprintf(stderr, "Mutation logic error 001.\n");
    exit(1);
  }
  qsFreeFullPathMatrix(fullpathmatrix);
}

void qsApplyMutation(struct QSTree *tree,
                        const uint16_t *fullpathmatrix,
                        uint64_t mutation_code_64) {
  uint16_t *utree = (uint16_t *) tree;
  uint16_t *mutation_code = (uint16_t *) &mutation_code_64;
  int mcode = mutation_code[0];
  if (mcode == 0) { // leaf swap
    int i = mutation_code[1];
    int j = mutation_code[2];
    int ni = utree[i];
    int nj = utree[j];
    QST_REMOVE_FROM_BOTH(uint16_t, utree, i, ni);
    QST_REMOVE_FROM_BOTH(uint16_t, utree, j, nj);
    QST_CONNECT_BOTH(uint16_t, utree, i, nj);
    QST_CONNECT_BOTH(uint16_t, utree, j, ni);
    qsNormalizeTree(tree);
    return;
  }
  if (mcode == 1) { // subtree transfer
    uint16_t path_buffer[16383];
    int k1 = mutation_code[1];
    int k2 = mutation_code[2];
    int m3 = mutation_code[3];
    qsPathFromTo(tree, fullpathmatrix, k1, k2, path_buffer);
    int i1 = path_buffer[1];
    int nlist = QST_NLIST_BASE(utree, i1);
    QST_REMOVE_FROM_BOTH(uint16_t, utree, k1, i1);
    int ms[3], mc=0, mo;
    for (mo = 0; mo < 3; mo++) {
      int mu = utree[nlist+mo];
      if (mu != 0xffff) {
        ms[mc++] = mu;
      }
    }
    int m1 = ms[0];
    int m2 = ms[1];
    QST_REMOVE_FROM_BOTH(uint16_t, utree, m1, i1);
    QST_REMOVE_FROM_BOTH(uint16_t, utree, m2, i1);
    QST_REMOVE_FROM_BOTH(uint16_t, utree, m3, k2);
    QST_CONNECT_BOTH(uint16_t, utree, m1, m2);
    QST_CONNECT_BOTH(uint16_t, utree, k2, i1);
    QST_CONNECT_BOTH(uint16_t, utree, m3, i1);
    QST_CONNECT_BOTH(uint16_t, utree, k1, i1);
    qsNormalizeTree(tree);
    return;
  }
  if (mcode == 2) { // subtree interchange
    uint16_t path_buffer[16383];
    int k1 = mutation_code[1];
    int k2 = mutation_code[2];
    int pathlen = qsPathFromTo(tree, fullpathmatrix, k1, k2, path_buffer);
    int n1 = path_buffer[1];
    int n2 = path_buffer[pathlen-2];
    QST_REMOVE_FROM_BOTH(uint16_t, utree, n1, k1);
    QST_REMOVE_FROM_BOTH(uint16_t, utree, n2, k2);
    QST_CONNECT_BOTH(uint16_t, utree, n1, k2);
    QST_CONNECT_BOTH(uint16_t, utree, n2, k1);
    qsNormalizeTree(tree);
    return;
  }
  fprintf(stderr, "Error, bad mutation code.\n");
  exit(1);
}

static __inline__ void fnv64Init(uint64_t *hval) {
  *hval = 14695981039346656037ULL;
}

static __inline__ void fnv64UpdateChar(uint64_t *hval, unsigned char ch) {
  *hval ^= ch;
  *hval *= 1099511628211U;
}

static __inline__ void fnv64UpdateBuffer(uint64_t *hval, const void *buf,
                                         uint64_t   len) {
  uint64_t i;
  for (i = 0; i < len; ++i) {
    unsigned char ch = ((unsigned char *) buf)[i];
    fnv64UpdateChar(hval, ch);
  }
}

static void genericResultHex64(char *result, int howBig, const uint64_t *hval) {
  int i;
  const static char *hex = "0123456789abcdef";
  uint64_t c = *hval;
  for (i = howBig-1; i >= 0; --i) {
    uint64_t n = c & 0xff;
    c >>= 8;
    result[2*i] = hex[n >> 4];
    result[2*i+1] = hex[n & 0x0f];
  }
  result[2*howBig] = '\0';
}


uint64_t qsTreeHash(const struct QSTree *tree) {
  uint16_t *tr = (uint16_t *) tree;
  int len = qsTreeAllocationSize(tr[-1]);
  tr -= 1;
  uint64_t hval;
  fnv64Init(&hval);
  fnv64UpdateBuffer(&hval, tr, len);
  return hval;
}

uint64_t qsTreeHashHex(const struct QSTree *tree, char hval[17]) {
  uint64_t result = qsTreeHash(tree);
  genericResultHex64(hval, 8, &result);
  return result;
}

struct QSTUInt64Node {
  struct QSTUInt64Node *next;
  uint64_t val;
};

struct QSTUInt64Table {
  uint32_t hashtab_size;
  struct QSTUInt64Node *hashtab;
};

struct QSTUInt64Table *qsNewUInt64Table(void) {
  struct QSTUInt64Table *hashtab = calloc(sizeof(struct QSTUInt64Table), 1);
  hashtab->hashtab_size = 111191;
  hashtab->hashtab = calloc(sizeof(struct QSTUInt64Node), hashtab->hashtab_size);
  return hashtab;
}

void qsAddUInt64ToTable(struct QSTUInt64Table *hashtab, uint64_t val) {
  uint64_t reduced = val % hashtab->hashtab_size;
  if (hashtab->hashtab[reduced].val == 0) {
    hashtab->hashtab[reduced].val = val;
    return;
  }
  struct QSTUInt64Node *cur = calloc(sizeof(struct QSTUInt64Node), 1);
  cur->val = val;
  cur->next = hashtab->hashtab[reduced].next;
  hashtab->hashtab[reduced].next = cur;
}

int qsIsUInt64InTable(const struct QSTUInt64Table *hashtab, uint64_t val) {
  uint64_t reduced = val % hashtab->hashtab_size;
  struct QSTUInt64Node *cur = &hashtab->hashtab[reduced];
  while (cur != NULL) {
    if (val == cur->val) {
      return 1;
    }
    cur = cur->next;
  }
  return 0;
}

void qsFreeUInt64Table(struct QSTUInt64Table *hashtab) {
  int i;
  for (i = 0; i < hashtab->hashtab_size; ++i) {
    struct QSTUInt64Node *cur = hashtab->hashtab[i].next;
    while (cur != NULL) {
      struct QSTUInt64Node *next_cur = cur->next;
      free(cur);
      cur = next_cur;
    }
  }
  free(hashtab->hashtab);
  free(hashtab);
}
