#ifndef __EXTMAT_H__
#define __EXTMAT_H__

#include <R.h>

/* External matrix structure */
typedef void (*mulfn) (double* out, const double* v, const void* matrix);
typedef unsigned (*infofn) (const void* matrix);

typedef struct {
  const char* type;
  void* matrix;
  mulfn mulfn;
  mulfn tmulfn;
  infofn ncol;
  infofn nrow;
} ext_matrix;

typedef SEXP (*extmat_fn_t)(SEXP);

SEXP is_extmat(SEXP ptr);

#endif /* __EXTMAT_H__ */
