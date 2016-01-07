#include <DipoleOperator.h>
#include <Utilities.h>
#include <stdlib.h>

void blDipoleOperatorDestroy(struct DipoleOperator *op) {
  op->destroy(op->ctx);
  free(op);
}

static void dipoleOperatorTLADestroy(void *ctx) {
  BL_UNUSED(ctx);
}

struct DipoleOperator *blDipoleOperatorTLACreate() {
  struct DipoleOperator *op = malloc(sizeof(*op));
  op->apply = 0;
  op->destroy = dipoleOperatorTLADestroy;
  op->ctx = 0;
  return op;
}
