#include <AtomFieldInteraction.h>
#include <stdlib.h>

struct BLAtomFieldInteraction* blAtomFieldInteractionCreate() {
  return 0;
}

void blAtomFieldInteractionDestroy(struct BLAtomFieldInteraction* atomFieldInteraction) {
  free(atomFieldInteraction);
}
