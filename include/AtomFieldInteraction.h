#ifndef ATOM_FIELD_INTERACTION_H
#define ATOM_FIELD_INTERACTION_H

#include <Ensemble.h>
#include <FieldState.h>
#include <ModeFunction.h>
#include <DipoleOperator.h>


struct BLAtomFieldInteraction;

struct BLAtomFieldInteraction* blAtomFieldInteractionCreate(int maxNumParticles,
    int internalStateSize, struct BLDipoleOperator *dipoleOperator,
    struct BLModeFunction *modeFunction);
void blAtomFieldInteractionDestroy(struct BLAtomFieldInteraction* atomFieldInteraction);
void blFieldAtomInteractionTakeStep(struct BLAtomFieldInteraction *atomFieldInteraction,
                                    double dt,
                                    struct BLFieldState *fieldState,
                                    struct BLEnsemble *ensemble);

#endif
