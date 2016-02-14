/*
Copyright 2014 Dominic Meiser

This file is part of BeamLaser.

BeamLaser is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

BeamLaser is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with BeamLaser.  If not, see <http://www.gnu.org/licenses/>.
*/
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
void blAtomFieldInteractionTakeStep(struct BLAtomFieldInteraction *atomFieldInteraction,
                                    double dt,
                                    struct BLFieldState *fieldState,
                                    struct BLEnsemble *ensemble);

#endif
