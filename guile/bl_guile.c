#include <libguile.h>
#include <Ensemble.h>
#include <math.h>

static scm_t_bits ensemble_tag;

struct ensemble {
  struct BLEnsemble *ensemble;
  SCM name;
};

static SCM
ensemble_create(SCM name, SCM max_num_ptcls, SCM internal_state_size)
{
  SCM smob;
  struct ensemble *ensemble = scm_gc_malloc(sizeof(*ensemble), "ensemble");
  ensemble->ensemble =
    scm_gc_malloc(sizeof(*ensemble->ensemble), "ensemble->ensemble");
  blEnsembleCreate(scm_to_int(max_num_ptcls), scm_to_int(internal_state_size),
      ensemble->ensemble);
  SCM_NEWSMOB(smob, ensemble_tag, ensemble);
  ensemble->name = name;
  return smob;
}

static SCM
ensemble_push(SCM dt, SCM ensemble_smob)
{
  struct ensemble *ensemble = (struct ensemble *) SCM_SMOB_DATA (ensemble_smob);
  blEnsemblePush(scm_to_double(dt), ensemble->ensemble);
  return ensemble_smob;
}

static SCM
ensemble_create_space(SCM numParticles, SCM ensemble_smob)
{
  struct ensemble *ensemble = (struct ensemble *) SCM_SMOB_DATA (ensemble_smob);
  blEnsemblePush(scm_to_int(numParticles), ensemble->ensemble);
  return ensemble_smob;
}

static SCM
ensemble_get_num_ptcls(SCM ensemble_smob)
{
  struct ensemble *ensemble = (struct ensemble *) SCM_SMOB_DATA (ensemble_smob);
  return scm_from_int(ensemble->ensemble->numPtcls);
}

static SCM
ensemble_set_num_ptcls(SCM n, SCM ensemble_smob)
{
  struct ensemble *ensemble = (struct ensemble *) SCM_SMOB_DATA (ensemble_smob);
  int numPtcls = scm_to_int(n);
  if (numPtcls < 0 || numPtcls >= ensemble->ensemble->maxNumPtcls) {
    scm_out_of_range("ensemble_set_num_ptcls", n);
  }
  ensemble->ensemble->numPtcls = scm_to_int(n);
  return ensemble_smob;
}

static SCM
ensemble_get_component(SCM component, SCM ensemble_smob, SCM array)
{
  struct ensemble *ensemble = (struct ensemble *) SCM_SMOB_DATA (ensemble_smob);

  int comp = scm_to_int(component);
  if (comp < 0 || comp > 5) {
    scm_out_of_range("ensemble_get_component", component);
  }
  scm_t_array_handle handle;
  scm_array_get_handle(array, &handle);
  double *a = scm_array_handle_f64_writable_elements(&handle);
  switch (comp) {
  case 0:
    memcpy(a, ensemble->ensemble->x, ensemble->ensemble->numPtcls * sizeof(*a));
    break;
  case 1:
    memcpy(a, ensemble->ensemble->y, ensemble->ensemble->numPtcls * sizeof(*a));
    break;
  case 2:
    memcpy(a, ensemble->ensemble->z, ensemble->ensemble->numPtcls * sizeof(*a));
    break;
  case 3:
    memcpy(a, ensemble->ensemble->vx, ensemble->ensemble->numPtcls * sizeof(*a));
    break;
  case 4:
    memcpy(a, ensemble->ensemble->vy, ensemble->ensemble->numPtcls * sizeof(*a));
    break;
  case 5:
    memcpy(a, ensemble->ensemble->vz, ensemble->ensemble->numPtcls * sizeof(*a));
    break;
  }
  scm_array_handle_release(&handle);
  return 0;
}

static SCM
ensemble_mark(SCM ensemble_smob)
{
  struct ensemble *ensemble = (struct ensemble *) SCM_SMOB_DATA (ensemble_smob);
  return ensemble->name;
}

static size_t
ensemble_free(SCM ensemble_smob)
{
  struct ensemble *ensemble = (struct ensemble *)SCM_SMOB_DATA(ensemble_smob);

  blEnsembleDestroy(ensemble->ensemble);
  scm_gc_free(ensemble->ensemble, sizeof(*ensemble->ensemble),
      "ensemble->ensemble");
  scm_gc_free(ensemble, sizeof(*ensemble), "ensemble");

  return 0;
}

static int
ensemble_print(SCM ensemble_smob, SCM port, scm_print_state *pstate)
{
  BL_UNUSED(pstate);
  struct ensemble *ensemble = (struct ensemble *) SCM_SMOB_DATA (ensemble_smob);
  scm_puts ("#<ensemble ", port);
  scm_display (ensemble->name, port);
  scm_puts (" ", port);
  char address_buffer[100];
  sprintf(address_buffer, "%p", ensemble->ensemble);
  scm_puts(address_buffer, port);
  scm_puts (" >", port);
  return 1;
}

void init_ensemble() {
  ensemble_tag = scm_make_smob_type("ensemble", sizeof(struct ensemble));
  scm_set_smob_print(ensemble_tag, ensemble_print);
  scm_set_smob_mark(ensemble_tag, ensemble_mark);
  scm_set_smob_free(ensemble_tag, ensemble_free);

  scm_c_define_gsubr("ensemble-create", 3, 0, 0, ensemble_create);
  scm_c_define_gsubr("ensemble-push", 2, 0, 0, ensemble_push);
  scm_c_define_gsubr("ensemble-create-space", 2, 0, 0, ensemble_create_space);
  scm_c_define_gsubr("ensemble-get-component", 2, 0, 0, ensemble_get_component);
  scm_c_define_gsubr("ensemble-get-num-ptcls", 1, 0, 0, ensemble_get_num_ptcls);
  scm_c_define_gsubr("ensemble-set-num-ptcls", 2, 0, 0, ensemble_set_num_ptcls);
}

