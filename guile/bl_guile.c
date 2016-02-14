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
ensemble_get_component(SCM component, SCM ensemble_smob)
{
  BL_UNUSED(component);
  BL_UNUSED(ensemble_smob);
  /*
   * This function requires a fair bit of code to create the arrays.
   * Would be better to create arrays on scheme side and to fill them on
   * the c side
   *
  struct ensemble *ensemble = (struct ensemble *) SCM_SMOB_DATA (ensemble_smob);

  SCM ra = scm_i_make_array(1);
  SCM_I_ARRAY_BASE(ra) = 0;
  scm_t_array_dim *s = SCM_I_ARRAY_DIMS (ra);
  s->lbnd = 0;
  s->ubnd = ensemble->ensemble->numPtcls;
  s->inc = 1;
  SCM a = scm_make_typed_array(scm_from_locale_string("f64"), 0, ra);
  */
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
}

