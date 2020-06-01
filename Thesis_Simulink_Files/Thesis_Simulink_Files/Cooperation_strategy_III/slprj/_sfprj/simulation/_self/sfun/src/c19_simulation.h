#ifndef __c19_simulation_h__
#define __c19_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc19_simulationInstanceStruct
#define typedef_SFc19_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c19_sfEvent;
  boolean_T c19_isStable;
  boolean_T c19_doneDoubleBufferReInit;
  uint8_T c19_is_active_c19_simulation;
  real_T c19_m_u;
  real_T c19_m_v;
  real_T c19_m_uv;
  real_T c19_J;
  real_T c19_l;
  real_T c19_Xu;
  real_T c19_Xuu;
  real_T c19_Yv;
  real_T c19_Yvv;
  real_T c19_Nr;
  real_T c19_Nrr;
} SFc19_simulationInstanceStruct;

#endif                                 /*typedef_SFc19_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c19_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c19_simulation_get_check_sum(mxArray *plhs[]);
extern void c19_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
