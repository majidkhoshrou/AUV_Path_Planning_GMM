#ifndef __c5_simulation_h__
#define __c5_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc5_simulationInstanceStruct
#define typedef_SFc5_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c5_sfEvent;
  boolean_T c5_isStable;
  boolean_T c5_doneDoubleBufferReInit;
  uint8_T c5_is_active_c5_simulation;
  real_T c5_k_u;
  real_T c5_k_r;
  real_T c5_m_u;
  real_T c5_m_v;
  real_T c5_m_uv;
  real_T c5_J;
  real_T c5_Xu;
  real_T c5_Xuu;
  real_T c5_Nr;
  real_T c5_Nrr;
} SFc5_simulationInstanceStruct;

#endif                                 /*typedef_SFc5_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c5_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c5_simulation_get_check_sum(mxArray *plhs[]);
extern void c5_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
