#ifndef __c16_simulation_h__
#define __c16_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc16_simulationInstanceStruct
#define typedef_SFc16_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c16_sfEvent;
  boolean_T c16_isStable;
  boolean_T c16_doneDoubleBufferReInit;
  uint8_T c16_is_active_c16_simulation;
  real_T c16_k_u;
  real_T c16_k_r;
  real_T c16_m_u;
  real_T c16_m_v;
  real_T c16_m_uv;
  real_T c16_J;
  real_T c16_Xu;
  real_T c16_Xuu;
  real_T c16_Nr;
  real_T c16_Nrr;
} SFc16_simulationInstanceStruct;

#endif                                 /*typedef_SFc16_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c16_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c16_simulation_get_check_sum(mxArray *plhs[]);
extern void c16_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
