#ifndef __c2_simulation_h__
#define __c2_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc2_simulationInstanceStruct
#define typedef_SFc2_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c2_sfEvent;
  boolean_T c2_isStable;
  boolean_T c2_doneDoubleBufferReInit;
  uint8_T c2_is_active_c2_simulation;
  real_T c2_m_u;
  real_T c2_m_v;
  real_T c2_m_uv;
  real_T c2_J;
  real_T c2_l;
  real_T c2_Xu;
  real_T c2_Xuu;
  real_T c2_Yv;
  real_T c2_Yvv;
  real_T c2_Nr;
  real_T c2_Nrr;
} SFc2_simulationInstanceStruct;

#endif                                 /*typedef_SFc2_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c2_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c2_simulation_get_check_sum(mxArray *plhs[]);
extern void c2_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
