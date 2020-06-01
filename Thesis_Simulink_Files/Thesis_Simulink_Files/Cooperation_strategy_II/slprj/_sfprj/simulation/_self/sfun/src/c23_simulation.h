#ifndef __c23_simulation_h__
#define __c23_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc23_simulationInstanceStruct
#define typedef_SFc23_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c23_sfEvent;
  boolean_T c23_isStable;
  boolean_T c23_doneDoubleBufferReInit;
  uint8_T c23_is_active_c23_simulation;
  real_T c23_n;
  real_T c23_k_xi;
  real_T c23_i;
} SFc23_simulationInstanceStruct;

#endif                                 /*typedef_SFc23_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c23_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c23_simulation_get_check_sum(mxArray *plhs[]);
extern void c23_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
