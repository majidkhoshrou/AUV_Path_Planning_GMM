#ifndef __c17_simulation_h__
#define __c17_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc17_simulationInstanceStruct
#define typedef_SFc17_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c17_sfEvent;
  boolean_T c17_isStable;
  boolean_T c17_doneDoubleBufferReInit;
  uint8_T c17_is_active_c17_simulation;
  real_T c17_n;
  real_T c17_k_xi;
  real_T c17_i;
} SFc17_simulationInstanceStruct;

#endif                                 /*typedef_SFc17_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c17_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c17_simulation_get_check_sum(mxArray *plhs[]);
extern void c17_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
