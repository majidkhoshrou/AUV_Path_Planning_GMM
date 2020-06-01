#ifndef __c15_simulation_h__
#define __c15_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc15_simulationInstanceStruct
#define typedef_SFc15_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c15_sfEvent;
  boolean_T c15_isStable;
  boolean_T c15_doneDoubleBufferReInit;
  uint8_T c15_is_active_c15_simulation;
  real_T c15_nsonar;
  real_T c15_i;
  real_T c15_kl[81];
} SFc15_simulationInstanceStruct;

#endif                                 /*typedef_SFc15_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c15_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c15_simulation_get_check_sum(mxArray *plhs[]);
extern void c15_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
