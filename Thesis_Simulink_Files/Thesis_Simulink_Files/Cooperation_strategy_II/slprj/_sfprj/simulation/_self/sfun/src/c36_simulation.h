#ifndef __c36_simulation_h__
#define __c36_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc36_simulationInstanceStruct
#define typedef_SFc36_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c36_sfEvent;
  boolean_T c36_isStable;
  boolean_T c36_doneDoubleBufferReInit;
  uint8_T c36_is_active_c36_simulation;
} SFc36_simulationInstanceStruct;

#endif                                 /*typedef_SFc36_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c36_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c36_simulation_get_check_sum(mxArray *plhs[]);
extern void c36_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
