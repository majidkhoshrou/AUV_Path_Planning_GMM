#ifndef __c22_simulation_h__
#define __c22_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc22_simulationInstanceStruct
#define typedef_SFc22_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c22_sfEvent;
  boolean_T c22_isStable;
  boolean_T c22_doneDoubleBufferReInit;
  uint8_T c22_is_active_c22_simulation;
} SFc22_simulationInstanceStruct;

#endif                                 /*typedef_SFc22_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c22_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c22_simulation_get_check_sum(mxArray *plhs[]);
extern void c22_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
