#ifndef __c27_simulation_h__
#define __c27_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc27_simulationInstanceStruct
#define typedef_SFc27_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c27_sfEvent;
  boolean_T c27_isStable;
  boolean_T c27_doneDoubleBufferReInit;
  uint8_T c27_is_active_c27_simulation;
} SFc27_simulationInstanceStruct;

#endif                                 /*typedef_SFc27_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c27_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c27_simulation_get_check_sum(mxArray *plhs[]);
extern void c27_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
