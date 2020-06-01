/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c30_simulation.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simulation_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c30_debug_family_names[9] = { "sonarmax", "i", "l", "nargin",
  "nargout", "sonar", "map", "n", "dist" };

/* Function Declarations */
static void initialize_c30_simulation(SFc30_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c30_simulation(SFc30_simulationInstanceStruct
  *chartInstance);
static void enable_c30_simulation(SFc30_simulationInstanceStruct *chartInstance);
static void disable_c30_simulation(SFc30_simulationInstanceStruct *chartInstance);
static void c30_update_debugger_state_c30_simulation
  (SFc30_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c30_simulation
  (SFc30_simulationInstanceStruct *chartInstance);
static void set_sim_state_c30_simulation(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_st);
static void finalize_c30_simulation(SFc30_simulationInstanceStruct
  *chartInstance);
static void sf_c30_simulation(SFc30_simulationInstanceStruct *chartInstance);
static void initSimStructsc30_simulation(SFc30_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c30_machineNumber, uint32_T
  c30_chartNumber);
static const mxArray *c30_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData);
static void c30_emlrt_marshallIn(SFc30_simulationInstanceStruct *chartInstance,
  const mxArray *c30_dist, const char_T *c30_identifier, real_T c30_y[4]);
static void c30_b_emlrt_marshallIn(SFc30_simulationInstanceStruct *chartInstance,
  const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId, real_T c30_y[4]);
static void c30_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c30_mxArrayInData, const char_T *c30_varName, void *c30_outData);
static const mxArray *c30_b_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData);
static real_T c30_c_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId);
static void c30_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c30_mxArrayInData, const char_T *c30_varName, void *c30_outData);
static const mxArray *c30_c_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData);
static const mxArray *c30_d_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData);
static const mxArray *c30_e_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData);
static boolean_T c30_d_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId);
static void c30_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c30_mxArrayInData, const char_T *c30_varName, void *c30_outData);
static void c30_info_helper(const mxArray **c30_info);
static const mxArray *c30_emlrt_marshallOut(char * c30_u);
static const mxArray *c30_b_emlrt_marshallOut(uint32_T c30_u);
static real_T c30_mpower(SFc30_simulationInstanceStruct *chartInstance, real_T
  c30_a);
static void c30_eml_scalar_eg(SFc30_simulationInstanceStruct *chartInstance);
static void c30_eml_error(SFc30_simulationInstanceStruct *chartInstance);
static const mxArray *c30_f_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData);
static int32_T c30_e_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId);
static void c30_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c30_mxArrayInData, const char_T *c30_varName, void *c30_outData);
static uint8_T c30_f_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_b_is_active_c30_simulation, const char_T
  *c30_identifier);
static uint8_T c30_g_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId);
static void init_dsm_address_info(SFc30_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c30_simulation(SFc30_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c30_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c30_is_active_c30_simulation = 0U;
}

static void initialize_params_c30_simulation(SFc30_simulationInstanceStruct
  *chartInstance)
{
  real_T c30_d0;
  sf_set_error_prefix_string(
    "Error evaluating data 'n' in the parent workspace.\n");
  sf_mex_import_named("n", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c30_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c30_n = c30_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c30_simulation(SFc30_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c30_simulation(SFc30_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c30_update_debugger_state_c30_simulation
  (SFc30_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c30_simulation
  (SFc30_simulationInstanceStruct *chartInstance)
{
  const mxArray *c30_st;
  const mxArray *c30_y = NULL;
  int32_T c30_i0;
  real_T c30_u[4];
  const mxArray *c30_b_y = NULL;
  uint8_T c30_hoistedGlobal;
  uint8_T c30_b_u;
  const mxArray *c30_c_y = NULL;
  real_T (*c30_dist)[4];
  c30_dist = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c30_st = NULL;
  c30_st = NULL;
  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_createcellarray(2), FALSE);
  for (c30_i0 = 0; c30_i0 < 4; c30_i0++) {
    c30_u[c30_i0] = (*c30_dist)[c30_i0];
  }

  c30_b_y = NULL;
  sf_mex_assign(&c30_b_y, sf_mex_create("y", c30_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_setcell(c30_y, 0, c30_b_y);
  c30_hoistedGlobal = chartInstance->c30_is_active_c30_simulation;
  c30_b_u = c30_hoistedGlobal;
  c30_c_y = NULL;
  sf_mex_assign(&c30_c_y, sf_mex_create("y", &c30_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c30_y, 1, c30_c_y);
  sf_mex_assign(&c30_st, c30_y, FALSE);
  return c30_st;
}

static void set_sim_state_c30_simulation(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_st)
{
  const mxArray *c30_u;
  real_T c30_dv0[4];
  int32_T c30_i1;
  real_T (*c30_dist)[4];
  c30_dist = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c30_doneDoubleBufferReInit = TRUE;
  c30_u = sf_mex_dup(c30_st);
  c30_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c30_u, 0)),
                       "dist", c30_dv0);
  for (c30_i1 = 0; c30_i1 < 4; c30_i1++) {
    (*c30_dist)[c30_i1] = c30_dv0[c30_i1];
  }

  chartInstance->c30_is_active_c30_simulation = c30_f_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c30_u, 1)),
     "is_active_c30_simulation");
  sf_mex_destroy(&c30_u);
  c30_update_debugger_state_c30_simulation(chartInstance);
  sf_mex_destroy(&c30_st);
}

static void finalize_c30_simulation(SFc30_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c30_simulation(SFc30_simulationInstanceStruct *chartInstance)
{
  int32_T c30_i2;
  int32_T c30_i3;
  int32_T c30_i4;
  real_T c30_hoistedGlobal;
  int32_T c30_i5;
  boolean_T c30_sonar[4];
  int32_T c30_i6;
  real_T c30_map[12];
  real_T c30_b_n;
  uint32_T c30_debug_family_var_map[9];
  boolean_T c30_sonarmax;
  real_T c30_i;
  real_T c30_l;
  real_T c30_nargin = 3.0;
  real_T c30_nargout = 1.0;
  real_T c30_dist[4];
  const mxArray *c30_y = NULL;
  int32_T c30_i7;
  boolean_T c30_varargin_1[4];
  boolean_T c30_mtmp;
  int32_T c30_itmp;
  int32_T c30_ix;
  int32_T c30_b_ix;
  boolean_T c30_a;
  boolean_T c30_b;
  boolean_T c30_p;
  boolean_T c30_b_mtmp;
  int32_T c30_b_itmp;
  boolean_T c30_extremum;
  int32_T c30_iindx;
  boolean_T c30_maxval;
  int32_T c30_b_iindx;
  real_T c30_indx;
  boolean_T c30_b_sonarmax;
  real_T c30_b_i;
  int32_T c30_i8;
  real_T c30_c_n;
  int32_T c30_i9;
  int32_T c30_b_l;
  real_T c30_x;
  real_T c30_b_x;
  int32_T c30_i10;
  real_T (*c30_b_dist)[4];
  real_T (*c30_b_map)[12];
  boolean_T (*c30_b_sonar)[4];
  c30_b_map = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 1);
  c30_b_dist = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
  c30_b_sonar = (boolean_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 29U, chartInstance->c30_sfEvent);
  for (c30_i2 = 0; c30_i2 < 4; c30_i2++) {
    _SFD_DATA_RANGE_CHECK((real_T)(*c30_b_sonar)[c30_i2], 0U);
  }

  for (c30_i3 = 0; c30_i3 < 4; c30_i3++) {
    _SFD_DATA_RANGE_CHECK((*c30_b_dist)[c30_i3], 1U);
  }

  for (c30_i4 = 0; c30_i4 < 12; c30_i4++) {
    _SFD_DATA_RANGE_CHECK((*c30_b_map)[c30_i4], 2U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c30_n, 3U);
  chartInstance->c30_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 29U, chartInstance->c30_sfEvent);
  c30_hoistedGlobal = chartInstance->c30_n;
  for (c30_i5 = 0; c30_i5 < 4; c30_i5++) {
    c30_sonar[c30_i5] = (*c30_b_sonar)[c30_i5];
  }

  for (c30_i6 = 0; c30_i6 < 12; c30_i6++) {
    c30_map[c30_i6] = (*c30_b_map)[c30_i6];
  }

  c30_b_n = c30_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 9U, 9U, c30_debug_family_names,
    c30_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c30_sonarmax, 0U, c30_e_sf_marshallOut,
    c30_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c30_i, 1U, c30_b_sf_marshallOut,
    c30_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c30_l, 2U, c30_b_sf_marshallOut,
    c30_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c30_nargin, 3U, c30_b_sf_marshallOut,
    c30_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c30_nargout, 4U, c30_b_sf_marshallOut,
    c30_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c30_sonar, 5U, c30_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c30_map, 6U, c30_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c30_b_n, 7U, c30_b_sf_marshallOut,
    c30_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c30_dist, 8U, c30_sf_marshallOut,
    c30_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c30_sfEvent, 3);
  if (c30_b_n < 20.0) {
  } else {
    c30_y = NULL;
    sf_mex_assign(&c30_y, sf_mex_create("y", "Assertion failed.", 15, 0U, 0U, 0U,
      2, 1, strlen("Assertion failed.")), FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, c30_y);
  }

  _SFD_EML_CALL(0U, chartInstance->c30_sfEvent, 4);
  for (c30_i7 = 0; c30_i7 < 4; c30_i7++) {
    c30_varargin_1[c30_i7] = c30_sonar[c30_i7];
  }

  c30_mtmp = c30_varargin_1[0];
  c30_itmp = 1;
  for (c30_ix = 2; c30_ix < 5; c30_ix++) {
    c30_b_ix = c30_ix;
    c30_a = c30_varargin_1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c30_b_ix), 1, 4, 1, 0) - 1];
    c30_b = c30_mtmp;
    c30_p = ((int32_T)c30_a > (int32_T)c30_b);
    if (c30_p) {
      c30_mtmp = c30_varargin_1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c30_b_ix), 1, 4, 1, 0) - 1];
      c30_itmp = c30_b_ix;
    }
  }

  c30_b_mtmp = c30_mtmp;
  c30_b_itmp = c30_itmp;
  c30_extremum = c30_b_mtmp;
  c30_iindx = c30_b_itmp;
  c30_maxval = c30_extremum;
  c30_b_iindx = c30_iindx;
  c30_indx = (real_T)c30_b_iindx;
  c30_b_sonarmax = c30_maxval;
  c30_b_i = c30_indx;
  c30_sonarmax = c30_b_sonarmax;
  c30_i = c30_b_i;
  _SFD_EML_CALL(0U, chartInstance->c30_sfEvent, 5);
  for (c30_i8 = 0; c30_i8 < 4; c30_i8++) {
    c30_dist[c30_i8] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c30_sfEvent, 6);
  if (CV_EML_IF(0, 1, 0, (real_T)c30_sonarmax > 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c30_sfEvent, 7);
    c30_c_n = c30_b_n;
    c30_i9 = (int32_T)c30_c_n;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c30_c_n, mxDOUBLE_CLASS, c30_i9);
    c30_l = 1.0;
    c30_b_l = 0;
    while (c30_b_l <= c30_i9 - 1) {
      c30_l = 1.0 + (real_T)c30_b_l;
      CV_EML_FOR(0, 1, 0, 1);
      _SFD_EML_CALL(0U, chartInstance->c30_sfEvent, 8);
      c30_x = c30_mpower(chartInstance, c30_map[3 * ((int32_T)(real_T)
        _SFD_EML_ARRAY_BOUNDS_CHECK("map", (int32_T)_SFD_INTEGER_CHECK("i",
        c30_i), 1, 4, 2, 0) - 1)] - c30_map[3 * ((int32_T)(real_T)
        _SFD_EML_ARRAY_BOUNDS_CHECK("map", (int32_T)_SFD_INTEGER_CHECK("l",
        c30_l), 1, 4, 2, 0) - 1)]) + c30_mpower(chartInstance, c30_map[1 + 3 *
        ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("map", (int32_T)
        _SFD_INTEGER_CHECK("i", c30_i), 1, 4, 2, 0) - 1)] - c30_map[1 + 3 *
        ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("map", (int32_T)
        _SFD_INTEGER_CHECK("l", c30_l), 1, 4, 2, 0) - 1)]);
      c30_b_x = c30_x;
      if (c30_b_x < 0.0) {
        c30_eml_error(chartInstance);
      }

      c30_b_x = muDoubleScalarSqrt(c30_b_x);
      c30_dist[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("dist", (int32_T)
        _SFD_INTEGER_CHECK("l", c30_l), 1, 4, 1, 0) - 1] = c30_b_x;
      c30_b_l++;
      _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
    }

    CV_EML_FOR(0, 1, 0, 0);
  }

  _SFD_EML_CALL(0U, chartInstance->c30_sfEvent, -8);
  _SFD_SYMBOL_SCOPE_POP();
  for (c30_i10 = 0; c30_i10 < 4; c30_i10++) {
    (*c30_b_dist)[c30_i10] = c30_dist[c30_i10];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 29U, chartInstance->c30_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc30_simulation(SFc30_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c30_machineNumber, uint32_T
  c30_chartNumber)
{
}

static const mxArray *c30_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData)
{
  const mxArray *c30_mxArrayOutData = NULL;
  int32_T c30_i11;
  real_T c30_b_inData[4];
  int32_T c30_i12;
  real_T c30_u[4];
  const mxArray *c30_y = NULL;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_mxArrayOutData = NULL;
  for (c30_i11 = 0; c30_i11 < 4; c30_i11++) {
    c30_b_inData[c30_i11] = (*(real_T (*)[4])c30_inData)[c30_i11];
  }

  for (c30_i12 = 0; c30_i12 < 4; c30_i12++) {
    c30_u[c30_i12] = c30_b_inData[c30_i12];
  }

  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", c30_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_assign(&c30_mxArrayOutData, c30_y, FALSE);
  return c30_mxArrayOutData;
}

static void c30_emlrt_marshallIn(SFc30_simulationInstanceStruct *chartInstance,
  const mxArray *c30_dist, const char_T *c30_identifier, real_T c30_y[4])
{
  emlrtMsgIdentifier c30_thisId;
  c30_thisId.fIdentifier = c30_identifier;
  c30_thisId.fParent = NULL;
  c30_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c30_dist), &c30_thisId, c30_y);
  sf_mex_destroy(&c30_dist);
}

static void c30_b_emlrt_marshallIn(SFc30_simulationInstanceStruct *chartInstance,
  const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId, real_T c30_y[4])
{
  real_T c30_dv1[4];
  int32_T c30_i13;
  sf_mex_import(c30_parentId, sf_mex_dup(c30_u), c30_dv1, 1, 0, 0U, 1, 0U, 1, 4);
  for (c30_i13 = 0; c30_i13 < 4; c30_i13++) {
    c30_y[c30_i13] = c30_dv1[c30_i13];
  }

  sf_mex_destroy(&c30_u);
}

static void c30_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c30_mxArrayInData, const char_T *c30_varName, void *c30_outData)
{
  const mxArray *c30_dist;
  const char_T *c30_identifier;
  emlrtMsgIdentifier c30_thisId;
  real_T c30_y[4];
  int32_T c30_i14;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_dist = sf_mex_dup(c30_mxArrayInData);
  c30_identifier = c30_varName;
  c30_thisId.fIdentifier = c30_identifier;
  c30_thisId.fParent = NULL;
  c30_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c30_dist), &c30_thisId, c30_y);
  sf_mex_destroy(&c30_dist);
  for (c30_i14 = 0; c30_i14 < 4; c30_i14++) {
    (*(real_T (*)[4])c30_outData)[c30_i14] = c30_y[c30_i14];
  }

  sf_mex_destroy(&c30_mxArrayInData);
}

static const mxArray *c30_b_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData)
{
  const mxArray *c30_mxArrayOutData = NULL;
  real_T c30_u;
  const mxArray *c30_y = NULL;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_mxArrayOutData = NULL;
  c30_u = *(real_T *)c30_inData;
  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", &c30_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c30_mxArrayOutData, c30_y, FALSE);
  return c30_mxArrayOutData;
}

static real_T c30_c_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId)
{
  real_T c30_y;
  real_T c30_d1;
  sf_mex_import(c30_parentId, sf_mex_dup(c30_u), &c30_d1, 1, 0, 0U, 0, 0U, 0);
  c30_y = c30_d1;
  sf_mex_destroy(&c30_u);
  return c30_y;
}

static void c30_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c30_mxArrayInData, const char_T *c30_varName, void *c30_outData)
{
  const mxArray *c30_b_n;
  const char_T *c30_identifier;
  emlrtMsgIdentifier c30_thisId;
  real_T c30_y;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_b_n = sf_mex_dup(c30_mxArrayInData);
  c30_identifier = c30_varName;
  c30_thisId.fIdentifier = c30_identifier;
  c30_thisId.fParent = NULL;
  c30_y = c30_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c30_b_n), &c30_thisId);
  sf_mex_destroy(&c30_b_n);
  *(real_T *)c30_outData = c30_y;
  sf_mex_destroy(&c30_mxArrayInData);
}

static const mxArray *c30_c_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData)
{
  const mxArray *c30_mxArrayOutData = NULL;
  int32_T c30_i15;
  int32_T c30_i16;
  int32_T c30_i17;
  real_T c30_b_inData[12];
  int32_T c30_i18;
  int32_T c30_i19;
  int32_T c30_i20;
  real_T c30_u[12];
  const mxArray *c30_y = NULL;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_mxArrayOutData = NULL;
  c30_i15 = 0;
  for (c30_i16 = 0; c30_i16 < 4; c30_i16++) {
    for (c30_i17 = 0; c30_i17 < 3; c30_i17++) {
      c30_b_inData[c30_i17 + c30_i15] = (*(real_T (*)[12])c30_inData)[c30_i17 +
        c30_i15];
    }

    c30_i15 += 3;
  }

  c30_i18 = 0;
  for (c30_i19 = 0; c30_i19 < 4; c30_i19++) {
    for (c30_i20 = 0; c30_i20 < 3; c30_i20++) {
      c30_u[c30_i20 + c30_i18] = c30_b_inData[c30_i20 + c30_i18];
    }

    c30_i18 += 3;
  }

  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", c30_u, 0, 0U, 1U, 0U, 2, 3, 4), FALSE);
  sf_mex_assign(&c30_mxArrayOutData, c30_y, FALSE);
  return c30_mxArrayOutData;
}

static const mxArray *c30_d_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData)
{
  const mxArray *c30_mxArrayOutData = NULL;
  int32_T c30_i21;
  boolean_T c30_b_inData[4];
  int32_T c30_i22;
  boolean_T c30_u[4];
  const mxArray *c30_y = NULL;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_mxArrayOutData = NULL;
  for (c30_i21 = 0; c30_i21 < 4; c30_i21++) {
    c30_b_inData[c30_i21] = (*(boolean_T (*)[4])c30_inData)[c30_i21];
  }

  for (c30_i22 = 0; c30_i22 < 4; c30_i22++) {
    c30_u[c30_i22] = c30_b_inData[c30_i22];
  }

  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", c30_u, 11, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_assign(&c30_mxArrayOutData, c30_y, FALSE);
  return c30_mxArrayOutData;
}

static const mxArray *c30_e_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData)
{
  const mxArray *c30_mxArrayOutData = NULL;
  boolean_T c30_u;
  const mxArray *c30_y = NULL;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_mxArrayOutData = NULL;
  c30_u = *(boolean_T *)c30_inData;
  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", &c30_u, 11, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c30_mxArrayOutData, c30_y, FALSE);
  return c30_mxArrayOutData;
}

static boolean_T c30_d_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId)
{
  boolean_T c30_y;
  boolean_T c30_b0;
  sf_mex_import(c30_parentId, sf_mex_dup(c30_u), &c30_b0, 1, 11, 0U, 0, 0U, 0);
  c30_y = c30_b0;
  sf_mex_destroy(&c30_u);
  return c30_y;
}

static void c30_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c30_mxArrayInData, const char_T *c30_varName, void *c30_outData)
{
  const mxArray *c30_sonarmax;
  const char_T *c30_identifier;
  emlrtMsgIdentifier c30_thisId;
  boolean_T c30_y;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_sonarmax = sf_mex_dup(c30_mxArrayInData);
  c30_identifier = c30_varName;
  c30_thisId.fIdentifier = c30_identifier;
  c30_thisId.fParent = NULL;
  c30_y = c30_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c30_sonarmax),
    &c30_thisId);
  sf_mex_destroy(&c30_sonarmax);
  *(boolean_T *)c30_outData = c30_y;
  sf_mex_destroy(&c30_mxArrayInData);
}

const mxArray *sf_c30_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c30_nameCaptureInfo = NULL;
  c30_nameCaptureInfo = NULL;
  sf_mex_assign(&c30_nameCaptureInfo, sf_mex_createstruct("structure", 2, 29, 1),
                FALSE);
  c30_info_helper(&c30_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c30_nameCaptureInfo);
  return c30_nameCaptureInfo;
}

static void c30_info_helper(const mxArray **c30_info)
{
  const mxArray *c30_rhs0 = NULL;
  const mxArray *c30_lhs0 = NULL;
  const mxArray *c30_rhs1 = NULL;
  const mxArray *c30_lhs1 = NULL;
  const mxArray *c30_rhs2 = NULL;
  const mxArray *c30_lhs2 = NULL;
  const mxArray *c30_rhs3 = NULL;
  const mxArray *c30_lhs3 = NULL;
  const mxArray *c30_rhs4 = NULL;
  const mxArray *c30_lhs4 = NULL;
  const mxArray *c30_rhs5 = NULL;
  const mxArray *c30_lhs5 = NULL;
  const mxArray *c30_rhs6 = NULL;
  const mxArray *c30_lhs6 = NULL;
  const mxArray *c30_rhs7 = NULL;
  const mxArray *c30_lhs7 = NULL;
  const mxArray *c30_rhs8 = NULL;
  const mxArray *c30_lhs8 = NULL;
  const mxArray *c30_rhs9 = NULL;
  const mxArray *c30_lhs9 = NULL;
  const mxArray *c30_rhs10 = NULL;
  const mxArray *c30_lhs10 = NULL;
  const mxArray *c30_rhs11 = NULL;
  const mxArray *c30_lhs11 = NULL;
  const mxArray *c30_rhs12 = NULL;
  const mxArray *c30_lhs12 = NULL;
  const mxArray *c30_rhs13 = NULL;
  const mxArray *c30_lhs13 = NULL;
  const mxArray *c30_rhs14 = NULL;
  const mxArray *c30_lhs14 = NULL;
  const mxArray *c30_rhs15 = NULL;
  const mxArray *c30_lhs15 = NULL;
  const mxArray *c30_rhs16 = NULL;
  const mxArray *c30_lhs16 = NULL;
  const mxArray *c30_rhs17 = NULL;
  const mxArray *c30_lhs17 = NULL;
  const mxArray *c30_rhs18 = NULL;
  const mxArray *c30_lhs18 = NULL;
  const mxArray *c30_rhs19 = NULL;
  const mxArray *c30_lhs19 = NULL;
  const mxArray *c30_rhs20 = NULL;
  const mxArray *c30_lhs20 = NULL;
  const mxArray *c30_rhs21 = NULL;
  const mxArray *c30_lhs21 = NULL;
  const mxArray *c30_rhs22 = NULL;
  const mxArray *c30_lhs22 = NULL;
  const mxArray *c30_rhs23 = NULL;
  const mxArray *c30_lhs23 = NULL;
  const mxArray *c30_rhs24 = NULL;
  const mxArray *c30_lhs24 = NULL;
  const mxArray *c30_rhs25 = NULL;
  const mxArray *c30_lhs25 = NULL;
  const mxArray *c30_rhs26 = NULL;
  const mxArray *c30_lhs26 = NULL;
  const mxArray *c30_rhs27 = NULL;
  const mxArray *c30_lhs27 = NULL;
  const mxArray *c30_rhs28 = NULL;
  const mxArray *c30_lhs28 = NULL;
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("max"), "name", "name", 0);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1311258916U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c30_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 1);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c30_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum"),
                  "context", "context", 2);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_const_nonsingleton_dim"),
                  "name", "name", 2);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1286822296U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c30_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum"),
                  "context", "context", 3);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 3);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c30_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum"),
                  "context", "context", 4);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 4);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c30_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 5);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 5);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c30_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 6);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("isnan"), "name", "name", 6);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c30_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 7);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c30_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 8);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 8);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 8);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c30_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 9);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c30_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 10);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 10);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c30_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 11);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("intmax"), "name", "name", 11);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c30_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 12);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_relop"), "name", "name",
                  12);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("function_handle"),
                  "dominantType", "dominantType", 12);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_relop.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1342454782U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c30_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(""), "context", "context", 13);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("mpower"), "name", "name", 13);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c30_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 14);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c30_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("ismatrix"), "name", "name",
                  15);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c30_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("power"), "name", "name", 16);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c30_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 17);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c30_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 18);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 18);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c30_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 19);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 19);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1358189740U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c30_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 20);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("floor"), "name", "name", 20);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c30_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 21);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c30_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 22);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1286822326U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c30_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 23);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 23);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c30_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 24);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("mtimes"), "name", "name", 24);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c30_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 25);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 25);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c30_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(""), "context", "context", 26);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("sqrt"), "name", "name", 26);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c30_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_error"), "name", "name",
                  27);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1343833958U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c30_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 28);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c30_info, c30_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(1286822338U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c30_info, c30_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c30_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c30_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c30_info, sf_mex_duplicatearraysafe(&c30_lhs28), "lhs", "lhs",
                  28);
  sf_mex_destroy(&c30_rhs0);
  sf_mex_destroy(&c30_lhs0);
  sf_mex_destroy(&c30_rhs1);
  sf_mex_destroy(&c30_lhs1);
  sf_mex_destroy(&c30_rhs2);
  sf_mex_destroy(&c30_lhs2);
  sf_mex_destroy(&c30_rhs3);
  sf_mex_destroy(&c30_lhs3);
  sf_mex_destroy(&c30_rhs4);
  sf_mex_destroy(&c30_lhs4);
  sf_mex_destroy(&c30_rhs5);
  sf_mex_destroy(&c30_lhs5);
  sf_mex_destroy(&c30_rhs6);
  sf_mex_destroy(&c30_lhs6);
  sf_mex_destroy(&c30_rhs7);
  sf_mex_destroy(&c30_lhs7);
  sf_mex_destroy(&c30_rhs8);
  sf_mex_destroy(&c30_lhs8);
  sf_mex_destroy(&c30_rhs9);
  sf_mex_destroy(&c30_lhs9);
  sf_mex_destroy(&c30_rhs10);
  sf_mex_destroy(&c30_lhs10);
  sf_mex_destroy(&c30_rhs11);
  sf_mex_destroy(&c30_lhs11);
  sf_mex_destroy(&c30_rhs12);
  sf_mex_destroy(&c30_lhs12);
  sf_mex_destroy(&c30_rhs13);
  sf_mex_destroy(&c30_lhs13);
  sf_mex_destroy(&c30_rhs14);
  sf_mex_destroy(&c30_lhs14);
  sf_mex_destroy(&c30_rhs15);
  sf_mex_destroy(&c30_lhs15);
  sf_mex_destroy(&c30_rhs16);
  sf_mex_destroy(&c30_lhs16);
  sf_mex_destroy(&c30_rhs17);
  sf_mex_destroy(&c30_lhs17);
  sf_mex_destroy(&c30_rhs18);
  sf_mex_destroy(&c30_lhs18);
  sf_mex_destroy(&c30_rhs19);
  sf_mex_destroy(&c30_lhs19);
  sf_mex_destroy(&c30_rhs20);
  sf_mex_destroy(&c30_lhs20);
  sf_mex_destroy(&c30_rhs21);
  sf_mex_destroy(&c30_lhs21);
  sf_mex_destroy(&c30_rhs22);
  sf_mex_destroy(&c30_lhs22);
  sf_mex_destroy(&c30_rhs23);
  sf_mex_destroy(&c30_lhs23);
  sf_mex_destroy(&c30_rhs24);
  sf_mex_destroy(&c30_lhs24);
  sf_mex_destroy(&c30_rhs25);
  sf_mex_destroy(&c30_lhs25);
  sf_mex_destroy(&c30_rhs26);
  sf_mex_destroy(&c30_lhs26);
  sf_mex_destroy(&c30_rhs27);
  sf_mex_destroy(&c30_lhs27);
  sf_mex_destroy(&c30_rhs28);
  sf_mex_destroy(&c30_lhs28);
}

static const mxArray *c30_emlrt_marshallOut(char * c30_u)
{
  const mxArray *c30_y = NULL;
  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", c30_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c30_u)), FALSE);
  return c30_y;
}

static const mxArray *c30_b_emlrt_marshallOut(uint32_T c30_u)
{
  const mxArray *c30_y = NULL;
  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", &c30_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c30_y;
}

static real_T c30_mpower(SFc30_simulationInstanceStruct *chartInstance, real_T
  c30_a)
{
  real_T c30_b_a;
  real_T c30_c_a;
  real_T c30_ak;
  real_T c30_d_a;
  real_T c30_e_a;
  real_T c30_b;
  c30_b_a = c30_a;
  c30_c_a = c30_b_a;
  c30_eml_scalar_eg(chartInstance);
  c30_ak = c30_c_a;
  c30_d_a = c30_ak;
  c30_eml_scalar_eg(chartInstance);
  c30_e_a = c30_d_a;
  c30_b = c30_d_a;
  return c30_e_a * c30_b;
}

static void c30_eml_scalar_eg(SFc30_simulationInstanceStruct *chartInstance)
{
}

static void c30_eml_error(SFc30_simulationInstanceStruct *chartInstance)
{
  int32_T c30_i23;
  static char_T c30_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c30_u[30];
  const mxArray *c30_y = NULL;
  int32_T c30_i24;
  static char_T c30_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c30_b_u[4];
  const mxArray *c30_b_y = NULL;
  for (c30_i23 = 0; c30_i23 < 30; c30_i23++) {
    c30_u[c30_i23] = c30_cv0[c30_i23];
  }

  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", c30_u, 10, 0U, 1U, 0U, 2, 1, 30),
                FALSE);
  for (c30_i24 = 0; c30_i24 < 4; c30_i24++) {
    c30_b_u[c30_i24] = c30_cv1[c30_i24];
  }

  c30_b_y = NULL;
  sf_mex_assign(&c30_b_y, sf_mex_create("y", c30_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U, 14,
    c30_y, 14, c30_b_y));
}

static const mxArray *c30_f_sf_marshallOut(void *chartInstanceVoid, void
  *c30_inData)
{
  const mxArray *c30_mxArrayOutData = NULL;
  int32_T c30_u;
  const mxArray *c30_y = NULL;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_mxArrayOutData = NULL;
  c30_u = *(int32_T *)c30_inData;
  c30_y = NULL;
  sf_mex_assign(&c30_y, sf_mex_create("y", &c30_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c30_mxArrayOutData, c30_y, FALSE);
  return c30_mxArrayOutData;
}

static int32_T c30_e_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId)
{
  int32_T c30_y;
  int32_T c30_i25;
  sf_mex_import(c30_parentId, sf_mex_dup(c30_u), &c30_i25, 1, 6, 0U, 0, 0U, 0);
  c30_y = c30_i25;
  sf_mex_destroy(&c30_u);
  return c30_y;
}

static void c30_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c30_mxArrayInData, const char_T *c30_varName, void *c30_outData)
{
  const mxArray *c30_b_sfEvent;
  const char_T *c30_identifier;
  emlrtMsgIdentifier c30_thisId;
  int32_T c30_y;
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)chartInstanceVoid;
  c30_b_sfEvent = sf_mex_dup(c30_mxArrayInData);
  c30_identifier = c30_varName;
  c30_thisId.fIdentifier = c30_identifier;
  c30_thisId.fParent = NULL;
  c30_y = c30_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c30_b_sfEvent),
    &c30_thisId);
  sf_mex_destroy(&c30_b_sfEvent);
  *(int32_T *)c30_outData = c30_y;
  sf_mex_destroy(&c30_mxArrayInData);
}

static uint8_T c30_f_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_b_is_active_c30_simulation, const char_T
  *c30_identifier)
{
  uint8_T c30_y;
  emlrtMsgIdentifier c30_thisId;
  c30_thisId.fIdentifier = c30_identifier;
  c30_thisId.fParent = NULL;
  c30_y = c30_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c30_b_is_active_c30_simulation), &c30_thisId);
  sf_mex_destroy(&c30_b_is_active_c30_simulation);
  return c30_y;
}

static uint8_T c30_g_emlrt_marshallIn(SFc30_simulationInstanceStruct
  *chartInstance, const mxArray *c30_u, const emlrtMsgIdentifier *c30_parentId)
{
  uint8_T c30_y;
  uint8_T c30_u0;
  sf_mex_import(c30_parentId, sf_mex_dup(c30_u), &c30_u0, 1, 3, 0U, 0, 0U, 0);
  c30_y = c30_u0;
  sf_mex_destroy(&c30_u);
  return c30_y;
}

static void init_dsm_address_info(SFc30_simulationInstanceStruct *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c30_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3483400384U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3438041391U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3347318242U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3466083732U);
}

mxArray *sf_c30_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("apEnm6j64xX2lsw7rrqo6D");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(1));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(4);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c30_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c30_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c30_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"dist\",},{M[8],M[0],T\"is_active_c30_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c30_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc30_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc30_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           30,
           1,
           1,
           4,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_simulationMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_simulationMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _simulationMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"sonar");
          _SFD_SET_DATA_PROPS(1,2,0,1,"dist");
          _SFD_SET_DATA_PROPS(2,1,1,0,"map");
          _SFD_SET_DATA_PROPS(3,10,0,0,"n");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,1,0,0,0,1,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,234);
        _SFD_CV_INIT_EML_IF(0,1,0,119,134,-1,234);
        _SFD_CV_INIT_EML_FOR(0,1,0,139,151,230);

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_UINT8,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c30_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c30_sf_marshallOut,(MexInFcnForType)
            c30_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c30_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c30_b_sf_marshallOut,(MexInFcnForType)
          c30_b_sf_marshallIn);

        {
          boolean_T (*c30_sonar)[4];
          real_T (*c30_dist)[4];
          real_T (*c30_map)[12];
          c30_map = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 1);
          c30_dist = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 1);
          c30_sonar = (boolean_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c30_sonar);
          _SFD_SET_DATA_VALUE_PTR(1U, *c30_dist);
          _SFD_SET_DATA_VALUE_PTR(2U, *c30_map);
          _SFD_SET_DATA_VALUE_PTR(3U, &chartInstance->c30_n);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _simulationMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "Zdcsl9NFBPhOh43b6o46NC";
}

static void sf_opaque_initialize_c30_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc30_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c30_simulation((SFc30_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c30_simulation((SFc30_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c30_simulation(void *chartInstanceVar)
{
  enable_c30_simulation((SFc30_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c30_simulation(void *chartInstanceVar)
{
  disable_c30_simulation((SFc30_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c30_simulation(void *chartInstanceVar)
{
  sf_c30_simulation((SFc30_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c30_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c30_simulation
    ((SFc30_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c30_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c30_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c30_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c30_simulation((SFc30_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c30_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c30_simulation(S);
}

static void sf_opaque_set_sim_state_c30_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c30_simulation(S, st);
}

static void sf_opaque_terminate_c30_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc30_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c30_simulation((SFc30_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc30_simulation((SFc30_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c30_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c30_simulation((SFc30_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c30_simulation(SimStruct *S)
{
  /* Actual parameters from chart:
     n
   */
  const char_T *rtParamNames[] = { "n" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for n*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      30);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,30,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,30,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,30);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,30,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,30,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,30);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3906513741U));
  ssSetChecksum1(S,(4028303289U));
  ssSetChecksum2(S,(450904539U));
  ssSetChecksum3(S,(2801608223U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c30_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c30_simulation(SimStruct *S)
{
  SFc30_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc30_simulationInstanceStruct *)utMalloc(sizeof
    (SFc30_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc30_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c30_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c30_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c30_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c30_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c30_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c30_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c30_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c30_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c30_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c30_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c30_simulation;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c30_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c30_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c30_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c30_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c30_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
