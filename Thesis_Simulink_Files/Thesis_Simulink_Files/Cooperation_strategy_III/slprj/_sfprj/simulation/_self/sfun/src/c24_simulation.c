/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c24_simulation.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simulation_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c24_debug_family_names[11] = { "pdataiStep", "nargin",
  "nargout", "pathdata", "gam", "id", "gamWP", "thetad", "P0", "v_L", "Circle" };

/* Function Declarations */
static void initialize_c24_simulation(SFc24_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c24_simulation(SFc24_simulationInstanceStruct
  *chartInstance);
static void enable_c24_simulation(SFc24_simulationInstanceStruct *chartInstance);
static void disable_c24_simulation(SFc24_simulationInstanceStruct *chartInstance);
static void c24_update_debugger_state_c24_simulation
  (SFc24_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c24_simulation
  (SFc24_simulationInstanceStruct *chartInstance);
static void set_sim_state_c24_simulation(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_st);
static void finalize_c24_simulation(SFc24_simulationInstanceStruct
  *chartInstance);
static void sf_c24_simulation(SFc24_simulationInstanceStruct *chartInstance);
static void initSimStructsc24_simulation(SFc24_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c24_machineNumber, uint32_T
  c24_chartNumber);
static const mxArray *c24_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData);
static real_T c24_emlrt_marshallIn(SFc24_simulationInstanceStruct *chartInstance,
  const mxArray *c24_Circle, const char_T *c24_identifier);
static real_T c24_b_emlrt_marshallIn(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId);
static void c24_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c24_mxArrayInData, const char_T *c24_varName, void *c24_outData);
static const mxArray *c24_b_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData);
static void c24_c_emlrt_marshallIn(SFc24_simulationInstanceStruct *chartInstance,
  const mxArray *c24_P0, const char_T *c24_identifier, real_T c24_y[2]);
static void c24_d_emlrt_marshallIn(SFc24_simulationInstanceStruct *chartInstance,
  const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId, real_T c24_y[2]);
static void c24_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c24_mxArrayInData, const char_T *c24_varName, void *c24_outData);
static const mxArray *c24_c_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData);
static const mxArray *c24_d_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData);
static void c24_e_emlrt_marshallIn(SFc24_simulationInstanceStruct *chartInstance,
  const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId, real_T c24_y[6]);
static void c24_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c24_mxArrayInData, const char_T *c24_varName, void *c24_outData);
static const mxArray *c24_e_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData);
static int32_T c24_f_emlrt_marshallIn(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId);
static void c24_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c24_mxArrayInData, const char_T *c24_varName, void *c24_outData);
static uint8_T c24_g_emlrt_marshallIn(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_b_is_active_c24_simulation, const char_T
  *c24_identifier);
static uint8_T c24_h_emlrt_marshallIn(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId);
static void init_dsm_address_info(SFc24_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c24_simulation(SFc24_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c24_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c24_is_active_c24_simulation = 0U;
}

static void initialize_params_c24_simulation(SFc24_simulationInstanceStruct
  *chartInstance)
{
}

static void enable_c24_simulation(SFc24_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c24_simulation(SFc24_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c24_update_debugger_state_c24_simulation
  (SFc24_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c24_simulation
  (SFc24_simulationInstanceStruct *chartInstance)
{
  const mxArray *c24_st;
  const mxArray *c24_y = NULL;
  real_T c24_hoistedGlobal;
  real_T c24_u;
  const mxArray *c24_b_y = NULL;
  int32_T c24_i0;
  real_T c24_b_u[2];
  const mxArray *c24_c_y = NULL;
  real_T c24_b_hoistedGlobal;
  real_T c24_c_u;
  const mxArray *c24_d_y = NULL;
  real_T c24_c_hoistedGlobal;
  real_T c24_d_u;
  const mxArray *c24_e_y = NULL;
  real_T c24_d_hoistedGlobal;
  real_T c24_e_u;
  const mxArray *c24_f_y = NULL;
  uint8_T c24_e_hoistedGlobal;
  uint8_T c24_f_u;
  const mxArray *c24_g_y = NULL;
  real_T *c24_Circle;
  real_T *c24_gamWP;
  real_T *c24_thetad;
  real_T *c24_v_L;
  real_T (*c24_P0)[2];
  c24_Circle = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c24_v_L = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c24_P0 = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c24_thetad = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c24_gamWP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c24_st = NULL;
  c24_st = NULL;
  c24_y = NULL;
  sf_mex_assign(&c24_y, sf_mex_createcellarray(6), FALSE);
  c24_hoistedGlobal = *c24_Circle;
  c24_u = c24_hoistedGlobal;
  c24_b_y = NULL;
  sf_mex_assign(&c24_b_y, sf_mex_create("y", &c24_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c24_y, 0, c24_b_y);
  for (c24_i0 = 0; c24_i0 < 2; c24_i0++) {
    c24_b_u[c24_i0] = (*c24_P0)[c24_i0];
  }

  c24_c_y = NULL;
  sf_mex_assign(&c24_c_y, sf_mex_create("y", c24_b_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c24_y, 1, c24_c_y);
  c24_b_hoistedGlobal = *c24_gamWP;
  c24_c_u = c24_b_hoistedGlobal;
  c24_d_y = NULL;
  sf_mex_assign(&c24_d_y, sf_mex_create("y", &c24_c_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c24_y, 2, c24_d_y);
  c24_c_hoistedGlobal = *c24_thetad;
  c24_d_u = c24_c_hoistedGlobal;
  c24_e_y = NULL;
  sf_mex_assign(&c24_e_y, sf_mex_create("y", &c24_d_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c24_y, 3, c24_e_y);
  c24_d_hoistedGlobal = *c24_v_L;
  c24_e_u = c24_d_hoistedGlobal;
  c24_f_y = NULL;
  sf_mex_assign(&c24_f_y, sf_mex_create("y", &c24_e_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c24_y, 4, c24_f_y);
  c24_e_hoistedGlobal = chartInstance->c24_is_active_c24_simulation;
  c24_f_u = c24_e_hoistedGlobal;
  c24_g_y = NULL;
  sf_mex_assign(&c24_g_y, sf_mex_create("y", &c24_f_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c24_y, 5, c24_g_y);
  sf_mex_assign(&c24_st, c24_y, FALSE);
  return c24_st;
}

static void set_sim_state_c24_simulation(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_st)
{
  const mxArray *c24_u;
  real_T c24_dv0[2];
  int32_T c24_i1;
  real_T *c24_Circle;
  real_T *c24_gamWP;
  real_T *c24_thetad;
  real_T *c24_v_L;
  real_T (*c24_P0)[2];
  c24_Circle = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c24_v_L = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c24_P0 = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c24_thetad = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c24_gamWP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c24_doneDoubleBufferReInit = TRUE;
  c24_u = sf_mex_dup(c24_st);
  *c24_Circle = c24_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c24_u, 0)), "Circle");
  c24_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c24_u, 1)),
    "P0", c24_dv0);
  for (c24_i1 = 0; c24_i1 < 2; c24_i1++) {
    (*c24_P0)[c24_i1] = c24_dv0[c24_i1];
  }

  *c24_gamWP = c24_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c24_u, 2)), "gamWP");
  *c24_thetad = c24_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c24_u, 3)), "thetad");
  *c24_v_L = c24_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c24_u,
    4)), "v_L");
  chartInstance->c24_is_active_c24_simulation = c24_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c24_u, 5)),
     "is_active_c24_simulation");
  sf_mex_destroy(&c24_u);
  c24_update_debugger_state_c24_simulation(chartInstance);
  sf_mex_destroy(&c24_st);
}

static void finalize_c24_simulation(SFc24_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c24_simulation(SFc24_simulationInstanceStruct *chartInstance)
{
  int32_T c24_i2;
  int32_T c24_i3;
  real_T c24_hoistedGlobal;
  real_T c24_b_hoistedGlobal;
  int32_T c24_i4;
  real_T c24_pathdata[144];
  real_T c24_gam;
  real_T c24_id;
  uint32_T c24_debug_family_var_map[11];
  real_T c24_pdataiStep[6];
  real_T c24_nargin = 3.0;
  real_T c24_nargout = 5.0;
  real_T c24_gamWP;
  real_T c24_thetad;
  real_T c24_P0[2];
  real_T c24_v_L;
  real_T c24_Circle;
  int32_T c24_b_id;
  int32_T c24_i5;
  int32_T c24_c_id;
  int32_T c24_i6;
  int32_T c24_d_id;
  int32_T c24_i7;
  int32_T c24_e_id;
  int32_T c24_i8;
  int32_T c24_f_id;
  int32_T c24_i9;
  int32_T c24_g_id;
  int32_T c24_i10;
  int32_T c24_i11;
  real_T *c24_b_gam;
  real_T *c24_b_gamWP;
  real_T *c24_b_thetad;
  real_T *c24_b_v_L;
  real_T *c24_h_id;
  real_T *c24_b_Circle;
  real_T (*c24_b_P0)[2];
  real_T (*c24_b_pathdata)[144];
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  boolean_T guard3 = FALSE;
  boolean_T guard4 = FALSE;
  boolean_T guard5 = FALSE;
  c24_b_Circle = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c24_h_id = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c24_b_v_L = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c24_b_P0 = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c24_b_thetad = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c24_b_gamWP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c24_b_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c24_b_pathdata = (real_T (*)[144])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 23U, chartInstance->c24_sfEvent);
  for (c24_i2 = 0; c24_i2 < 144; c24_i2++) {
    _SFD_DATA_RANGE_CHECK((*c24_b_pathdata)[c24_i2], 0U);
  }

  _SFD_DATA_RANGE_CHECK(*c24_b_gam, 1U);
  _SFD_DATA_RANGE_CHECK(*c24_b_gamWP, 2U);
  _SFD_DATA_RANGE_CHECK(*c24_b_thetad, 3U);
  for (c24_i3 = 0; c24_i3 < 2; c24_i3++) {
    _SFD_DATA_RANGE_CHECK((*c24_b_P0)[c24_i3], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*c24_b_v_L, 5U);
  _SFD_DATA_RANGE_CHECK(*c24_h_id, 6U);
  _SFD_DATA_RANGE_CHECK(*c24_b_Circle, 7U);
  chartInstance->c24_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 23U, chartInstance->c24_sfEvent);
  c24_hoistedGlobal = *c24_b_gam;
  c24_b_hoistedGlobal = *c24_h_id;
  for (c24_i4 = 0; c24_i4 < 144; c24_i4++) {
    c24_pathdata[c24_i4] = (*c24_b_pathdata)[c24_i4];
  }

  c24_gam = c24_hoistedGlobal;
  c24_id = c24_b_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 11U, 11U, c24_debug_family_names,
    c24_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c24_pdataiStep, 0U, c24_d_sf_marshallOut,
    c24_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c24_nargin, 1U, c24_sf_marshallOut,
    c24_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c24_nargout, 2U, c24_sf_marshallOut,
    c24_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c24_pathdata, 3U, c24_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c24_gam, 4U, c24_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c24_id, 5U, c24_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c24_gamWP, 6U, c24_sf_marshallOut,
    c24_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c24_thetad, 7U, c24_sf_marshallOut,
    c24_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c24_P0, 8U, c24_b_sf_marshallOut,
    c24_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c24_v_L, 9U, c24_sf_marshallOut,
    c24_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c24_Circle, 10U, c24_sf_marshallOut,
    c24_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 3);
  guard1 = FALSE;
  if (CV_EML_COND(0, 1, 0, c24_gam >= 0.0)) {
    if (CV_EML_COND(0, 1, 1, c24_gam < c24_pathdata[6 + 36 * ((int32_T)(real_T)
          _SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)_SFD_INTEGER_CHECK(
            "id", c24_id), 1, 4, 3, 0) - 1)])) {
      CV_EML_MCDC(0, 1, 0, TRUE);
      CV_EML_IF(0, 1, 0, TRUE);
      _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 4);
      c24_b_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
        (int32_T)_SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1;
      for (c24_i5 = 0; c24_i5 < 6; c24_i5++) {
        c24_pdataiStep[c24_i5] = c24_pathdata[c24_i5 + 36 * c24_b_id];
      }
    } else {
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    CV_EML_MCDC(0, 1, 0, FALSE);
    CV_EML_IF(0, 1, 0, FALSE);
    _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 5);
    guard2 = FALSE;
    if (CV_EML_COND(0, 1, 2, c24_gam >= c24_pathdata[6 + 36 * ((int32_T)(real_T)
          _SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)_SFD_INTEGER_CHECK(
            "id", c24_id), 1, 4, 3, 0) - 1)])) {
      if (CV_EML_COND(0, 1, 3, c24_gam < c24_pathdata[12 + 36 * ((int32_T)
            (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
             _SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1)])) {
        CV_EML_MCDC(0, 1, 1, TRUE);
        CV_EML_IF(0, 1, 1, TRUE);
        _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 6);
        c24_c_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
          (int32_T)_SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1;
        for (c24_i6 = 0; c24_i6 < 6; c24_i6++) {
          c24_pdataiStep[c24_i6] = c24_pathdata[6 + (c24_i6 + 36 * c24_c_id)];
        }
      } else {
        guard2 = TRUE;
      }
    } else {
      guard2 = TRUE;
    }

    if (guard2 == TRUE) {
      CV_EML_MCDC(0, 1, 1, FALSE);
      CV_EML_IF(0, 1, 1, FALSE);
      _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 7);
      guard3 = FALSE;
      if (CV_EML_COND(0, 1, 4, c24_gam >= c24_pathdata[12 + 36 * ((int32_T)
            (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
             _SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1)])) {
        if (CV_EML_COND(0, 1, 5, c24_gam < c24_pathdata[18 + 36 * ((int32_T)
              (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
               _SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1)])) {
          CV_EML_MCDC(0, 1, 2, TRUE);
          CV_EML_IF(0, 1, 2, TRUE);
          _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 8);
          c24_d_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
            (int32_T)_SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1;
          for (c24_i7 = 0; c24_i7 < 6; c24_i7++) {
            c24_pdataiStep[c24_i7] = c24_pathdata[12 + (c24_i7 + 36 * c24_d_id)];
          }
        } else {
          guard3 = TRUE;
        }
      } else {
        guard3 = TRUE;
      }

      if (guard3 == TRUE) {
        CV_EML_MCDC(0, 1, 2, FALSE);
        CV_EML_IF(0, 1, 2, FALSE);
        _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 9);
        guard4 = FALSE;
        if (CV_EML_COND(0, 1, 6, c24_gam >= c24_pathdata[18 + 36 * ((int32_T)
              (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
               _SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1)])) {
          if (CV_EML_COND(0, 1, 7, c24_gam < c24_pathdata[24 + 36 * ((int32_T)
                (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
                 _SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1)])) {
            CV_EML_MCDC(0, 1, 3, TRUE);
            CV_EML_IF(0, 1, 3, TRUE);
            _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 10);
            c24_e_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
              (int32_T)_SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1;
            for (c24_i8 = 0; c24_i8 < 6; c24_i8++) {
              c24_pdataiStep[c24_i8] = c24_pathdata[18 + (c24_i8 + 36 * c24_e_id)];
            }
          } else {
            guard4 = TRUE;
          }
        } else {
          guard4 = TRUE;
        }

        if (guard4 == TRUE) {
          CV_EML_MCDC(0, 1, 3, FALSE);
          CV_EML_IF(0, 1, 3, FALSE);
          _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 11);
          guard5 = FALSE;
          if (CV_EML_COND(0, 1, 8, c24_gam >= c24_pathdata[24 + 36 * ((int32_T)
                (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
                 _SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1)])) {
            if (CV_EML_COND(0, 1, 9, c24_gam < c24_pathdata[30 + 36 * ((int32_T)
                  (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
                   _SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1)])) {
              CV_EML_MCDC(0, 1, 4, TRUE);
              CV_EML_IF(0, 1, 4, TRUE);
              _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 12);
              c24_f_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
                (int32_T)_SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1;
              for (c24_i9 = 0; c24_i9 < 6; c24_i9++) {
                c24_pdataiStep[c24_i9] = c24_pathdata[24 + (c24_i9 + 36 *
                  c24_f_id)];
              }
            } else {
              guard5 = TRUE;
            }
          } else {
            guard5 = TRUE;
          }

          if (guard5 == TRUE) {
            CV_EML_MCDC(0, 1, 4, FALSE);
            CV_EML_IF(0, 1, 4, FALSE);
            _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 14);
            c24_g_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
              (int32_T)_SFD_INTEGER_CHECK("id", c24_id), 1, 4, 3, 0) - 1;
            for (c24_i10 = 0; c24_i10 < 6; c24_i10++) {
              c24_pdataiStep[c24_i10] = c24_pathdata[30 + (c24_i10 + 36 *
                c24_g_id)];
            }
          }
        }
      }
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 17);
  c24_gamWP = c24_pdataiStep[0];
  _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 18);
  c24_thetad = c24_pdataiStep[1];
  _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 19);
  c24_P0[0] = c24_pdataiStep[2];
  c24_P0[1] = c24_pdataiStep[3];
  _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 20);
  c24_v_L = c24_pdataiStep[4];
  _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, 21);
  c24_Circle = c24_pdataiStep[5];
  _SFD_EML_CALL(0U, chartInstance->c24_sfEvent, -21);
  _SFD_SYMBOL_SCOPE_POP();
  *c24_b_gamWP = c24_gamWP;
  *c24_b_thetad = c24_thetad;
  for (c24_i11 = 0; c24_i11 < 2; c24_i11++) {
    (*c24_b_P0)[c24_i11] = c24_P0[c24_i11];
  }

  *c24_b_v_L = c24_v_L;
  *c24_b_Circle = c24_Circle;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 23U, chartInstance->c24_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc24_simulation(SFc24_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c24_machineNumber, uint32_T
  c24_chartNumber)
{
}

static const mxArray *c24_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData)
{
  const mxArray *c24_mxArrayOutData = NULL;
  real_T c24_u;
  const mxArray *c24_y = NULL;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_mxArrayOutData = NULL;
  c24_u = *(real_T *)c24_inData;
  c24_y = NULL;
  sf_mex_assign(&c24_y, sf_mex_create("y", &c24_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c24_mxArrayOutData, c24_y, FALSE);
  return c24_mxArrayOutData;
}

static real_T c24_emlrt_marshallIn(SFc24_simulationInstanceStruct *chartInstance,
  const mxArray *c24_Circle, const char_T *c24_identifier)
{
  real_T c24_y;
  emlrtMsgIdentifier c24_thisId;
  c24_thisId.fIdentifier = c24_identifier;
  c24_thisId.fParent = NULL;
  c24_y = c24_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c24_Circle),
    &c24_thisId);
  sf_mex_destroy(&c24_Circle);
  return c24_y;
}

static real_T c24_b_emlrt_marshallIn(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId)
{
  real_T c24_y;
  real_T c24_d0;
  sf_mex_import(c24_parentId, sf_mex_dup(c24_u), &c24_d0, 1, 0, 0U, 0, 0U, 0);
  c24_y = c24_d0;
  sf_mex_destroy(&c24_u);
  return c24_y;
}

static void c24_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c24_mxArrayInData, const char_T *c24_varName, void *c24_outData)
{
  const mxArray *c24_Circle;
  const char_T *c24_identifier;
  emlrtMsgIdentifier c24_thisId;
  real_T c24_y;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_Circle = sf_mex_dup(c24_mxArrayInData);
  c24_identifier = c24_varName;
  c24_thisId.fIdentifier = c24_identifier;
  c24_thisId.fParent = NULL;
  c24_y = c24_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c24_Circle),
    &c24_thisId);
  sf_mex_destroy(&c24_Circle);
  *(real_T *)c24_outData = c24_y;
  sf_mex_destroy(&c24_mxArrayInData);
}

static const mxArray *c24_b_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData)
{
  const mxArray *c24_mxArrayOutData = NULL;
  int32_T c24_i12;
  real_T c24_b_inData[2];
  int32_T c24_i13;
  real_T c24_u[2];
  const mxArray *c24_y = NULL;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_mxArrayOutData = NULL;
  for (c24_i12 = 0; c24_i12 < 2; c24_i12++) {
    c24_b_inData[c24_i12] = (*(real_T (*)[2])c24_inData)[c24_i12];
  }

  for (c24_i13 = 0; c24_i13 < 2; c24_i13++) {
    c24_u[c24_i13] = c24_b_inData[c24_i13];
  }

  c24_y = NULL;
  sf_mex_assign(&c24_y, sf_mex_create("y", c24_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c24_mxArrayOutData, c24_y, FALSE);
  return c24_mxArrayOutData;
}

static void c24_c_emlrt_marshallIn(SFc24_simulationInstanceStruct *chartInstance,
  const mxArray *c24_P0, const char_T *c24_identifier, real_T c24_y[2])
{
  emlrtMsgIdentifier c24_thisId;
  c24_thisId.fIdentifier = c24_identifier;
  c24_thisId.fParent = NULL;
  c24_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c24_P0), &c24_thisId, c24_y);
  sf_mex_destroy(&c24_P0);
}

static void c24_d_emlrt_marshallIn(SFc24_simulationInstanceStruct *chartInstance,
  const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId, real_T c24_y[2])
{
  real_T c24_dv1[2];
  int32_T c24_i14;
  sf_mex_import(c24_parentId, sf_mex_dup(c24_u), c24_dv1, 1, 0, 0U, 1, 0U, 1, 2);
  for (c24_i14 = 0; c24_i14 < 2; c24_i14++) {
    c24_y[c24_i14] = c24_dv1[c24_i14];
  }

  sf_mex_destroy(&c24_u);
}

static void c24_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c24_mxArrayInData, const char_T *c24_varName, void *c24_outData)
{
  const mxArray *c24_P0;
  const char_T *c24_identifier;
  emlrtMsgIdentifier c24_thisId;
  real_T c24_y[2];
  int32_T c24_i15;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_P0 = sf_mex_dup(c24_mxArrayInData);
  c24_identifier = c24_varName;
  c24_thisId.fIdentifier = c24_identifier;
  c24_thisId.fParent = NULL;
  c24_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c24_P0), &c24_thisId, c24_y);
  sf_mex_destroy(&c24_P0);
  for (c24_i15 = 0; c24_i15 < 2; c24_i15++) {
    (*(real_T (*)[2])c24_outData)[c24_i15] = c24_y[c24_i15];
  }

  sf_mex_destroy(&c24_mxArrayInData);
}

static const mxArray *c24_c_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData)
{
  const mxArray *c24_mxArrayOutData = NULL;
  int32_T c24_i16;
  int32_T c24_i17;
  int32_T c24_i18;
  int32_T c24_i19;
  int32_T c24_i20;
  real_T c24_b_inData[144];
  int32_T c24_i21;
  int32_T c24_i22;
  int32_T c24_i23;
  int32_T c24_i24;
  int32_T c24_i25;
  real_T c24_u[144];
  const mxArray *c24_y = NULL;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_mxArrayOutData = NULL;
  c24_i16 = 0;
  for (c24_i17 = 0; c24_i17 < 4; c24_i17++) {
    c24_i18 = 0;
    for (c24_i19 = 0; c24_i19 < 6; c24_i19++) {
      for (c24_i20 = 0; c24_i20 < 6; c24_i20++) {
        c24_b_inData[(c24_i20 + c24_i18) + c24_i16] = (*(real_T (*)[144])
          c24_inData)[(c24_i20 + c24_i18) + c24_i16];
      }

      c24_i18 += 6;
    }

    c24_i16 += 36;
  }

  c24_i21 = 0;
  for (c24_i22 = 0; c24_i22 < 4; c24_i22++) {
    c24_i23 = 0;
    for (c24_i24 = 0; c24_i24 < 6; c24_i24++) {
      for (c24_i25 = 0; c24_i25 < 6; c24_i25++) {
        c24_u[(c24_i25 + c24_i23) + c24_i21] = c24_b_inData[(c24_i25 + c24_i23)
          + c24_i21];
      }

      c24_i23 += 6;
    }

    c24_i21 += 36;
  }

  c24_y = NULL;
  sf_mex_assign(&c24_y, sf_mex_create("y", c24_u, 0, 0U, 1U, 0U, 3, 6, 6, 4),
                FALSE);
  sf_mex_assign(&c24_mxArrayOutData, c24_y, FALSE);
  return c24_mxArrayOutData;
}

static const mxArray *c24_d_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData)
{
  const mxArray *c24_mxArrayOutData = NULL;
  int32_T c24_i26;
  real_T c24_b_inData[6];
  int32_T c24_i27;
  real_T c24_u[6];
  const mxArray *c24_y = NULL;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_mxArrayOutData = NULL;
  for (c24_i26 = 0; c24_i26 < 6; c24_i26++) {
    c24_b_inData[c24_i26] = (*(real_T (*)[6])c24_inData)[c24_i26];
  }

  for (c24_i27 = 0; c24_i27 < 6; c24_i27++) {
    c24_u[c24_i27] = c24_b_inData[c24_i27];
  }

  c24_y = NULL;
  sf_mex_assign(&c24_y, sf_mex_create("y", c24_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_assign(&c24_mxArrayOutData, c24_y, FALSE);
  return c24_mxArrayOutData;
}

static void c24_e_emlrt_marshallIn(SFc24_simulationInstanceStruct *chartInstance,
  const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId, real_T c24_y[6])
{
  real_T c24_dv2[6];
  int32_T c24_i28;
  sf_mex_import(c24_parentId, sf_mex_dup(c24_u), c24_dv2, 1, 0, 0U, 1, 0U, 1, 6);
  for (c24_i28 = 0; c24_i28 < 6; c24_i28++) {
    c24_y[c24_i28] = c24_dv2[c24_i28];
  }

  sf_mex_destroy(&c24_u);
}

static void c24_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c24_mxArrayInData, const char_T *c24_varName, void *c24_outData)
{
  const mxArray *c24_pdataiStep;
  const char_T *c24_identifier;
  emlrtMsgIdentifier c24_thisId;
  real_T c24_y[6];
  int32_T c24_i29;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_pdataiStep = sf_mex_dup(c24_mxArrayInData);
  c24_identifier = c24_varName;
  c24_thisId.fIdentifier = c24_identifier;
  c24_thisId.fParent = NULL;
  c24_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c24_pdataiStep), &c24_thisId,
    c24_y);
  sf_mex_destroy(&c24_pdataiStep);
  for (c24_i29 = 0; c24_i29 < 6; c24_i29++) {
    (*(real_T (*)[6])c24_outData)[c24_i29] = c24_y[c24_i29];
  }

  sf_mex_destroy(&c24_mxArrayInData);
}

const mxArray *sf_c24_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c24_nameCaptureInfo = NULL;
  c24_nameCaptureInfo = NULL;
  sf_mex_assign(&c24_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), FALSE);
  return c24_nameCaptureInfo;
}

static const mxArray *c24_e_sf_marshallOut(void *chartInstanceVoid, void
  *c24_inData)
{
  const mxArray *c24_mxArrayOutData = NULL;
  int32_T c24_u;
  const mxArray *c24_y = NULL;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_mxArrayOutData = NULL;
  c24_u = *(int32_T *)c24_inData;
  c24_y = NULL;
  sf_mex_assign(&c24_y, sf_mex_create("y", &c24_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c24_mxArrayOutData, c24_y, FALSE);
  return c24_mxArrayOutData;
}

static int32_T c24_f_emlrt_marshallIn(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId)
{
  int32_T c24_y;
  int32_T c24_i30;
  sf_mex_import(c24_parentId, sf_mex_dup(c24_u), &c24_i30, 1, 6, 0U, 0, 0U, 0);
  c24_y = c24_i30;
  sf_mex_destroy(&c24_u);
  return c24_y;
}

static void c24_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c24_mxArrayInData, const char_T *c24_varName, void *c24_outData)
{
  const mxArray *c24_b_sfEvent;
  const char_T *c24_identifier;
  emlrtMsgIdentifier c24_thisId;
  int32_T c24_y;
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)chartInstanceVoid;
  c24_b_sfEvent = sf_mex_dup(c24_mxArrayInData);
  c24_identifier = c24_varName;
  c24_thisId.fIdentifier = c24_identifier;
  c24_thisId.fParent = NULL;
  c24_y = c24_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c24_b_sfEvent),
    &c24_thisId);
  sf_mex_destroy(&c24_b_sfEvent);
  *(int32_T *)c24_outData = c24_y;
  sf_mex_destroy(&c24_mxArrayInData);
}

static uint8_T c24_g_emlrt_marshallIn(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_b_is_active_c24_simulation, const char_T
  *c24_identifier)
{
  uint8_T c24_y;
  emlrtMsgIdentifier c24_thisId;
  c24_thisId.fIdentifier = c24_identifier;
  c24_thisId.fParent = NULL;
  c24_y = c24_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c24_b_is_active_c24_simulation), &c24_thisId);
  sf_mex_destroy(&c24_b_is_active_c24_simulation);
  return c24_y;
}

static uint8_T c24_h_emlrt_marshallIn(SFc24_simulationInstanceStruct
  *chartInstance, const mxArray *c24_u, const emlrtMsgIdentifier *c24_parentId)
{
  uint8_T c24_y;
  uint8_T c24_u0;
  sf_mex_import(c24_parentId, sf_mex_dup(c24_u), &c24_u0, 1, 3, 0U, 0, 0U, 0);
  c24_y = c24_u0;
  sf_mex_destroy(&c24_u);
  return c24_y;
}

static void init_dsm_address_info(SFc24_simulationInstanceStruct *chartInstance)
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

void sf_c24_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1460442087U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1535620015U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2703913299U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(714419128U);
}

mxArray *sf_c24_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("lqOmPKQg6fA36l4NJ2CljF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,3,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(6);
      pr[2] = (double)(4);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c24_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c24_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c24_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x6'type','srcId','name','auxInfo'{{M[1],M[23],T\"Circle\",},{M[1],M[17],T\"P0\",},{M[1],M[9],T\"gamWP\",},{M[1],M[16],T\"thetad\",},{M[1],M[18],T\"v_L\",},{M[8],M[0],T\"is_active_c24_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 6, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c24_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc24_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc24_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           24,
           1,
           1,
           8,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"pathdata");
          _SFD_SET_DATA_PROPS(1,1,1,0,"gam");
          _SFD_SET_DATA_PROPS(2,2,0,1,"gamWP");
          _SFD_SET_DATA_PROPS(3,2,0,1,"thetad");
          _SFD_SET_DATA_PROPS(4,2,0,1,"P0");
          _SFD_SET_DATA_PROPS(5,2,0,1,"v_L");
          _SFD_SET_DATA_PROPS(6,1,1,0,"id");
          _SFD_SET_DATA_PROPS(7,2,0,1,"Circle");
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
        _SFD_CV_INIT_EML(0,1,1,5,0,0,0,0,0,10,5);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,741);
        _SFD_CV_INIT_EML_IF(0,1,0,64,105,145,572);
        _SFD_CV_INIT_EML_IF(0,1,1,145,201,241,572);
        _SFD_CV_INIT_EML_IF(0,1,2,241,297,337,572);
        _SFD_CV_INIT_EML_IF(0,1,3,337,393,433,572);
        _SFD_CV_INIT_EML_IF(0,1,4,433,489,529,572);

        {
          static int condStart[] = { 72, 83 };

          static int condEnd[] = { 79, 105 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,0,72,105,2,0,&(condStart[0]),&(condEnd[0]),3,
                                &(pfixExpr[0]));
        }

        {
          static int condStart[] = { 153, 179 };

          static int condEnd[] = { 175, 201 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,1,153,201,2,2,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 249, 275 };

          static int condEnd[] = { 271, 297 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,2,249,297,2,4,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 345, 371 };

          static int condEnd[] = { 367, 393 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,3,345,393,2,6,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 441, 467 };

          static int condEnd[] = { 463, 489 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,4,441,489,2,8,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          unsigned int dimVector[3];
          dimVector[0]= 6;
          dimVector[1]= 6;
          dimVector[2]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,3,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c24_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c24_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c24_sf_marshallOut,(MexInFcnForType)c24_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c24_sf_marshallOut,(MexInFcnForType)c24_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c24_b_sf_marshallOut,(MexInFcnForType)
            c24_b_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c24_sf_marshallOut,(MexInFcnForType)c24_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c24_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c24_sf_marshallOut,(MexInFcnForType)c24_sf_marshallIn);

        {
          real_T *c24_gam;
          real_T *c24_gamWP;
          real_T *c24_thetad;
          real_T *c24_v_L;
          real_T *c24_id;
          real_T *c24_Circle;
          real_T (*c24_pathdata)[144];
          real_T (*c24_P0)[2];
          c24_Circle = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
          c24_id = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c24_v_L = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
          c24_P0 = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
          c24_thetad = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c24_gamWP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          c24_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c24_pathdata = (real_T (*)[144])ssGetInputPortSignal(chartInstance->S,
            0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c24_pathdata);
          _SFD_SET_DATA_VALUE_PTR(1U, c24_gam);
          _SFD_SET_DATA_VALUE_PTR(2U, c24_gamWP);
          _SFD_SET_DATA_VALUE_PTR(3U, c24_thetad);
          _SFD_SET_DATA_VALUE_PTR(4U, *c24_P0);
          _SFD_SET_DATA_VALUE_PTR(5U, c24_v_L);
          _SFD_SET_DATA_VALUE_PTR(6U, c24_id);
          _SFD_SET_DATA_VALUE_PTR(7U, c24_Circle);
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
  return "K4ZOrZMb09CExf3z87c9rF";
}

static void sf_opaque_initialize_c24_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc24_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c24_simulation((SFc24_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c24_simulation((SFc24_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c24_simulation(void *chartInstanceVar)
{
  enable_c24_simulation((SFc24_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c24_simulation(void *chartInstanceVar)
{
  disable_c24_simulation((SFc24_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c24_simulation(void *chartInstanceVar)
{
  sf_c24_simulation((SFc24_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c24_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c24_simulation
    ((SFc24_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c24_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c24_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c24_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c24_simulation((SFc24_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c24_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c24_simulation(S);
}

static void sf_opaque_set_sim_state_c24_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c24_simulation(S, st);
}

static void sf_opaque_terminate_c24_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc24_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c24_simulation((SFc24_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc24_simulation((SFc24_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c24_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c24_simulation((SFc24_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c24_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      24);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,24,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,24,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,24);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,24,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,24,5);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=5; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,24);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3441174101U));
  ssSetChecksum1(S,(1112925735U));
  ssSetChecksum2(S,(3434207747U));
  ssSetChecksum3(S,(2156099875U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c24_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c24_simulation(SimStruct *S)
{
  SFc24_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc24_simulationInstanceStruct *)utMalloc(sizeof
    (SFc24_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc24_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c24_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c24_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c24_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c24_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c24_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c24_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c24_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c24_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c24_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c24_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c24_simulation;
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

void c24_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c24_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c24_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c24_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c24_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
