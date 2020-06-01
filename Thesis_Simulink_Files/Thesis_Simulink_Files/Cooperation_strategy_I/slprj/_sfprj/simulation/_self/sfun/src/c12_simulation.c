/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c12_simulation.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simulation_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c12_debug_family_names[11] = { "pdataiStep", "nargin",
  "nargout", "pathdata", "gam", "id", "gamWP", "thetad", "P0", "v_L", "Circle" };

/* Function Declarations */
static void initialize_c12_simulation(SFc12_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c12_simulation(SFc12_simulationInstanceStruct
  *chartInstance);
static void enable_c12_simulation(SFc12_simulationInstanceStruct *chartInstance);
static void disable_c12_simulation(SFc12_simulationInstanceStruct *chartInstance);
static void c12_update_debugger_state_c12_simulation
  (SFc12_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c12_simulation
  (SFc12_simulationInstanceStruct *chartInstance);
static void set_sim_state_c12_simulation(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_st);
static void finalize_c12_simulation(SFc12_simulationInstanceStruct
  *chartInstance);
static void sf_c12_simulation(SFc12_simulationInstanceStruct *chartInstance);
static void initSimStructsc12_simulation(SFc12_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c12_machineNumber, uint32_T
  c12_chartNumber);
static const mxArray *c12_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static real_T c12_emlrt_marshallIn(SFc12_simulationInstanceStruct *chartInstance,
  const mxArray *c12_Circle, const char_T *c12_identifier);
static real_T c12_b_emlrt_marshallIn(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId);
static void c12_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static const mxArray *c12_b_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static void c12_c_emlrt_marshallIn(SFc12_simulationInstanceStruct *chartInstance,
  const mxArray *c12_P0, const char_T *c12_identifier, real_T c12_y[2]);
static void c12_d_emlrt_marshallIn(SFc12_simulationInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[2]);
static void c12_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static const mxArray *c12_c_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static const mxArray *c12_d_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static void c12_e_emlrt_marshallIn(SFc12_simulationInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[6]);
static void c12_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static const mxArray *c12_e_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static int32_T c12_f_emlrt_marshallIn(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId);
static void c12_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static uint8_T c12_g_emlrt_marshallIn(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_b_is_active_c12_simulation, const char_T
  *c12_identifier);
static uint8_T c12_h_emlrt_marshallIn(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId);
static void init_dsm_address_info(SFc12_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c12_simulation(SFc12_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c12_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c12_is_active_c12_simulation = 0U;
}

static void initialize_params_c12_simulation(SFc12_simulationInstanceStruct
  *chartInstance)
{
}

static void enable_c12_simulation(SFc12_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c12_simulation(SFc12_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c12_update_debugger_state_c12_simulation
  (SFc12_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c12_simulation
  (SFc12_simulationInstanceStruct *chartInstance)
{
  const mxArray *c12_st;
  const mxArray *c12_y = NULL;
  real_T c12_hoistedGlobal;
  real_T c12_u;
  const mxArray *c12_b_y = NULL;
  int32_T c12_i0;
  real_T c12_b_u[2];
  const mxArray *c12_c_y = NULL;
  real_T c12_b_hoistedGlobal;
  real_T c12_c_u;
  const mxArray *c12_d_y = NULL;
  real_T c12_c_hoistedGlobal;
  real_T c12_d_u;
  const mxArray *c12_e_y = NULL;
  real_T c12_d_hoistedGlobal;
  real_T c12_e_u;
  const mxArray *c12_f_y = NULL;
  uint8_T c12_e_hoistedGlobal;
  uint8_T c12_f_u;
  const mxArray *c12_g_y = NULL;
  real_T *c12_Circle;
  real_T *c12_gamWP;
  real_T *c12_thetad;
  real_T *c12_v_L;
  real_T (*c12_P0)[2];
  c12_Circle = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c12_v_L = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c12_P0 = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c12_thetad = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c12_gamWP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c12_st = NULL;
  c12_st = NULL;
  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_createcellarray(6), FALSE);
  c12_hoistedGlobal = *c12_Circle;
  c12_u = c12_hoistedGlobal;
  c12_b_y = NULL;
  sf_mex_assign(&c12_b_y, sf_mex_create("y", &c12_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c12_y, 0, c12_b_y);
  for (c12_i0 = 0; c12_i0 < 2; c12_i0++) {
    c12_b_u[c12_i0] = (*c12_P0)[c12_i0];
  }

  c12_c_y = NULL;
  sf_mex_assign(&c12_c_y, sf_mex_create("y", c12_b_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c12_y, 1, c12_c_y);
  c12_b_hoistedGlobal = *c12_gamWP;
  c12_c_u = c12_b_hoistedGlobal;
  c12_d_y = NULL;
  sf_mex_assign(&c12_d_y, sf_mex_create("y", &c12_c_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c12_y, 2, c12_d_y);
  c12_c_hoistedGlobal = *c12_thetad;
  c12_d_u = c12_c_hoistedGlobal;
  c12_e_y = NULL;
  sf_mex_assign(&c12_e_y, sf_mex_create("y", &c12_d_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c12_y, 3, c12_e_y);
  c12_d_hoistedGlobal = *c12_v_L;
  c12_e_u = c12_d_hoistedGlobal;
  c12_f_y = NULL;
  sf_mex_assign(&c12_f_y, sf_mex_create("y", &c12_e_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c12_y, 4, c12_f_y);
  c12_e_hoistedGlobal = chartInstance->c12_is_active_c12_simulation;
  c12_f_u = c12_e_hoistedGlobal;
  c12_g_y = NULL;
  sf_mex_assign(&c12_g_y, sf_mex_create("y", &c12_f_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c12_y, 5, c12_g_y);
  sf_mex_assign(&c12_st, c12_y, FALSE);
  return c12_st;
}

static void set_sim_state_c12_simulation(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_st)
{
  const mxArray *c12_u;
  real_T c12_dv0[2];
  int32_T c12_i1;
  real_T *c12_Circle;
  real_T *c12_gamWP;
  real_T *c12_thetad;
  real_T *c12_v_L;
  real_T (*c12_P0)[2];
  c12_Circle = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c12_v_L = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c12_P0 = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c12_thetad = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c12_gamWP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c12_doneDoubleBufferReInit = TRUE;
  c12_u = sf_mex_dup(c12_st);
  *c12_Circle = c12_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c12_u, 0)), "Circle");
  c12_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c12_u, 1)),
    "P0", c12_dv0);
  for (c12_i1 = 0; c12_i1 < 2; c12_i1++) {
    (*c12_P0)[c12_i1] = c12_dv0[c12_i1];
  }

  *c12_gamWP = c12_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c12_u, 2)), "gamWP");
  *c12_thetad = c12_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c12_u, 3)), "thetad");
  *c12_v_L = c12_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c12_u,
    4)), "v_L");
  chartInstance->c12_is_active_c12_simulation = c12_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c12_u, 5)),
     "is_active_c12_simulation");
  sf_mex_destroy(&c12_u);
  c12_update_debugger_state_c12_simulation(chartInstance);
  sf_mex_destroy(&c12_st);
}

static void finalize_c12_simulation(SFc12_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c12_simulation(SFc12_simulationInstanceStruct *chartInstance)
{
  int32_T c12_i2;
  int32_T c12_i3;
  real_T c12_hoistedGlobal;
  real_T c12_b_hoistedGlobal;
  int32_T c12_i4;
  real_T c12_pathdata[144];
  real_T c12_gam;
  real_T c12_id;
  uint32_T c12_debug_family_var_map[11];
  real_T c12_pdataiStep[6];
  real_T c12_nargin = 3.0;
  real_T c12_nargout = 5.0;
  real_T c12_gamWP;
  real_T c12_thetad;
  real_T c12_P0[2];
  real_T c12_v_L;
  real_T c12_Circle;
  int32_T c12_b_id;
  int32_T c12_i5;
  int32_T c12_c_id;
  int32_T c12_i6;
  int32_T c12_d_id;
  int32_T c12_i7;
  int32_T c12_e_id;
  int32_T c12_i8;
  int32_T c12_f_id;
  int32_T c12_i9;
  int32_T c12_g_id;
  int32_T c12_i10;
  int32_T c12_i11;
  real_T *c12_b_gam;
  real_T *c12_b_gamWP;
  real_T *c12_b_thetad;
  real_T *c12_b_v_L;
  real_T *c12_h_id;
  real_T *c12_b_Circle;
  real_T (*c12_b_P0)[2];
  real_T (*c12_b_pathdata)[144];
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  boolean_T guard3 = FALSE;
  boolean_T guard4 = FALSE;
  boolean_T guard5 = FALSE;
  c12_b_Circle = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c12_h_id = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c12_b_v_L = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c12_b_P0 = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c12_b_thetad = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c12_b_gamWP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c12_b_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c12_b_pathdata = (real_T (*)[144])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 11U, chartInstance->c12_sfEvent);
  for (c12_i2 = 0; c12_i2 < 144; c12_i2++) {
    _SFD_DATA_RANGE_CHECK((*c12_b_pathdata)[c12_i2], 0U);
  }

  _SFD_DATA_RANGE_CHECK(*c12_b_gam, 1U);
  _SFD_DATA_RANGE_CHECK(*c12_b_gamWP, 2U);
  _SFD_DATA_RANGE_CHECK(*c12_b_thetad, 3U);
  for (c12_i3 = 0; c12_i3 < 2; c12_i3++) {
    _SFD_DATA_RANGE_CHECK((*c12_b_P0)[c12_i3], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*c12_b_v_L, 5U);
  _SFD_DATA_RANGE_CHECK(*c12_h_id, 6U);
  _SFD_DATA_RANGE_CHECK(*c12_b_Circle, 7U);
  chartInstance->c12_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 11U, chartInstance->c12_sfEvent);
  c12_hoistedGlobal = *c12_b_gam;
  c12_b_hoistedGlobal = *c12_h_id;
  for (c12_i4 = 0; c12_i4 < 144; c12_i4++) {
    c12_pathdata[c12_i4] = (*c12_b_pathdata)[c12_i4];
  }

  c12_gam = c12_hoistedGlobal;
  c12_id = c12_b_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 11U, 11U, c12_debug_family_names,
    c12_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c12_pdataiStep, 0U, c12_d_sf_marshallOut,
    c12_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_nargin, 1U, c12_sf_marshallOut,
    c12_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_nargout, 2U, c12_sf_marshallOut,
    c12_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c12_pathdata, 3U, c12_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c12_gam, 4U, c12_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c12_id, 5U, c12_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_gamWP, 6U, c12_sf_marshallOut,
    c12_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_thetad, 7U, c12_sf_marshallOut,
    c12_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c12_P0, 8U, c12_b_sf_marshallOut,
    c12_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_v_L, 9U, c12_sf_marshallOut,
    c12_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_Circle, 10U, c12_sf_marshallOut,
    c12_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 3);
  guard1 = FALSE;
  if (CV_EML_COND(0, 1, 0, c12_gam >= 0.0)) {
    if (CV_EML_COND(0, 1, 1, c12_gam < c12_pathdata[6 + 36 * ((int32_T)(real_T)
          _SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)_SFD_INTEGER_CHECK(
            "id", c12_id), 1, 4, 3, 0) - 1)])) {
      CV_EML_MCDC(0, 1, 0, TRUE);
      CV_EML_IF(0, 1, 0, TRUE);
      _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 4);
      c12_b_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
        (int32_T)_SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1;
      for (c12_i5 = 0; c12_i5 < 6; c12_i5++) {
        c12_pdataiStep[c12_i5] = c12_pathdata[c12_i5 + 36 * c12_b_id];
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
    _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 5);
    guard2 = FALSE;
    if (CV_EML_COND(0, 1, 2, c12_gam >= c12_pathdata[6 + 36 * ((int32_T)(real_T)
          _SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)_SFD_INTEGER_CHECK(
            "id", c12_id), 1, 4, 3, 0) - 1)])) {
      if (CV_EML_COND(0, 1, 3, c12_gam < c12_pathdata[12 + 36 * ((int32_T)
            (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
             _SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1)])) {
        CV_EML_MCDC(0, 1, 1, TRUE);
        CV_EML_IF(0, 1, 1, TRUE);
        _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 6);
        c12_c_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
          (int32_T)_SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1;
        for (c12_i6 = 0; c12_i6 < 6; c12_i6++) {
          c12_pdataiStep[c12_i6] = c12_pathdata[6 + (c12_i6 + 36 * c12_c_id)];
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
      _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 7);
      guard3 = FALSE;
      if (CV_EML_COND(0, 1, 4, c12_gam >= c12_pathdata[12 + 36 * ((int32_T)
            (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
             _SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1)])) {
        if (CV_EML_COND(0, 1, 5, c12_gam < c12_pathdata[18 + 36 * ((int32_T)
              (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
               _SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1)])) {
          CV_EML_MCDC(0, 1, 2, TRUE);
          CV_EML_IF(0, 1, 2, TRUE);
          _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 8);
          c12_d_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
            (int32_T)_SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1;
          for (c12_i7 = 0; c12_i7 < 6; c12_i7++) {
            c12_pdataiStep[c12_i7] = c12_pathdata[12 + (c12_i7 + 36 * c12_d_id)];
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
        _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 9);
        guard4 = FALSE;
        if (CV_EML_COND(0, 1, 6, c12_gam >= c12_pathdata[18 + 36 * ((int32_T)
              (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
               _SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1)])) {
          if (CV_EML_COND(0, 1, 7, c12_gam < c12_pathdata[24 + 36 * ((int32_T)
                (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
                 _SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1)])) {
            CV_EML_MCDC(0, 1, 3, TRUE);
            CV_EML_IF(0, 1, 3, TRUE);
            _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 10);
            c12_e_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
              (int32_T)_SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1;
            for (c12_i8 = 0; c12_i8 < 6; c12_i8++) {
              c12_pdataiStep[c12_i8] = c12_pathdata[18 + (c12_i8 + 36 * c12_e_id)];
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
          _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 11);
          guard5 = FALSE;
          if (CV_EML_COND(0, 1, 8, c12_gam >= c12_pathdata[24 + 36 * ((int32_T)
                (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
                 _SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1)])) {
            if (CV_EML_COND(0, 1, 9, c12_gam < c12_pathdata[30 + 36 * ((int32_T)
                  (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata", (int32_T)
                   _SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1)])) {
              CV_EML_MCDC(0, 1, 4, TRUE);
              CV_EML_IF(0, 1, 4, TRUE);
              _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 12);
              c12_f_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
                (int32_T)_SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1;
              for (c12_i9 = 0; c12_i9 < 6; c12_i9++) {
                c12_pdataiStep[c12_i9] = c12_pathdata[24 + (c12_i9 + 36 *
                  c12_f_id)];
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
            _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 14);
            c12_g_id = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("pathdata",
              (int32_T)_SFD_INTEGER_CHECK("id", c12_id), 1, 4, 3, 0) - 1;
            for (c12_i10 = 0; c12_i10 < 6; c12_i10++) {
              c12_pdataiStep[c12_i10] = c12_pathdata[30 + (c12_i10 + 36 *
                c12_g_id)];
            }
          }
        }
      }
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 17);
  c12_gamWP = c12_pdataiStep[0];
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 18);
  c12_thetad = c12_pdataiStep[1];
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 19);
  c12_P0[0] = c12_pdataiStep[2];
  c12_P0[1] = c12_pdataiStep[3];
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 20);
  c12_v_L = c12_pdataiStep[4];
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 21);
  c12_Circle = c12_pdataiStep[5];
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, -21);
  _SFD_SYMBOL_SCOPE_POP();
  *c12_b_gamWP = c12_gamWP;
  *c12_b_thetad = c12_thetad;
  for (c12_i11 = 0; c12_i11 < 2; c12_i11++) {
    (*c12_b_P0)[c12_i11] = c12_P0[c12_i11];
  }

  *c12_b_v_L = c12_v_L;
  *c12_b_Circle = c12_Circle;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 11U, chartInstance->c12_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc12_simulation(SFc12_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c12_machineNumber, uint32_T
  c12_chartNumber)
{
}

static const mxArray *c12_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  real_T c12_u;
  const mxArray *c12_y = NULL;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  c12_u = *(real_T *)c12_inData;
  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", &c12_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static real_T c12_emlrt_marshallIn(SFc12_simulationInstanceStruct *chartInstance,
  const mxArray *c12_Circle, const char_T *c12_identifier)
{
  real_T c12_y;
  emlrtMsgIdentifier c12_thisId;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_y = c12_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_Circle),
    &c12_thisId);
  sf_mex_destroy(&c12_Circle);
  return c12_y;
}

static real_T c12_b_emlrt_marshallIn(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId)
{
  real_T c12_y;
  real_T c12_d0;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), &c12_d0, 1, 0, 0U, 0, 0U, 0);
  c12_y = c12_d0;
  sf_mex_destroy(&c12_u);
  return c12_y;
}

static void c12_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_Circle;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  real_T c12_y;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_Circle = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_y = c12_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_Circle),
    &c12_thisId);
  sf_mex_destroy(&c12_Circle);
  *(real_T *)c12_outData = c12_y;
  sf_mex_destroy(&c12_mxArrayInData);
}

static const mxArray *c12_b_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_i12;
  real_T c12_b_inData[2];
  int32_T c12_i13;
  real_T c12_u[2];
  const mxArray *c12_y = NULL;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  for (c12_i12 = 0; c12_i12 < 2; c12_i12++) {
    c12_b_inData[c12_i12] = (*(real_T (*)[2])c12_inData)[c12_i12];
  }

  for (c12_i13 = 0; c12_i13 < 2; c12_i13++) {
    c12_u[c12_i13] = c12_b_inData[c12_i13];
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static void c12_c_emlrt_marshallIn(SFc12_simulationInstanceStruct *chartInstance,
  const mxArray *c12_P0, const char_T *c12_identifier, real_T c12_y[2])
{
  emlrtMsgIdentifier c12_thisId;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_P0), &c12_thisId, c12_y);
  sf_mex_destroy(&c12_P0);
}

static void c12_d_emlrt_marshallIn(SFc12_simulationInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[2])
{
  real_T c12_dv1[2];
  int32_T c12_i14;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), c12_dv1, 1, 0, 0U, 1, 0U, 1, 2);
  for (c12_i14 = 0; c12_i14 < 2; c12_i14++) {
    c12_y[c12_i14] = c12_dv1[c12_i14];
  }

  sf_mex_destroy(&c12_u);
}

static void c12_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_P0;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  real_T c12_y[2];
  int32_T c12_i15;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_P0 = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_P0), &c12_thisId, c12_y);
  sf_mex_destroy(&c12_P0);
  for (c12_i15 = 0; c12_i15 < 2; c12_i15++) {
    (*(real_T (*)[2])c12_outData)[c12_i15] = c12_y[c12_i15];
  }

  sf_mex_destroy(&c12_mxArrayInData);
}

static const mxArray *c12_c_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_i16;
  int32_T c12_i17;
  int32_T c12_i18;
  int32_T c12_i19;
  int32_T c12_i20;
  real_T c12_b_inData[144];
  int32_T c12_i21;
  int32_T c12_i22;
  int32_T c12_i23;
  int32_T c12_i24;
  int32_T c12_i25;
  real_T c12_u[144];
  const mxArray *c12_y = NULL;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  c12_i16 = 0;
  for (c12_i17 = 0; c12_i17 < 4; c12_i17++) {
    c12_i18 = 0;
    for (c12_i19 = 0; c12_i19 < 6; c12_i19++) {
      for (c12_i20 = 0; c12_i20 < 6; c12_i20++) {
        c12_b_inData[(c12_i20 + c12_i18) + c12_i16] = (*(real_T (*)[144])
          c12_inData)[(c12_i20 + c12_i18) + c12_i16];
      }

      c12_i18 += 6;
    }

    c12_i16 += 36;
  }

  c12_i21 = 0;
  for (c12_i22 = 0; c12_i22 < 4; c12_i22++) {
    c12_i23 = 0;
    for (c12_i24 = 0; c12_i24 < 6; c12_i24++) {
      for (c12_i25 = 0; c12_i25 < 6; c12_i25++) {
        c12_u[(c12_i25 + c12_i23) + c12_i21] = c12_b_inData[(c12_i25 + c12_i23)
          + c12_i21];
      }

      c12_i23 += 6;
    }

    c12_i21 += 36;
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 0, 0U, 1U, 0U, 3, 6, 6, 4),
                FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static const mxArray *c12_d_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_i26;
  real_T c12_b_inData[6];
  int32_T c12_i27;
  real_T c12_u[6];
  const mxArray *c12_y = NULL;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  for (c12_i26 = 0; c12_i26 < 6; c12_i26++) {
    c12_b_inData[c12_i26] = (*(real_T (*)[6])c12_inData)[c12_i26];
  }

  for (c12_i27 = 0; c12_i27 < 6; c12_i27++) {
    c12_u[c12_i27] = c12_b_inData[c12_i27];
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static void c12_e_emlrt_marshallIn(SFc12_simulationInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[6])
{
  real_T c12_dv2[6];
  int32_T c12_i28;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), c12_dv2, 1, 0, 0U, 1, 0U, 1, 6);
  for (c12_i28 = 0; c12_i28 < 6; c12_i28++) {
    c12_y[c12_i28] = c12_dv2[c12_i28];
  }

  sf_mex_destroy(&c12_u);
}

static void c12_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_pdataiStep;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  real_T c12_y[6];
  int32_T c12_i29;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_pdataiStep = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_pdataiStep), &c12_thisId,
    c12_y);
  sf_mex_destroy(&c12_pdataiStep);
  for (c12_i29 = 0; c12_i29 < 6; c12_i29++) {
    (*(real_T (*)[6])c12_outData)[c12_i29] = c12_y[c12_i29];
  }

  sf_mex_destroy(&c12_mxArrayInData);
}

const mxArray *sf_c12_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c12_nameCaptureInfo = NULL;
  c12_nameCaptureInfo = NULL;
  sf_mex_assign(&c12_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), FALSE);
  return c12_nameCaptureInfo;
}

static const mxArray *c12_e_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_u;
  const mxArray *c12_y = NULL;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  c12_u = *(int32_T *)c12_inData;
  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", &c12_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static int32_T c12_f_emlrt_marshallIn(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId)
{
  int32_T c12_y;
  int32_T c12_i30;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), &c12_i30, 1, 6, 0U, 0, 0U, 0);
  c12_y = c12_i30;
  sf_mex_destroy(&c12_u);
  return c12_y;
}

static void c12_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_b_sfEvent;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  int32_T c12_y;
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)chartInstanceVoid;
  c12_b_sfEvent = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_y = c12_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_b_sfEvent),
    &c12_thisId);
  sf_mex_destroy(&c12_b_sfEvent);
  *(int32_T *)c12_outData = c12_y;
  sf_mex_destroy(&c12_mxArrayInData);
}

static uint8_T c12_g_emlrt_marshallIn(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_b_is_active_c12_simulation, const char_T
  *c12_identifier)
{
  uint8_T c12_y;
  emlrtMsgIdentifier c12_thisId;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_y = c12_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c12_b_is_active_c12_simulation), &c12_thisId);
  sf_mex_destroy(&c12_b_is_active_c12_simulation);
  return c12_y;
}

static uint8_T c12_h_emlrt_marshallIn(SFc12_simulationInstanceStruct
  *chartInstance, const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId)
{
  uint8_T c12_y;
  uint8_T c12_u0;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), &c12_u0, 1, 3, 0U, 0, 0U, 0);
  c12_y = c12_u0;
  sf_mex_destroy(&c12_u);
  return c12_y;
}

static void init_dsm_address_info(SFc12_simulationInstanceStruct *chartInstance)
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

void sf_c12_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1460442087U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1535620015U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2703913299U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(714419128U);
}

mxArray *sf_c12_simulation_get_autoinheritance_info(void)
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

mxArray *sf_c12_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c12_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c12_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x6'type','srcId','name','auxInfo'{{M[1],M[23],T\"Circle\",},{M[1],M[17],T\"P0\",},{M[1],M[9],T\"gamWP\",},{M[1],M[16],T\"thetad\",},{M[1],M[18],T\"v_L\",},{M[8],M[0],T\"is_active_c12_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 6, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c12_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc12_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc12_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           12,
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
            1.0,0,0,(MexFcnForType)c12_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c12_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c12_sf_marshallOut,(MexInFcnForType)c12_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c12_sf_marshallOut,(MexInFcnForType)c12_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c12_b_sf_marshallOut,(MexInFcnForType)
            c12_b_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c12_sf_marshallOut,(MexInFcnForType)c12_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c12_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c12_sf_marshallOut,(MexInFcnForType)c12_sf_marshallIn);

        {
          real_T *c12_gam;
          real_T *c12_gamWP;
          real_T *c12_thetad;
          real_T *c12_v_L;
          real_T *c12_id;
          real_T *c12_Circle;
          real_T (*c12_pathdata)[144];
          real_T (*c12_P0)[2];
          c12_Circle = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
          c12_id = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c12_v_L = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
          c12_P0 = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
          c12_thetad = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c12_gamWP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          c12_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c12_pathdata = (real_T (*)[144])ssGetInputPortSignal(chartInstance->S,
            0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c12_pathdata);
          _SFD_SET_DATA_VALUE_PTR(1U, c12_gam);
          _SFD_SET_DATA_VALUE_PTR(2U, c12_gamWP);
          _SFD_SET_DATA_VALUE_PTR(3U, c12_thetad);
          _SFD_SET_DATA_VALUE_PTR(4U, *c12_P0);
          _SFD_SET_DATA_VALUE_PTR(5U, c12_v_L);
          _SFD_SET_DATA_VALUE_PTR(6U, c12_id);
          _SFD_SET_DATA_VALUE_PTR(7U, c12_Circle);
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

static void sf_opaque_initialize_c12_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc12_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c12_simulation((SFc12_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c12_simulation((SFc12_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c12_simulation(void *chartInstanceVar)
{
  enable_c12_simulation((SFc12_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c12_simulation(void *chartInstanceVar)
{
  disable_c12_simulation((SFc12_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c12_simulation(void *chartInstanceVar)
{
  sf_c12_simulation((SFc12_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c12_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c12_simulation
    ((SFc12_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c12_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c12_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c12_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c12_simulation((SFc12_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c12_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c12_simulation(S);
}

static void sf_opaque_set_sim_state_c12_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c12_simulation(S, st);
}

static void sf_opaque_terminate_c12_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc12_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c12_simulation((SFc12_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc12_simulation((SFc12_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c12_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c12_simulation((SFc12_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c12_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      12);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,12,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,12,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,12);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,12,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,12,5);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,12);
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

static void mdlRTW_c12_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c12_simulation(SimStruct *S)
{
  SFc12_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc12_simulationInstanceStruct *)utMalloc(sizeof
    (SFc12_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc12_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c12_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c12_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c12_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c12_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c12_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c12_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c12_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c12_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c12_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c12_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c12_simulation;
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

void c12_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c12_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c12_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c12_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c12_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
