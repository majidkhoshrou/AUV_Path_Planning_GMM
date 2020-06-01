/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c21_simulation.h"
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
static const char * c21_debug_family_names[13] = { "R0", "R", "alpha0", "nargin",
  "nargout", "gam", "gamWP", "thetad", "P0", "v_L", "Radius", "Pd", "dPd" };

/* Function Declarations */
static void initialize_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance);
static void enable_c21_simulation(SFc21_simulationInstanceStruct *chartInstance);
static void disable_c21_simulation(SFc21_simulationInstanceStruct *chartInstance);
static void c21_update_debugger_state_c21_simulation
  (SFc21_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c21_simulation
  (SFc21_simulationInstanceStruct *chartInstance);
static void set_sim_state_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_st);
static void finalize_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance);
static void sf_c21_simulation(SFc21_simulationInstanceStruct *chartInstance);
static void c21_chartstep_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc21_simulation(SFc21_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c21_machineNumber, uint32_T
  c21_chartNumber);
static const mxArray *c21_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static void c21_emlrt_marshallIn(SFc21_simulationInstanceStruct *chartInstance,
  const mxArray *c21_dPd, const char_T *c21_identifier, real_T c21_y[2]);
static void c21_b_emlrt_marshallIn(SFc21_simulationInstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[2]);
static void c21_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static const mxArray *c21_b_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static real_T c21_c_emlrt_marshallIn(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId);
static void c21_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static void c21_info_helper(const mxArray **c21_info);
static const mxArray *c21_emlrt_marshallOut(char * c21_u);
static const mxArray *c21_b_emlrt_marshallOut(uint32_T c21_u);
static const mxArray *c21_c_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData);
static int32_T c21_d_emlrt_marshallIn(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId);
static void c21_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData);
static uint8_T c21_e_emlrt_marshallIn(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_b_is_active_c21_simulation, const char_T
  *c21_identifier);
static uint8_T c21_f_emlrt_marshallIn(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId);
static void init_dsm_address_info(SFc21_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c21_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c21_is_active_c21_simulation = 0U;
}

static void initialize_params_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance)
{
}

static void enable_c21_simulation(SFc21_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c21_simulation(SFc21_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c21_update_debugger_state_c21_simulation
  (SFc21_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c21_simulation
  (SFc21_simulationInstanceStruct *chartInstance)
{
  const mxArray *c21_st;
  const mxArray *c21_y = NULL;
  int32_T c21_i0;
  real_T c21_u[2];
  const mxArray *c21_b_y = NULL;
  int32_T c21_i1;
  real_T c21_b_u[2];
  const mxArray *c21_c_y = NULL;
  uint8_T c21_hoistedGlobal;
  uint8_T c21_c_u;
  const mxArray *c21_d_y = NULL;
  real_T (*c21_dPd)[2];
  real_T (*c21_Pd)[2];
  c21_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c21_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c21_st = NULL;
  c21_st = NULL;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_createcellarray(3), FALSE);
  for (c21_i0 = 0; c21_i0 < 2; c21_i0++) {
    c21_u[c21_i0] = (*c21_Pd)[c21_i0];
  }

  c21_b_y = NULL;
  sf_mex_assign(&c21_b_y, sf_mex_create("y", c21_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_setcell(c21_y, 0, c21_b_y);
  for (c21_i1 = 0; c21_i1 < 2; c21_i1++) {
    c21_b_u[c21_i1] = (*c21_dPd)[c21_i1];
  }

  c21_c_y = NULL;
  sf_mex_assign(&c21_c_y, sf_mex_create("y", c21_b_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c21_y, 1, c21_c_y);
  c21_hoistedGlobal = chartInstance->c21_is_active_c21_simulation;
  c21_c_u = c21_hoistedGlobal;
  c21_d_y = NULL;
  sf_mex_assign(&c21_d_y, sf_mex_create("y", &c21_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c21_y, 2, c21_d_y);
  sf_mex_assign(&c21_st, c21_y, FALSE);
  return c21_st;
}

static void set_sim_state_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_st)
{
  const mxArray *c21_u;
  real_T c21_dv0[2];
  int32_T c21_i2;
  real_T c21_dv1[2];
  int32_T c21_i3;
  real_T (*c21_Pd)[2];
  real_T (*c21_dPd)[2];
  c21_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c21_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c21_doneDoubleBufferReInit = TRUE;
  c21_u = sf_mex_dup(c21_st);
  c21_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c21_u, 0)), "Pd",
                       c21_dv0);
  for (c21_i2 = 0; c21_i2 < 2; c21_i2++) {
    (*c21_Pd)[c21_i2] = c21_dv0[c21_i2];
  }

  c21_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c21_u, 1)),
                       "dPd", c21_dv1);
  for (c21_i3 = 0; c21_i3 < 2; c21_i3++) {
    (*c21_dPd)[c21_i3] = c21_dv1[c21_i3];
  }

  chartInstance->c21_is_active_c21_simulation = c21_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c21_u, 2)),
     "is_active_c21_simulation");
  sf_mex_destroy(&c21_u);
  c21_update_debugger_state_c21_simulation(chartInstance);
  sf_mex_destroy(&c21_st);
}

static void finalize_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c21_simulation(SFc21_simulationInstanceStruct *chartInstance)
{
  int32_T c21_i4;
  int32_T c21_i5;
  int32_T c21_i6;
  real_T *c21_gam;
  real_T *c21_gamWP;
  real_T *c21_thetad;
  real_T *c21_v_L;
  real_T *c21_Radius;
  real_T (*c21_P0)[2];
  real_T (*c21_dPd)[2];
  real_T (*c21_Pd)[2];
  c21_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c21_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c21_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c21_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c21_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c21_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c21_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c21_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 20U, chartInstance->c21_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c21_gam, 0U);
  for (c21_i4 = 0; c21_i4 < 2; c21_i4++) {
    _SFD_DATA_RANGE_CHECK((*c21_Pd)[c21_i4], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c21_gamWP, 2U);
  for (c21_i5 = 0; c21_i5 < 2; c21_i5++) {
    _SFD_DATA_RANGE_CHECK((*c21_dPd)[c21_i5], 3U);
  }

  _SFD_DATA_RANGE_CHECK(*c21_thetad, 4U);
  for (c21_i6 = 0; c21_i6 < 2; c21_i6++) {
    _SFD_DATA_RANGE_CHECK((*c21_P0)[c21_i6], 5U);
  }

  _SFD_DATA_RANGE_CHECK(*c21_v_L, 6U);
  _SFD_DATA_RANGE_CHECK(*c21_Radius, 7U);
  chartInstance->c21_sfEvent = CALL_EVENT;
  c21_chartstep_c21_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c21_chartstep_c21_simulation(SFc21_simulationInstanceStruct
  *chartInstance)
{
  real_T c21_hoistedGlobal;
  real_T c21_b_hoistedGlobal;
  real_T c21_c_hoistedGlobal;
  real_T c21_d_hoistedGlobal;
  real_T c21_e_hoistedGlobal;
  real_T c21_gam;
  real_T c21_gamWP;
  real_T c21_thetad;
  int32_T c21_i7;
  real_T c21_P0[2];
  real_T c21_v_L;
  real_T c21_Radius;
  uint32_T c21_debug_family_var_map[13];
  real_T c21_R0[2];
  real_T c21_R;
  real_T c21_alpha0;
  real_T c21_nargin = 6.0;
  real_T c21_nargout = 2.0;
  real_T c21_Pd[2];
  real_T c21_dPd[2];
  int32_T c21_i8;
  int32_T c21_i9;
  real_T c21_a;
  real_T c21_b;
  real_T c21_y;
  real_T c21_A;
  real_T c21_B;
  real_T c21_x;
  real_T c21_b_y;
  real_T c21_b_x;
  real_T c21_c_y;
  real_T c21_d_y;
  real_T c21_c_x;
  real_T c21_d_x;
  real_T c21_b_a;
  real_T c21_b_b;
  real_T c21_e_y;
  real_T c21_c_a;
  real_T c21_c_b;
  real_T c21_f_y;
  real_T c21_b_A;
  real_T c21_b_B;
  real_T c21_e_x;
  real_T c21_g_y;
  real_T c21_f_x;
  real_T c21_h_y;
  real_T c21_i_y;
  real_T c21_g_x;
  real_T c21_h_x;
  real_T c21_d_a;
  real_T c21_d_b;
  real_T c21_j_y;
  real_T c21_e_a;
  real_T c21_e_b;
  real_T c21_k_y;
  real_T c21_c_A;
  real_T c21_c_B;
  real_T c21_i_x;
  real_T c21_l_y;
  real_T c21_j_x;
  real_T c21_m_y;
  real_T c21_n_y;
  real_T c21_k_x;
  real_T c21_l_x;
  real_T c21_f_a;
  real_T c21_f_b;
  real_T c21_o_y;
  real_T c21_g_a;
  real_T c21_g_b;
  real_T c21_p_y;
  real_T c21_d_A;
  real_T c21_d_B;
  real_T c21_m_x;
  real_T c21_q_y;
  real_T c21_n_x;
  real_T c21_r_y;
  real_T c21_s_y;
  real_T c21_o_x;
  real_T c21_p_x;
  real_T c21_h_a;
  real_T c21_h_b;
  real_T c21_t_y;
  real_T c21_u_y[2];
  int32_T c21_i10;
  int32_T c21_i11;
  real_T c21_q_x;
  real_T c21_r_x;
  real_T c21_i_a;
  real_T c21_i_b;
  real_T c21_v_y;
  real_T c21_s_x;
  real_T c21_t_x;
  real_T c21_j_a;
  real_T c21_j_b;
  real_T c21_w_y;
  int32_T c21_i12;
  real_T c21_k_a[2];
  real_T c21_k_b;
  int32_T c21_i13;
  int32_T c21_i14;
  int32_T c21_i15;
  int32_T c21_i16;
  real_T *c21_b_Radius;
  real_T *c21_b_v_L;
  real_T *c21_b_thetad;
  real_T *c21_b_gamWP;
  real_T *c21_b_gam;
  real_T (*c21_b_Pd)[2];
  real_T (*c21_b_dPd)[2];
  real_T (*c21_b_P0)[2];
  c21_b_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c21_b_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c21_b_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c21_b_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c21_b_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c21_b_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c21_b_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c21_b_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 20U, chartInstance->c21_sfEvent);
  c21_hoistedGlobal = *c21_b_gam;
  c21_b_hoistedGlobal = *c21_b_gamWP;
  c21_c_hoistedGlobal = *c21_b_thetad;
  c21_d_hoistedGlobal = *c21_b_v_L;
  c21_e_hoistedGlobal = *c21_b_Radius;
  c21_gam = c21_hoistedGlobal;
  c21_gamWP = c21_b_hoistedGlobal;
  c21_thetad = c21_c_hoistedGlobal;
  for (c21_i7 = 0; c21_i7 < 2; c21_i7++) {
    c21_P0[c21_i7] = (*c21_b_P0)[c21_i7];
  }

  c21_v_L = c21_d_hoistedGlobal;
  c21_Radius = c21_e_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 13U, 13U, c21_debug_family_names,
    c21_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_R0, 0U, c21_sf_marshallOut,
    c21_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_R, 1U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_alpha0, 2U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_nargin, 3U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c21_nargout, 4U, c21_b_sf_marshallOut,
    c21_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c21_gam, 5U, c21_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c21_gamWP, 6U, c21_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c21_thetad, 7U, c21_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c21_P0, 8U, c21_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c21_v_L, 9U, c21_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c21_Radius, 10U, c21_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_Pd, 11U, c21_sf_marshallOut,
    c21_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c21_dPd, 12U, c21_sf_marshallOut,
    c21_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 3);
  if (CV_EML_IF(0, 1, 0, c21_Radius != 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 4);
    for (c21_i8 = 0; c21_i8 < 2; c21_i8++) {
      c21_R0[c21_i8] = c21_P0[c21_i8];
    }

    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 5);
    c21_R = c21_Radius;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 6);
    c21_alpha0 = c21_thetad;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 7);
    for (c21_i9 = 0; c21_i9 < 2; c21_i9++) {
      c21_dPd[c21_i9] = 0.0;
    }

    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 8);
    c21_a = c21_gam - c21_gamWP;
    c21_b = c21_v_L;
    c21_y = c21_a * c21_b;
    c21_A = c21_y;
    c21_B = c21_R;
    c21_x = c21_A;
    c21_b_y = c21_B;
    c21_b_x = c21_x;
    c21_c_y = c21_b_y;
    c21_d_y = c21_b_x / c21_c_y;
    c21_c_x = c21_alpha0 + c21_d_y;
    c21_d_x = c21_c_x;
    c21_d_x = muDoubleScalarSin(c21_d_x);
    c21_b_a = -c21_d_x;
    c21_b_b = c21_v_L;
    c21_e_y = c21_b_a * c21_b_b;
    c21_dPd[0] = c21_e_y;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 9);
    c21_c_a = c21_gam - c21_gamWP;
    c21_c_b = c21_v_L;
    c21_f_y = c21_c_a * c21_c_b;
    c21_b_A = c21_f_y;
    c21_b_B = c21_R;
    c21_e_x = c21_b_A;
    c21_g_y = c21_b_B;
    c21_f_x = c21_e_x;
    c21_h_y = c21_g_y;
    c21_i_y = c21_f_x / c21_h_y;
    c21_g_x = c21_alpha0 + c21_i_y;
    c21_h_x = c21_g_x;
    c21_h_x = muDoubleScalarCos(c21_h_x);
    c21_d_a = c21_h_x;
    c21_d_b = c21_v_L;
    c21_j_y = c21_d_a * c21_d_b;
    c21_dPd[1] = c21_j_y;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 10);
    c21_e_a = c21_gam - c21_gamWP;
    c21_e_b = c21_v_L;
    c21_k_y = c21_e_a * c21_e_b;
    c21_c_A = c21_k_y;
    c21_c_B = c21_R;
    c21_i_x = c21_c_A;
    c21_l_y = c21_c_B;
    c21_j_x = c21_i_x;
    c21_m_y = c21_l_y;
    c21_n_y = c21_j_x / c21_m_y;
    c21_k_x = c21_alpha0 + c21_n_y;
    c21_l_x = c21_k_x;
    c21_l_x = muDoubleScalarCos(c21_l_x);
    c21_f_a = c21_l_x;
    c21_f_b = c21_R;
    c21_o_y = c21_f_a * c21_f_b;
    c21_g_a = c21_gam - c21_gamWP;
    c21_g_b = c21_v_L;
    c21_p_y = c21_g_a * c21_g_b;
    c21_d_A = c21_p_y;
    c21_d_B = c21_R;
    c21_m_x = c21_d_A;
    c21_q_y = c21_d_B;
    c21_n_x = c21_m_x;
    c21_r_y = c21_q_y;
    c21_s_y = c21_n_x / c21_r_y;
    c21_o_x = c21_alpha0 + c21_s_y;
    c21_p_x = c21_o_x;
    c21_p_x = muDoubleScalarSin(c21_p_x);
    c21_h_a = c21_p_x;
    c21_h_b = c21_R;
    c21_t_y = c21_h_a * c21_h_b;
    c21_u_y[0] = c21_o_y;
    c21_u_y[1] = c21_t_y;
    for (c21_i10 = 0; c21_i10 < 2; c21_i10++) {
      c21_Pd[c21_i10] = c21_R0[c21_i10] + c21_u_y[c21_i10];
    }
  } else {
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 13);
    for (c21_i11 = 0; c21_i11 < 2; c21_i11++) {
      c21_dPd[c21_i11] = 0.0;
    }

    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 14);
    c21_q_x = c21_thetad;
    c21_r_x = c21_q_x;
    c21_r_x = muDoubleScalarCos(c21_r_x);
    c21_i_a = c21_r_x;
    c21_i_b = c21_v_L;
    c21_v_y = c21_i_a * c21_i_b;
    c21_dPd[0] = c21_v_y;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 15);
    c21_s_x = c21_thetad;
    c21_t_x = c21_s_x;
    c21_t_x = muDoubleScalarSin(c21_t_x);
    c21_j_a = c21_t_x;
    c21_j_b = c21_v_L;
    c21_w_y = c21_j_a * c21_j_b;
    c21_dPd[1] = c21_w_y;
    _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, 16);
    for (c21_i12 = 0; c21_i12 < 2; c21_i12++) {
      c21_k_a[c21_i12] = c21_dPd[c21_i12];
    }

    c21_k_b = c21_gam - c21_gamWP;
    for (c21_i13 = 0; c21_i13 < 2; c21_i13++) {
      c21_k_a[c21_i13] *= c21_k_b;
    }

    for (c21_i14 = 0; c21_i14 < 2; c21_i14++) {
      c21_Pd[c21_i14] = c21_P0[c21_i14] + c21_k_a[c21_i14];
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c21_sfEvent, -16);
  _SFD_SYMBOL_SCOPE_POP();
  for (c21_i15 = 0; c21_i15 < 2; c21_i15++) {
    (*c21_b_Pd)[c21_i15] = c21_Pd[c21_i15];
  }

  for (c21_i16 = 0; c21_i16 < 2; c21_i16++) {
    (*c21_b_dPd)[c21_i16] = c21_dPd[c21_i16];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 20U, chartInstance->c21_sfEvent);
}

static void initSimStructsc21_simulation(SFc21_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c21_machineNumber, uint32_T
  c21_chartNumber)
{
}

static const mxArray *c21_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  int32_T c21_i17;
  real_T c21_b_inData[2];
  int32_T c21_i18;
  real_T c21_u[2];
  const mxArray *c21_y = NULL;
  SFc21_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc21_simulationInstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  for (c21_i17 = 0; c21_i17 < 2; c21_i17++) {
    c21_b_inData[c21_i17] = (*(real_T (*)[2])c21_inData)[c21_i17];
  }

  for (c21_i18 = 0; c21_i18 < 2; c21_i18++) {
    c21_u[c21_i18] = c21_b_inData[c21_i18];
  }

  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", c21_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, FALSE);
  return c21_mxArrayOutData;
}

static void c21_emlrt_marshallIn(SFc21_simulationInstanceStruct *chartInstance,
  const mxArray *c21_dPd, const char_T *c21_identifier, real_T c21_y[2])
{
  emlrtMsgIdentifier c21_thisId;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_dPd), &c21_thisId, c21_y);
  sf_mex_destroy(&c21_dPd);
}

static void c21_b_emlrt_marshallIn(SFc21_simulationInstanceStruct *chartInstance,
  const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId, real_T c21_y[2])
{
  real_T c21_dv2[2];
  int32_T c21_i19;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), c21_dv2, 1, 0, 0U, 1, 0U, 1, 2);
  for (c21_i19 = 0; c21_i19 < 2; c21_i19++) {
    c21_y[c21_i19] = c21_dv2[c21_i19];
  }

  sf_mex_destroy(&c21_u);
}

static void c21_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_dPd;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  real_T c21_y[2];
  int32_T c21_i20;
  SFc21_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc21_simulationInstanceStruct *)chartInstanceVoid;
  c21_dPd = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_dPd), &c21_thisId, c21_y);
  sf_mex_destroy(&c21_dPd);
  for (c21_i20 = 0; c21_i20 < 2; c21_i20++) {
    (*(real_T (*)[2])c21_outData)[c21_i20] = c21_y[c21_i20];
  }

  sf_mex_destroy(&c21_mxArrayInData);
}

static const mxArray *c21_b_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  real_T c21_u;
  const mxArray *c21_y = NULL;
  SFc21_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc21_simulationInstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  c21_u = *(real_T *)c21_inData;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", &c21_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, FALSE);
  return c21_mxArrayOutData;
}

static real_T c21_c_emlrt_marshallIn(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId)
{
  real_T c21_y;
  real_T c21_d0;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), &c21_d0, 1, 0, 0U, 0, 0U, 0);
  c21_y = c21_d0;
  sf_mex_destroy(&c21_u);
  return c21_y;
}

static void c21_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_nargout;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  real_T c21_y;
  SFc21_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc21_simulationInstanceStruct *)chartInstanceVoid;
  c21_nargout = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_y = c21_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_nargout),
    &c21_thisId);
  sf_mex_destroy(&c21_nargout);
  *(real_T *)c21_outData = c21_y;
  sf_mex_destroy(&c21_mxArrayInData);
}

const mxArray *sf_c21_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c21_nameCaptureInfo = NULL;
  c21_nameCaptureInfo = NULL;
  sf_mex_assign(&c21_nameCaptureInfo, sf_mex_createstruct("structure", 2, 11, 1),
                FALSE);
  c21_info_helper(&c21_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c21_nameCaptureInfo);
  return c21_nameCaptureInfo;
}

static void c21_info_helper(const mxArray **c21_info)
{
  const mxArray *c21_rhs0 = NULL;
  const mxArray *c21_lhs0 = NULL;
  const mxArray *c21_rhs1 = NULL;
  const mxArray *c21_lhs1 = NULL;
  const mxArray *c21_rhs2 = NULL;
  const mxArray *c21_lhs2 = NULL;
  const mxArray *c21_rhs3 = NULL;
  const mxArray *c21_lhs3 = NULL;
  const mxArray *c21_rhs4 = NULL;
  const mxArray *c21_lhs4 = NULL;
  const mxArray *c21_rhs5 = NULL;
  const mxArray *c21_lhs5 = NULL;
  const mxArray *c21_rhs6 = NULL;
  const mxArray *c21_lhs6 = NULL;
  const mxArray *c21_rhs7 = NULL;
  const mxArray *c21_lhs7 = NULL;
  const mxArray *c21_rhs8 = NULL;
  const mxArray *c21_lhs8 = NULL;
  const mxArray *c21_rhs9 = NULL;
  const mxArray *c21_lhs9 = NULL;
  const mxArray *c21_rhs10 = NULL;
  const mxArray *c21_lhs10 = NULL;
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("mtimes"), "name", "name", 0);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c21_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 1);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c21_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("mrdivide"), "name", "name",
                  2);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c21_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 3);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("rdivide"), "name", "name", 3);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c21_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c21_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 5);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c21_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_div"), "name", "name", 6);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c21_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 7);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("sin"), "name", "name", 7);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c21_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 8);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c21_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(""), "context", "context", 9);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("cos"), "name", "name", 9);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c21_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 10);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c21_info, c21_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c21_info, c21_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c21_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c21_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c21_info, sf_mex_duplicatearraysafe(&c21_lhs10), "lhs", "lhs",
                  10);
  sf_mex_destroy(&c21_rhs0);
  sf_mex_destroy(&c21_lhs0);
  sf_mex_destroy(&c21_rhs1);
  sf_mex_destroy(&c21_lhs1);
  sf_mex_destroy(&c21_rhs2);
  sf_mex_destroy(&c21_lhs2);
  sf_mex_destroy(&c21_rhs3);
  sf_mex_destroy(&c21_lhs3);
  sf_mex_destroy(&c21_rhs4);
  sf_mex_destroy(&c21_lhs4);
  sf_mex_destroy(&c21_rhs5);
  sf_mex_destroy(&c21_lhs5);
  sf_mex_destroy(&c21_rhs6);
  sf_mex_destroy(&c21_lhs6);
  sf_mex_destroy(&c21_rhs7);
  sf_mex_destroy(&c21_lhs7);
  sf_mex_destroy(&c21_rhs8);
  sf_mex_destroy(&c21_lhs8);
  sf_mex_destroy(&c21_rhs9);
  sf_mex_destroy(&c21_lhs9);
  sf_mex_destroy(&c21_rhs10);
  sf_mex_destroy(&c21_lhs10);
}

static const mxArray *c21_emlrt_marshallOut(char * c21_u)
{
  const mxArray *c21_y = NULL;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", c21_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c21_u)), FALSE);
  return c21_y;
}

static const mxArray *c21_b_emlrt_marshallOut(uint32_T c21_u)
{
  const mxArray *c21_y = NULL;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", &c21_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c21_y;
}

static const mxArray *c21_c_sf_marshallOut(void *chartInstanceVoid, void
  *c21_inData)
{
  const mxArray *c21_mxArrayOutData = NULL;
  int32_T c21_u;
  const mxArray *c21_y = NULL;
  SFc21_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc21_simulationInstanceStruct *)chartInstanceVoid;
  c21_mxArrayOutData = NULL;
  c21_u = *(int32_T *)c21_inData;
  c21_y = NULL;
  sf_mex_assign(&c21_y, sf_mex_create("y", &c21_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c21_mxArrayOutData, c21_y, FALSE);
  return c21_mxArrayOutData;
}

static int32_T c21_d_emlrt_marshallIn(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId)
{
  int32_T c21_y;
  int32_T c21_i21;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), &c21_i21, 1, 6, 0U, 0, 0U, 0);
  c21_y = c21_i21;
  sf_mex_destroy(&c21_u);
  return c21_y;
}

static void c21_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c21_mxArrayInData, const char_T *c21_varName, void *c21_outData)
{
  const mxArray *c21_b_sfEvent;
  const char_T *c21_identifier;
  emlrtMsgIdentifier c21_thisId;
  int32_T c21_y;
  SFc21_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc21_simulationInstanceStruct *)chartInstanceVoid;
  c21_b_sfEvent = sf_mex_dup(c21_mxArrayInData);
  c21_identifier = c21_varName;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_y = c21_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c21_b_sfEvent),
    &c21_thisId);
  sf_mex_destroy(&c21_b_sfEvent);
  *(int32_T *)c21_outData = c21_y;
  sf_mex_destroy(&c21_mxArrayInData);
}

static uint8_T c21_e_emlrt_marshallIn(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_b_is_active_c21_simulation, const char_T
  *c21_identifier)
{
  uint8_T c21_y;
  emlrtMsgIdentifier c21_thisId;
  c21_thisId.fIdentifier = c21_identifier;
  c21_thisId.fParent = NULL;
  c21_y = c21_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c21_b_is_active_c21_simulation), &c21_thisId);
  sf_mex_destroy(&c21_b_is_active_c21_simulation);
  return c21_y;
}

static uint8_T c21_f_emlrt_marshallIn(SFc21_simulationInstanceStruct
  *chartInstance, const mxArray *c21_u, const emlrtMsgIdentifier *c21_parentId)
{
  uint8_T c21_y;
  uint8_T c21_u0;
  sf_mex_import(c21_parentId, sf_mex_dup(c21_u), &c21_u0, 1, 3, 0U, 0, 0U, 0);
  c21_y = c21_u0;
  sf_mex_destroy(&c21_u);
  return c21_y;
}

static void init_dsm_address_info(SFc21_simulationInstanceStruct *chartInstance)
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

void sf_c21_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1094584075U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3420162766U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(887642851U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1149413169U);
}

mxArray *sf_c21_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("1i0keOOofqyCVWTwJwtZXH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
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
      pr[0] = (double)(2);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c21_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c21_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c21_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[7],T\"Pd\",},{M[1],M[5],T\"dPd\",},{M[8],M[0],T\"is_active_c21_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c21_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc21_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc21_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           21,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"gam");
          _SFD_SET_DATA_PROPS(1,2,0,1,"Pd");
          _SFD_SET_DATA_PROPS(2,1,1,0,"gamWP");
          _SFD_SET_DATA_PROPS(3,2,0,1,"dPd");
          _SFD_SET_DATA_PROPS(4,1,1,0,"thetad");
          _SFD_SET_DATA_PROPS(5,1,1,0,"P0");
          _SFD_SET_DATA_PROPS(6,1,1,0,"v_L");
          _SFD_SET_DATA_PROPS(7,1,1,0,"Radius");
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
        _SFD_CV_INIT_EML(0,1,1,1,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,474);
        _SFD_CV_INIT_EML_IF(0,1,0,67,81,356,472);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c21_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c21_sf_marshallOut,(MexInFcnForType)
            c21_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c21_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c21_sf_marshallOut,(MexInFcnForType)
            c21_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c21_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c21_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c21_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c21_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c21_gam;
          real_T *c21_gamWP;
          real_T *c21_thetad;
          real_T *c21_v_L;
          real_T *c21_Radius;
          real_T (*c21_Pd)[2];
          real_T (*c21_dPd)[2];
          real_T (*c21_P0)[2];
          c21_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
          c21_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c21_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
          c21_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c21_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
          c21_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c21_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
          c21_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c21_gam);
          _SFD_SET_DATA_VALUE_PTR(1U, *c21_Pd);
          _SFD_SET_DATA_VALUE_PTR(2U, c21_gamWP);
          _SFD_SET_DATA_VALUE_PTR(3U, *c21_dPd);
          _SFD_SET_DATA_VALUE_PTR(4U, c21_thetad);
          _SFD_SET_DATA_VALUE_PTR(5U, *c21_P0);
          _SFD_SET_DATA_VALUE_PTR(6U, c21_v_L);
          _SFD_SET_DATA_VALUE_PTR(7U, c21_Radius);
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
  return "MVWvJkFnvz9zAFQJSl6x4E";
}

static void sf_opaque_initialize_c21_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc21_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c21_simulation((SFc21_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c21_simulation((SFc21_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c21_simulation(void *chartInstanceVar)
{
  enable_c21_simulation((SFc21_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c21_simulation(void *chartInstanceVar)
{
  disable_c21_simulation((SFc21_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c21_simulation(void *chartInstanceVar)
{
  sf_c21_simulation((SFc21_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c21_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c21_simulation
    ((SFc21_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c21_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c21_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c21_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c21_simulation((SFc21_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c21_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c21_simulation(S);
}

static void sf_opaque_set_sim_state_c21_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c21_simulation(S, st);
}

static void sf_opaque_terminate_c21_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc21_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c21_simulation((SFc21_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc21_simulation((SFc21_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c21_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c21_simulation((SFc21_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c21_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      21);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,21,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,21,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,21);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,21,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,21,2);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=2; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 6; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,21);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(271783144U));
  ssSetChecksum1(S,(1705197971U));
  ssSetChecksum2(S,(2628351788U));
  ssSetChecksum3(S,(2027562728U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c21_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c21_simulation(SimStruct *S)
{
  SFc21_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc21_simulationInstanceStruct *)utMalloc(sizeof
    (SFc21_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc21_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c21_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c21_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c21_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c21_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c21_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c21_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c21_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c21_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c21_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c21_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c21_simulation;
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

void c21_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c21_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c21_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c21_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c21_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
