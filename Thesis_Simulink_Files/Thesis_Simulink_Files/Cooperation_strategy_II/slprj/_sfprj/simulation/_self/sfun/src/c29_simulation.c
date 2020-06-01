/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c29_simulation.h"
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
static const char * c29_debug_family_names[13] = { "R0", "R", "alpha0", "nargin",
  "nargout", "gam", "gamWP", "thetad", "P0", "v_L", "Radius", "Pd", "dPd" };

/* Function Declarations */
static void initialize_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance);
static void enable_c29_simulation(SFc29_simulationInstanceStruct *chartInstance);
static void disable_c29_simulation(SFc29_simulationInstanceStruct *chartInstance);
static void c29_update_debugger_state_c29_simulation
  (SFc29_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c29_simulation
  (SFc29_simulationInstanceStruct *chartInstance);
static void set_sim_state_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_st);
static void finalize_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance);
static void sf_c29_simulation(SFc29_simulationInstanceStruct *chartInstance);
static void c29_chartstep_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc29_simulation(SFc29_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c29_machineNumber, uint32_T
  c29_chartNumber);
static const mxArray *c29_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData);
static void c29_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_dPd, const char_T *c29_identifier, real_T c29_y[2]);
static void c29_b_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId, real_T c29_y[2]);
static void c29_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData);
static const mxArray *c29_b_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData);
static real_T c29_c_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId);
static void c29_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData);
static void c29_info_helper(const mxArray **c29_info);
static const mxArray *c29_emlrt_marshallOut(char * c29_u);
static const mxArray *c29_b_emlrt_marshallOut(uint32_T c29_u);
static const mxArray *c29_c_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData);
static int32_T c29_d_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId);
static void c29_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData);
static uint8_T c29_e_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_b_is_active_c29_simulation, const char_T
  *c29_identifier);
static uint8_T c29_f_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId);
static void init_dsm_address_info(SFc29_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c29_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c29_is_active_c29_simulation = 0U;
}

static void initialize_params_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance)
{
}

static void enable_c29_simulation(SFc29_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c29_simulation(SFc29_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c29_update_debugger_state_c29_simulation
  (SFc29_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c29_simulation
  (SFc29_simulationInstanceStruct *chartInstance)
{
  const mxArray *c29_st;
  const mxArray *c29_y = NULL;
  int32_T c29_i0;
  real_T c29_u[2];
  const mxArray *c29_b_y = NULL;
  int32_T c29_i1;
  real_T c29_b_u[2];
  const mxArray *c29_c_y = NULL;
  uint8_T c29_hoistedGlobal;
  uint8_T c29_c_u;
  const mxArray *c29_d_y = NULL;
  real_T (*c29_dPd)[2];
  real_T (*c29_Pd)[2];
  c29_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c29_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c29_st = NULL;
  c29_st = NULL;
  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_createcellarray(3), FALSE);
  for (c29_i0 = 0; c29_i0 < 2; c29_i0++) {
    c29_u[c29_i0] = (*c29_Pd)[c29_i0];
  }

  c29_b_y = NULL;
  sf_mex_assign(&c29_b_y, sf_mex_create("y", c29_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_setcell(c29_y, 0, c29_b_y);
  for (c29_i1 = 0; c29_i1 < 2; c29_i1++) {
    c29_b_u[c29_i1] = (*c29_dPd)[c29_i1];
  }

  c29_c_y = NULL;
  sf_mex_assign(&c29_c_y, sf_mex_create("y", c29_b_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c29_y, 1, c29_c_y);
  c29_hoistedGlobal = chartInstance->c29_is_active_c29_simulation;
  c29_c_u = c29_hoistedGlobal;
  c29_d_y = NULL;
  sf_mex_assign(&c29_d_y, sf_mex_create("y", &c29_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c29_y, 2, c29_d_y);
  sf_mex_assign(&c29_st, c29_y, FALSE);
  return c29_st;
}

static void set_sim_state_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_st)
{
  const mxArray *c29_u;
  real_T c29_dv0[2];
  int32_T c29_i2;
  real_T c29_dv1[2];
  int32_T c29_i3;
  real_T (*c29_Pd)[2];
  real_T (*c29_dPd)[2];
  c29_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c29_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c29_doneDoubleBufferReInit = TRUE;
  c29_u = sf_mex_dup(c29_st);
  c29_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c29_u, 0)), "Pd",
                       c29_dv0);
  for (c29_i2 = 0; c29_i2 < 2; c29_i2++) {
    (*c29_Pd)[c29_i2] = c29_dv0[c29_i2];
  }

  c29_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c29_u, 1)),
                       "dPd", c29_dv1);
  for (c29_i3 = 0; c29_i3 < 2; c29_i3++) {
    (*c29_dPd)[c29_i3] = c29_dv1[c29_i3];
  }

  chartInstance->c29_is_active_c29_simulation = c29_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c29_u, 2)),
     "is_active_c29_simulation");
  sf_mex_destroy(&c29_u);
  c29_update_debugger_state_c29_simulation(chartInstance);
  sf_mex_destroy(&c29_st);
}

static void finalize_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c29_simulation(SFc29_simulationInstanceStruct *chartInstance)
{
  int32_T c29_i4;
  int32_T c29_i5;
  int32_T c29_i6;
  real_T *c29_gam;
  real_T *c29_gamWP;
  real_T *c29_thetad;
  real_T *c29_v_L;
  real_T *c29_Radius;
  real_T (*c29_P0)[2];
  real_T (*c29_dPd)[2];
  real_T (*c29_Pd)[2];
  c29_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c29_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c29_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c29_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c29_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c29_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c29_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c29_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 28U, chartInstance->c29_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c29_gam, 0U);
  for (c29_i4 = 0; c29_i4 < 2; c29_i4++) {
    _SFD_DATA_RANGE_CHECK((*c29_Pd)[c29_i4], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c29_gamWP, 2U);
  for (c29_i5 = 0; c29_i5 < 2; c29_i5++) {
    _SFD_DATA_RANGE_CHECK((*c29_dPd)[c29_i5], 3U);
  }

  _SFD_DATA_RANGE_CHECK(*c29_thetad, 4U);
  for (c29_i6 = 0; c29_i6 < 2; c29_i6++) {
    _SFD_DATA_RANGE_CHECK((*c29_P0)[c29_i6], 5U);
  }

  _SFD_DATA_RANGE_CHECK(*c29_v_L, 6U);
  _SFD_DATA_RANGE_CHECK(*c29_Radius, 7U);
  chartInstance->c29_sfEvent = CALL_EVENT;
  c29_chartstep_c29_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c29_chartstep_c29_simulation(SFc29_simulationInstanceStruct
  *chartInstance)
{
  real_T c29_hoistedGlobal;
  real_T c29_b_hoistedGlobal;
  real_T c29_c_hoistedGlobal;
  real_T c29_d_hoistedGlobal;
  real_T c29_e_hoistedGlobal;
  real_T c29_gam;
  real_T c29_gamWP;
  real_T c29_thetad;
  int32_T c29_i7;
  real_T c29_P0[2];
  real_T c29_v_L;
  real_T c29_Radius;
  uint32_T c29_debug_family_var_map[13];
  real_T c29_R0[2];
  real_T c29_R;
  real_T c29_alpha0;
  real_T c29_nargin = 6.0;
  real_T c29_nargout = 2.0;
  real_T c29_Pd[2];
  real_T c29_dPd[2];
  int32_T c29_i8;
  int32_T c29_i9;
  real_T c29_a;
  real_T c29_b;
  real_T c29_y;
  real_T c29_A;
  real_T c29_B;
  real_T c29_x;
  real_T c29_b_y;
  real_T c29_b_x;
  real_T c29_c_y;
  real_T c29_d_y;
  real_T c29_c_x;
  real_T c29_d_x;
  real_T c29_b_a;
  real_T c29_b_b;
  real_T c29_e_y;
  real_T c29_c_a;
  real_T c29_c_b;
  real_T c29_f_y;
  real_T c29_b_A;
  real_T c29_b_B;
  real_T c29_e_x;
  real_T c29_g_y;
  real_T c29_f_x;
  real_T c29_h_y;
  real_T c29_i_y;
  real_T c29_g_x;
  real_T c29_h_x;
  real_T c29_d_a;
  real_T c29_d_b;
  real_T c29_j_y;
  real_T c29_e_a;
  real_T c29_e_b;
  real_T c29_k_y;
  real_T c29_c_A;
  real_T c29_c_B;
  real_T c29_i_x;
  real_T c29_l_y;
  real_T c29_j_x;
  real_T c29_m_y;
  real_T c29_n_y;
  real_T c29_k_x;
  real_T c29_l_x;
  real_T c29_f_a;
  real_T c29_f_b;
  real_T c29_o_y;
  real_T c29_g_a;
  real_T c29_g_b;
  real_T c29_p_y;
  real_T c29_d_A;
  real_T c29_d_B;
  real_T c29_m_x;
  real_T c29_q_y;
  real_T c29_n_x;
  real_T c29_r_y;
  real_T c29_s_y;
  real_T c29_o_x;
  real_T c29_p_x;
  real_T c29_h_a;
  real_T c29_h_b;
  real_T c29_t_y;
  real_T c29_u_y[2];
  int32_T c29_i10;
  int32_T c29_i11;
  real_T c29_q_x;
  real_T c29_r_x;
  real_T c29_i_a;
  real_T c29_i_b;
  real_T c29_v_y;
  real_T c29_s_x;
  real_T c29_t_x;
  real_T c29_j_a;
  real_T c29_j_b;
  real_T c29_w_y;
  int32_T c29_i12;
  real_T c29_k_a[2];
  real_T c29_k_b;
  int32_T c29_i13;
  int32_T c29_i14;
  int32_T c29_i15;
  int32_T c29_i16;
  real_T *c29_b_Radius;
  real_T *c29_b_v_L;
  real_T *c29_b_thetad;
  real_T *c29_b_gamWP;
  real_T *c29_b_gam;
  real_T (*c29_b_Pd)[2];
  real_T (*c29_b_dPd)[2];
  real_T (*c29_b_P0)[2];
  c29_b_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c29_b_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c29_b_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c29_b_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c29_b_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c29_b_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c29_b_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c29_b_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 28U, chartInstance->c29_sfEvent);
  c29_hoistedGlobal = *c29_b_gam;
  c29_b_hoistedGlobal = *c29_b_gamWP;
  c29_c_hoistedGlobal = *c29_b_thetad;
  c29_d_hoistedGlobal = *c29_b_v_L;
  c29_e_hoistedGlobal = *c29_b_Radius;
  c29_gam = c29_hoistedGlobal;
  c29_gamWP = c29_b_hoistedGlobal;
  c29_thetad = c29_c_hoistedGlobal;
  for (c29_i7 = 0; c29_i7 < 2; c29_i7++) {
    c29_P0[c29_i7] = (*c29_b_P0)[c29_i7];
  }

  c29_v_L = c29_d_hoistedGlobal;
  c29_Radius = c29_e_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 13U, 13U, c29_debug_family_names,
    c29_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c29_R0, 0U, c29_sf_marshallOut,
    c29_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_R, 1U, c29_b_sf_marshallOut,
    c29_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_alpha0, 2U, c29_b_sf_marshallOut,
    c29_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_nargin, 3U, c29_b_sf_marshallOut,
    c29_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_nargout, 4U, c29_b_sf_marshallOut,
    c29_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c29_gam, 5U, c29_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c29_gamWP, 6U, c29_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c29_thetad, 7U, c29_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c29_P0, 8U, c29_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c29_v_L, 9U, c29_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c29_Radius, 10U, c29_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c29_Pd, 11U, c29_sf_marshallOut,
    c29_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c29_dPd, 12U, c29_sf_marshallOut,
    c29_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 3);
  if (CV_EML_IF(0, 1, 0, c29_Radius != 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 4);
    for (c29_i8 = 0; c29_i8 < 2; c29_i8++) {
      c29_R0[c29_i8] = c29_P0[c29_i8];
    }

    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 5);
    c29_R = c29_Radius;
    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 6);
    c29_alpha0 = c29_thetad;
    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 7);
    for (c29_i9 = 0; c29_i9 < 2; c29_i9++) {
      c29_dPd[c29_i9] = 0.0;
    }

    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 8);
    c29_a = c29_gam - c29_gamWP;
    c29_b = c29_v_L;
    c29_y = c29_a * c29_b;
    c29_A = c29_y;
    c29_B = c29_R;
    c29_x = c29_A;
    c29_b_y = c29_B;
    c29_b_x = c29_x;
    c29_c_y = c29_b_y;
    c29_d_y = c29_b_x / c29_c_y;
    c29_c_x = c29_alpha0 + c29_d_y;
    c29_d_x = c29_c_x;
    c29_d_x = muDoubleScalarSin(c29_d_x);
    c29_b_a = -c29_d_x;
    c29_b_b = c29_v_L;
    c29_e_y = c29_b_a * c29_b_b;
    c29_dPd[0] = c29_e_y;
    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 9);
    c29_c_a = c29_gam - c29_gamWP;
    c29_c_b = c29_v_L;
    c29_f_y = c29_c_a * c29_c_b;
    c29_b_A = c29_f_y;
    c29_b_B = c29_R;
    c29_e_x = c29_b_A;
    c29_g_y = c29_b_B;
    c29_f_x = c29_e_x;
    c29_h_y = c29_g_y;
    c29_i_y = c29_f_x / c29_h_y;
    c29_g_x = c29_alpha0 + c29_i_y;
    c29_h_x = c29_g_x;
    c29_h_x = muDoubleScalarCos(c29_h_x);
    c29_d_a = c29_h_x;
    c29_d_b = c29_v_L;
    c29_j_y = c29_d_a * c29_d_b;
    c29_dPd[1] = c29_j_y;
    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 10);
    c29_e_a = c29_gam - c29_gamWP;
    c29_e_b = c29_v_L;
    c29_k_y = c29_e_a * c29_e_b;
    c29_c_A = c29_k_y;
    c29_c_B = c29_R;
    c29_i_x = c29_c_A;
    c29_l_y = c29_c_B;
    c29_j_x = c29_i_x;
    c29_m_y = c29_l_y;
    c29_n_y = c29_j_x / c29_m_y;
    c29_k_x = c29_alpha0 + c29_n_y;
    c29_l_x = c29_k_x;
    c29_l_x = muDoubleScalarCos(c29_l_x);
    c29_f_a = c29_l_x;
    c29_f_b = c29_R;
    c29_o_y = c29_f_a * c29_f_b;
    c29_g_a = c29_gam - c29_gamWP;
    c29_g_b = c29_v_L;
    c29_p_y = c29_g_a * c29_g_b;
    c29_d_A = c29_p_y;
    c29_d_B = c29_R;
    c29_m_x = c29_d_A;
    c29_q_y = c29_d_B;
    c29_n_x = c29_m_x;
    c29_r_y = c29_q_y;
    c29_s_y = c29_n_x / c29_r_y;
    c29_o_x = c29_alpha0 + c29_s_y;
    c29_p_x = c29_o_x;
    c29_p_x = muDoubleScalarSin(c29_p_x);
    c29_h_a = c29_p_x;
    c29_h_b = c29_R;
    c29_t_y = c29_h_a * c29_h_b;
    c29_u_y[0] = c29_o_y;
    c29_u_y[1] = c29_t_y;
    for (c29_i10 = 0; c29_i10 < 2; c29_i10++) {
      c29_Pd[c29_i10] = c29_R0[c29_i10] + c29_u_y[c29_i10];
    }
  } else {
    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 13);
    for (c29_i11 = 0; c29_i11 < 2; c29_i11++) {
      c29_dPd[c29_i11] = 0.0;
    }

    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 14);
    c29_q_x = c29_thetad;
    c29_r_x = c29_q_x;
    c29_r_x = muDoubleScalarCos(c29_r_x);
    c29_i_a = c29_r_x;
    c29_i_b = c29_v_L;
    c29_v_y = c29_i_a * c29_i_b;
    c29_dPd[0] = c29_v_y;
    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 15);
    c29_s_x = c29_thetad;
    c29_t_x = c29_s_x;
    c29_t_x = muDoubleScalarSin(c29_t_x);
    c29_j_a = c29_t_x;
    c29_j_b = c29_v_L;
    c29_w_y = c29_j_a * c29_j_b;
    c29_dPd[1] = c29_w_y;
    _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 16);
    for (c29_i12 = 0; c29_i12 < 2; c29_i12++) {
      c29_k_a[c29_i12] = c29_dPd[c29_i12];
    }

    c29_k_b = c29_gam - c29_gamWP;
    for (c29_i13 = 0; c29_i13 < 2; c29_i13++) {
      c29_k_a[c29_i13] *= c29_k_b;
    }

    for (c29_i14 = 0; c29_i14 < 2; c29_i14++) {
      c29_Pd[c29_i14] = c29_P0[c29_i14] + c29_k_a[c29_i14];
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, -16);
  _SFD_SYMBOL_SCOPE_POP();
  for (c29_i15 = 0; c29_i15 < 2; c29_i15++) {
    (*c29_b_Pd)[c29_i15] = c29_Pd[c29_i15];
  }

  for (c29_i16 = 0; c29_i16 < 2; c29_i16++) {
    (*c29_b_dPd)[c29_i16] = c29_dPd[c29_i16];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 28U, chartInstance->c29_sfEvent);
}

static void initSimStructsc29_simulation(SFc29_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c29_machineNumber, uint32_T
  c29_chartNumber)
{
}

static const mxArray *c29_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData)
{
  const mxArray *c29_mxArrayOutData = NULL;
  int32_T c29_i17;
  real_T c29_b_inData[2];
  int32_T c29_i18;
  real_T c29_u[2];
  const mxArray *c29_y = NULL;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_mxArrayOutData = NULL;
  for (c29_i17 = 0; c29_i17 < 2; c29_i17++) {
    c29_b_inData[c29_i17] = (*(real_T (*)[2])c29_inData)[c29_i17];
  }

  for (c29_i18 = 0; c29_i18 < 2; c29_i18++) {
    c29_u[c29_i18] = c29_b_inData[c29_i18];
  }

  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", c29_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c29_mxArrayOutData, c29_y, FALSE);
  return c29_mxArrayOutData;
}

static void c29_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_dPd, const char_T *c29_identifier, real_T c29_y[2])
{
  emlrtMsgIdentifier c29_thisId;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_dPd), &c29_thisId, c29_y);
  sf_mex_destroy(&c29_dPd);
}

static void c29_b_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId, real_T c29_y[2])
{
  real_T c29_dv2[2];
  int32_T c29_i19;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), c29_dv2, 1, 0, 0U, 1, 0U, 1, 2);
  for (c29_i19 = 0; c29_i19 < 2; c29_i19++) {
    c29_y[c29_i19] = c29_dv2[c29_i19];
  }

  sf_mex_destroy(&c29_u);
}

static void c29_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData)
{
  const mxArray *c29_dPd;
  const char_T *c29_identifier;
  emlrtMsgIdentifier c29_thisId;
  real_T c29_y[2];
  int32_T c29_i20;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_dPd = sf_mex_dup(c29_mxArrayInData);
  c29_identifier = c29_varName;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_dPd), &c29_thisId, c29_y);
  sf_mex_destroy(&c29_dPd);
  for (c29_i20 = 0; c29_i20 < 2; c29_i20++) {
    (*(real_T (*)[2])c29_outData)[c29_i20] = c29_y[c29_i20];
  }

  sf_mex_destroy(&c29_mxArrayInData);
}

static const mxArray *c29_b_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData)
{
  const mxArray *c29_mxArrayOutData = NULL;
  real_T c29_u;
  const mxArray *c29_y = NULL;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_mxArrayOutData = NULL;
  c29_u = *(real_T *)c29_inData;
  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", &c29_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c29_mxArrayOutData, c29_y, FALSE);
  return c29_mxArrayOutData;
}

static real_T c29_c_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId)
{
  real_T c29_y;
  real_T c29_d0;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), &c29_d0, 1, 0, 0U, 0, 0U, 0);
  c29_y = c29_d0;
  sf_mex_destroy(&c29_u);
  return c29_y;
}

static void c29_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData)
{
  const mxArray *c29_nargout;
  const char_T *c29_identifier;
  emlrtMsgIdentifier c29_thisId;
  real_T c29_y;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_nargout = sf_mex_dup(c29_mxArrayInData);
  c29_identifier = c29_varName;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_y = c29_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_nargout),
    &c29_thisId);
  sf_mex_destroy(&c29_nargout);
  *(real_T *)c29_outData = c29_y;
  sf_mex_destroy(&c29_mxArrayInData);
}

const mxArray *sf_c29_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c29_nameCaptureInfo = NULL;
  c29_nameCaptureInfo = NULL;
  sf_mex_assign(&c29_nameCaptureInfo, sf_mex_createstruct("structure", 2, 11, 1),
                FALSE);
  c29_info_helper(&c29_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c29_nameCaptureInfo);
  return c29_nameCaptureInfo;
}

static void c29_info_helper(const mxArray **c29_info)
{
  const mxArray *c29_rhs0 = NULL;
  const mxArray *c29_lhs0 = NULL;
  const mxArray *c29_rhs1 = NULL;
  const mxArray *c29_lhs1 = NULL;
  const mxArray *c29_rhs2 = NULL;
  const mxArray *c29_lhs2 = NULL;
  const mxArray *c29_rhs3 = NULL;
  const mxArray *c29_lhs3 = NULL;
  const mxArray *c29_rhs4 = NULL;
  const mxArray *c29_lhs4 = NULL;
  const mxArray *c29_rhs5 = NULL;
  const mxArray *c29_lhs5 = NULL;
  const mxArray *c29_rhs6 = NULL;
  const mxArray *c29_lhs6 = NULL;
  const mxArray *c29_rhs7 = NULL;
  const mxArray *c29_lhs7 = NULL;
  const mxArray *c29_rhs8 = NULL;
  const mxArray *c29_lhs8 = NULL;
  const mxArray *c29_rhs9 = NULL;
  const mxArray *c29_lhs9 = NULL;
  const mxArray *c29_rhs10 = NULL;
  const mxArray *c29_lhs10 = NULL;
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("mtimes"), "name", "name", 0);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c29_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 1);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c29_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("mrdivide"), "name", "name",
                  2);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c29_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 3);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("rdivide"), "name", "name", 3);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c29_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c29_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 5);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c29_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_div"), "name", "name", 6);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c29_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "context", "context", 7);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("sin"), "name", "name", 7);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c29_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 8);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c29_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "context", "context", 9);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("cos"), "name", "name", 9);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c29_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 10);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c29_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs10), "lhs", "lhs",
                  10);
  sf_mex_destroy(&c29_rhs0);
  sf_mex_destroy(&c29_lhs0);
  sf_mex_destroy(&c29_rhs1);
  sf_mex_destroy(&c29_lhs1);
  sf_mex_destroy(&c29_rhs2);
  sf_mex_destroy(&c29_lhs2);
  sf_mex_destroy(&c29_rhs3);
  sf_mex_destroy(&c29_lhs3);
  sf_mex_destroy(&c29_rhs4);
  sf_mex_destroy(&c29_lhs4);
  sf_mex_destroy(&c29_rhs5);
  sf_mex_destroy(&c29_lhs5);
  sf_mex_destroy(&c29_rhs6);
  sf_mex_destroy(&c29_lhs6);
  sf_mex_destroy(&c29_rhs7);
  sf_mex_destroy(&c29_lhs7);
  sf_mex_destroy(&c29_rhs8);
  sf_mex_destroy(&c29_lhs8);
  sf_mex_destroy(&c29_rhs9);
  sf_mex_destroy(&c29_lhs9);
  sf_mex_destroy(&c29_rhs10);
  sf_mex_destroy(&c29_lhs10);
}

static const mxArray *c29_emlrt_marshallOut(char * c29_u)
{
  const mxArray *c29_y = NULL;
  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", c29_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c29_u)), FALSE);
  return c29_y;
}

static const mxArray *c29_b_emlrt_marshallOut(uint32_T c29_u)
{
  const mxArray *c29_y = NULL;
  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", &c29_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c29_y;
}

static const mxArray *c29_c_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData)
{
  const mxArray *c29_mxArrayOutData = NULL;
  int32_T c29_u;
  const mxArray *c29_y = NULL;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_mxArrayOutData = NULL;
  c29_u = *(int32_T *)c29_inData;
  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", &c29_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c29_mxArrayOutData, c29_y, FALSE);
  return c29_mxArrayOutData;
}

static int32_T c29_d_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId)
{
  int32_T c29_y;
  int32_T c29_i21;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), &c29_i21, 1, 6, 0U, 0, 0U, 0);
  c29_y = c29_i21;
  sf_mex_destroy(&c29_u);
  return c29_y;
}

static void c29_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData)
{
  const mxArray *c29_b_sfEvent;
  const char_T *c29_identifier;
  emlrtMsgIdentifier c29_thisId;
  int32_T c29_y;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_b_sfEvent = sf_mex_dup(c29_mxArrayInData);
  c29_identifier = c29_varName;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_y = c29_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_b_sfEvent),
    &c29_thisId);
  sf_mex_destroy(&c29_b_sfEvent);
  *(int32_T *)c29_outData = c29_y;
  sf_mex_destroy(&c29_mxArrayInData);
}

static uint8_T c29_e_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_b_is_active_c29_simulation, const char_T
  *c29_identifier)
{
  uint8_T c29_y;
  emlrtMsgIdentifier c29_thisId;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_y = c29_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c29_b_is_active_c29_simulation), &c29_thisId);
  sf_mex_destroy(&c29_b_is_active_c29_simulation);
  return c29_y;
}

static uint8_T c29_f_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId)
{
  uint8_T c29_y;
  uint8_T c29_u0;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), &c29_u0, 1, 3, 0U, 0, 0U, 0);
  c29_y = c29_u0;
  sf_mex_destroy(&c29_u);
  return c29_y;
}

static void init_dsm_address_info(SFc29_simulationInstanceStruct *chartInstance)
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

void sf_c29_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1094584075U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3420162766U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(887642851U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1149413169U);
}

mxArray *sf_c29_simulation_get_autoinheritance_info(void)
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

mxArray *sf_c29_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c29_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c29_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[7],T\"Pd\",},{M[1],M[5],T\"dPd\",},{M[8],M[0],T\"is_active_c29_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c29_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc29_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc29_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           29,
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
          (MexFcnForType)c29_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_sf_marshallOut,(MexInFcnForType)
            c29_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c29_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_sf_marshallOut,(MexInFcnForType)
            c29_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c29_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c29_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c29_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c29_gam;
          real_T *c29_gamWP;
          real_T *c29_thetad;
          real_T *c29_v_L;
          real_T *c29_Radius;
          real_T (*c29_Pd)[2];
          real_T (*c29_dPd)[2];
          real_T (*c29_P0)[2];
          c29_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
          c29_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c29_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
          c29_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c29_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
          c29_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c29_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
          c29_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c29_gam);
          _SFD_SET_DATA_VALUE_PTR(1U, *c29_Pd);
          _SFD_SET_DATA_VALUE_PTR(2U, c29_gamWP);
          _SFD_SET_DATA_VALUE_PTR(3U, *c29_dPd);
          _SFD_SET_DATA_VALUE_PTR(4U, c29_thetad);
          _SFD_SET_DATA_VALUE_PTR(5U, *c29_P0);
          _SFD_SET_DATA_VALUE_PTR(6U, c29_v_L);
          _SFD_SET_DATA_VALUE_PTR(7U, c29_Radius);
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

static void sf_opaque_initialize_c29_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc29_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c29_simulation((SFc29_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c29_simulation((SFc29_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c29_simulation(void *chartInstanceVar)
{
  enable_c29_simulation((SFc29_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c29_simulation(void *chartInstanceVar)
{
  disable_c29_simulation((SFc29_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c29_simulation(void *chartInstanceVar)
{
  sf_c29_simulation((SFc29_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c29_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c29_simulation
    ((SFc29_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c29_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c29_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c29_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c29_simulation((SFc29_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c29_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c29_simulation(S);
}

static void sf_opaque_set_sim_state_c29_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c29_simulation(S, st);
}

static void sf_opaque_terminate_c29_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc29_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c29_simulation((SFc29_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc29_simulation((SFc29_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c29_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c29_simulation((SFc29_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c29_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      29);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,29,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,29,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,29);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,29,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,29,2);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,29);
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

static void mdlRTW_c29_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c29_simulation(SimStruct *S)
{
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)utMalloc(sizeof
    (SFc29_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc29_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c29_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c29_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c29_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c29_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c29_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c29_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c29_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c29_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c29_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c29_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c29_simulation;
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

void c29_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c29_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c29_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c29_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c29_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
