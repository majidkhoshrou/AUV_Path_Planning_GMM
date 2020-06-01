/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c11_simulation.h"
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
static const char * c11_debug_family_names[13] = { "R0", "R", "alpha0", "nargin",
  "nargout", "gam", "gamWP", "thetad", "P0", "v_L", "Radius", "Pd", "dPd" };

/* Function Declarations */
static void initialize_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance);
static void enable_c11_simulation(SFc11_simulationInstanceStruct *chartInstance);
static void disable_c11_simulation(SFc11_simulationInstanceStruct *chartInstance);
static void c11_update_debugger_state_c11_simulation
  (SFc11_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c11_simulation
  (SFc11_simulationInstanceStruct *chartInstance);
static void set_sim_state_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_st);
static void finalize_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance);
static void sf_c11_simulation(SFc11_simulationInstanceStruct *chartInstance);
static void c11_chartstep_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc11_simulation(SFc11_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c11_machineNumber, uint32_T
  c11_chartNumber);
static const mxArray *c11_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static void c11_emlrt_marshallIn(SFc11_simulationInstanceStruct *chartInstance,
  const mxArray *c11_dPd, const char_T *c11_identifier, real_T c11_y[2]);
static void c11_b_emlrt_marshallIn(SFc11_simulationInstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId, real_T c11_y[2]);
static void c11_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData);
static const mxArray *c11_b_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static real_T c11_c_emlrt_marshallIn(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId);
static void c11_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData);
static void c11_info_helper(const mxArray **c11_info);
static const mxArray *c11_emlrt_marshallOut(char * c11_u);
static const mxArray *c11_b_emlrt_marshallOut(uint32_T c11_u);
static const mxArray *c11_c_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData);
static int32_T c11_d_emlrt_marshallIn(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId);
static void c11_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData);
static uint8_T c11_e_emlrt_marshallIn(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_b_is_active_c11_simulation, const char_T
  *c11_identifier);
static uint8_T c11_f_emlrt_marshallIn(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId);
static void init_dsm_address_info(SFc11_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c11_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c11_is_active_c11_simulation = 0U;
}

static void initialize_params_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance)
{
}

static void enable_c11_simulation(SFc11_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c11_simulation(SFc11_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c11_update_debugger_state_c11_simulation
  (SFc11_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c11_simulation
  (SFc11_simulationInstanceStruct *chartInstance)
{
  const mxArray *c11_st;
  const mxArray *c11_y = NULL;
  int32_T c11_i0;
  real_T c11_u[2];
  const mxArray *c11_b_y = NULL;
  int32_T c11_i1;
  real_T c11_b_u[2];
  const mxArray *c11_c_y = NULL;
  uint8_T c11_hoistedGlobal;
  uint8_T c11_c_u;
  const mxArray *c11_d_y = NULL;
  real_T (*c11_dPd)[2];
  real_T (*c11_Pd)[2];
  c11_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c11_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c11_st = NULL;
  c11_st = NULL;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_createcellarray(3), FALSE);
  for (c11_i0 = 0; c11_i0 < 2; c11_i0++) {
    c11_u[c11_i0] = (*c11_Pd)[c11_i0];
  }

  c11_b_y = NULL;
  sf_mex_assign(&c11_b_y, sf_mex_create("y", c11_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_setcell(c11_y, 0, c11_b_y);
  for (c11_i1 = 0; c11_i1 < 2; c11_i1++) {
    c11_b_u[c11_i1] = (*c11_dPd)[c11_i1];
  }

  c11_c_y = NULL;
  sf_mex_assign(&c11_c_y, sf_mex_create("y", c11_b_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c11_y, 1, c11_c_y);
  c11_hoistedGlobal = chartInstance->c11_is_active_c11_simulation;
  c11_c_u = c11_hoistedGlobal;
  c11_d_y = NULL;
  sf_mex_assign(&c11_d_y, sf_mex_create("y", &c11_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c11_y, 2, c11_d_y);
  sf_mex_assign(&c11_st, c11_y, FALSE);
  return c11_st;
}

static void set_sim_state_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_st)
{
  const mxArray *c11_u;
  real_T c11_dv0[2];
  int32_T c11_i2;
  real_T c11_dv1[2];
  int32_T c11_i3;
  real_T (*c11_Pd)[2];
  real_T (*c11_dPd)[2];
  c11_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c11_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c11_doneDoubleBufferReInit = TRUE;
  c11_u = sf_mex_dup(c11_st);
  c11_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c11_u, 0)), "Pd",
                       c11_dv0);
  for (c11_i2 = 0; c11_i2 < 2; c11_i2++) {
    (*c11_Pd)[c11_i2] = c11_dv0[c11_i2];
  }

  c11_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c11_u, 1)),
                       "dPd", c11_dv1);
  for (c11_i3 = 0; c11_i3 < 2; c11_i3++) {
    (*c11_dPd)[c11_i3] = c11_dv1[c11_i3];
  }

  chartInstance->c11_is_active_c11_simulation = c11_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c11_u, 2)),
     "is_active_c11_simulation");
  sf_mex_destroy(&c11_u);
  c11_update_debugger_state_c11_simulation(chartInstance);
  sf_mex_destroy(&c11_st);
}

static void finalize_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c11_simulation(SFc11_simulationInstanceStruct *chartInstance)
{
  int32_T c11_i4;
  int32_T c11_i5;
  int32_T c11_i6;
  real_T *c11_gam;
  real_T *c11_gamWP;
  real_T *c11_thetad;
  real_T *c11_v_L;
  real_T *c11_Radius;
  real_T (*c11_P0)[2];
  real_T (*c11_dPd)[2];
  real_T (*c11_Pd)[2];
  c11_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c11_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c11_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c11_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c11_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c11_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c11_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c11_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 10U, chartInstance->c11_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c11_gam, 0U);
  for (c11_i4 = 0; c11_i4 < 2; c11_i4++) {
    _SFD_DATA_RANGE_CHECK((*c11_Pd)[c11_i4], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c11_gamWP, 2U);
  for (c11_i5 = 0; c11_i5 < 2; c11_i5++) {
    _SFD_DATA_RANGE_CHECK((*c11_dPd)[c11_i5], 3U);
  }

  _SFD_DATA_RANGE_CHECK(*c11_thetad, 4U);
  for (c11_i6 = 0; c11_i6 < 2; c11_i6++) {
    _SFD_DATA_RANGE_CHECK((*c11_P0)[c11_i6], 5U);
  }

  _SFD_DATA_RANGE_CHECK(*c11_v_L, 6U);
  _SFD_DATA_RANGE_CHECK(*c11_Radius, 7U);
  chartInstance->c11_sfEvent = CALL_EVENT;
  c11_chartstep_c11_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c11_chartstep_c11_simulation(SFc11_simulationInstanceStruct
  *chartInstance)
{
  real_T c11_hoistedGlobal;
  real_T c11_b_hoistedGlobal;
  real_T c11_c_hoistedGlobal;
  real_T c11_d_hoistedGlobal;
  real_T c11_e_hoistedGlobal;
  real_T c11_gam;
  real_T c11_gamWP;
  real_T c11_thetad;
  int32_T c11_i7;
  real_T c11_P0[2];
  real_T c11_v_L;
  real_T c11_Radius;
  uint32_T c11_debug_family_var_map[13];
  real_T c11_R0[2];
  real_T c11_R;
  real_T c11_alpha0;
  real_T c11_nargin = 6.0;
  real_T c11_nargout = 2.0;
  real_T c11_Pd[2];
  real_T c11_dPd[2];
  int32_T c11_i8;
  int32_T c11_i9;
  real_T c11_a;
  real_T c11_b;
  real_T c11_y;
  real_T c11_A;
  real_T c11_B;
  real_T c11_x;
  real_T c11_b_y;
  real_T c11_b_x;
  real_T c11_c_y;
  real_T c11_d_y;
  real_T c11_c_x;
  real_T c11_d_x;
  real_T c11_b_a;
  real_T c11_b_b;
  real_T c11_e_y;
  real_T c11_c_a;
  real_T c11_c_b;
  real_T c11_f_y;
  real_T c11_b_A;
  real_T c11_b_B;
  real_T c11_e_x;
  real_T c11_g_y;
  real_T c11_f_x;
  real_T c11_h_y;
  real_T c11_i_y;
  real_T c11_g_x;
  real_T c11_h_x;
  real_T c11_d_a;
  real_T c11_d_b;
  real_T c11_j_y;
  real_T c11_e_a;
  real_T c11_e_b;
  real_T c11_k_y;
  real_T c11_c_A;
  real_T c11_c_B;
  real_T c11_i_x;
  real_T c11_l_y;
  real_T c11_j_x;
  real_T c11_m_y;
  real_T c11_n_y;
  real_T c11_k_x;
  real_T c11_l_x;
  real_T c11_f_a;
  real_T c11_f_b;
  real_T c11_o_y;
  real_T c11_g_a;
  real_T c11_g_b;
  real_T c11_p_y;
  real_T c11_d_A;
  real_T c11_d_B;
  real_T c11_m_x;
  real_T c11_q_y;
  real_T c11_n_x;
  real_T c11_r_y;
  real_T c11_s_y;
  real_T c11_o_x;
  real_T c11_p_x;
  real_T c11_h_a;
  real_T c11_h_b;
  real_T c11_t_y;
  real_T c11_u_y[2];
  int32_T c11_i10;
  int32_T c11_i11;
  real_T c11_q_x;
  real_T c11_r_x;
  real_T c11_i_a;
  real_T c11_i_b;
  real_T c11_v_y;
  real_T c11_s_x;
  real_T c11_t_x;
  real_T c11_j_a;
  real_T c11_j_b;
  real_T c11_w_y;
  int32_T c11_i12;
  real_T c11_k_a[2];
  real_T c11_k_b;
  int32_T c11_i13;
  int32_T c11_i14;
  int32_T c11_i15;
  int32_T c11_i16;
  real_T *c11_b_Radius;
  real_T *c11_b_v_L;
  real_T *c11_b_thetad;
  real_T *c11_b_gamWP;
  real_T *c11_b_gam;
  real_T (*c11_b_Pd)[2];
  real_T (*c11_b_dPd)[2];
  real_T (*c11_b_P0)[2];
  c11_b_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c11_b_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c11_b_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c11_b_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c11_b_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c11_b_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c11_b_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c11_b_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 10U, chartInstance->c11_sfEvent);
  c11_hoistedGlobal = *c11_b_gam;
  c11_b_hoistedGlobal = *c11_b_gamWP;
  c11_c_hoistedGlobal = *c11_b_thetad;
  c11_d_hoistedGlobal = *c11_b_v_L;
  c11_e_hoistedGlobal = *c11_b_Radius;
  c11_gam = c11_hoistedGlobal;
  c11_gamWP = c11_b_hoistedGlobal;
  c11_thetad = c11_c_hoistedGlobal;
  for (c11_i7 = 0; c11_i7 < 2; c11_i7++) {
    c11_P0[c11_i7] = (*c11_b_P0)[c11_i7];
  }

  c11_v_L = c11_d_hoistedGlobal;
  c11_Radius = c11_e_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 13U, 13U, c11_debug_family_names,
    c11_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_R0, 0U, c11_sf_marshallOut,
    c11_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_R, 1U, c11_b_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_alpha0, 2U, c11_b_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_nargin, 3U, c11_b_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c11_nargout, 4U, c11_b_sf_marshallOut,
    c11_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c11_gam, 5U, c11_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c11_gamWP, 6U, c11_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c11_thetad, 7U, c11_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c11_P0, 8U, c11_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c11_v_L, 9U, c11_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c11_Radius, 10U, c11_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_Pd, 11U, c11_sf_marshallOut,
    c11_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c11_dPd, 12U, c11_sf_marshallOut,
    c11_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 3);
  if (CV_EML_IF(0, 1, 0, c11_Radius != 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 4);
    for (c11_i8 = 0; c11_i8 < 2; c11_i8++) {
      c11_R0[c11_i8] = c11_P0[c11_i8];
    }

    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 5);
    c11_R = c11_Radius;
    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 6);
    c11_alpha0 = c11_thetad;
    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 7);
    for (c11_i9 = 0; c11_i9 < 2; c11_i9++) {
      c11_dPd[c11_i9] = 0.0;
    }

    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 8);
    c11_a = c11_gam - c11_gamWP;
    c11_b = c11_v_L;
    c11_y = c11_a * c11_b;
    c11_A = c11_y;
    c11_B = c11_R;
    c11_x = c11_A;
    c11_b_y = c11_B;
    c11_b_x = c11_x;
    c11_c_y = c11_b_y;
    c11_d_y = c11_b_x / c11_c_y;
    c11_c_x = c11_alpha0 + c11_d_y;
    c11_d_x = c11_c_x;
    c11_d_x = muDoubleScalarSin(c11_d_x);
    c11_b_a = -c11_d_x;
    c11_b_b = c11_v_L;
    c11_e_y = c11_b_a * c11_b_b;
    c11_dPd[0] = c11_e_y;
    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 9);
    c11_c_a = c11_gam - c11_gamWP;
    c11_c_b = c11_v_L;
    c11_f_y = c11_c_a * c11_c_b;
    c11_b_A = c11_f_y;
    c11_b_B = c11_R;
    c11_e_x = c11_b_A;
    c11_g_y = c11_b_B;
    c11_f_x = c11_e_x;
    c11_h_y = c11_g_y;
    c11_i_y = c11_f_x / c11_h_y;
    c11_g_x = c11_alpha0 + c11_i_y;
    c11_h_x = c11_g_x;
    c11_h_x = muDoubleScalarCos(c11_h_x);
    c11_d_a = c11_h_x;
    c11_d_b = c11_v_L;
    c11_j_y = c11_d_a * c11_d_b;
    c11_dPd[1] = c11_j_y;
    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 10);
    c11_e_a = c11_gam - c11_gamWP;
    c11_e_b = c11_v_L;
    c11_k_y = c11_e_a * c11_e_b;
    c11_c_A = c11_k_y;
    c11_c_B = c11_R;
    c11_i_x = c11_c_A;
    c11_l_y = c11_c_B;
    c11_j_x = c11_i_x;
    c11_m_y = c11_l_y;
    c11_n_y = c11_j_x / c11_m_y;
    c11_k_x = c11_alpha0 + c11_n_y;
    c11_l_x = c11_k_x;
    c11_l_x = muDoubleScalarCos(c11_l_x);
    c11_f_a = c11_l_x;
    c11_f_b = c11_R;
    c11_o_y = c11_f_a * c11_f_b;
    c11_g_a = c11_gam - c11_gamWP;
    c11_g_b = c11_v_L;
    c11_p_y = c11_g_a * c11_g_b;
    c11_d_A = c11_p_y;
    c11_d_B = c11_R;
    c11_m_x = c11_d_A;
    c11_q_y = c11_d_B;
    c11_n_x = c11_m_x;
    c11_r_y = c11_q_y;
    c11_s_y = c11_n_x / c11_r_y;
    c11_o_x = c11_alpha0 + c11_s_y;
    c11_p_x = c11_o_x;
    c11_p_x = muDoubleScalarSin(c11_p_x);
    c11_h_a = c11_p_x;
    c11_h_b = c11_R;
    c11_t_y = c11_h_a * c11_h_b;
    c11_u_y[0] = c11_o_y;
    c11_u_y[1] = c11_t_y;
    for (c11_i10 = 0; c11_i10 < 2; c11_i10++) {
      c11_Pd[c11_i10] = c11_R0[c11_i10] + c11_u_y[c11_i10];
    }
  } else {
    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 13);
    for (c11_i11 = 0; c11_i11 < 2; c11_i11++) {
      c11_dPd[c11_i11] = 0.0;
    }

    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 14);
    c11_q_x = c11_thetad;
    c11_r_x = c11_q_x;
    c11_r_x = muDoubleScalarCos(c11_r_x);
    c11_i_a = c11_r_x;
    c11_i_b = c11_v_L;
    c11_v_y = c11_i_a * c11_i_b;
    c11_dPd[0] = c11_v_y;
    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 15);
    c11_s_x = c11_thetad;
    c11_t_x = c11_s_x;
    c11_t_x = muDoubleScalarSin(c11_t_x);
    c11_j_a = c11_t_x;
    c11_j_b = c11_v_L;
    c11_w_y = c11_j_a * c11_j_b;
    c11_dPd[1] = c11_w_y;
    _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, 16);
    for (c11_i12 = 0; c11_i12 < 2; c11_i12++) {
      c11_k_a[c11_i12] = c11_dPd[c11_i12];
    }

    c11_k_b = c11_gam - c11_gamWP;
    for (c11_i13 = 0; c11_i13 < 2; c11_i13++) {
      c11_k_a[c11_i13] *= c11_k_b;
    }

    for (c11_i14 = 0; c11_i14 < 2; c11_i14++) {
      c11_Pd[c11_i14] = c11_P0[c11_i14] + c11_k_a[c11_i14];
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c11_sfEvent, -16);
  _SFD_SYMBOL_SCOPE_POP();
  for (c11_i15 = 0; c11_i15 < 2; c11_i15++) {
    (*c11_b_Pd)[c11_i15] = c11_Pd[c11_i15];
  }

  for (c11_i16 = 0; c11_i16 < 2; c11_i16++) {
    (*c11_b_dPd)[c11_i16] = c11_dPd[c11_i16];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 10U, chartInstance->c11_sfEvent);
}

static void initSimStructsc11_simulation(SFc11_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c11_machineNumber, uint32_T
  c11_chartNumber)
{
}

static const mxArray *c11_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  int32_T c11_i17;
  real_T c11_b_inData[2];
  int32_T c11_i18;
  real_T c11_u[2];
  const mxArray *c11_y = NULL;
  SFc11_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc11_simulationInstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  for (c11_i17 = 0; c11_i17 < 2; c11_i17++) {
    c11_b_inData[c11_i17] = (*(real_T (*)[2])c11_inData)[c11_i17];
  }

  for (c11_i18 = 0; c11_i18 < 2; c11_i18++) {
    c11_u[c11_i18] = c11_b_inData[c11_i18];
  }

  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", c11_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, FALSE);
  return c11_mxArrayOutData;
}

static void c11_emlrt_marshallIn(SFc11_simulationInstanceStruct *chartInstance,
  const mxArray *c11_dPd, const char_T *c11_identifier, real_T c11_y[2])
{
  emlrtMsgIdentifier c11_thisId;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_dPd), &c11_thisId, c11_y);
  sf_mex_destroy(&c11_dPd);
}

static void c11_b_emlrt_marshallIn(SFc11_simulationInstanceStruct *chartInstance,
  const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId, real_T c11_y[2])
{
  real_T c11_dv2[2];
  int32_T c11_i19;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), c11_dv2, 1, 0, 0U, 1, 0U, 1, 2);
  for (c11_i19 = 0; c11_i19 < 2; c11_i19++) {
    c11_y[c11_i19] = c11_dv2[c11_i19];
  }

  sf_mex_destroy(&c11_u);
}

static void c11_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData)
{
  const mxArray *c11_dPd;
  const char_T *c11_identifier;
  emlrtMsgIdentifier c11_thisId;
  real_T c11_y[2];
  int32_T c11_i20;
  SFc11_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc11_simulationInstanceStruct *)chartInstanceVoid;
  c11_dPd = sf_mex_dup(c11_mxArrayInData);
  c11_identifier = c11_varName;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_dPd), &c11_thisId, c11_y);
  sf_mex_destroy(&c11_dPd);
  for (c11_i20 = 0; c11_i20 < 2; c11_i20++) {
    (*(real_T (*)[2])c11_outData)[c11_i20] = c11_y[c11_i20];
  }

  sf_mex_destroy(&c11_mxArrayInData);
}

static const mxArray *c11_b_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  real_T c11_u;
  const mxArray *c11_y = NULL;
  SFc11_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc11_simulationInstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  c11_u = *(real_T *)c11_inData;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", &c11_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, FALSE);
  return c11_mxArrayOutData;
}

static real_T c11_c_emlrt_marshallIn(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId)
{
  real_T c11_y;
  real_T c11_d0;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), &c11_d0, 1, 0, 0U, 0, 0U, 0);
  c11_y = c11_d0;
  sf_mex_destroy(&c11_u);
  return c11_y;
}

static void c11_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData)
{
  const mxArray *c11_nargout;
  const char_T *c11_identifier;
  emlrtMsgIdentifier c11_thisId;
  real_T c11_y;
  SFc11_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc11_simulationInstanceStruct *)chartInstanceVoid;
  c11_nargout = sf_mex_dup(c11_mxArrayInData);
  c11_identifier = c11_varName;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_y = c11_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_nargout),
    &c11_thisId);
  sf_mex_destroy(&c11_nargout);
  *(real_T *)c11_outData = c11_y;
  sf_mex_destroy(&c11_mxArrayInData);
}

const mxArray *sf_c11_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c11_nameCaptureInfo = NULL;
  c11_nameCaptureInfo = NULL;
  sf_mex_assign(&c11_nameCaptureInfo, sf_mex_createstruct("structure", 2, 11, 1),
                FALSE);
  c11_info_helper(&c11_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c11_nameCaptureInfo);
  return c11_nameCaptureInfo;
}

static void c11_info_helper(const mxArray **c11_info)
{
  const mxArray *c11_rhs0 = NULL;
  const mxArray *c11_lhs0 = NULL;
  const mxArray *c11_rhs1 = NULL;
  const mxArray *c11_lhs1 = NULL;
  const mxArray *c11_rhs2 = NULL;
  const mxArray *c11_lhs2 = NULL;
  const mxArray *c11_rhs3 = NULL;
  const mxArray *c11_lhs3 = NULL;
  const mxArray *c11_rhs4 = NULL;
  const mxArray *c11_lhs4 = NULL;
  const mxArray *c11_rhs5 = NULL;
  const mxArray *c11_lhs5 = NULL;
  const mxArray *c11_rhs6 = NULL;
  const mxArray *c11_lhs6 = NULL;
  const mxArray *c11_rhs7 = NULL;
  const mxArray *c11_lhs7 = NULL;
  const mxArray *c11_rhs8 = NULL;
  const mxArray *c11_lhs8 = NULL;
  const mxArray *c11_rhs9 = NULL;
  const mxArray *c11_lhs9 = NULL;
  const mxArray *c11_rhs10 = NULL;
  const mxArray *c11_lhs10 = NULL;
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("mtimes"), "name", "name", 0);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c11_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 1);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c11_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("mrdivide"), "name", "name",
                  2);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c11_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 3);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("rdivide"), "name", "name", 3);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c11_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c11_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 5);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c11_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_div"), "name", "name", 6);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c11_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "context", "context", 7);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("sin"), "name", "name", 7);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c11_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 8);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c11_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(""), "context", "context", 9);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("cos"), "name", "name", 9);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c11_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 10);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c11_info, c11_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c11_info, c11_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c11_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c11_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c11_info, sf_mex_duplicatearraysafe(&c11_lhs10), "lhs", "lhs",
                  10);
  sf_mex_destroy(&c11_rhs0);
  sf_mex_destroy(&c11_lhs0);
  sf_mex_destroy(&c11_rhs1);
  sf_mex_destroy(&c11_lhs1);
  sf_mex_destroy(&c11_rhs2);
  sf_mex_destroy(&c11_lhs2);
  sf_mex_destroy(&c11_rhs3);
  sf_mex_destroy(&c11_lhs3);
  sf_mex_destroy(&c11_rhs4);
  sf_mex_destroy(&c11_lhs4);
  sf_mex_destroy(&c11_rhs5);
  sf_mex_destroy(&c11_lhs5);
  sf_mex_destroy(&c11_rhs6);
  sf_mex_destroy(&c11_lhs6);
  sf_mex_destroy(&c11_rhs7);
  sf_mex_destroy(&c11_lhs7);
  sf_mex_destroy(&c11_rhs8);
  sf_mex_destroy(&c11_lhs8);
  sf_mex_destroy(&c11_rhs9);
  sf_mex_destroy(&c11_lhs9);
  sf_mex_destroy(&c11_rhs10);
  sf_mex_destroy(&c11_lhs10);
}

static const mxArray *c11_emlrt_marshallOut(char * c11_u)
{
  const mxArray *c11_y = NULL;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", c11_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c11_u)), FALSE);
  return c11_y;
}

static const mxArray *c11_b_emlrt_marshallOut(uint32_T c11_u)
{
  const mxArray *c11_y = NULL;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", &c11_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c11_y;
}

static const mxArray *c11_c_sf_marshallOut(void *chartInstanceVoid, void
  *c11_inData)
{
  const mxArray *c11_mxArrayOutData = NULL;
  int32_T c11_u;
  const mxArray *c11_y = NULL;
  SFc11_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc11_simulationInstanceStruct *)chartInstanceVoid;
  c11_mxArrayOutData = NULL;
  c11_u = *(int32_T *)c11_inData;
  c11_y = NULL;
  sf_mex_assign(&c11_y, sf_mex_create("y", &c11_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c11_mxArrayOutData, c11_y, FALSE);
  return c11_mxArrayOutData;
}

static int32_T c11_d_emlrt_marshallIn(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId)
{
  int32_T c11_y;
  int32_T c11_i21;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), &c11_i21, 1, 6, 0U, 0, 0U, 0);
  c11_y = c11_i21;
  sf_mex_destroy(&c11_u);
  return c11_y;
}

static void c11_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c11_mxArrayInData, const char_T *c11_varName, void *c11_outData)
{
  const mxArray *c11_b_sfEvent;
  const char_T *c11_identifier;
  emlrtMsgIdentifier c11_thisId;
  int32_T c11_y;
  SFc11_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc11_simulationInstanceStruct *)chartInstanceVoid;
  c11_b_sfEvent = sf_mex_dup(c11_mxArrayInData);
  c11_identifier = c11_varName;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_y = c11_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c11_b_sfEvent),
    &c11_thisId);
  sf_mex_destroy(&c11_b_sfEvent);
  *(int32_T *)c11_outData = c11_y;
  sf_mex_destroy(&c11_mxArrayInData);
}

static uint8_T c11_e_emlrt_marshallIn(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_b_is_active_c11_simulation, const char_T
  *c11_identifier)
{
  uint8_T c11_y;
  emlrtMsgIdentifier c11_thisId;
  c11_thisId.fIdentifier = c11_identifier;
  c11_thisId.fParent = NULL;
  c11_y = c11_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c11_b_is_active_c11_simulation), &c11_thisId);
  sf_mex_destroy(&c11_b_is_active_c11_simulation);
  return c11_y;
}

static uint8_T c11_f_emlrt_marshallIn(SFc11_simulationInstanceStruct
  *chartInstance, const mxArray *c11_u, const emlrtMsgIdentifier *c11_parentId)
{
  uint8_T c11_y;
  uint8_T c11_u0;
  sf_mex_import(c11_parentId, sf_mex_dup(c11_u), &c11_u0, 1, 3, 0U, 0, 0U, 0);
  c11_y = c11_u0;
  sf_mex_destroy(&c11_u);
  return c11_y;
}

static void init_dsm_address_info(SFc11_simulationInstanceStruct *chartInstance)
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

void sf_c11_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1094584075U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3420162766U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(887642851U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1149413169U);
}

mxArray *sf_c11_simulation_get_autoinheritance_info(void)
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

mxArray *sf_c11_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c11_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c11_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[7],T\"Pd\",},{M[1],M[5],T\"dPd\",},{M[8],M[0],T\"is_active_c11_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c11_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc11_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc11_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           11,
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
          (MexFcnForType)c11_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c11_sf_marshallOut,(MexInFcnForType)
            c11_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c11_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c11_sf_marshallOut,(MexInFcnForType)
            c11_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c11_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c11_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c11_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c11_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c11_gam;
          real_T *c11_gamWP;
          real_T *c11_thetad;
          real_T *c11_v_L;
          real_T *c11_Radius;
          real_T (*c11_Pd)[2];
          real_T (*c11_dPd)[2];
          real_T (*c11_P0)[2];
          c11_Radius = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
          c11_v_L = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c11_P0 = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
          c11_thetad = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c11_dPd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
          c11_gamWP = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c11_Pd = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
          c11_gam = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c11_gam);
          _SFD_SET_DATA_VALUE_PTR(1U, *c11_Pd);
          _SFD_SET_DATA_VALUE_PTR(2U, c11_gamWP);
          _SFD_SET_DATA_VALUE_PTR(3U, *c11_dPd);
          _SFD_SET_DATA_VALUE_PTR(4U, c11_thetad);
          _SFD_SET_DATA_VALUE_PTR(5U, *c11_P0);
          _SFD_SET_DATA_VALUE_PTR(6U, c11_v_L);
          _SFD_SET_DATA_VALUE_PTR(7U, c11_Radius);
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

static void sf_opaque_initialize_c11_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc11_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c11_simulation((SFc11_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c11_simulation((SFc11_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c11_simulation(void *chartInstanceVar)
{
  enable_c11_simulation((SFc11_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c11_simulation(void *chartInstanceVar)
{
  disable_c11_simulation((SFc11_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c11_simulation(void *chartInstanceVar)
{
  sf_c11_simulation((SFc11_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c11_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c11_simulation
    ((SFc11_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c11_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c11_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c11_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c11_simulation((SFc11_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c11_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c11_simulation(S);
}

static void sf_opaque_set_sim_state_c11_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c11_simulation(S, st);
}

static void sf_opaque_terminate_c11_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc11_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c11_simulation((SFc11_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc11_simulation((SFc11_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c11_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c11_simulation((SFc11_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c11_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      11);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,11,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,11,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,11);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,11,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,11,2);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,11);
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

static void mdlRTW_c11_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c11_simulation(SimStruct *S)
{
  SFc11_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc11_simulationInstanceStruct *)utMalloc(sizeof
    (SFc11_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc11_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c11_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c11_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c11_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c11_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c11_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c11_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c11_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c11_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c11_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c11_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c11_simulation;
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

void c11_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c11_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c11_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c11_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c11_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
