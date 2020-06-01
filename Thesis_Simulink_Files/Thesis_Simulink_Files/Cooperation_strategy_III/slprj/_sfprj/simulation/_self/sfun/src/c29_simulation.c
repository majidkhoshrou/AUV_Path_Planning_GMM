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
static const char * c29_debug_family_names[15] = { "u", "v", "F_km1", "U_km1",
  "nargin", "nargout", "Sigma_km1", "State_km1", "T_est", "Vc", "Q", "Psi",
  "Usens", "State_predict", "Sigma_predict" };

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
  const mxArray *c29_Sigma_predict, const char_T *c29_identifier, real_T c29_y
  [81]);
static void c29_b_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId, real_T c29_y[81]);
static void c29_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData);
static const mxArray *c29_b_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData);
static void c29_c_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_State_predict, const char_T *c29_identifier, real_T c29_y[9]);
static void c29_d_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId, real_T c29_y[9]);
static void c29_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData);
static const mxArray *c29_c_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData);
static const mxArray *c29_d_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData);
static const mxArray *c29_e_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData);
static real_T c29_e_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId);
static void c29_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData);
static void c29_info_helper(const mxArray **c29_info);
static const mxArray *c29_emlrt_marshallOut(char * c29_u);
static const mxArray *c29_b_emlrt_marshallOut(uint32_T c29_u);
static void c29_eml_scalar_eg(SFc29_simulationInstanceStruct *chartInstance);
static void c29_b_eml_scalar_eg(SFc29_simulationInstanceStruct *chartInstance);
static void c29_eml_xgemm(SFc29_simulationInstanceStruct *chartInstance, real_T
  c29_A[81], real_T c29_B[81], real_T c29_C[81], real_T c29_b_C[81]);
static const mxArray *c29_f_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData);
static int32_T c29_f_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId);
static void c29_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData);
static uint8_T c29_g_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_b_is_active_c29_simulation, const char_T
  *c29_identifier);
static uint8_T c29_h_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId);
static void c29_b_eml_xgemm(SFc29_simulationInstanceStruct *chartInstance,
  real_T c29_A[81], real_T c29_B[81], real_T c29_C[81]);
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
  real_T c29_d0;
  sf_set_error_prefix_string(
    "Error evaluating data 'T_est' in the parent workspace.\n");
  sf_mex_import_named("T_est", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c29_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c29_T_est = c29_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
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
  real_T c29_u[81];
  const mxArray *c29_b_y = NULL;
  int32_T c29_i1;
  real_T c29_b_u[9];
  const mxArray *c29_c_y = NULL;
  uint8_T c29_hoistedGlobal;
  uint8_T c29_c_u;
  const mxArray *c29_d_y = NULL;
  real_T (*c29_State_predict)[9];
  real_T (*c29_Sigma_predict)[81];
  c29_Sigma_predict = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S, 2);
  c29_State_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c29_st = NULL;
  c29_st = NULL;
  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_createcellarray(3), FALSE);
  for (c29_i0 = 0; c29_i0 < 81; c29_i0++) {
    c29_u[c29_i0] = (*c29_Sigma_predict)[c29_i0];
  }

  c29_b_y = NULL;
  sf_mex_assign(&c29_b_y, sf_mex_create("y", c29_u, 0, 0U, 1U, 0U, 2, 9, 9),
                FALSE);
  sf_mex_setcell(c29_y, 0, c29_b_y);
  for (c29_i1 = 0; c29_i1 < 9; c29_i1++) {
    c29_b_u[c29_i1] = (*c29_State_predict)[c29_i1];
  }

  c29_c_y = NULL;
  sf_mex_assign(&c29_c_y, sf_mex_create("y", c29_b_u, 0, 0U, 1U, 0U, 1, 9),
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
  real_T c29_dv0[81];
  int32_T c29_i2;
  real_T c29_dv1[9];
  int32_T c29_i3;
  real_T (*c29_Sigma_predict)[81];
  real_T (*c29_State_predict)[9];
  c29_Sigma_predict = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S, 2);
  c29_State_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c29_doneDoubleBufferReInit = TRUE;
  c29_u = sf_mex_dup(c29_st);
  c29_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c29_u, 0)),
                       "Sigma_predict", c29_dv0);
  for (c29_i2 = 0; c29_i2 < 81; c29_i2++) {
    (*c29_Sigma_predict)[c29_i2] = c29_dv0[c29_i2];
  }

  c29_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c29_u, 1)),
    "State_predict", c29_dv1);
  for (c29_i3 = 0; c29_i3 < 9; c29_i3++) {
    (*c29_State_predict)[c29_i3] = c29_dv1[c29_i3];
  }

  chartInstance->c29_is_active_c29_simulation = c29_g_emlrt_marshallIn
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
  int32_T c29_i7;
  int32_T c29_i8;
  int32_T c29_i9;
  int32_T c29_i10;
  real_T *c29_Psi;
  real_T (*c29_Sigma_predict)[81];
  real_T (*c29_State_predict)[9];
  real_T (*c29_Usens)[3];
  real_T (*c29_Q)[81];
  real_T (*c29_Vc)[2];
  real_T (*c29_State_km1)[9];
  real_T (*c29_Sigma_km1)[81];
  c29_Sigma_predict = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S, 2);
  c29_State_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c29_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c29_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c29_Q = (real_T (*)[81])ssGetInputPortSignal(chartInstance->S, 3);
  c29_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
  c29_State_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 1);
  c29_Sigma_km1 = (real_T (*)[81])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 28U, chartInstance->c29_sfEvent);
  for (c29_i4 = 0; c29_i4 < 81; c29_i4++) {
    _SFD_DATA_RANGE_CHECK((*c29_Sigma_km1)[c29_i4], 0U);
  }

  for (c29_i5 = 0; c29_i5 < 9; c29_i5++) {
    _SFD_DATA_RANGE_CHECK((*c29_State_km1)[c29_i5], 1U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c29_T_est, 2U);
  for (c29_i6 = 0; c29_i6 < 2; c29_i6++) {
    _SFD_DATA_RANGE_CHECK((*c29_Vc)[c29_i6], 3U);
  }

  for (c29_i7 = 0; c29_i7 < 81; c29_i7++) {
    _SFD_DATA_RANGE_CHECK((*c29_Q)[c29_i7], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*c29_Psi, 5U);
  for (c29_i8 = 0; c29_i8 < 3; c29_i8++) {
    _SFD_DATA_RANGE_CHECK((*c29_Usens)[c29_i8], 6U);
  }

  for (c29_i9 = 0; c29_i9 < 9; c29_i9++) {
    _SFD_DATA_RANGE_CHECK((*c29_State_predict)[c29_i9], 7U);
  }

  for (c29_i10 = 0; c29_i10 < 81; c29_i10++) {
    _SFD_DATA_RANGE_CHECK((*c29_Sigma_predict)[c29_i10], 8U);
  }

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
  int32_T c29_i11;
  real_T c29_Sigma_km1[81];
  int32_T c29_i12;
  real_T c29_State_km1[9];
  real_T c29_b_T_est;
  int32_T c29_i13;
  real_T c29_Vc[2];
  int32_T c29_i14;
  real_T c29_Q[81];
  real_T c29_Psi;
  int32_T c29_i15;
  real_T c29_Usens[3];
  uint32_T c29_debug_family_var_map[15];
  real_T c29_u;
  real_T c29_v;
  real_T c29_F_km1[81];
  real_T c29_U_km1[81];
  real_T c29_nargin = 7.0;
  real_T c29_nargout = 2.0;
  real_T c29_State_predict[9];
  real_T c29_Sigma_predict[81];
  int32_T c29_i16;
  static real_T c29_a[81] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c29_x;
  real_T c29_b_x;
  real_T c29_b_a;
  real_T c29_b;
  real_T c29_y;
  real_T c29_c_x;
  real_T c29_d_x;
  real_T c29_c_a;
  real_T c29_b_b;
  real_T c29_b_y;
  real_T c29_e_x;
  real_T c29_f_x;
  real_T c29_d_a;
  real_T c29_c_b;
  real_T c29_c_y;
  real_T c29_g_x;
  real_T c29_h_x;
  real_T c29_e_a;
  real_T c29_d_b;
  real_T c29_d_y;
  real_T c29_f_a;
  real_T c29_e_b;
  real_T c29_e_y;
  real_T c29_i_x;
  real_T c29_j_x;
  real_T c29_g_a;
  real_T c29_f_b;
  real_T c29_f_y;
  real_T c29_k_x;
  real_T c29_l_x;
  real_T c29_h_a;
  real_T c29_g_b;
  real_T c29_g_y;
  real_T c29_m_x;
  real_T c29_n_x;
  real_T c29_i_a;
  real_T c29_h_b;
  real_T c29_h_y;
  real_T c29_o_x;
  real_T c29_p_x;
  real_T c29_j_a;
  real_T c29_i_b;
  real_T c29_i_y;
  real_T c29_k_a;
  real_T c29_j_b;
  real_T c29_j_y;
  int32_T c29_i17;
  int32_T c29_i18;
  int32_T c29_i19;
  int32_T c29_i20;
  int32_T c29_i21;
  int32_T c29_i22;
  static real_T c29_dv2[9] = { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  int32_T c29_i23;
  real_T c29_k_b[9];
  int32_T c29_i24;
  real_T c29_k_y[9];
  int32_T c29_i25;
  int32_T c29_i26;
  real_T c29_l_a;
  real_T c29_l_b;
  real_T c29_l_y;
  real_T c29_m_a;
  real_T c29_m_b;
  real_T c29_m_y;
  real_T c29_q_x;
  real_T c29_r_x;
  real_T c29_n_a;
  real_T c29_n_b;
  real_T c29_n_y;
  real_T c29_s_x;
  real_T c29_t_x;
  real_T c29_o_a;
  real_T c29_o_b;
  real_T c29_o_y;
  real_T c29_p_a;
  real_T c29_p_b;
  real_T c29_p_y;
  real_T c29_u_x;
  real_T c29_v_x;
  real_T c29_q_a;
  real_T c29_q_b;
  real_T c29_q_y;
  real_T c29_w_x;
  real_T c29_x_x;
  real_T c29_r_a;
  real_T c29_r_b;
  real_T c29_r_y;
  real_T c29_s_a;
  real_T c29_s_b;
  real_T c29_s_y;
  real_T c29_t_y[9];
  int32_T c29_i27;
  real_T c29_u_y[9];
  int32_T c29_i28;
  int32_T c29_i29;
  int32_T c29_i30;
  real_T c29_t_b[81];
  int32_T c29_i31;
  real_T c29_v_y[81];
  int32_T c29_i32;
  real_T c29_t_a[81];
  int32_T c29_i33;
  real_T c29_u_b[81];
  int32_T c29_i34;
  real_T c29_w_y[81];
  int32_T c29_i35;
  real_T c29_x_y[81];
  int32_T c29_i36;
  real_T c29_u_a[81];
  int32_T c29_i37;
  real_T c29_v_a[81];
  int32_T c29_i38;
  int32_T c29_i39;
  int32_T c29_i40;
  real_T c29_w_a[81];
  int32_T c29_i41;
  real_T c29_v_b[81];
  int32_T c29_i42;
  int32_T c29_i43;
  int32_T c29_i44;
  int32_T c29_i45;
  int32_T c29_i46;
  int32_T c29_i47;
  real_T c29_y_y[81];
  int32_T c29_i48;
  real_T c29_w_b[81];
  int32_T c29_i49;
  int32_T c29_i50;
  int32_T c29_i51;
  real_T *c29_b_Psi;
  real_T (*c29_b_State_predict)[9];
  real_T (*c29_b_Sigma_predict)[81];
  real_T (*c29_b_Usens)[3];
  real_T (*c29_b_Q)[81];
  real_T (*c29_b_Vc)[2];
  real_T (*c29_b_State_km1)[9];
  real_T (*c29_b_Sigma_km1)[81];
  c29_b_Sigma_predict = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S,
    2);
  c29_b_State_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c29_b_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c29_b_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c29_b_Q = (real_T (*)[81])ssGetInputPortSignal(chartInstance->S, 3);
  c29_b_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
  c29_b_State_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 1);
  c29_b_Sigma_km1 = (real_T (*)[81])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 28U, chartInstance->c29_sfEvent);
  c29_hoistedGlobal = chartInstance->c29_T_est;
  c29_b_hoistedGlobal = *c29_b_Psi;
  for (c29_i11 = 0; c29_i11 < 81; c29_i11++) {
    c29_Sigma_km1[c29_i11] = (*c29_b_Sigma_km1)[c29_i11];
  }

  for (c29_i12 = 0; c29_i12 < 9; c29_i12++) {
    c29_State_km1[c29_i12] = (*c29_b_State_km1)[c29_i12];
  }

  c29_b_T_est = c29_hoistedGlobal;
  for (c29_i13 = 0; c29_i13 < 2; c29_i13++) {
    c29_Vc[c29_i13] = (*c29_b_Vc)[c29_i13];
  }

  for (c29_i14 = 0; c29_i14 < 81; c29_i14++) {
    c29_Q[c29_i14] = (*c29_b_Q)[c29_i14];
  }

  c29_Psi = c29_b_hoistedGlobal;
  for (c29_i15 = 0; c29_i15 < 3; c29_i15++) {
    c29_Usens[c29_i15] = (*c29_b_Usens)[c29_i15];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 15U, 15U, c29_debug_family_names,
    c29_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_u, 0U, c29_d_sf_marshallOut,
    c29_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_v, 1U, c29_d_sf_marshallOut,
    c29_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c29_F_km1, 2U, c29_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c29_U_km1, 3U, c29_sf_marshallOut,
    c29_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_nargin, 4U, c29_d_sf_marshallOut,
    c29_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_nargout, 5U, c29_d_sf_marshallOut,
    c29_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c29_Sigma_km1, 6U, c29_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c29_State_km1, 7U, c29_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c29_b_T_est, 8U, c29_d_sf_marshallOut,
    c29_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c29_Vc, 9U, c29_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c29_Q, 10U, c29_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c29_Psi, 11U, c29_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c29_Usens, 12U, c29_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c29_State_predict, 13U,
    c29_b_sf_marshallOut, c29_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c29_Sigma_predict, 14U,
    c29_sf_marshallOut, c29_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 4);
  c29_u = c29_Usens[0];
  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 5);
  c29_v = c29_Usens[1];
  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 8);
  for (c29_i16 = 0; c29_i16 < 81; c29_i16++) {
    c29_F_km1[c29_i16] = c29_a[c29_i16];
  }

  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 18);
  c29_x = c29_Psi;
  c29_b_x = c29_x;
  c29_b_x = muDoubleScalarCos(c29_b_x);
  c29_b_a = c29_b_x;
  c29_b = c29_b_T_est;
  c29_y = c29_b_a * c29_b;
  c29_c_x = c29_Psi;
  c29_d_x = c29_c_x;
  c29_d_x = muDoubleScalarSin(c29_d_x);
  c29_c_a = -c29_d_x;
  c29_b_b = c29_b_T_est;
  c29_b_y = c29_c_a * c29_b_b;
  c29_e_x = c29_Psi;
  c29_f_x = c29_e_x;
  c29_f_x = muDoubleScalarSin(c29_f_x);
  c29_d_a = -c29_u;
  c29_c_b = c29_f_x;
  c29_c_y = c29_d_a * c29_c_b;
  c29_g_x = c29_Psi;
  c29_h_x = c29_g_x;
  c29_h_x = muDoubleScalarCos(c29_h_x);
  c29_e_a = c29_v;
  c29_d_b = c29_h_x;
  c29_d_y = c29_e_a * c29_d_b;
  c29_f_a = c29_c_y - c29_d_y;
  c29_e_b = c29_b_T_est;
  c29_e_y = c29_f_a * c29_e_b;
  c29_i_x = c29_Psi;
  c29_j_x = c29_i_x;
  c29_j_x = muDoubleScalarSin(c29_j_x);
  c29_g_a = c29_j_x;
  c29_f_b = c29_b_T_est;
  c29_f_y = c29_g_a * c29_f_b;
  c29_k_x = c29_Psi;
  c29_l_x = c29_k_x;
  c29_l_x = muDoubleScalarCos(c29_l_x);
  c29_h_a = c29_l_x;
  c29_g_b = c29_b_T_est;
  c29_g_y = c29_h_a * c29_g_b;
  c29_m_x = c29_Psi;
  c29_n_x = c29_m_x;
  c29_n_x = muDoubleScalarCos(c29_n_x);
  c29_i_a = c29_u;
  c29_h_b = c29_n_x;
  c29_h_y = c29_i_a * c29_h_b;
  c29_o_x = c29_Psi;
  c29_p_x = c29_o_x;
  c29_p_x = muDoubleScalarSin(c29_p_x);
  c29_j_a = c29_v;
  c29_i_b = c29_p_x;
  c29_i_y = c29_j_a * c29_i_b;
  c29_k_a = c29_h_y - c29_i_y;
  c29_j_b = c29_b_T_est;
  c29_j_y = c29_k_a * c29_j_b;
  c29_U_km1[0] = c29_y;
  c29_U_km1[9] = c29_b_y;
  c29_U_km1[18] = c29_e_y;
  c29_i17 = 0;
  for (c29_i18 = 0; c29_i18 < 6; c29_i18++) {
    c29_U_km1[c29_i17 + 27] = 0.0;
    c29_i17 += 9;
  }

  c29_U_km1[1] = c29_f_y;
  c29_U_km1[10] = c29_g_y;
  c29_U_km1[19] = c29_j_y;
  c29_i19 = 0;
  for (c29_i20 = 0; c29_i20 < 6; c29_i20++) {
    c29_U_km1[c29_i19 + 28] = 0.0;
    c29_i19 += 9;
  }

  c29_i21 = 0;
  for (c29_i22 = 0; c29_i22 < 9; c29_i22++) {
    c29_U_km1[c29_i21 + 2] = c29_dv2[c29_i22];
    c29_i21 += 9;
  }

  c29_U_km1[3] = 0.0;
  c29_U_km1[12] = 0.0;
  c29_U_km1[21] = 0.0;
  c29_U_km1[30] = c29_b_T_est;
  c29_U_km1[39] = 0.0;
  c29_U_km1[48] = 0.0;
  c29_U_km1[57] = 0.0;
  c29_U_km1[66] = 0.0;
  c29_U_km1[75] = 0.0;
  c29_U_km1[4] = 0.0;
  c29_U_km1[13] = 0.0;
  c29_U_km1[22] = 0.0;
  c29_U_km1[31] = 0.0;
  c29_U_km1[40] = c29_b_T_est;
  c29_U_km1[49] = 0.0;
  c29_U_km1[58] = 0.0;
  c29_U_km1[67] = 0.0;
  c29_U_km1[76] = 0.0;
  c29_U_km1[5] = 0.0;
  c29_U_km1[14] = 0.0;
  c29_U_km1[23] = 0.0;
  c29_U_km1[32] = 0.0;
  c29_U_km1[41] = 0.0;
  c29_U_km1[50] = c29_b_T_est;
  c29_U_km1[59] = 0.0;
  c29_U_km1[68] = 0.0;
  c29_U_km1[77] = 0.0;
  c29_U_km1[6] = 0.0;
  c29_U_km1[15] = 0.0;
  c29_U_km1[24] = 0.0;
  c29_U_km1[33] = 0.0;
  c29_U_km1[42] = 0.0;
  c29_U_km1[51] = 0.0;
  c29_U_km1[60] = c29_b_T_est;
  c29_U_km1[69] = 0.0;
  c29_U_km1[78] = 0.0;
  c29_U_km1[7] = 0.0;
  c29_U_km1[16] = 0.0;
  c29_U_km1[25] = 0.0;
  c29_U_km1[34] = 0.0;
  c29_U_km1[43] = 0.0;
  c29_U_km1[52] = 0.0;
  c29_U_km1[61] = 0.0;
  c29_U_km1[70] = c29_b_T_est;
  c29_U_km1[79] = 0.0;
  c29_U_km1[8] = 0.0;
  c29_U_km1[17] = 0.0;
  c29_U_km1[26] = 0.0;
  c29_U_km1[35] = 0.0;
  c29_U_km1[44] = 0.0;
  c29_U_km1[53] = 0.0;
  c29_U_km1[62] = 0.0;
  c29_U_km1[71] = 0.0;
  c29_U_km1[80] = c29_b_T_est;
  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 31);
  for (c29_i23 = 0; c29_i23 < 9; c29_i23++) {
    c29_k_b[c29_i23] = c29_State_km1[c29_i23];
  }

  c29_eml_scalar_eg(chartInstance);
  c29_eml_scalar_eg(chartInstance);
  for (c29_i24 = 0; c29_i24 < 9; c29_i24++) {
    c29_k_y[c29_i24] = 0.0;
    c29_i25 = 0;
    for (c29_i26 = 0; c29_i26 < 9; c29_i26++) {
      c29_k_y[c29_i24] += c29_a[c29_i25 + c29_i24] * c29_k_b[c29_i26];
      c29_i25 += 9;
    }
  }

  c29_l_a = c29_Vc[0];
  c29_l_b = c29_b_T_est;
  c29_l_y = c29_l_a * c29_l_b;
  c29_m_a = c29_Vc[1];
  c29_m_b = c29_b_T_est;
  c29_m_y = c29_m_a * c29_m_b;
  c29_q_x = c29_Psi;
  c29_r_x = c29_q_x;
  c29_r_x = muDoubleScalarCos(c29_r_x);
  c29_n_a = c29_u;
  c29_n_b = c29_r_x;
  c29_n_y = c29_n_a * c29_n_b;
  c29_s_x = c29_Psi;
  c29_t_x = c29_s_x;
  c29_t_x = muDoubleScalarSin(c29_t_x);
  c29_o_a = c29_v;
  c29_o_b = c29_t_x;
  c29_o_y = c29_o_a * c29_o_b;
  c29_p_a = c29_n_y - c29_o_y;
  c29_p_b = c29_b_T_est;
  c29_p_y = c29_p_a * c29_p_b;
  c29_u_x = c29_Psi;
  c29_v_x = c29_u_x;
  c29_v_x = muDoubleScalarSin(c29_v_x);
  c29_q_a = c29_u;
  c29_q_b = c29_v_x;
  c29_q_y = c29_q_a * c29_q_b;
  c29_w_x = c29_Psi;
  c29_x_x = c29_w_x;
  c29_x_x = muDoubleScalarCos(c29_x_x);
  c29_r_a = c29_v;
  c29_r_b = c29_x_x;
  c29_r_y = c29_r_a * c29_r_b;
  c29_s_a = c29_q_y + c29_r_y;
  c29_s_b = c29_b_T_est;
  c29_s_y = c29_s_a * c29_s_b;
  c29_t_y[0] = c29_l_y;
  c29_t_y[1] = c29_m_y;
  c29_t_y[2] = 0.0;
  for (c29_i27 = 0; c29_i27 < 6; c29_i27++) {
    c29_t_y[c29_i27 + 3] = 0.0;
  }

  c29_u_y[0] = c29_p_y;
  c29_u_y[1] = c29_s_y;
  c29_u_y[2] = c29_Psi;
  for (c29_i28 = 0; c29_i28 < 6; c29_i28++) {
    c29_u_y[c29_i28 + 3] = 0.0;
  }

  for (c29_i29 = 0; c29_i29 < 9; c29_i29++) {
    c29_State_predict[c29_i29] = (c29_k_y[c29_i29] + c29_t_y[c29_i29]) +
      c29_u_y[c29_i29];
  }

  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, 36);
  for (c29_i30 = 0; c29_i30 < 81; c29_i30++) {
    c29_t_b[c29_i30] = c29_Sigma_km1[c29_i30];
  }

  c29_b_eml_scalar_eg(chartInstance);
  c29_b_eml_scalar_eg(chartInstance);
  for (c29_i31 = 0; c29_i31 < 81; c29_i31++) {
    c29_v_y[c29_i31] = 0.0;
  }

  for (c29_i32 = 0; c29_i32 < 81; c29_i32++) {
    c29_t_a[c29_i32] = c29_a[c29_i32];
  }

  for (c29_i33 = 0; c29_i33 < 81; c29_i33++) {
    c29_u_b[c29_i33] = c29_t_b[c29_i33];
  }

  c29_b_eml_xgemm(chartInstance, c29_t_a, c29_u_b, c29_v_y);
  c29_b_eml_scalar_eg(chartInstance);
  c29_b_eml_scalar_eg(chartInstance);
  for (c29_i34 = 0; c29_i34 < 81; c29_i34++) {
    c29_w_y[c29_i34] = 0.0;
  }

  for (c29_i35 = 0; c29_i35 < 81; c29_i35++) {
    c29_x_y[c29_i35] = c29_v_y[c29_i35];
  }

  for (c29_i36 = 0; c29_i36 < 81; c29_i36++) {
    c29_u_a[c29_i36] = c29_a[c29_i36];
  }

  c29_b_eml_xgemm(chartInstance, c29_x_y, c29_u_a, c29_w_y);
  for (c29_i37 = 0; c29_i37 < 81; c29_i37++) {
    c29_v_a[c29_i37] = c29_U_km1[c29_i37];
  }

  for (c29_i38 = 0; c29_i38 < 81; c29_i38++) {
    c29_t_b[c29_i38] = c29_Q[c29_i38];
  }

  c29_b_eml_scalar_eg(chartInstance);
  c29_b_eml_scalar_eg(chartInstance);
  for (c29_i39 = 0; c29_i39 < 81; c29_i39++) {
    c29_v_y[c29_i39] = 0.0;
  }

  for (c29_i40 = 0; c29_i40 < 81; c29_i40++) {
    c29_w_a[c29_i40] = c29_v_a[c29_i40];
  }

  for (c29_i41 = 0; c29_i41 < 81; c29_i41++) {
    c29_v_b[c29_i41] = c29_t_b[c29_i41];
  }

  c29_b_eml_xgemm(chartInstance, c29_w_a, c29_v_b, c29_v_y);
  c29_i42 = 0;
  for (c29_i43 = 0; c29_i43 < 9; c29_i43++) {
    c29_i44 = 0;
    for (c29_i45 = 0; c29_i45 < 9; c29_i45++) {
      c29_t_b[c29_i45 + c29_i42] = c29_U_km1[c29_i44 + c29_i43];
      c29_i44 += 9;
    }

    c29_i42 += 9;
  }

  c29_b_eml_scalar_eg(chartInstance);
  c29_b_eml_scalar_eg(chartInstance);
  for (c29_i46 = 0; c29_i46 < 81; c29_i46++) {
    c29_v_a[c29_i46] = 0.0;
  }

  for (c29_i47 = 0; c29_i47 < 81; c29_i47++) {
    c29_y_y[c29_i47] = c29_v_y[c29_i47];
  }

  for (c29_i48 = 0; c29_i48 < 81; c29_i48++) {
    c29_w_b[c29_i48] = c29_t_b[c29_i48];
  }

  c29_b_eml_xgemm(chartInstance, c29_y_y, c29_w_b, c29_v_a);
  for (c29_i49 = 0; c29_i49 < 81; c29_i49++) {
    c29_Sigma_predict[c29_i49] = c29_w_y[c29_i49] + c29_v_a[c29_i49];
  }

  _SFD_EML_CALL(0U, chartInstance->c29_sfEvent, -36);
  _SFD_SYMBOL_SCOPE_POP();
  for (c29_i50 = 0; c29_i50 < 9; c29_i50++) {
    (*c29_b_State_predict)[c29_i50] = c29_State_predict[c29_i50];
  }

  for (c29_i51 = 0; c29_i51 < 81; c29_i51++) {
    (*c29_b_Sigma_predict)[c29_i51] = c29_Sigma_predict[c29_i51];
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
  int32_T c29_i52;
  int32_T c29_i53;
  int32_T c29_i54;
  real_T c29_b_inData[81];
  int32_T c29_i55;
  int32_T c29_i56;
  int32_T c29_i57;
  real_T c29_u[81];
  const mxArray *c29_y = NULL;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_mxArrayOutData = NULL;
  c29_i52 = 0;
  for (c29_i53 = 0; c29_i53 < 9; c29_i53++) {
    for (c29_i54 = 0; c29_i54 < 9; c29_i54++) {
      c29_b_inData[c29_i54 + c29_i52] = (*(real_T (*)[81])c29_inData)[c29_i54 +
        c29_i52];
    }

    c29_i52 += 9;
  }

  c29_i55 = 0;
  for (c29_i56 = 0; c29_i56 < 9; c29_i56++) {
    for (c29_i57 = 0; c29_i57 < 9; c29_i57++) {
      c29_u[c29_i57 + c29_i55] = c29_b_inData[c29_i57 + c29_i55];
    }

    c29_i55 += 9;
  }

  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", c29_u, 0, 0U, 1U, 0U, 2, 9, 9), FALSE);
  sf_mex_assign(&c29_mxArrayOutData, c29_y, FALSE);
  return c29_mxArrayOutData;
}

static void c29_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_Sigma_predict, const char_T *c29_identifier, real_T c29_y
  [81])
{
  emlrtMsgIdentifier c29_thisId;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_Sigma_predict),
    &c29_thisId, c29_y);
  sf_mex_destroy(&c29_Sigma_predict);
}

static void c29_b_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId, real_T c29_y[81])
{
  real_T c29_dv3[81];
  int32_T c29_i58;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), c29_dv3, 1, 0, 0U, 1, 0U, 2, 9,
                9);
  for (c29_i58 = 0; c29_i58 < 81; c29_i58++) {
    c29_y[c29_i58] = c29_dv3[c29_i58];
  }

  sf_mex_destroy(&c29_u);
}

static void c29_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData)
{
  const mxArray *c29_Sigma_predict;
  const char_T *c29_identifier;
  emlrtMsgIdentifier c29_thisId;
  real_T c29_y[81];
  int32_T c29_i59;
  int32_T c29_i60;
  int32_T c29_i61;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_Sigma_predict = sf_mex_dup(c29_mxArrayInData);
  c29_identifier = c29_varName;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_Sigma_predict),
    &c29_thisId, c29_y);
  sf_mex_destroy(&c29_Sigma_predict);
  c29_i59 = 0;
  for (c29_i60 = 0; c29_i60 < 9; c29_i60++) {
    for (c29_i61 = 0; c29_i61 < 9; c29_i61++) {
      (*(real_T (*)[81])c29_outData)[c29_i61 + c29_i59] = c29_y[c29_i61 +
        c29_i59];
    }

    c29_i59 += 9;
  }

  sf_mex_destroy(&c29_mxArrayInData);
}

static const mxArray *c29_b_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData)
{
  const mxArray *c29_mxArrayOutData = NULL;
  int32_T c29_i62;
  real_T c29_b_inData[9];
  int32_T c29_i63;
  real_T c29_u[9];
  const mxArray *c29_y = NULL;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_mxArrayOutData = NULL;
  for (c29_i62 = 0; c29_i62 < 9; c29_i62++) {
    c29_b_inData[c29_i62] = (*(real_T (*)[9])c29_inData)[c29_i62];
  }

  for (c29_i63 = 0; c29_i63 < 9; c29_i63++) {
    c29_u[c29_i63] = c29_b_inData[c29_i63];
  }

  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", c29_u, 0, 0U, 1U, 0U, 1, 9), FALSE);
  sf_mex_assign(&c29_mxArrayOutData, c29_y, FALSE);
  return c29_mxArrayOutData;
}

static void c29_c_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_State_predict, const char_T *c29_identifier, real_T c29_y[9])
{
  emlrtMsgIdentifier c29_thisId;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_State_predict),
    &c29_thisId, c29_y);
  sf_mex_destroy(&c29_State_predict);
}

static void c29_d_emlrt_marshallIn(SFc29_simulationInstanceStruct *chartInstance,
  const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId, real_T c29_y[9])
{
  real_T c29_dv4[9];
  int32_T c29_i64;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), c29_dv4, 1, 0, 0U, 1, 0U, 1, 9);
  for (c29_i64 = 0; c29_i64 < 9; c29_i64++) {
    c29_y[c29_i64] = c29_dv4[c29_i64];
  }

  sf_mex_destroy(&c29_u);
}

static void c29_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData)
{
  const mxArray *c29_State_predict;
  const char_T *c29_identifier;
  emlrtMsgIdentifier c29_thisId;
  real_T c29_y[9];
  int32_T c29_i65;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_State_predict = sf_mex_dup(c29_mxArrayInData);
  c29_identifier = c29_varName;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_State_predict),
    &c29_thisId, c29_y);
  sf_mex_destroy(&c29_State_predict);
  for (c29_i65 = 0; c29_i65 < 9; c29_i65++) {
    (*(real_T (*)[9])c29_outData)[c29_i65] = c29_y[c29_i65];
  }

  sf_mex_destroy(&c29_mxArrayInData);
}

static const mxArray *c29_c_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData)
{
  const mxArray *c29_mxArrayOutData = NULL;
  int32_T c29_i66;
  real_T c29_b_inData[3];
  int32_T c29_i67;
  real_T c29_u[3];
  const mxArray *c29_y = NULL;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_mxArrayOutData = NULL;
  for (c29_i66 = 0; c29_i66 < 3; c29_i66++) {
    c29_b_inData[c29_i66] = (*(real_T (*)[3])c29_inData)[c29_i66];
  }

  for (c29_i67 = 0; c29_i67 < 3; c29_i67++) {
    c29_u[c29_i67] = c29_b_inData[c29_i67];
  }

  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", c29_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c29_mxArrayOutData, c29_y, FALSE);
  return c29_mxArrayOutData;
}

static const mxArray *c29_d_sf_marshallOut(void *chartInstanceVoid, void
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

static const mxArray *c29_e_sf_marshallOut(void *chartInstanceVoid, void
  *c29_inData)
{
  const mxArray *c29_mxArrayOutData = NULL;
  int32_T c29_i68;
  real_T c29_b_inData[2];
  int32_T c29_i69;
  real_T c29_u[2];
  const mxArray *c29_y = NULL;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_mxArrayOutData = NULL;
  for (c29_i68 = 0; c29_i68 < 2; c29_i68++) {
    c29_b_inData[c29_i68] = (*(real_T (*)[2])c29_inData)[c29_i68];
  }

  for (c29_i69 = 0; c29_i69 < 2; c29_i69++) {
    c29_u[c29_i69] = c29_b_inData[c29_i69];
  }

  c29_y = NULL;
  sf_mex_assign(&c29_y, sf_mex_create("y", c29_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c29_mxArrayOutData, c29_y, FALSE);
  return c29_mxArrayOutData;
}

static real_T c29_e_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId)
{
  real_T c29_y;
  real_T c29_d1;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), &c29_d1, 1, 0, 0U, 0, 0U, 0);
  c29_y = c29_d1;
  sf_mex_destroy(&c29_u);
  return c29_y;
}

static void c29_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c29_mxArrayInData, const char_T *c29_varName, void *c29_outData)
{
  const mxArray *c29_b_T_est;
  const char_T *c29_identifier;
  emlrtMsgIdentifier c29_thisId;
  real_T c29_y;
  SFc29_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc29_simulationInstanceStruct *)chartInstanceVoid;
  c29_b_T_est = sf_mex_dup(c29_mxArrayInData);
  c29_identifier = c29_varName;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_y = c29_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_b_T_est),
    &c29_thisId);
  sf_mex_destroy(&c29_b_T_est);
  *(real_T *)c29_outData = c29_y;
  sf_mex_destroy(&c29_mxArrayInData);
}

const mxArray *sf_c29_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c29_nameCaptureInfo = NULL;
  c29_nameCaptureInfo = NULL;
  sf_mex_assign(&c29_nameCaptureInfo, sf_mex_createstruct("structure", 2, 14, 1),
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
  const mxArray *c29_rhs11 = NULL;
  const mxArray *c29_lhs11 = NULL;
  const mxArray *c29_rhs12 = NULL;
  const mxArray *c29_lhs12 = NULL;
  const mxArray *c29_rhs13 = NULL;
  const mxArray *c29_lhs13 = NULL;
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
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
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("mtimes"), "name", "name", 2);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 3);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
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
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("sin"), "name", "name", 4);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 5);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 6);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
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
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  8);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
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
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 9);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 10);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("mtimes"), "name", "name", 10);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
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
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 11);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 11);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c29_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c29_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 13);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c29_info, c29_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c29_info, c29_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c29_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c29_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c29_info, sf_mex_duplicatearraysafe(&c29_lhs13), "lhs", "lhs",
                  13);
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
  sf_mex_destroy(&c29_rhs11);
  sf_mex_destroy(&c29_lhs11);
  sf_mex_destroy(&c29_rhs12);
  sf_mex_destroy(&c29_lhs12);
  sf_mex_destroy(&c29_rhs13);
  sf_mex_destroy(&c29_lhs13);
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

static void c29_eml_scalar_eg(SFc29_simulationInstanceStruct *chartInstance)
{
}

static void c29_b_eml_scalar_eg(SFc29_simulationInstanceStruct *chartInstance)
{
}

static void c29_eml_xgemm(SFc29_simulationInstanceStruct *chartInstance, real_T
  c29_A[81], real_T c29_B[81], real_T c29_C[81], real_T c29_b_C[81])
{
  int32_T c29_i70;
  int32_T c29_i71;
  real_T c29_b_A[81];
  int32_T c29_i72;
  real_T c29_b_B[81];
  for (c29_i70 = 0; c29_i70 < 81; c29_i70++) {
    c29_b_C[c29_i70] = c29_C[c29_i70];
  }

  for (c29_i71 = 0; c29_i71 < 81; c29_i71++) {
    c29_b_A[c29_i71] = c29_A[c29_i71];
  }

  for (c29_i72 = 0; c29_i72 < 81; c29_i72++) {
    c29_b_B[c29_i72] = c29_B[c29_i72];
  }

  c29_b_eml_xgemm(chartInstance, c29_b_A, c29_b_B, c29_b_C);
}

static const mxArray *c29_f_sf_marshallOut(void *chartInstanceVoid, void
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

static int32_T c29_f_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId)
{
  int32_T c29_y;
  int32_T c29_i73;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), &c29_i73, 1, 6, 0U, 0, 0U, 0);
  c29_y = c29_i73;
  sf_mex_destroy(&c29_u);
  return c29_y;
}

static void c29_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
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
  c29_y = c29_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c29_b_sfEvent),
    &c29_thisId);
  sf_mex_destroy(&c29_b_sfEvent);
  *(int32_T *)c29_outData = c29_y;
  sf_mex_destroy(&c29_mxArrayInData);
}

static uint8_T c29_g_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_b_is_active_c29_simulation, const char_T
  *c29_identifier)
{
  uint8_T c29_y;
  emlrtMsgIdentifier c29_thisId;
  c29_thisId.fIdentifier = c29_identifier;
  c29_thisId.fParent = NULL;
  c29_y = c29_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c29_b_is_active_c29_simulation), &c29_thisId);
  sf_mex_destroy(&c29_b_is_active_c29_simulation);
  return c29_y;
}

static uint8_T c29_h_emlrt_marshallIn(SFc29_simulationInstanceStruct
  *chartInstance, const mxArray *c29_u, const emlrtMsgIdentifier *c29_parentId)
{
  uint8_T c29_y;
  uint8_T c29_u0;
  sf_mex_import(c29_parentId, sf_mex_dup(c29_u), &c29_u0, 1, 3, 0U, 0, 0U, 0);
  c29_y = c29_u0;
  sf_mex_destroy(&c29_u);
  return c29_y;
}

static void c29_b_eml_xgemm(SFc29_simulationInstanceStruct *chartInstance,
  real_T c29_A[81], real_T c29_B[81], real_T c29_C[81])
{
  real_T c29_alpha1;
  real_T c29_beta1;
  char_T c29_TRANSB;
  char_T c29_TRANSA;
  ptrdiff_t c29_m_t;
  ptrdiff_t c29_n_t;
  ptrdiff_t c29_k_t;
  ptrdiff_t c29_lda_t;
  ptrdiff_t c29_ldb_t;
  ptrdiff_t c29_ldc_t;
  double * c29_alpha1_t;
  double * c29_Aia0_t;
  double * c29_Bib0_t;
  double * c29_beta1_t;
  double * c29_Cic0_t;
  c29_alpha1 = 1.0;
  c29_beta1 = 0.0;
  c29_TRANSB = 'N';
  c29_TRANSA = 'N';
  c29_m_t = (ptrdiff_t)(9);
  c29_n_t = (ptrdiff_t)(9);
  c29_k_t = (ptrdiff_t)(9);
  c29_lda_t = (ptrdiff_t)(9);
  c29_ldb_t = (ptrdiff_t)(9);
  c29_ldc_t = (ptrdiff_t)(9);
  c29_alpha1_t = (double *)(&c29_alpha1);
  c29_Aia0_t = (double *)(&c29_A[0]);
  c29_Bib0_t = (double *)(&c29_B[0]);
  c29_beta1_t = (double *)(&c29_beta1);
  c29_Cic0_t = (double *)(&c29_C[0]);
  dgemm(&c29_TRANSA, &c29_TRANSB, &c29_m_t, &c29_n_t, &c29_k_t, c29_alpha1_t,
        c29_Aia0_t, &c29_lda_t, c29_Bib0_t, &c29_ldb_t, c29_beta1_t, c29_Cic0_t,
        &c29_ldc_t);
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
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4111358928U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(201005745U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2882677057U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3352912400U);
}

mxArray *sf_c29_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("uG2trjyup9dMsOszMva5p");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(9);
      pr[1] = (double)(9);
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
      pr[0] = (double)(9);
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
      pr[0] = (double)(9);
      pr[1] = (double)(9);
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
      pr[0] = (double)(3);
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

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(9);
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
      pr[0] = (double)(9);
      pr[1] = (double)(9);
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
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[8],T\"Sigma_predict\",},{M[1],M[5],T\"State_predict\",},{M[8],M[0],T\"is_active_c29_simulation\",}}"
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
           9,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"Sigma_km1");
          _SFD_SET_DATA_PROPS(1,1,1,0,"State_km1");
          _SFD_SET_DATA_PROPS(2,10,0,0,"T_est");
          _SFD_SET_DATA_PROPS(3,1,1,0,"Vc");
          _SFD_SET_DATA_PROPS(4,1,1,0,"Q");
          _SFD_SET_DATA_PROPS(5,1,1,0,"Psi");
          _SFD_SET_DATA_PROPS(6,1,1,0,"Usens");
          _SFD_SET_DATA_PROPS(7,2,0,1,"State_predict");
          _SFD_SET_DATA_PROPS(8,2,0,1,"Sigma_predict");
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
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,1646);

        {
          unsigned int dimVector[2];
          dimVector[0]= 9;
          dimVector[1]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c29_d_sf_marshallOut,(MexInFcnForType)
          c29_c_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 9;
          dimVector[1]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c29_d_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_b_sf_marshallOut,(MexInFcnForType)
            c29_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 9;
          dimVector[1]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c29_sf_marshallOut,(MexInFcnForType)
            c29_sf_marshallIn);
        }

        {
          real_T *c29_Psi;
          real_T (*c29_Sigma_km1)[81];
          real_T (*c29_State_km1)[9];
          real_T (*c29_Vc)[2];
          real_T (*c29_Q)[81];
          real_T (*c29_Usens)[3];
          real_T (*c29_State_predict)[9];
          real_T (*c29_Sigma_predict)[81];
          c29_Sigma_predict = (real_T (*)[81])ssGetOutputPortSignal
            (chartInstance->S, 2);
          c29_State_predict = (real_T (*)[9])ssGetOutputPortSignal
            (chartInstance->S, 1);
          c29_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
          c29_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c29_Q = (real_T (*)[81])ssGetInputPortSignal(chartInstance->S, 3);
          c29_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
          c29_State_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S,
            1);
          c29_Sigma_km1 = (real_T (*)[81])ssGetInputPortSignal(chartInstance->S,
            0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c29_Sigma_km1);
          _SFD_SET_DATA_VALUE_PTR(1U, *c29_State_km1);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c29_T_est);
          _SFD_SET_DATA_VALUE_PTR(3U, *c29_Vc);
          _SFD_SET_DATA_VALUE_PTR(4U, *c29_Q);
          _SFD_SET_DATA_VALUE_PTR(5U, c29_Psi);
          _SFD_SET_DATA_VALUE_PTR(6U, *c29_Usens);
          _SFD_SET_DATA_VALUE_PTR(7U, *c29_State_predict);
          _SFD_SET_DATA_VALUE_PTR(8U, *c29_Sigma_predict);
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
  return "VRvN14rWB8XykPyyzWPveD";
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
  /* Actual parameters from chart:
     T_est
   */
  const char_T *rtParamNames[] = { "T_est" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for T_est*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);
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
  ssSetChecksum0(S,(1782325652U));
  ssSetChecksum1(S,(1991135352U));
  ssSetChecksum2(S,(2774158480U));
  ssSetChecksum3(S,(739624091U));
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
