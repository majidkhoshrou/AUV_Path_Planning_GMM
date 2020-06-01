/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c19_simulation.h"
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
static const char * c19_debug_family_names[25] = { "d_u", "d_v", "d_r", "du",
  "dv", "dr", "nargin", "nargout", "u", "v", "r", "Fs", "Fp", "m_u", "m_v",
  "m_uv", "J", "l", "Xu", "Xuu", "Yv", "Yvv", "Nr", "Nrr", "dU" };

/* Function Declarations */
static void initialize_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance);
static void enable_c19_simulation(SFc19_simulationInstanceStruct *chartInstance);
static void disable_c19_simulation(SFc19_simulationInstanceStruct *chartInstance);
static void c19_update_debugger_state_c19_simulation
  (SFc19_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c19_simulation
  (SFc19_simulationInstanceStruct *chartInstance);
static void set_sim_state_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_st);
static void finalize_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance);
static void sf_c19_simulation(SFc19_simulationInstanceStruct *chartInstance);
static void c19_chartstep_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc19_simulation(SFc19_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c19_machineNumber, uint32_T
  c19_chartNumber);
static const mxArray *c19_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData);
static void c19_emlrt_marshallIn(SFc19_simulationInstanceStruct *chartInstance,
  const mxArray *c19_dU, const char_T *c19_identifier, real_T c19_y[3]);
static void c19_b_emlrt_marshallIn(SFc19_simulationInstanceStruct *chartInstance,
  const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId, real_T c19_y[3]);
static void c19_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData);
static const mxArray *c19_b_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData);
static real_T c19_c_emlrt_marshallIn(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId);
static void c19_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData);
static void c19_info_helper(const mxArray **c19_info);
static const mxArray *c19_emlrt_marshallOut(char * c19_u);
static const mxArray *c19_b_emlrt_marshallOut(uint32_T c19_u);
static const mxArray *c19_c_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData);
static int32_T c19_d_emlrt_marshallIn(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId);
static void c19_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData);
static uint8_T c19_e_emlrt_marshallIn(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_b_is_active_c19_simulation, const char_T
  *c19_identifier);
static uint8_T c19_f_emlrt_marshallIn(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId);
static void init_dsm_address_info(SFc19_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c19_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c19_is_active_c19_simulation = 0U;
}

static void initialize_params_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance)
{
  real_T c19_d0;
  real_T c19_d1;
  real_T c19_d2;
  real_T c19_d3;
  real_T c19_d4;
  real_T c19_d5;
  real_T c19_d6;
  real_T c19_d7;
  real_T c19_d8;
  real_T c19_d9;
  real_T c19_d10;
  sf_set_error_prefix_string(
    "Error evaluating data 'm_u' in the parent workspace.\n");
  sf_mex_import_named("m_u", sf_mex_get_sfun_param(chartInstance->S, 8, 0),
                      &c19_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_m_u = c19_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'm_v' in the parent workspace.\n");
  sf_mex_import_named("m_v", sf_mex_get_sfun_param(chartInstance->S, 10, 0),
                      &c19_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_m_v = c19_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'm_uv' in the parent workspace.\n");
  sf_mex_import_named("m_uv", sf_mex_get_sfun_param(chartInstance->S, 9, 0),
                      &c19_d2, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_m_uv = c19_d2;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'J' in the parent workspace.\n");
  sf_mex_import_named("J", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c19_d3, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_J = c19_d3;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'l' in the parent workspace.\n");
  sf_mex_import_named("l", sf_mex_get_sfun_param(chartInstance->S, 7, 0),
                      &c19_d4, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_l = c19_d4;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Xu' in the parent workspace.\n");
  sf_mex_import_named("Xu", sf_mex_get_sfun_param(chartInstance->S, 3, 0),
                      &c19_d5, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_Xu = c19_d5;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Xuu' in the parent workspace.\n");
  sf_mex_import_named("Xuu", sf_mex_get_sfun_param(chartInstance->S, 4, 0),
                      &c19_d6, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_Xuu = c19_d6;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Yv' in the parent workspace.\n");
  sf_mex_import_named("Yv", sf_mex_get_sfun_param(chartInstance->S, 5, 0),
                      &c19_d7, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_Yv = c19_d7;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Yvv' in the parent workspace.\n");
  sf_mex_import_named("Yvv", sf_mex_get_sfun_param(chartInstance->S, 6, 0),
                      &c19_d8, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_Yvv = c19_d8;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Nr' in the parent workspace.\n");
  sf_mex_import_named("Nr", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      &c19_d9, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_Nr = c19_d9;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Nrr' in the parent workspace.\n");
  sf_mex_import_named("Nrr", sf_mex_get_sfun_param(chartInstance->S, 2, 0),
                      &c19_d10, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c19_Nrr = c19_d10;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c19_simulation(SFc19_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c19_simulation(SFc19_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c19_update_debugger_state_c19_simulation
  (SFc19_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c19_simulation
  (SFc19_simulationInstanceStruct *chartInstance)
{
  const mxArray *c19_st;
  const mxArray *c19_y = NULL;
  int32_T c19_i0;
  real_T c19_u[3];
  const mxArray *c19_b_y = NULL;
  uint8_T c19_hoistedGlobal;
  uint8_T c19_b_u;
  const mxArray *c19_c_y = NULL;
  real_T (*c19_dU)[3];
  c19_dU = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c19_st = NULL;
  c19_st = NULL;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_createcellarray(2), FALSE);
  for (c19_i0 = 0; c19_i0 < 3; c19_i0++) {
    c19_u[c19_i0] = (*c19_dU)[c19_i0];
  }

  c19_b_y = NULL;
  sf_mex_assign(&c19_b_y, sf_mex_create("y", c19_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_setcell(c19_y, 0, c19_b_y);
  c19_hoistedGlobal = chartInstance->c19_is_active_c19_simulation;
  c19_b_u = c19_hoistedGlobal;
  c19_c_y = NULL;
  sf_mex_assign(&c19_c_y, sf_mex_create("y", &c19_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c19_y, 1, c19_c_y);
  sf_mex_assign(&c19_st, c19_y, FALSE);
  return c19_st;
}

static void set_sim_state_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_st)
{
  const mxArray *c19_u;
  real_T c19_dv0[3];
  int32_T c19_i1;
  real_T (*c19_dU)[3];
  c19_dU = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c19_doneDoubleBufferReInit = TRUE;
  c19_u = sf_mex_dup(c19_st);
  c19_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c19_u, 0)), "dU",
                       c19_dv0);
  for (c19_i1 = 0; c19_i1 < 3; c19_i1++) {
    (*c19_dU)[c19_i1] = c19_dv0[c19_i1];
  }

  chartInstance->c19_is_active_c19_simulation = c19_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c19_u, 1)),
     "is_active_c19_simulation");
  sf_mex_destroy(&c19_u);
  c19_update_debugger_state_c19_simulation(chartInstance);
  sf_mex_destroy(&c19_st);
}

static void finalize_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c19_simulation(SFc19_simulationInstanceStruct *chartInstance)
{
  int32_T c19_i2;
  real_T *c19_u;
  real_T *c19_v;
  real_T *c19_r;
  real_T *c19_Fs;
  real_T *c19_Fp;
  real_T (*c19_dU)[3];
  c19_Fp = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c19_Fs = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c19_r = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c19_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c19_dU = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c19_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 18U, chartInstance->c19_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c19_u, 0U);
  for (c19_i2 = 0; c19_i2 < 3; c19_i2++) {
    _SFD_DATA_RANGE_CHECK((*c19_dU)[c19_i2], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c19_v, 2U);
  _SFD_DATA_RANGE_CHECK(*c19_r, 3U);
  _SFD_DATA_RANGE_CHECK(*c19_Fs, 4U);
  _SFD_DATA_RANGE_CHECK(*c19_Fp, 5U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_m_u, 6U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_m_v, 7U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_m_uv, 8U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_J, 9U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_l, 10U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_Xu, 11U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_Xuu, 12U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_Yv, 13U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_Yvv, 14U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_Nr, 15U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c19_Nrr, 16U);
  chartInstance->c19_sfEvent = CALL_EVENT;
  c19_chartstep_c19_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c19_chartstep_c19_simulation(SFc19_simulationInstanceStruct
  *chartInstance)
{
  real_T c19_hoistedGlobal;
  real_T c19_b_hoistedGlobal;
  real_T c19_c_hoistedGlobal;
  real_T c19_d_hoistedGlobal;
  real_T c19_e_hoistedGlobal;
  real_T c19_f_hoistedGlobal;
  real_T c19_g_hoistedGlobal;
  real_T c19_h_hoistedGlobal;
  real_T c19_i_hoistedGlobal;
  real_T c19_j_hoistedGlobal;
  real_T c19_k_hoistedGlobal;
  real_T c19_l_hoistedGlobal;
  real_T c19_m_hoistedGlobal;
  real_T c19_n_hoistedGlobal;
  real_T c19_o_hoistedGlobal;
  real_T c19_p_hoistedGlobal;
  real_T c19_u;
  real_T c19_v;
  real_T c19_r;
  real_T c19_Fs;
  real_T c19_Fp;
  real_T c19_b_m_u;
  real_T c19_b_m_v;
  real_T c19_b_m_uv;
  real_T c19_b_J;
  real_T c19_b_l;
  real_T c19_b_Xu;
  real_T c19_b_Xuu;
  real_T c19_b_Yv;
  real_T c19_b_Yvv;
  real_T c19_b_Nr;
  real_T c19_b_Nrr;
  uint32_T c19_debug_family_var_map[25];
  real_T c19_d_u;
  real_T c19_d_v;
  real_T c19_d_r;
  real_T c19_du;
  real_T c19_dv;
  real_T c19_dr;
  real_T c19_nargin = 16.0;
  real_T c19_nargout = 1.0;
  real_T c19_dU[3];
  real_T c19_x;
  real_T c19_b_x;
  real_T c19_y;
  real_T c19_a;
  real_T c19_b;
  real_T c19_b_y;
  real_T c19_c_x;
  real_T c19_d_x;
  real_T c19_c_y;
  real_T c19_b_a;
  real_T c19_b_b;
  real_T c19_d_y;
  real_T c19_e_x;
  real_T c19_f_x;
  real_T c19_e_y;
  real_T c19_c_a;
  real_T c19_c_b;
  real_T c19_f_y;
  real_T c19_d_a;
  real_T c19_d_b;
  real_T c19_g_y;
  real_T c19_e_a;
  real_T c19_e_b;
  real_T c19_h_y;
  real_T c19_f_a;
  real_T c19_f_b;
  real_T c19_i_y;
  real_T c19_A;
  real_T c19_B;
  real_T c19_g_x;
  real_T c19_j_y;
  real_T c19_h_x;
  real_T c19_k_y;
  real_T c19_g_a;
  real_T c19_g_b;
  real_T c19_l_y;
  real_T c19_h_a;
  real_T c19_h_b;
  real_T c19_m_y;
  real_T c19_i_a;
  real_T c19_i_b;
  real_T c19_n_y;
  real_T c19_b_A;
  real_T c19_b_B;
  real_T c19_i_x;
  real_T c19_o_y;
  real_T c19_j_x;
  real_T c19_p_y;
  real_T c19_j_a;
  real_T c19_j_b;
  real_T c19_q_y;
  real_T c19_k_a;
  real_T c19_k_b;
  real_T c19_r_y;
  real_T c19_l_a;
  real_T c19_l_b;
  real_T c19_s_y;
  real_T c19_m_a;
  real_T c19_m_b;
  real_T c19_t_y;
  real_T c19_c_A;
  real_T c19_c_B;
  real_T c19_k_x;
  real_T c19_u_y;
  real_T c19_l_x;
  real_T c19_v_y;
  real_T c19_b_du[3];
  int32_T c19_i3;
  int32_T c19_i4;
  real_T *c19_b_Fp;
  real_T *c19_b_Fs;
  real_T *c19_b_r;
  real_T *c19_b_v;
  real_T *c19_b_u;
  real_T (*c19_b_dU)[3];
  c19_b_Fp = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c19_b_Fs = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c19_b_r = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c19_b_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c19_b_dU = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c19_b_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 18U, chartInstance->c19_sfEvent);
  c19_hoistedGlobal = *c19_b_u;
  c19_b_hoistedGlobal = *c19_b_v;
  c19_c_hoistedGlobal = *c19_b_r;
  c19_d_hoistedGlobal = *c19_b_Fs;
  c19_e_hoistedGlobal = *c19_b_Fp;
  c19_f_hoistedGlobal = chartInstance->c19_m_u;
  c19_g_hoistedGlobal = chartInstance->c19_m_v;
  c19_h_hoistedGlobal = chartInstance->c19_m_uv;
  c19_i_hoistedGlobal = chartInstance->c19_J;
  c19_j_hoistedGlobal = chartInstance->c19_l;
  c19_k_hoistedGlobal = chartInstance->c19_Xu;
  c19_l_hoistedGlobal = chartInstance->c19_Xuu;
  c19_m_hoistedGlobal = chartInstance->c19_Yv;
  c19_n_hoistedGlobal = chartInstance->c19_Yvv;
  c19_o_hoistedGlobal = chartInstance->c19_Nr;
  c19_p_hoistedGlobal = chartInstance->c19_Nrr;
  c19_u = c19_hoistedGlobal;
  c19_v = c19_b_hoistedGlobal;
  c19_r = c19_c_hoistedGlobal;
  c19_Fs = c19_d_hoistedGlobal;
  c19_Fp = c19_e_hoistedGlobal;
  c19_b_m_u = c19_f_hoistedGlobal;
  c19_b_m_v = c19_g_hoistedGlobal;
  c19_b_m_uv = c19_h_hoistedGlobal;
  c19_b_J = c19_i_hoistedGlobal;
  c19_b_l = c19_j_hoistedGlobal;
  c19_b_Xu = c19_k_hoistedGlobal;
  c19_b_Xuu = c19_l_hoistedGlobal;
  c19_b_Yv = c19_m_hoistedGlobal;
  c19_b_Yvv = c19_n_hoistedGlobal;
  c19_b_Nr = c19_o_hoistedGlobal;
  c19_b_Nrr = c19_p_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 25U, 25U, c19_debug_family_names,
    c19_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_d_u, 0U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_d_v, 1U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_d_r, 2U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_du, 3U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_dv, 4U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_dr, 5U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_nargin, 6U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_nargout, 7U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c19_u, 8U, c19_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c19_v, 9U, c19_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c19_r, 10U, c19_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c19_Fs, 11U, c19_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c19_Fp, 12U, c19_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_m_u, 13U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_m_v, 14U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_m_uv, 15U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_J, 16U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_l, 17U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_Xu, 18U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_Xuu, 19U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_Yv, 20U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_Yvv, 21U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_Nr, 22U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c19_b_Nrr, 23U, c19_b_sf_marshallOut,
    c19_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c19_dU, 24U, c19_sf_marshallOut,
    c19_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c19_sfEvent, 3);
  c19_x = c19_u;
  c19_b_x = c19_x;
  c19_y = muDoubleScalarAbs(c19_b_x);
  c19_a = c19_b_Xuu;
  c19_b = c19_y;
  c19_b_y = c19_a * c19_b;
  c19_d_u = -c19_b_Xu - c19_b_y;
  _SFD_EML_CALL(0U, chartInstance->c19_sfEvent, 4);
  c19_c_x = c19_v;
  c19_d_x = c19_c_x;
  c19_c_y = muDoubleScalarAbs(c19_d_x);
  c19_b_a = c19_b_Yvv;
  c19_b_b = c19_c_y;
  c19_d_y = c19_b_a * c19_b_b;
  c19_d_v = -c19_b_Yv - c19_d_y;
  _SFD_EML_CALL(0U, chartInstance->c19_sfEvent, 5);
  c19_e_x = c19_r;
  c19_f_x = c19_e_x;
  c19_e_y = muDoubleScalarAbs(c19_f_x);
  c19_c_a = c19_b_Nrr;
  c19_c_b = c19_e_y;
  c19_f_y = c19_c_a * c19_c_b;
  c19_d_r = -c19_b_Nr - c19_f_y;
  _SFD_EML_CALL(0U, chartInstance->c19_sfEvent, 7);
  c19_d_a = c19_b_m_v;
  c19_d_b = c19_v;
  c19_g_y = c19_d_a * c19_d_b;
  c19_e_a = c19_g_y;
  c19_e_b = c19_r;
  c19_h_y = c19_e_a * c19_e_b;
  c19_f_a = c19_d_u;
  c19_f_b = c19_u;
  c19_i_y = c19_f_a * c19_f_b;
  c19_A = ((c19_h_y - c19_i_y) + c19_Fs) + c19_Fp;
  c19_B = c19_b_m_u;
  c19_g_x = c19_A;
  c19_j_y = c19_B;
  c19_h_x = c19_g_x;
  c19_k_y = c19_j_y;
  c19_du = c19_h_x / c19_k_y;
  _SFD_EML_CALL(0U, chartInstance->c19_sfEvent, 8);
  c19_g_a = -c19_b_m_u;
  c19_g_b = c19_u;
  c19_l_y = c19_g_a * c19_g_b;
  c19_h_a = c19_l_y;
  c19_h_b = c19_r;
  c19_m_y = c19_h_a * c19_h_b;
  c19_i_a = c19_d_v;
  c19_i_b = c19_v;
  c19_n_y = c19_i_a * c19_i_b;
  c19_b_A = c19_m_y - c19_n_y;
  c19_b_B = c19_b_m_v;
  c19_i_x = c19_b_A;
  c19_o_y = c19_b_B;
  c19_j_x = c19_i_x;
  c19_p_y = c19_o_y;
  c19_dv = c19_j_x / c19_p_y;
  _SFD_EML_CALL(0U, chartInstance->c19_sfEvent, 9);
  c19_j_a = c19_b_m_uv;
  c19_j_b = c19_u;
  c19_q_y = c19_j_a * c19_j_b;
  c19_k_a = c19_q_y;
  c19_k_b = c19_v;
  c19_r_y = c19_k_a * c19_k_b;
  c19_l_a = c19_d_r;
  c19_l_b = c19_r;
  c19_s_y = c19_l_a * c19_l_b;
  c19_m_a = c19_b_l;
  c19_m_b = c19_Fs - c19_Fp;
  c19_t_y = c19_m_a * c19_m_b;
  c19_c_A = (c19_r_y - c19_s_y) + c19_t_y;
  c19_c_B = c19_b_J;
  c19_k_x = c19_c_A;
  c19_u_y = c19_c_B;
  c19_l_x = c19_k_x;
  c19_v_y = c19_u_y;
  c19_dr = c19_l_x / c19_v_y;
  _SFD_EML_CALL(0U, chartInstance->c19_sfEvent, 11);
  c19_b_du[0] = c19_du;
  c19_b_du[1] = c19_dv;
  c19_b_du[2] = c19_dr;
  for (c19_i3 = 0; c19_i3 < 3; c19_i3++) {
    c19_dU[c19_i3] = c19_b_du[c19_i3];
  }

  _SFD_EML_CALL(0U, chartInstance->c19_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  for (c19_i4 = 0; c19_i4 < 3; c19_i4++) {
    (*c19_b_dU)[c19_i4] = c19_dU[c19_i4];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 18U, chartInstance->c19_sfEvent);
}

static void initSimStructsc19_simulation(SFc19_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c19_machineNumber, uint32_T
  c19_chartNumber)
{
}

static const mxArray *c19_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData)
{
  const mxArray *c19_mxArrayOutData = NULL;
  int32_T c19_i5;
  real_T c19_b_inData[3];
  int32_T c19_i6;
  real_T c19_u[3];
  const mxArray *c19_y = NULL;
  SFc19_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc19_simulationInstanceStruct *)chartInstanceVoid;
  c19_mxArrayOutData = NULL;
  for (c19_i5 = 0; c19_i5 < 3; c19_i5++) {
    c19_b_inData[c19_i5] = (*(real_T (*)[3])c19_inData)[c19_i5];
  }

  for (c19_i6 = 0; c19_i6 < 3; c19_i6++) {
    c19_u[c19_i6] = c19_b_inData[c19_i6];
  }

  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", c19_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c19_mxArrayOutData, c19_y, FALSE);
  return c19_mxArrayOutData;
}

static void c19_emlrt_marshallIn(SFc19_simulationInstanceStruct *chartInstance,
  const mxArray *c19_dU, const char_T *c19_identifier, real_T c19_y[3])
{
  emlrtMsgIdentifier c19_thisId;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c19_dU), &c19_thisId, c19_y);
  sf_mex_destroy(&c19_dU);
}

static void c19_b_emlrt_marshallIn(SFc19_simulationInstanceStruct *chartInstance,
  const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId, real_T c19_y[3])
{
  real_T c19_dv1[3];
  int32_T c19_i7;
  sf_mex_import(c19_parentId, sf_mex_dup(c19_u), c19_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c19_i7 = 0; c19_i7 < 3; c19_i7++) {
    c19_y[c19_i7] = c19_dv1[c19_i7];
  }

  sf_mex_destroy(&c19_u);
}

static void c19_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData)
{
  const mxArray *c19_dU;
  const char_T *c19_identifier;
  emlrtMsgIdentifier c19_thisId;
  real_T c19_y[3];
  int32_T c19_i8;
  SFc19_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc19_simulationInstanceStruct *)chartInstanceVoid;
  c19_dU = sf_mex_dup(c19_mxArrayInData);
  c19_identifier = c19_varName;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c19_dU), &c19_thisId, c19_y);
  sf_mex_destroy(&c19_dU);
  for (c19_i8 = 0; c19_i8 < 3; c19_i8++) {
    (*(real_T (*)[3])c19_outData)[c19_i8] = c19_y[c19_i8];
  }

  sf_mex_destroy(&c19_mxArrayInData);
}

static const mxArray *c19_b_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData)
{
  const mxArray *c19_mxArrayOutData = NULL;
  real_T c19_u;
  const mxArray *c19_y = NULL;
  SFc19_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc19_simulationInstanceStruct *)chartInstanceVoid;
  c19_mxArrayOutData = NULL;
  c19_u = *(real_T *)c19_inData;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", &c19_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c19_mxArrayOutData, c19_y, FALSE);
  return c19_mxArrayOutData;
}

static real_T c19_c_emlrt_marshallIn(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId)
{
  real_T c19_y;
  real_T c19_d11;
  sf_mex_import(c19_parentId, sf_mex_dup(c19_u), &c19_d11, 1, 0, 0U, 0, 0U, 0);
  c19_y = c19_d11;
  sf_mex_destroy(&c19_u);
  return c19_y;
}

static void c19_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData)
{
  const mxArray *c19_b_Nrr;
  const char_T *c19_identifier;
  emlrtMsgIdentifier c19_thisId;
  real_T c19_y;
  SFc19_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc19_simulationInstanceStruct *)chartInstanceVoid;
  c19_b_Nrr = sf_mex_dup(c19_mxArrayInData);
  c19_identifier = c19_varName;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_y = c19_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c19_b_Nrr),
    &c19_thisId);
  sf_mex_destroy(&c19_b_Nrr);
  *(real_T *)c19_outData = c19_y;
  sf_mex_destroy(&c19_mxArrayInData);
}

const mxArray *sf_c19_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c19_nameCaptureInfo = NULL;
  c19_nameCaptureInfo = NULL;
  sf_mex_assign(&c19_nameCaptureInfo, sf_mex_createstruct("structure", 2, 10, 1),
                FALSE);
  c19_info_helper(&c19_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c19_nameCaptureInfo);
  return c19_nameCaptureInfo;
}

static void c19_info_helper(const mxArray **c19_info)
{
  const mxArray *c19_rhs0 = NULL;
  const mxArray *c19_lhs0 = NULL;
  const mxArray *c19_rhs1 = NULL;
  const mxArray *c19_lhs1 = NULL;
  const mxArray *c19_rhs2 = NULL;
  const mxArray *c19_lhs2 = NULL;
  const mxArray *c19_rhs3 = NULL;
  const mxArray *c19_lhs3 = NULL;
  const mxArray *c19_rhs4 = NULL;
  const mxArray *c19_lhs4 = NULL;
  const mxArray *c19_rhs5 = NULL;
  const mxArray *c19_lhs5 = NULL;
  const mxArray *c19_rhs6 = NULL;
  const mxArray *c19_lhs6 = NULL;
  const mxArray *c19_rhs7 = NULL;
  const mxArray *c19_lhs7 = NULL;
  const mxArray *c19_rhs8 = NULL;
  const mxArray *c19_lhs8 = NULL;
  const mxArray *c19_rhs9 = NULL;
  const mxArray *c19_lhs9 = NULL;
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("abs"), "name", "name", 0);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c19_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c19_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 2);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1286822312U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c19_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(""), "context", "context", 3);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("mtimes"), "name", "name", 3);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c19_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 4);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c19_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(""), "context", "context", 5);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("mrdivide"), "name", "name",
                  5);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c19_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 6);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("rdivide"), "name", "name", 6);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c19_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 7);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c19_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 8);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c19_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("eml_div"), "name", "name", 9);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c19_info, c19_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c19_info, c19_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c19_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c19_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c19_info, sf_mex_duplicatearraysafe(&c19_lhs9), "lhs", "lhs",
                  9);
  sf_mex_destroy(&c19_rhs0);
  sf_mex_destroy(&c19_lhs0);
  sf_mex_destroy(&c19_rhs1);
  sf_mex_destroy(&c19_lhs1);
  sf_mex_destroy(&c19_rhs2);
  sf_mex_destroy(&c19_lhs2);
  sf_mex_destroy(&c19_rhs3);
  sf_mex_destroy(&c19_lhs3);
  sf_mex_destroy(&c19_rhs4);
  sf_mex_destroy(&c19_lhs4);
  sf_mex_destroy(&c19_rhs5);
  sf_mex_destroy(&c19_lhs5);
  sf_mex_destroy(&c19_rhs6);
  sf_mex_destroy(&c19_lhs6);
  sf_mex_destroy(&c19_rhs7);
  sf_mex_destroy(&c19_lhs7);
  sf_mex_destroy(&c19_rhs8);
  sf_mex_destroy(&c19_lhs8);
  sf_mex_destroy(&c19_rhs9);
  sf_mex_destroy(&c19_lhs9);
}

static const mxArray *c19_emlrt_marshallOut(char * c19_u)
{
  const mxArray *c19_y = NULL;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", c19_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c19_u)), FALSE);
  return c19_y;
}

static const mxArray *c19_b_emlrt_marshallOut(uint32_T c19_u)
{
  const mxArray *c19_y = NULL;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", &c19_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c19_y;
}

static const mxArray *c19_c_sf_marshallOut(void *chartInstanceVoid, void
  *c19_inData)
{
  const mxArray *c19_mxArrayOutData = NULL;
  int32_T c19_u;
  const mxArray *c19_y = NULL;
  SFc19_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc19_simulationInstanceStruct *)chartInstanceVoid;
  c19_mxArrayOutData = NULL;
  c19_u = *(int32_T *)c19_inData;
  c19_y = NULL;
  sf_mex_assign(&c19_y, sf_mex_create("y", &c19_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c19_mxArrayOutData, c19_y, FALSE);
  return c19_mxArrayOutData;
}

static int32_T c19_d_emlrt_marshallIn(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId)
{
  int32_T c19_y;
  int32_T c19_i9;
  sf_mex_import(c19_parentId, sf_mex_dup(c19_u), &c19_i9, 1, 6, 0U, 0, 0U, 0);
  c19_y = c19_i9;
  sf_mex_destroy(&c19_u);
  return c19_y;
}

static void c19_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c19_mxArrayInData, const char_T *c19_varName, void *c19_outData)
{
  const mxArray *c19_b_sfEvent;
  const char_T *c19_identifier;
  emlrtMsgIdentifier c19_thisId;
  int32_T c19_y;
  SFc19_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc19_simulationInstanceStruct *)chartInstanceVoid;
  c19_b_sfEvent = sf_mex_dup(c19_mxArrayInData);
  c19_identifier = c19_varName;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_y = c19_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c19_b_sfEvent),
    &c19_thisId);
  sf_mex_destroy(&c19_b_sfEvent);
  *(int32_T *)c19_outData = c19_y;
  sf_mex_destroy(&c19_mxArrayInData);
}

static uint8_T c19_e_emlrt_marshallIn(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_b_is_active_c19_simulation, const char_T
  *c19_identifier)
{
  uint8_T c19_y;
  emlrtMsgIdentifier c19_thisId;
  c19_thisId.fIdentifier = c19_identifier;
  c19_thisId.fParent = NULL;
  c19_y = c19_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c19_b_is_active_c19_simulation), &c19_thisId);
  sf_mex_destroy(&c19_b_is_active_c19_simulation);
  return c19_y;
}

static uint8_T c19_f_emlrt_marshallIn(SFc19_simulationInstanceStruct
  *chartInstance, const mxArray *c19_u, const emlrtMsgIdentifier *c19_parentId)
{
  uint8_T c19_y;
  uint8_T c19_u0;
  sf_mex_import(c19_parentId, sf_mex_dup(c19_u), &c19_u0, 1, 3, 0U, 0, 0U, 0);
  c19_y = c19_u0;
  sf_mex_destroy(&c19_u);
  return c19_y;
}

static void init_dsm_address_info(SFc19_simulationInstanceStruct *chartInstance)
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

void sf_c19_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(585514460U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2558334882U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1301590099U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3651518508U);
}

mxArray *sf_c19_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("L1U9KWlF16TwgOJsxCOdKC");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,11,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,6,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,6,"type",mxType);
    }

    mxSetField(mxData,6,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,7,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,7,"type",mxType);
    }

    mxSetField(mxData,7,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,8,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,8,"type",mxType);
    }

    mxSetField(mxData,8,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,9,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,9,"type",mxType);
    }

    mxSetField(mxData,9,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,10,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,10,"type",mxType);
    }

    mxSetField(mxData,10,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
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

mxArray *sf_c19_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c19_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c19_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"dU\",},{M[8],M[0],T\"is_active_c19_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c19_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc19_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc19_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           19,
           1,
           1,
           17,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"u");
          _SFD_SET_DATA_PROPS(1,2,0,1,"dU");
          _SFD_SET_DATA_PROPS(2,1,1,0,"v");
          _SFD_SET_DATA_PROPS(3,1,1,0,"r");
          _SFD_SET_DATA_PROPS(4,1,1,0,"Fs");
          _SFD_SET_DATA_PROPS(5,1,1,0,"Fp");
          _SFD_SET_DATA_PROPS(6,10,0,0,"m_u");
          _SFD_SET_DATA_PROPS(7,10,0,0,"m_v");
          _SFD_SET_DATA_PROPS(8,10,0,0,"m_uv");
          _SFD_SET_DATA_PROPS(9,10,0,0,"J");
          _SFD_SET_DATA_PROPS(10,10,0,0,"l");
          _SFD_SET_DATA_PROPS(11,10,0,0,"Xu");
          _SFD_SET_DATA_PROPS(12,10,0,0,"Xuu");
          _SFD_SET_DATA_PROPS(13,10,0,0,"Yv");
          _SFD_SET_DATA_PROPS(14,10,0,0,"Yvv");
          _SFD_SET_DATA_PROPS(15,10,0,0,"Nr");
          _SFD_SET_DATA_PROPS(16,10,0,0,"Nrr");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,270);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c19_sf_marshallOut,(MexInFcnForType)
            c19_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(13,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(14,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(15,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(16,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c19_b_sf_marshallOut,(MexInFcnForType)
          c19_b_sf_marshallIn);

        {
          real_T *c19_u;
          real_T *c19_v;
          real_T *c19_r;
          real_T *c19_Fs;
          real_T *c19_Fp;
          real_T (*c19_dU)[3];
          c19_Fp = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c19_Fs = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c19_r = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c19_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c19_dU = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
          c19_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c19_u);
          _SFD_SET_DATA_VALUE_PTR(1U, *c19_dU);
          _SFD_SET_DATA_VALUE_PTR(2U, c19_v);
          _SFD_SET_DATA_VALUE_PTR(3U, c19_r);
          _SFD_SET_DATA_VALUE_PTR(4U, c19_Fs);
          _SFD_SET_DATA_VALUE_PTR(5U, c19_Fp);
          _SFD_SET_DATA_VALUE_PTR(6U, &chartInstance->c19_m_u);
          _SFD_SET_DATA_VALUE_PTR(7U, &chartInstance->c19_m_v);
          _SFD_SET_DATA_VALUE_PTR(8U, &chartInstance->c19_m_uv);
          _SFD_SET_DATA_VALUE_PTR(9U, &chartInstance->c19_J);
          _SFD_SET_DATA_VALUE_PTR(10U, &chartInstance->c19_l);
          _SFD_SET_DATA_VALUE_PTR(11U, &chartInstance->c19_Xu);
          _SFD_SET_DATA_VALUE_PTR(12U, &chartInstance->c19_Xuu);
          _SFD_SET_DATA_VALUE_PTR(13U, &chartInstance->c19_Yv);
          _SFD_SET_DATA_VALUE_PTR(14U, &chartInstance->c19_Yvv);
          _SFD_SET_DATA_VALUE_PTR(15U, &chartInstance->c19_Nr);
          _SFD_SET_DATA_VALUE_PTR(16U, &chartInstance->c19_Nrr);
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
  return "fBk9zgNY5ZyqKmoxPyrvRF";
}

static void sf_opaque_initialize_c19_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc19_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c19_simulation((SFc19_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c19_simulation((SFc19_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c19_simulation(void *chartInstanceVar)
{
  enable_c19_simulation((SFc19_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c19_simulation(void *chartInstanceVar)
{
  disable_c19_simulation((SFc19_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c19_simulation(void *chartInstanceVar)
{
  sf_c19_simulation((SFc19_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c19_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c19_simulation
    ((SFc19_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c19_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c19_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c19_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c19_simulation((SFc19_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c19_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c19_simulation(S);
}

static void sf_opaque_set_sim_state_c19_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c19_simulation(S, st);
}

static void sf_opaque_terminate_c19_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc19_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c19_simulation((SFc19_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc19_simulation((SFc19_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c19_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c19_simulation((SFc19_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c19_simulation(SimStruct *S)
{
  /* Actual parameters from chart:
     J Nr Nrr Xu Xuu Yv Yvv l m_u m_uv m_v
   */
  const char_T *rtParamNames[] = { "J", "Nr", "Nrr", "Xu", "Xuu", "Yv", "Yvv",
    "l", "m_u", "m_uv", "m_v" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for J*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for Nr*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);

  /* registration for Nrr*/
  ssRegDlgParamAsRunTimeParam(S, 2, 2, rtParamNames[2], SS_DOUBLE);

  /* registration for Xu*/
  ssRegDlgParamAsRunTimeParam(S, 3, 3, rtParamNames[3], SS_DOUBLE);

  /* registration for Xuu*/
  ssRegDlgParamAsRunTimeParam(S, 4, 4, rtParamNames[4], SS_DOUBLE);

  /* registration for Yv*/
  ssRegDlgParamAsRunTimeParam(S, 5, 5, rtParamNames[5], SS_DOUBLE);

  /* registration for Yvv*/
  ssRegDlgParamAsRunTimeParam(S, 6, 6, rtParamNames[6], SS_DOUBLE);

  /* registration for l*/
  ssRegDlgParamAsRunTimeParam(S, 7, 7, rtParamNames[7], SS_DOUBLE);

  /* registration for m_u*/
  ssRegDlgParamAsRunTimeParam(S, 8, 8, rtParamNames[8], SS_DOUBLE);

  /* registration for m_uv*/
  ssRegDlgParamAsRunTimeParam(S, 9, 9, rtParamNames[9], SS_DOUBLE);

  /* registration for m_v*/
  ssRegDlgParamAsRunTimeParam(S, 10, 10, rtParamNames[10], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      19);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,19,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,19,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,19);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,19,5);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,19,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 5; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,19);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(4257671422U));
  ssSetChecksum1(S,(507799627U));
  ssSetChecksum2(S,(1027883903U));
  ssSetChecksum3(S,(3805618979U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c19_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c19_simulation(SimStruct *S)
{
  SFc19_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc19_simulationInstanceStruct *)utMalloc(sizeof
    (SFc19_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc19_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c19_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c19_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c19_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c19_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c19_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c19_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c19_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c19_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c19_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c19_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c19_simulation;
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

void c19_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c19_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c19_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c19_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c19_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
