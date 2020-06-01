/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c9_simulation.h"
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
static const char * c9_debug_family_names[11] = { "D", "A", "L_d", "V_tilde",
  "nargin", "nargout", "Gam_est", "n", "k_xi", "i", "v_tilde" };

/* Function Declarations */
static void initialize_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance);
static void enable_c9_simulation(SFc9_simulationInstanceStruct *chartInstance);
static void disable_c9_simulation(SFc9_simulationInstanceStruct *chartInstance);
static void c9_update_debugger_state_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c9_simulation(SFc9_simulationInstanceStruct *
  chartInstance);
static void set_sim_state_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance, const mxArray *c9_st);
static void finalize_c9_simulation(SFc9_simulationInstanceStruct *chartInstance);
static void sf_c9_simulation(SFc9_simulationInstanceStruct *chartInstance);
static void c9_chartstep_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc9_simulation(SFc9_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c9_machineNumber, uint32_T
  c9_chartNumber);
static const mxArray *c9_sf_marshallOut(void *chartInstanceVoid, void *c9_inData);
static real_T c9_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_v_tilde, const char_T *c9_identifier);
static real_T c9_b_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId);
static void c9_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData);
static const mxArray *c9_b_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData);
static const mxArray *c9_c_sf_marshallOut(void *chartInstanceVoid, real_T
  c9_inData_data[20], int32_T c9_inData_sizes[1]);
static void c9_c_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, real_T c9_y_data
  [20], int32_T c9_y_sizes[1]);
static void c9_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, real_T c9_outData_data[20],
  int32_T c9_outData_sizes[1]);
static const mxArray *c9_d_sf_marshallOut(void *chartInstanceVoid, real_T
  c9_inData_data[400], int32_T c9_inData_sizes[2]);
static void c9_d_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, real_T c9_y_data
  [400], int32_T c9_y_sizes[2]);
static void c9_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, real_T c9_outData_data[400],
  int32_T c9_outData_sizes[2]);
static void c9_info_helper(const mxArray **c9_info);
static const mxArray *c9_emlrt_marshallOut(char * c9_u);
static const mxArray *c9_b_emlrt_marshallOut(uint32_T c9_u);
static void c9_b_info_helper(const mxArray **c9_info);
static void c9_c_info_helper(const mxArray **c9_info);
static void c9_d_info_helper(const mxArray **c9_info);
static void c9_diag(SFc9_simulationInstanceStruct *chartInstance, real_T
                    c9_v_data[20], int32_T c9_v_sizes[2], real_T c9_d_data[400],
                    int32_T c9_d_sizes[2]);
static void c9_check_forloop_overflow_error(SFc9_simulationInstanceStruct
  *chartInstance, boolean_T c9_overflow);
static void c9_eye(SFc9_simulationInstanceStruct *chartInstance, real_T
                   c9_varargin_1, real_T c9_I_data[400], int32_T c9_I_sizes[2]);
static void c9_eml_assert_valid_size_arg(SFc9_simulationInstanceStruct
  *chartInstance, real_T c9_varargin_1);
static void c9_mldivide(SFc9_simulationInstanceStruct *chartInstance, real_T
  c9_A_data[400], int32_T c9_A_sizes[2], real_T c9_B_data[400], int32_T
  c9_B_sizes[2], real_T c9_Y_data[400], int32_T c9_Y_sizes[2]);
static void c9_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance);
static void c9_eml_lusolve(SFc9_simulationInstanceStruct *chartInstance, real_T
  c9_A_data[400], int32_T c9_A_sizes[2], real_T c9_B_data[400], int32_T
  c9_B_sizes[2], real_T c9_X_data[400], int32_T c9_X_sizes[2]);
static void c9_realmin(SFc9_simulationInstanceStruct *chartInstance);
static void c9_eps(SFc9_simulationInstanceStruct *chartInstance);
static void c9_b_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance);
static void c9_eml_signed_integer_colon(SFc9_simulationInstanceStruct
  *chartInstance, int32_T c9_b, int32_T c9_y_data[20], int32_T c9_y_sizes[2]);
static boolean_T c9_eml_use_refblas(SFc9_simulationInstanceStruct *chartInstance);
static real_T c9_abs(SFc9_simulationInstanceStruct *chartInstance, real_T c9_x);
static void c9_eml_xswap(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0, int32_T
  c9_incx, int32_T c9_iy0, int32_T c9_incy, real_T c9_b_x_data[400], int32_T
  c9_b_x_sizes[2]);
static void c9_eml_xger(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_alpha1, int32_T c9_ix0, int32_T c9_iy0,
  int32_T c9_incy, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T c9_ia0,
  int32_T c9_lda, real_T c9_b_A_data[400], int32_T c9_b_A_sizes[2]);
static void c9_c_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance);
static void c9_eml_warning(SFc9_simulationInstanceStruct *chartInstance);
static void c9_eml_xtrsm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B_data[400], int32_T c9_B_sizes[2], int32_T c9_ldb, real_T
  c9_b_B_data[400], int32_T c9_b_B_sizes[2]);
static void c9_d_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance);
static void c9_b_eml_xtrsm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B_data[400], int32_T c9_B_sizes[2], int32_T c9_ldb, real_T
  c9_b_B_data[400], int32_T c9_b_B_sizes[2]);
static void c9_eml_qrsolve(SFc9_simulationInstanceStruct *chartInstance, real_T
  c9_A_data[400], int32_T c9_A_sizes[2], real_T c9_B_data[400], int32_T
  c9_B_sizes[2], real_T c9_Y_data[400], int32_T c9_Y_sizes[2]);
static real_T c9_sqrt(SFc9_simulationInstanceStruct *chartInstance, real_T c9_x);
static real_T c9_eml_xnrm2(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0);
static int32_T c9_eml_ixamax(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_b_n, real_T c9_x_data[20], int32_T c9_x_sizes[1], int32_T c9_ix0);
static void c9_b_eml_xswap(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0, int32_T
  c9_iy0, real_T c9_b_x_data[400], int32_T c9_b_x_sizes[2]);
static void c9_eml_matlab_zlarfg(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_b_n, real_T c9_alpha1, real_T c9_x_data[400], int32_T c9_x_sizes[2],
  int32_T c9_ix0, real_T *c9_b_alpha1, real_T c9_b_x_data[400], int32_T
  c9_b_x_sizes[2], real_T *c9_tau);
static void c9_eml_xscal(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_a, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T
  c9_ix0, real_T c9_b_x_data[400], int32_T c9_b_x_sizes[2]);
static void c9_b_eml_matlab_zlarfg(SFc9_simulationInstanceStruct *chartInstance,
  real_T c9_alpha1, real_T c9_x, real_T *c9_b_alpha1, real_T *c9_b_x, real_T
  *c9_tau);
static void c9_b_eml_xnrm2(SFc9_simulationInstanceStruct *chartInstance);
static real_T c9_b_eml_xscal(SFc9_simulationInstanceStruct *chartInstance,
  real_T c9_x);
static void c9_e_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance);
static void c9_eml_matlab_zlarf(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_m, int32_T c9_b_n, int32_T c9_iv0, real_T c9_tau, real_T c9_C_data
  [400], int32_T c9_C_sizes[2], int32_T c9_ic0, int32_T c9_ldc, real_T
  c9_work_data[20], int32_T c9_work_sizes[1], real_T c9_b_C_data[400], int32_T
  c9_b_C_sizes[2], real_T c9_b_work_data[20], int32_T c9_b_work_sizes[1]);
static void c9_eml_xgemv(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_ia0, int32_T c9_lda, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T
  c9_ix0, real_T c9_y_data[20], int32_T c9_y_sizes[1], real_T c9_b_y_data[20],
  int32_T c9_b_y_sizes[1]);
static void c9_b_eml_xger(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_alpha1, int32_T c9_ix0, real_T c9_y_data[20],
  int32_T c9_y_sizes[1], real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_ia0, int32_T c9_lda, real_T c9_b_A_data[400], int32_T c9_b_A_sizes[2]);
static real_T c9_c_eml_xnrm2(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0);
static void c9_eml_flt2str(SFc9_simulationInstanceStruct *chartInstance, real_T
  c9_x, char_T c9_str[14]);
static void c9_b_eml_warning(SFc9_simulationInstanceStruct *chartInstance,
  real_T c9_varargin_2, char_T c9_varargin_3[14]);
static void c9_f_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance);
static void c9_eml_xgemm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_k, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B[4], int32_T c9_ldb, real_T c9_C_data[20], int32_T
  c9_C_sizes[1], int32_T c9_ldc, real_T c9_b_C_data[20], int32_T c9_b_C_sizes[1]);
static void c9_below_threshold(SFc9_simulationInstanceStruct *chartInstance);
static void c9_tanh(SFc9_simulationInstanceStruct *chartInstance, real_T
                    c9_x_data[20], int32_T c9_x_sizes[1], real_T c9_b_x_data[20],
                    int32_T c9_b_x_sizes[1]);
static void c9_e_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_sprintf, const char_T *c9_identifier, char_T c9_y[14]);
static void c9_f_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, char_T c9_y[14]);
static const mxArray *c9_e_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData);
static int32_T c9_g_emlrt_marshallIn(SFc9_simulationInstanceStruct
  *chartInstance, const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId);
static void c9_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData);
static uint8_T c9_h_emlrt_marshallIn(SFc9_simulationInstanceStruct
  *chartInstance, const mxArray *c9_b_is_active_c9_simulation, const char_T
  *c9_identifier);
static uint8_T c9_i_emlrt_marshallIn(SFc9_simulationInstanceStruct
  *chartInstance, const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId);
static void c9_c_eml_xswap(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0, int32_T
  c9_incx, int32_T c9_iy0, int32_T c9_incy);
static void c9_c_eml_xger(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_alpha1, int32_T c9_ix0, int32_T c9_iy0,
  int32_T c9_incy, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T c9_ia0,
  int32_T c9_lda);
static void c9_c_eml_xtrsm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B_data[400], int32_T c9_B_sizes[2], int32_T c9_ldb);
static void c9_d_eml_xtrsm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B_data[400], int32_T c9_B_sizes[2], int32_T c9_ldb);
static void c9_b_sqrt(SFc9_simulationInstanceStruct *chartInstance, real_T *c9_x);
static void c9_d_eml_xswap(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0, int32_T
  c9_iy0);
static real_T c9_c_eml_matlab_zlarfg(SFc9_simulationInstanceStruct
  *chartInstance, int32_T c9_b_n, real_T *c9_alpha1, real_T c9_x_data[400],
  int32_T c9_x_sizes[2], int32_T c9_ix0);
static void c9_c_eml_xscal(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_a, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T
  c9_ix0);
static real_T c9_d_eml_matlab_zlarfg(SFc9_simulationInstanceStruct
  *chartInstance, real_T *c9_alpha1, real_T *c9_x);
static void c9_d_eml_xscal(SFc9_simulationInstanceStruct *chartInstance, real_T *
  c9_x);
static void c9_b_eml_matlab_zlarf(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_m, int32_T c9_b_n, int32_T c9_iv0, real_T c9_tau, real_T c9_C_data
  [400], int32_T c9_C_sizes[2], int32_T c9_ic0, int32_T c9_ldc, real_T
  c9_work_data[20], int32_T c9_work_sizes[1]);
static void c9_b_eml_xgemv(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_ia0, int32_T c9_lda, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T
  c9_ix0, real_T c9_y_data[20], int32_T c9_y_sizes[1]);
static void c9_d_eml_xger(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_alpha1, int32_T c9_ix0, real_T c9_y_data[20],
  int32_T c9_y_sizes[1], real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_ia0, int32_T c9_lda);
static void c9_b_eml_xgemm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_k, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B[4], int32_T c9_ldb, real_T c9_C_data[20], int32_T
  c9_C_sizes[1], int32_T c9_ldc);
static void c9_b_tanh(SFc9_simulationInstanceStruct *chartInstance, real_T
                      c9_x_data[20], int32_T c9_x_sizes[1]);
static void init_dsm_address_info(SFc9_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c9_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c9_is_active_c9_simulation = 0U;
}

static void initialize_params_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance)
{
  real_T c9_d0;
  real_T c9_d1;
  real_T c9_d2;
  sf_set_error_prefix_string(
    "Error evaluating data 'n' in the parent workspace.\n");
  sf_mex_import_named("n", sf_mex_get_sfun_param(chartInstance->S, 2, 0), &c9_d0,
                      0, 0, 0U, 0, 0U, 0);
  chartInstance->c9_n = c9_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'k_xi' in the parent workspace.\n");
  sf_mex_import_named("k_xi", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      &c9_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c9_k_xi = c9_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'i' in the parent workspace.\n");
  sf_mex_import_named("i", sf_mex_get_sfun_param(chartInstance->S, 0, 0), &c9_d2,
                      0, 0, 0U, 0, 0U, 0);
  chartInstance->c9_i = c9_d2;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c9_simulation(SFc9_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c9_simulation(SFc9_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c9_update_debugger_state_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c9_simulation(SFc9_simulationInstanceStruct *
  chartInstance)
{
  const mxArray *c9_st;
  const mxArray *c9_y = NULL;
  real_T c9_hoistedGlobal;
  real_T c9_u;
  const mxArray *c9_b_y = NULL;
  uint8_T c9_b_hoistedGlobal;
  uint8_T c9_b_u;
  const mxArray *c9_c_y = NULL;
  real_T *c9_v_tilde;
  c9_v_tilde = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c9_st = NULL;
  c9_st = NULL;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_createcellarray(2), FALSE);
  c9_hoistedGlobal = *c9_v_tilde;
  c9_u = c9_hoistedGlobal;
  c9_b_y = NULL;
  sf_mex_assign(&c9_b_y, sf_mex_create("y", &c9_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c9_y, 0, c9_b_y);
  c9_b_hoistedGlobal = chartInstance->c9_is_active_c9_simulation;
  c9_b_u = c9_b_hoistedGlobal;
  c9_c_y = NULL;
  sf_mex_assign(&c9_c_y, sf_mex_create("y", &c9_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c9_y, 1, c9_c_y);
  sf_mex_assign(&c9_st, c9_y, FALSE);
  return c9_st;
}

static void set_sim_state_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance, const mxArray *c9_st)
{
  const mxArray *c9_u;
  real_T *c9_v_tilde;
  c9_v_tilde = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c9_doneDoubleBufferReInit = TRUE;
  c9_u = sf_mex_dup(c9_st);
  *c9_v_tilde = c9_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c9_u, 0)), "v_tilde");
  chartInstance->c9_is_active_c9_simulation = c9_h_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c9_u, 1)),
     "is_active_c9_simulation");
  sf_mex_destroy(&c9_u);
  c9_update_debugger_state_c9_simulation(chartInstance);
  sf_mex_destroy(&c9_st);
}

static void finalize_c9_simulation(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void sf_c9_simulation(SFc9_simulationInstanceStruct *chartInstance)
{
  int32_T c9_i0;
  real_T *c9_v_tilde;
  real_T (*c9_Gam_est)[4];
  c9_v_tilde = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c9_Gam_est = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 8U, chartInstance->c9_sfEvent);
  for (c9_i0 = 0; c9_i0 < 4; c9_i0++) {
    _SFD_DATA_RANGE_CHECK((*c9_Gam_est)[c9_i0], 0U);
  }

  _SFD_DATA_RANGE_CHECK(*c9_v_tilde, 1U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c9_n, 2U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c9_k_xi, 3U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c9_i, 4U);
  chartInstance->c9_sfEvent = CALL_EVENT;
  c9_chartstep_c9_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c9_chartstep_c9_simulation(SFc9_simulationInstanceStruct
  *chartInstance)
{
  real_T c9_hoistedGlobal;
  real_T c9_b_hoistedGlobal;
  real_T c9_c_hoistedGlobal;
  int32_T c9_i1;
  real_T c9_Gam_est[4];
  real_T c9_b_n;
  real_T c9_b_k_xi;
  real_T c9_b_i;
  uint32_T c9_debug_family_var_map[11];
  int32_T c9_D_sizes[2];
  real_T c9_D_data[400];
  int32_T c9_A_sizes[2];
  real_T c9_A_data[400];
  int32_T c9_L_d_sizes[2];
  real_T c9_L_d_data[400];
  int32_T c9_V_tilde_sizes;
  real_T c9_V_tilde_data[20];
  real_T c9_nargin = 4.0;
  real_T c9_nargout = 1.0;
  real_T c9_v_tilde;
  const mxArray *c9_y = NULL;
  real_T c9_a;
  int32_T c9_b_sizes[2];
  int32_T c9_iv0[2];
  int32_T c9_b;
  int32_T c9_b_b;
  int32_T c9_loop_ub;
  int32_T c9_i2;
  real_T c9_b_data[20];
  int32_T c9_c_b;
  int32_T c9_d_b;
  int32_T c9_e_b;
  int32_T c9_f_b;
  int32_T c9_b_loop_ub;
  int32_T c9_i3;
  int32_T c9_b_b_sizes[2];
  int32_T c9_g_b;
  int32_T c9_h_b;
  int32_T c9_c_loop_ub;
  int32_T c9_i4;
  real_T c9_b_b_data[20];
  int32_T c9_tmp_sizes[2];
  real_T c9_tmp_data[400];
  int32_T c9_D;
  int32_T c9_b_D;
  int32_T c9_d_loop_ub;
  int32_T c9_i5;
  int32_T c9_i6;
  int32_T c9_iv1[2];
  int32_T c9_a_sizes[2];
  int32_T c9_iv2[2];
  int32_T c9_b_a;
  int32_T c9_c_a;
  int32_T c9_e_loop_ub;
  int32_T c9_i7;
  real_T c9_a_data[400];
  int32_T c9_b_tmp_sizes[2];
  real_T c9_b_tmp_data[400];
  int32_T c9_i8;
  int32_T c9_d_a[2];
  int32_T c9_i9;
  int32_T c9_iv3[2];
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_f_loop_ub;
  int32_T c9_i10;
  int32_T c9_i11;
  int32_T c9_c_D[2];
  int32_T c9_i12;
  int32_T c9_c_A[2];
  int32_T c9_b_D_sizes[2];
  int32_T c9_d_D;
  int32_T c9_e_D;
  int32_T c9_g_loop_ub;
  int32_T c9_i13;
  real_T c9_b_D_data[400];
  int32_T c9_c_D_sizes[2];
  int32_T c9_f_D;
  int32_T c9_g_D;
  int32_T c9_h_loop_ub;
  int32_T c9_i14;
  real_T c9_c_D_data[400];
  int32_T c9_c_tmp_sizes[2];
  real_T c9_c_tmp_data[400];
  int32_T c9_L_d;
  int32_T c9_b_L_d;
  int32_T c9_i_loop_ub;
  int32_T c9_i15;
  int32_T c9_e_a;
  int32_T c9_f_a;
  int32_T c9_j_loop_ub;
  int32_T c9_i16;
  int32_T c9_i17;
  real_T c9_i_b[4];
  boolean_T c9_innerDimOk;
  boolean_T c9_b0;
  boolean_T c9_b1;
  int32_T c9_i18;
  static char_T c9_cv0[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D', 'y',
    'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p', 'a',
    'n', 's', 'i', 'o', 'n' };

  char_T c9_u[45];
  const mxArray *c9_b_y = NULL;
  int32_T c9_i19;
  static char_T c9_cv1[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  char_T c9_b_u[21];
  const mxArray *c9_c_y = NULL;
  int32_T c9_y_sizes;
  int32_T c9_k_loop_ub;
  int32_T c9_i20;
  real_T c9_y_data[20];
  int32_T c9_l_loop_ub;
  int32_T c9_i21;
  int32_T c9_k;
  real_T c9_dv0[2];
  int32_T c9_d_tmp_sizes;
  int32_T c9_m_loop_ub;
  int32_T c9_i22;
  real_T c9_d_tmp_data[20];
  int32_T c9_m;
  int32_T c9_d_y;
  int32_T c9_n_loop_ub;
  int32_T c9_i23;
  int32_T c9_b_a_sizes[2];
  int32_T c9_g_a;
  int32_T c9_h_a;
  int32_T c9_o_loop_ub;
  int32_T c9_i24;
  real_T c9_b_a_data[400];
  int32_T c9_i25;
  real_T c9_j_b[4];
  real_T c9_i_a;
  int32_T c9_p_loop_ub;
  int32_T c9_i26;
  real_T *c9_b_v_tilde;
  real_T (*c9_b_Gam_est)[4];
  c9_b_v_tilde = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c9_b_Gam_est = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 8U, chartInstance->c9_sfEvent);
  c9_hoistedGlobal = chartInstance->c9_n;
  c9_b_hoistedGlobal = chartInstance->c9_k_xi;
  c9_c_hoistedGlobal = chartInstance->c9_i;
  for (c9_i1 = 0; c9_i1 < 4; c9_i1++) {
    c9_Gam_est[c9_i1] = (*c9_b_Gam_est)[c9_i1];
  }

  c9_b_n = c9_hoistedGlobal;
  c9_b_k_xi = c9_b_hoistedGlobal;
  c9_b_i = c9_c_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 11U, 11U, c9_debug_family_names,
    c9_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_IMPORTABLE(c9_D_data, (const int32_T *)
    &c9_D_sizes, NULL, 0, 0, (void *)c9_d_sf_marshallOut, (void *)
    c9_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_IMPORTABLE(c9_A_data, (const int32_T *)
    &c9_A_sizes, NULL, 0, 1, (void *)c9_d_sf_marshallOut, (void *)
    c9_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_IMPORTABLE(c9_L_d_data, (const int32_T *)
    &c9_L_d_sizes, NULL, 0, 2, (void *)c9_d_sf_marshallOut, (void *)
    c9_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_IMPORTABLE(c9_V_tilde_data, (const int32_T *)
    &c9_V_tilde_sizes, NULL, 0, 3, (void *)c9_c_sf_marshallOut, (void *)
    c9_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_nargin, 4U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_nargout, 5U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c9_Gam_est, 6U, c9_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_b_n, 7U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_b_k_xi, 8U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_b_i, 9U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c9_v_tilde, 10U, c9_sf_marshallOut,
    c9_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 2);
  if (c9_b_n < 20.0) {
  } else {
    c9_y = NULL;
    sf_mex_assign(&c9_y, sf_mex_create("y", "Assertion failed.", 15, 0U, 0U, 0U,
      2, 1, strlen("Assertion failed.")), FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, c9_y);
  }

  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 3);
  c9_a = c9_b_n - 1.0;
  c9_b_sizes[0] = 1;
  c9_iv0[0] = 1;
  c9_iv0[1] = (int32_T)_SFD_INTEGER_CHECK("n", _SFD_NON_NEGATIVE_CHECK("n",
    c9_b_n));
  c9_b_sizes[1] = c9_iv0[1];
  c9_b = c9_b_sizes[0];
  c9_b_b = c9_b_sizes[1];
  c9_loop_ub = (int32_T)_SFD_INTEGER_CHECK("n", _SFD_NON_NEGATIVE_CHECK("n",
    c9_b_n)) - 1;
  for (c9_i2 = 0; c9_i2 <= c9_loop_ub; c9_i2++) {
    c9_b_data[c9_i2] = 1.0;
  }

  c9_b_sizes[0] = 1;
  c9_b_sizes[1];
  c9_c_b = c9_b_sizes[0];
  c9_d_b = c9_b_sizes[1];
  c9_e_b = c9_b_sizes[0];
  c9_f_b = c9_b_sizes[1];
  c9_b_loop_ub = c9_e_b * c9_f_b - 1;
  for (c9_i3 = 0; c9_i3 <= c9_b_loop_ub; c9_i3++) {
    c9_b_data[c9_i3] *= c9_a;
  }

  c9_b_b_sizes[0] = 1;
  c9_b_b_sizes[1] = c9_b_sizes[1];
  c9_g_b = c9_b_b_sizes[0];
  c9_h_b = c9_b_b_sizes[1];
  c9_c_loop_ub = c9_b_sizes[0] * c9_b_sizes[1] - 1;
  for (c9_i4 = 0; c9_i4 <= c9_c_loop_ub; c9_i4++) {
    c9_b_b_data[c9_i4] = c9_b_data[c9_i4];
  }

  c9_diag(chartInstance, c9_b_b_data, c9_b_b_sizes, c9_tmp_data, c9_tmp_sizes);
  c9_D_sizes[0] = c9_tmp_sizes[0];
  c9_D_sizes[1] = c9_tmp_sizes[1];
  c9_D = c9_D_sizes[0];
  c9_b_D = c9_D_sizes[1];
  c9_d_loop_ub = c9_tmp_sizes[0] * c9_tmp_sizes[1] - 1;
  for (c9_i5 = 0; c9_i5 <= c9_d_loop_ub; c9_i5++) {
    c9_D_data[c9_i5] = c9_tmp_data[c9_i5];
  }

  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 4);
  c9_i6 = (int32_T)c9_b_n;
  c9_iv1[0] = c9_i6;
  c9_iv1[1] = c9_i6;
  c9_a_sizes[0] = c9_iv1[0];
  c9_iv2[0] = c9_i6;
  c9_iv2[1] = c9_i6;
  c9_a_sizes[1] = c9_iv2[1];
  c9_b_a = c9_a_sizes[0];
  c9_c_a = c9_a_sizes[1];
  c9_e_loop_ub = c9_i6 * c9_i6 - 1;
  for (c9_i7 = 0; c9_i7 <= c9_e_loop_ub; c9_i7++) {
    c9_a_data[c9_i7] = 1.0;
  }

  c9_eye(chartInstance, c9_b_n, c9_b_tmp_data, c9_b_tmp_sizes);
  for (c9_i8 = 0; c9_i8 < 2; c9_i8++) {
    c9_d_a[c9_i8] = c9_a_sizes[c9_i8];
  }

  for (c9_i9 = 0; c9_i9 < 2; c9_i9++) {
    c9_iv3[c9_i9] = c9_b_tmp_sizes[c9_i9];
  }

  _SFD_SIZE_EQ_CHECK_ND(c9_d_a, c9_iv3, 2);
  c9_A_sizes[0] = c9_a_sizes[0];
  c9_A_sizes[1] = c9_a_sizes[1];
  c9_A = c9_A_sizes[0];
  c9_b_A = c9_A_sizes[1];
  c9_f_loop_ub = c9_a_sizes[0] * c9_a_sizes[1] - 1;
  for (c9_i10 = 0; c9_i10 <= c9_f_loop_ub; c9_i10++) {
    c9_A_data[c9_i10] = c9_a_data[c9_i10] - c9_b_tmp_data[c9_i10];
  }

  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 5);
  for (c9_i11 = 0; c9_i11 < 2; c9_i11++) {
    c9_c_D[c9_i11] = c9_D_sizes[c9_i11];
  }

  for (c9_i12 = 0; c9_i12 < 2; c9_i12++) {
    c9_c_A[c9_i12] = c9_A_sizes[c9_i12];
  }

  _SFD_SIZE_EQ_CHECK_ND(c9_c_D, c9_c_A, 2);
  c9_b_D_sizes[0] = c9_D_sizes[0];
  c9_b_D_sizes[1] = c9_D_sizes[1];
  c9_d_D = c9_b_D_sizes[0];
  c9_e_D = c9_b_D_sizes[1];
  c9_g_loop_ub = c9_D_sizes[0] * c9_D_sizes[1] - 1;
  for (c9_i13 = 0; c9_i13 <= c9_g_loop_ub; c9_i13++) {
    c9_b_D_data[c9_i13] = c9_D_data[c9_i13];
  }

  c9_c_D_sizes[0] = c9_D_sizes[0];
  c9_c_D_sizes[1] = c9_D_sizes[1];
  c9_f_D = c9_c_D_sizes[0];
  c9_g_D = c9_c_D_sizes[1];
  c9_h_loop_ub = c9_D_sizes[0] * c9_D_sizes[1] - 1;
  for (c9_i14 = 0; c9_i14 <= c9_h_loop_ub; c9_i14++) {
    c9_c_D_data[c9_i14] = c9_D_data[c9_i14] - c9_A_data[c9_i14];
  }

  c9_mldivide(chartInstance, c9_b_D_data, c9_b_D_sizes, c9_c_D_data,
              c9_c_D_sizes, c9_c_tmp_data, c9_c_tmp_sizes);
  c9_L_d_sizes[0] = c9_c_tmp_sizes[0];
  c9_L_d_sizes[1] = c9_c_tmp_sizes[1];
  c9_L_d = c9_L_d_sizes[0];
  c9_b_L_d = c9_L_d_sizes[1];
  c9_i_loop_ub = c9_c_tmp_sizes[0] * c9_c_tmp_sizes[1] - 1;
  for (c9_i15 = 0; c9_i15 <= c9_i_loop_ub; c9_i15++) {
    c9_L_d_data[c9_i15] = c9_c_tmp_data[c9_i15];
  }

  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 6);
  c9_a_sizes[0] = c9_L_d_sizes[0];
  c9_a_sizes[1] = c9_L_d_sizes[1];
  c9_e_a = c9_a_sizes[0];
  c9_f_a = c9_a_sizes[1];
  c9_j_loop_ub = c9_L_d_sizes[0] * c9_L_d_sizes[1] - 1;
  for (c9_i16 = 0; c9_i16 <= c9_j_loop_ub; c9_i16++) {
    c9_a_data[c9_i16] = c9_L_d_data[c9_i16];
  }

  for (c9_i17 = 0; c9_i17 < 4; c9_i17++) {
    c9_i_b[c9_i17] = c9_Gam_est[c9_i17];
  }

  c9_innerDimOk = ((real_T)c9_a_sizes[1] == 4.0);
  if (!c9_innerDimOk) {
    c9_b0 = (c9_a_sizes[0] == 1);
    c9_b1 = (c9_a_sizes[1] == 1);
    if (c9_b0 && c9_b1) {
      for (c9_i18 = 0; c9_i18 < 45; c9_i18++) {
        c9_u[c9_i18] = c9_cv0[c9_i18];
      }

      c9_b_y = NULL;
      sf_mex_assign(&c9_b_y, sf_mex_create("y", c9_u, 10, 0U, 1U, 0U, 2, 1, 45),
                    FALSE);
      sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
        14, c9_b_y));
    } else {
      for (c9_i19 = 0; c9_i19 < 21; c9_i19++) {
        c9_b_u[c9_i19] = c9_cv1[c9_i19];
      }

      c9_c_y = NULL;
      sf_mex_assign(&c9_c_y, sf_mex_create("y", c9_b_u, 10, 0U, 1U, 0U, 2, 1, 21),
                    FALSE);
      sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
        14, c9_c_y));
    }
  }

  if ((real_T)c9_a_sizes[1] == 1.0) {
    c9_y_sizes = c9_a_sizes[0];
    c9_k_loop_ub = c9_a_sizes[0] - 1;
    for (c9_i20 = 0; c9_i20 <= c9_k_loop_ub; c9_i20++) {
      c9_y_data[c9_i20] = 0.0;
      c9_l_loop_ub = c9_a_sizes[1] - 1;
      for (c9_i21 = 0; c9_i21 <= c9_l_loop_ub; c9_i21++) {
        c9_y_data[c9_i20] += c9_a_data[c9_i20 + c9_a_sizes[0] * c9_i21] *
          c9_i_b[c9_i21];
      }
    }
  } else {
    c9_k = (int32_T)(real_T)c9_a_sizes[1];
    c9_f_eml_scalar_eg(chartInstance);
    c9_dv0[0] = (real_T)c9_a_sizes[0];
    c9_dv0[1] = 1.0;
    c9_d_tmp_sizes = (int32_T)c9_dv0[0];
    c9_m_loop_ub = (int32_T)c9_dv0[0] - 1;
    for (c9_i22 = 0; c9_i22 <= c9_m_loop_ub; c9_i22++) {
      c9_d_tmp_data[c9_i22] = 0.0;
    }

    c9_y_sizes = c9_d_tmp_sizes;
    c9_m = (int32_T)(real_T)c9_a_sizes[0];
    c9_f_eml_scalar_eg(chartInstance);
    c9_d_y = c9_y_sizes;
    c9_y_sizes = c9_d_y;
    c9_n_loop_ub = c9_d_y - 1;
    for (c9_i23 = 0; c9_i23 <= c9_n_loop_ub; c9_i23++) {
      c9_y_data[c9_i23] = 0.0;
    }

    c9_b_a_sizes[0] = c9_a_sizes[0];
    c9_b_a_sizes[1] = c9_a_sizes[1];
    c9_g_a = c9_b_a_sizes[0];
    c9_h_a = c9_b_a_sizes[1];
    c9_o_loop_ub = c9_a_sizes[0] * c9_a_sizes[1] - 1;
    for (c9_i24 = 0; c9_i24 <= c9_o_loop_ub; c9_i24++) {
      c9_b_a_data[c9_i24] = c9_a_data[c9_i24];
    }

    for (c9_i25 = 0; c9_i25 < 4; c9_i25++) {
      c9_j_b[c9_i25] = c9_i_b[c9_i25];
    }

    c9_b_eml_xgemm(chartInstance, c9_m, c9_k, c9_b_a_data, c9_b_a_sizes, c9_m,
                   c9_j_b, c9_k, c9_y_data, &c9_y_sizes, c9_m);
  }

  c9_i_a = -c9_b_k_xi;
  c9_b_tanh(chartInstance, c9_y_data, &c9_y_sizes);
  c9_V_tilde_sizes = c9_y_sizes;
  c9_p_loop_ub = c9_y_sizes - 1;
  for (c9_i26 = 0; c9_i26 <= c9_p_loop_ub; c9_i26++) {
    c9_V_tilde_data[c9_i26] = c9_i_a * c9_y_data[c9_i26];
  }

  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, 7);
  c9_v_tilde = c9_V_tilde_data[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
    "V_tilde", (int32_T)_SFD_INTEGER_CHECK("i", c9_b_i), 1, c9_V_tilde_sizes, 1,
    0) - 1];
  _SFD_EML_CALL(0U, chartInstance->c9_sfEvent, -7);
  _SFD_SYMBOL_SCOPE_POP();
  *c9_b_v_tilde = c9_v_tilde;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 8U, chartInstance->c9_sfEvent);
}

static void initSimStructsc9_simulation(SFc9_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c9_machineNumber, uint32_T
  c9_chartNumber)
{
}

static const mxArray *c9_sf_marshallOut(void *chartInstanceVoid, void *c9_inData)
{
  const mxArray *c9_mxArrayOutData = NULL;
  real_T c9_u;
  const mxArray *c9_y = NULL;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  c9_u = *(real_T *)c9_inData;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", &c9_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, FALSE);
  return c9_mxArrayOutData;
}

static real_T c9_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_v_tilde, const char_T *c9_identifier)
{
  real_T c9_y;
  emlrtMsgIdentifier c9_thisId;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_y = c9_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_v_tilde), &c9_thisId);
  sf_mex_destroy(&c9_v_tilde);
  return c9_y;
}

static real_T c9_b_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId)
{
  real_T c9_y;
  real_T c9_d3;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), &c9_d3, 1, 0, 0U, 0, 0U, 0);
  c9_y = c9_d3;
  sf_mex_destroy(&c9_u);
  return c9_y;
}

static void c9_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData)
{
  const mxArray *c9_v_tilde;
  const char_T *c9_identifier;
  emlrtMsgIdentifier c9_thisId;
  real_T c9_y;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_v_tilde = sf_mex_dup(c9_mxArrayInData);
  c9_identifier = c9_varName;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_y = c9_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_v_tilde), &c9_thisId);
  sf_mex_destroy(&c9_v_tilde);
  *(real_T *)c9_outData = c9_y;
  sf_mex_destroy(&c9_mxArrayInData);
}

static const mxArray *c9_b_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData)
{
  const mxArray *c9_mxArrayOutData = NULL;
  int32_T c9_i27;
  real_T c9_b_inData[4];
  int32_T c9_i28;
  real_T c9_u[4];
  const mxArray *c9_y = NULL;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  for (c9_i27 = 0; c9_i27 < 4; c9_i27++) {
    c9_b_inData[c9_i27] = (*(real_T (*)[4])c9_inData)[c9_i27];
  }

  for (c9_i28 = 0; c9_i28 < 4; c9_i28++) {
    c9_u[c9_i28] = c9_b_inData[c9_i28];
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, FALSE);
  return c9_mxArrayOutData;
}

static const mxArray *c9_c_sf_marshallOut(void *chartInstanceVoid, real_T
  c9_inData_data[20], int32_T c9_inData_sizes[1])
{
  const mxArray *c9_mxArrayOutData = NULL;
  int32_T c9_b_inData_sizes;
  int32_T c9_loop_ub;
  int32_T c9_i29;
  real_T c9_b_inData_data[20];
  int32_T c9_u_sizes;
  int32_T c9_b_loop_ub;
  int32_T c9_i30;
  real_T c9_u_data[20];
  const mxArray *c9_y = NULL;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  c9_b_inData_sizes = c9_inData_sizes[0];
  c9_loop_ub = c9_inData_sizes[0] - 1;
  for (c9_i29 = 0; c9_i29 <= c9_loop_ub; c9_i29++) {
    c9_b_inData_data[c9_i29] = c9_inData_data[c9_i29];
  }

  c9_u_sizes = c9_b_inData_sizes;
  c9_b_loop_ub = c9_b_inData_sizes - 1;
  for (c9_i30 = 0; c9_i30 <= c9_b_loop_ub; c9_i30++) {
    c9_u_data[c9_i30] = c9_b_inData_data[c9_i30];
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u_data, 0, 0U, 1U, 0U, 1,
    c9_u_sizes), FALSE);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, FALSE);
  return c9_mxArrayOutData;
}

static void c9_c_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, real_T c9_y_data
  [20], int32_T c9_y_sizes[1])
{
  static uint32_T c9_uv0[1] = { 20U };

  uint32_T c9_uv1[1];
  static boolean_T c9_bv0[1] = { TRUE };

  boolean_T c9_bv1[1];
  int32_T c9_tmp_sizes;
  real_T c9_tmp_data[20];
  int32_T c9_loop_ub;
  int32_T c9_i31;
  c9_uv1[0] = c9_uv0[0];
  c9_bv1[0] = c9_bv0[0];
  sf_mex_import_vs(c9_parentId, sf_mex_dup(c9_u), c9_tmp_data, 1, 0, 0U, 1, 0U,
                   1, c9_bv1, c9_uv1, &c9_tmp_sizes);
  c9_y_sizes[0] = c9_tmp_sizes;
  c9_loop_ub = c9_tmp_sizes - 1;
  for (c9_i31 = 0; c9_i31 <= c9_loop_ub; c9_i31++) {
    c9_y_data[c9_i31] = c9_tmp_data[c9_i31];
  }

  sf_mex_destroy(&c9_u);
}

static void c9_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, real_T c9_outData_data[20],
  int32_T c9_outData_sizes[1])
{
  const mxArray *c9_V_tilde;
  const char_T *c9_identifier;
  emlrtMsgIdentifier c9_thisId;
  int32_T c9_y_sizes;
  real_T c9_y_data[20];
  int32_T c9_loop_ub;
  int32_T c9_i32;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_V_tilde = sf_mex_dup(c9_mxArrayInData);
  c9_identifier = c9_varName;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_V_tilde), &c9_thisId,
                        c9_y_data, *(int32_T (*)[1])&c9_y_sizes);
  sf_mex_destroy(&c9_V_tilde);
  c9_outData_sizes[0] = c9_y_sizes;
  c9_loop_ub = c9_y_sizes - 1;
  for (c9_i32 = 0; c9_i32 <= c9_loop_ub; c9_i32++) {
    c9_outData_data[c9_i32] = c9_y_data[c9_i32];
  }

  sf_mex_destroy(&c9_mxArrayInData);
}

static const mxArray *c9_d_sf_marshallOut(void *chartInstanceVoid, real_T
  c9_inData_data[400], int32_T c9_inData_sizes[2])
{
  const mxArray *c9_mxArrayOutData = NULL;
  int32_T c9_b_inData_sizes[2];
  int32_T c9_loop_ub;
  int32_T c9_i33;
  int32_T c9_b_loop_ub;
  int32_T c9_i34;
  real_T c9_b_inData_data[400];
  int32_T c9_u_sizes[2];
  int32_T c9_c_loop_ub;
  int32_T c9_i35;
  int32_T c9_d_loop_ub;
  int32_T c9_i36;
  real_T c9_u_data[400];
  const mxArray *c9_y = NULL;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  c9_b_inData_sizes[0] = c9_inData_sizes[0];
  c9_b_inData_sizes[1] = c9_inData_sizes[1];
  c9_loop_ub = c9_inData_sizes[1] - 1;
  for (c9_i33 = 0; c9_i33 <= c9_loop_ub; c9_i33++) {
    c9_b_loop_ub = c9_inData_sizes[0] - 1;
    for (c9_i34 = 0; c9_i34 <= c9_b_loop_ub; c9_i34++) {
      c9_b_inData_data[c9_i34 + c9_b_inData_sizes[0] * c9_i33] =
        c9_inData_data[c9_i34 + c9_inData_sizes[0] * c9_i33];
    }
  }

  c9_u_sizes[0] = c9_b_inData_sizes[0];
  c9_u_sizes[1] = c9_b_inData_sizes[1];
  c9_c_loop_ub = c9_b_inData_sizes[1] - 1;
  for (c9_i35 = 0; c9_i35 <= c9_c_loop_ub; c9_i35++) {
    c9_d_loop_ub = c9_b_inData_sizes[0] - 1;
    for (c9_i36 = 0; c9_i36 <= c9_d_loop_ub; c9_i36++) {
      c9_u_data[c9_i36 + c9_u_sizes[0] * c9_i35] = c9_b_inData_data[c9_i36 +
        c9_b_inData_sizes[0] * c9_i35];
    }
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u_data, 0, 0U, 1U, 0U, 2,
    c9_u_sizes[0], c9_u_sizes[1]), FALSE);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, FALSE);
  return c9_mxArrayOutData;
}

static void c9_d_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, real_T c9_y_data
  [400], int32_T c9_y_sizes[2])
{
  int32_T c9_i37;
  uint32_T c9_uv2[2];
  int32_T c9_i38;
  boolean_T c9_bv2[2];
  int32_T c9_tmp_sizes[2];
  real_T c9_tmp_data[400];
  int32_T c9_y;
  int32_T c9_b_y;
  int32_T c9_loop_ub;
  int32_T c9_i39;
  for (c9_i37 = 0; c9_i37 < 2; c9_i37++) {
    c9_uv2[c9_i37] = 20U;
  }

  for (c9_i38 = 0; c9_i38 < 2; c9_i38++) {
    c9_bv2[c9_i38] = TRUE;
  }

  sf_mex_import_vs(c9_parentId, sf_mex_dup(c9_u), c9_tmp_data, 1, 0, 0U, 1, 0U,
                   2, c9_bv2, c9_uv2, c9_tmp_sizes);
  c9_y_sizes[0] = c9_tmp_sizes[0];
  c9_y_sizes[1] = c9_tmp_sizes[1];
  c9_y = c9_y_sizes[0];
  c9_b_y = c9_y_sizes[1];
  c9_loop_ub = c9_tmp_sizes[0] * c9_tmp_sizes[1] - 1;
  for (c9_i39 = 0; c9_i39 <= c9_loop_ub; c9_i39++) {
    c9_y_data[c9_i39] = c9_tmp_data[c9_i39];
  }

  sf_mex_destroy(&c9_u);
}

static void c9_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, real_T c9_outData_data[400],
  int32_T c9_outData_sizes[2])
{
  const mxArray *c9_L_d;
  const char_T *c9_identifier;
  emlrtMsgIdentifier c9_thisId;
  int32_T c9_y_sizes[2];
  real_T c9_y_data[400];
  int32_T c9_loop_ub;
  int32_T c9_i40;
  int32_T c9_b_loop_ub;
  int32_T c9_i41;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_L_d = sf_mex_dup(c9_mxArrayInData);
  c9_identifier = c9_varName;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_L_d), &c9_thisId, c9_y_data,
                        c9_y_sizes);
  sf_mex_destroy(&c9_L_d);
  c9_outData_sizes[0] = c9_y_sizes[0];
  c9_outData_sizes[1] = c9_y_sizes[1];
  c9_loop_ub = c9_y_sizes[1] - 1;
  for (c9_i40 = 0; c9_i40 <= c9_loop_ub; c9_i40++) {
    c9_b_loop_ub = c9_y_sizes[0] - 1;
    for (c9_i41 = 0; c9_i41 <= c9_b_loop_ub; c9_i41++) {
      c9_outData_data[c9_i41 + c9_outData_sizes[0] * c9_i40] = c9_y_data[c9_i41
        + c9_y_sizes[0] * c9_i40];
    }
  }

  sf_mex_destroy(&c9_mxArrayInData);
}

const mxArray *sf_c9_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c9_nameCaptureInfo = NULL;
  c9_nameCaptureInfo = NULL;
  sf_mex_assign(&c9_nameCaptureInfo, sf_mex_createstruct("structure", 2, 250, 1),
                FALSE);
  c9_info_helper(&c9_nameCaptureInfo);
  c9_b_info_helper(&c9_nameCaptureInfo);
  c9_c_info_helper(&c9_nameCaptureInfo);
  c9_d_info_helper(&c9_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c9_nameCaptureInfo);
  return c9_nameCaptureInfo;
}

static void c9_info_helper(const mxArray **c9_info)
{
  const mxArray *c9_rhs0 = NULL;
  const mxArray *c9_lhs0 = NULL;
  const mxArray *c9_rhs1 = NULL;
  const mxArray *c9_lhs1 = NULL;
  const mxArray *c9_rhs2 = NULL;
  const mxArray *c9_lhs2 = NULL;
  const mxArray *c9_rhs3 = NULL;
  const mxArray *c9_lhs3 = NULL;
  const mxArray *c9_rhs4 = NULL;
  const mxArray *c9_lhs4 = NULL;
  const mxArray *c9_rhs5 = NULL;
  const mxArray *c9_lhs5 = NULL;
  const mxArray *c9_rhs6 = NULL;
  const mxArray *c9_lhs6 = NULL;
  const mxArray *c9_rhs7 = NULL;
  const mxArray *c9_lhs7 = NULL;
  const mxArray *c9_rhs8 = NULL;
  const mxArray *c9_lhs8 = NULL;
  const mxArray *c9_rhs9 = NULL;
  const mxArray *c9_lhs9 = NULL;
  const mxArray *c9_rhs10 = NULL;
  const mxArray *c9_lhs10 = NULL;
  const mxArray *c9_rhs11 = NULL;
  const mxArray *c9_lhs11 = NULL;
  const mxArray *c9_rhs12 = NULL;
  const mxArray *c9_lhs12 = NULL;
  const mxArray *c9_rhs13 = NULL;
  const mxArray *c9_lhs13 = NULL;
  const mxArray *c9_rhs14 = NULL;
  const mxArray *c9_lhs14 = NULL;
  const mxArray *c9_rhs15 = NULL;
  const mxArray *c9_lhs15 = NULL;
  const mxArray *c9_rhs16 = NULL;
  const mxArray *c9_lhs16 = NULL;
  const mxArray *c9_rhs17 = NULL;
  const mxArray *c9_lhs17 = NULL;
  const mxArray *c9_rhs18 = NULL;
  const mxArray *c9_lhs18 = NULL;
  const mxArray *c9_rhs19 = NULL;
  const mxArray *c9_lhs19 = NULL;
  const mxArray *c9_rhs20 = NULL;
  const mxArray *c9_lhs20 = NULL;
  const mxArray *c9_rhs21 = NULL;
  const mxArray *c9_lhs21 = NULL;
  const mxArray *c9_rhs22 = NULL;
  const mxArray *c9_lhs22 = NULL;
  const mxArray *c9_rhs23 = NULL;
  const mxArray *c9_lhs23 = NULL;
  const mxArray *c9_rhs24 = NULL;
  const mxArray *c9_lhs24 = NULL;
  const mxArray *c9_rhs25 = NULL;
  const mxArray *c9_lhs25 = NULL;
  const mxArray *c9_rhs26 = NULL;
  const mxArray *c9_lhs26 = NULL;
  const mxArray *c9_rhs27 = NULL;
  const mxArray *c9_lhs27 = NULL;
  const mxArray *c9_rhs28 = NULL;
  const mxArray *c9_lhs28 = NULL;
  const mxArray *c9_rhs29 = NULL;
  const mxArray *c9_lhs29 = NULL;
  const mxArray *c9_rhs30 = NULL;
  const mxArray *c9_lhs30 = NULL;
  const mxArray *c9_rhs31 = NULL;
  const mxArray *c9_lhs31 = NULL;
  const mxArray *c9_rhs32 = NULL;
  const mxArray *c9_lhs32 = NULL;
  const mxArray *c9_rhs33 = NULL;
  const mxArray *c9_lhs33 = NULL;
  const mxArray *c9_rhs34 = NULL;
  const mxArray *c9_lhs34 = NULL;
  const mxArray *c9_rhs35 = NULL;
  const mxArray *c9_lhs35 = NULL;
  const mxArray *c9_rhs36 = NULL;
  const mxArray *c9_lhs36 = NULL;
  const mxArray *c9_rhs37 = NULL;
  const mxArray *c9_lhs37 = NULL;
  const mxArray *c9_rhs38 = NULL;
  const mxArray *c9_lhs38 = NULL;
  const mxArray *c9_rhs39 = NULL;
  const mxArray *c9_lhs39 = NULL;
  const mxArray *c9_rhs40 = NULL;
  const mxArray *c9_lhs40 = NULL;
  const mxArray *c9_rhs41 = NULL;
  const mxArray *c9_lhs41 = NULL;
  const mxArray *c9_rhs42 = NULL;
  const mxArray *c9_lhs42 = NULL;
  const mxArray *c9_rhs43 = NULL;
  const mxArray *c9_lhs43 = NULL;
  const mxArray *c9_rhs44 = NULL;
  const mxArray *c9_lhs44 = NULL;
  const mxArray *c9_rhs45 = NULL;
  const mxArray *c9_lhs45 = NULL;
  const mxArray *c9_rhs46 = NULL;
  const mxArray *c9_lhs46 = NULL;
  const mxArray *c9_rhs47 = NULL;
  const mxArray *c9_lhs47 = NULL;
  const mxArray *c9_rhs48 = NULL;
  const mxArray *c9_lhs48 = NULL;
  const mxArray *c9_rhs49 = NULL;
  const mxArray *c9_lhs49 = NULL;
  const mxArray *c9_rhs50 = NULL;
  const mxArray *c9_lhs50 = NULL;
  const mxArray *c9_rhs51 = NULL;
  const mxArray *c9_lhs51 = NULL;
  const mxArray *c9_rhs52 = NULL;
  const mxArray *c9_lhs52 = NULL;
  const mxArray *c9_rhs53 = NULL;
  const mxArray *c9_lhs53 = NULL;
  const mxArray *c9_rhs54 = NULL;
  const mxArray *c9_lhs54 = NULL;
  const mxArray *c9_rhs55 = NULL;
  const mxArray *c9_lhs55 = NULL;
  const mxArray *c9_rhs56 = NULL;
  const mxArray *c9_lhs56 = NULL;
  const mxArray *c9_rhs57 = NULL;
  const mxArray *c9_lhs57 = NULL;
  const mxArray *c9_rhs58 = NULL;
  const mxArray *c9_lhs58 = NULL;
  const mxArray *c9_rhs59 = NULL;
  const mxArray *c9_lhs59 = NULL;
  const mxArray *c9_rhs60 = NULL;
  const mxArray *c9_lhs60 = NULL;
  const mxArray *c9_rhs61 = NULL;
  const mxArray *c9_lhs61 = NULL;
  const mxArray *c9_rhs62 = NULL;
  const mxArray *c9_lhs62 = NULL;
  const mxArray *c9_rhs63 = NULL;
  const mxArray *c9_lhs63 = NULL;
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mtimes"), "name", "name", 0);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c9_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 1);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c9_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("diag"), "name", "name", 2);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c9_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("ismatrix"), "name", "name", 3);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c9_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 4);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c9_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 5);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c9_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 6);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c9_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 7);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmax"), "name", "name", 7);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c9_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 8);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eye"), "name", "name", 8);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1368186630U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c9_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 9);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1368186630U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c9_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 10);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 10);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c9_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 11);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("isinf"), "name", "name", 11);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c9_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 12);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c9_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 13);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 13);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822382U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c9_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 14);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmax"), "name", "name", 14);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c9_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 15);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmin"), "name", "name", 15);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c9_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 16);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 16);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c9_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 17);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 17);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c9_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 18);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmin"), "name", "name", 18);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c9_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 19);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmin"), "name", "name", 19);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c9_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 20);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmax"), "name", "name", 20);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c9_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!numel_for_size"),
                  "context", "context", 21);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mtimes"), "name", "name", 21);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c9_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 22);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 22);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c9_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 23);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 23);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c9_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 24);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c9_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 25);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mldivide"), "name", "name", 25);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c9_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "context",
                  "context", 26);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 26);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c9_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "context",
                  "context", 27);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_lusolve"), "name", "name",
                  27);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1309454796U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c9_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 28);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c9_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 29);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 29);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c9_rhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 30);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xgetrf"), "name", "name",
                  30);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822406U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c9_rhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "context", "context", 31);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_lapack_xgetrf"), "name",
                  "name", 31);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822410U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c9_rhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "context", "context", 32);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_matlab_zgetrf"), "name",
                  "name", 32);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1302692594U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c9_rhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 33);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("realmin"), "name", "name", 33);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1307654842U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c9_rhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 34);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_realmin"), "name", "name",
                  34);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 34);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1307654844U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c9_rhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 35);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c9_rhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eps"), "name", "name", 36);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 36);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c9_rhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 37);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 37);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822382U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c9_rhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_eps"), "name", "name", 38);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 38);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c9_rhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 39);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 39);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c9_rhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 40);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("min"), "name", "name", 40);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 40);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1311258918U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c9_rhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 41);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c9_rhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 42);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 42);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 42);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c9_rhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 43);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 43);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 43);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1358189740U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c9_rhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 44);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 44);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c9_rhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 45);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 45);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 45);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 45);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c9_rhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 46);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 46);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 46);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c9_rhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("colon"), "name", "name", 47);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1366165842U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c9_rhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 48);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("colon"), "name", "name", 48);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 48);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1366165842U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c9_rhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 49);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 49);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c9_rhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 50);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 50);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 50);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c9_rhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 51);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("floor"), "name", "name", 51);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c9_rhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 52);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 52);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c9_rhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 53);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 53);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822326U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c9_rhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 54);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmin"), "name", "name", 54);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 54);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c9_rhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 55);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmax"), "name", "name", 55);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 55);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c9_rhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 56);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmin"), "name", "name", 56);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 56);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c9_rhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 57);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmax"), "name", "name", 57);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 57);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c9_rhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 58);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_isa_uint"), "name", "name",
                  58);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 58);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 58);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822384U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c9_rhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 59);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_unsigned_class"), "name",
                  "name", 59);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174180U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c9_rhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "context", "context", 60);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 60);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c9_rhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 61);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 61);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c9_rhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 62);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmax"), "name", "name", 62);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 62);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c9_rhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 63);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_isa_uint"), "name", "name",
                  63);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 63);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 63);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822384U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c9_rhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c9_rhs0);
  sf_mex_destroy(&c9_lhs0);
  sf_mex_destroy(&c9_rhs1);
  sf_mex_destroy(&c9_lhs1);
  sf_mex_destroy(&c9_rhs2);
  sf_mex_destroy(&c9_lhs2);
  sf_mex_destroy(&c9_rhs3);
  sf_mex_destroy(&c9_lhs3);
  sf_mex_destroy(&c9_rhs4);
  sf_mex_destroy(&c9_lhs4);
  sf_mex_destroy(&c9_rhs5);
  sf_mex_destroy(&c9_lhs5);
  sf_mex_destroy(&c9_rhs6);
  sf_mex_destroy(&c9_lhs6);
  sf_mex_destroy(&c9_rhs7);
  sf_mex_destroy(&c9_lhs7);
  sf_mex_destroy(&c9_rhs8);
  sf_mex_destroy(&c9_lhs8);
  sf_mex_destroy(&c9_rhs9);
  sf_mex_destroy(&c9_lhs9);
  sf_mex_destroy(&c9_rhs10);
  sf_mex_destroy(&c9_lhs10);
  sf_mex_destroy(&c9_rhs11);
  sf_mex_destroy(&c9_lhs11);
  sf_mex_destroy(&c9_rhs12);
  sf_mex_destroy(&c9_lhs12);
  sf_mex_destroy(&c9_rhs13);
  sf_mex_destroy(&c9_lhs13);
  sf_mex_destroy(&c9_rhs14);
  sf_mex_destroy(&c9_lhs14);
  sf_mex_destroy(&c9_rhs15);
  sf_mex_destroy(&c9_lhs15);
  sf_mex_destroy(&c9_rhs16);
  sf_mex_destroy(&c9_lhs16);
  sf_mex_destroy(&c9_rhs17);
  sf_mex_destroy(&c9_lhs17);
  sf_mex_destroy(&c9_rhs18);
  sf_mex_destroy(&c9_lhs18);
  sf_mex_destroy(&c9_rhs19);
  sf_mex_destroy(&c9_lhs19);
  sf_mex_destroy(&c9_rhs20);
  sf_mex_destroy(&c9_lhs20);
  sf_mex_destroy(&c9_rhs21);
  sf_mex_destroy(&c9_lhs21);
  sf_mex_destroy(&c9_rhs22);
  sf_mex_destroy(&c9_lhs22);
  sf_mex_destroy(&c9_rhs23);
  sf_mex_destroy(&c9_lhs23);
  sf_mex_destroy(&c9_rhs24);
  sf_mex_destroy(&c9_lhs24);
  sf_mex_destroy(&c9_rhs25);
  sf_mex_destroy(&c9_lhs25);
  sf_mex_destroy(&c9_rhs26);
  sf_mex_destroy(&c9_lhs26);
  sf_mex_destroy(&c9_rhs27);
  sf_mex_destroy(&c9_lhs27);
  sf_mex_destroy(&c9_rhs28);
  sf_mex_destroy(&c9_lhs28);
  sf_mex_destroy(&c9_rhs29);
  sf_mex_destroy(&c9_lhs29);
  sf_mex_destroy(&c9_rhs30);
  sf_mex_destroy(&c9_lhs30);
  sf_mex_destroy(&c9_rhs31);
  sf_mex_destroy(&c9_lhs31);
  sf_mex_destroy(&c9_rhs32);
  sf_mex_destroy(&c9_lhs32);
  sf_mex_destroy(&c9_rhs33);
  sf_mex_destroy(&c9_lhs33);
  sf_mex_destroy(&c9_rhs34);
  sf_mex_destroy(&c9_lhs34);
  sf_mex_destroy(&c9_rhs35);
  sf_mex_destroy(&c9_lhs35);
  sf_mex_destroy(&c9_rhs36);
  sf_mex_destroy(&c9_lhs36);
  sf_mex_destroy(&c9_rhs37);
  sf_mex_destroy(&c9_lhs37);
  sf_mex_destroy(&c9_rhs38);
  sf_mex_destroy(&c9_lhs38);
  sf_mex_destroy(&c9_rhs39);
  sf_mex_destroy(&c9_lhs39);
  sf_mex_destroy(&c9_rhs40);
  sf_mex_destroy(&c9_lhs40);
  sf_mex_destroy(&c9_rhs41);
  sf_mex_destroy(&c9_lhs41);
  sf_mex_destroy(&c9_rhs42);
  sf_mex_destroy(&c9_lhs42);
  sf_mex_destroy(&c9_rhs43);
  sf_mex_destroy(&c9_lhs43);
  sf_mex_destroy(&c9_rhs44);
  sf_mex_destroy(&c9_lhs44);
  sf_mex_destroy(&c9_rhs45);
  sf_mex_destroy(&c9_lhs45);
  sf_mex_destroy(&c9_rhs46);
  sf_mex_destroy(&c9_lhs46);
  sf_mex_destroy(&c9_rhs47);
  sf_mex_destroy(&c9_lhs47);
  sf_mex_destroy(&c9_rhs48);
  sf_mex_destroy(&c9_lhs48);
  sf_mex_destroy(&c9_rhs49);
  sf_mex_destroy(&c9_lhs49);
  sf_mex_destroy(&c9_rhs50);
  sf_mex_destroy(&c9_lhs50);
  sf_mex_destroy(&c9_rhs51);
  sf_mex_destroy(&c9_lhs51);
  sf_mex_destroy(&c9_rhs52);
  sf_mex_destroy(&c9_lhs52);
  sf_mex_destroy(&c9_rhs53);
  sf_mex_destroy(&c9_lhs53);
  sf_mex_destroy(&c9_rhs54);
  sf_mex_destroy(&c9_lhs54);
  sf_mex_destroy(&c9_rhs55);
  sf_mex_destroy(&c9_lhs55);
  sf_mex_destroy(&c9_rhs56);
  sf_mex_destroy(&c9_lhs56);
  sf_mex_destroy(&c9_rhs57);
  sf_mex_destroy(&c9_lhs57);
  sf_mex_destroy(&c9_rhs58);
  sf_mex_destroy(&c9_lhs58);
  sf_mex_destroy(&c9_rhs59);
  sf_mex_destroy(&c9_lhs59);
  sf_mex_destroy(&c9_rhs60);
  sf_mex_destroy(&c9_lhs60);
  sf_mex_destroy(&c9_rhs61);
  sf_mex_destroy(&c9_lhs61);
  sf_mex_destroy(&c9_rhs62);
  sf_mex_destroy(&c9_lhs62);
  sf_mex_destroy(&c9_rhs63);
  sf_mex_destroy(&c9_lhs63);
}

static const mxArray *c9_emlrt_marshallOut(char * c9_u)
{
  const mxArray *c9_y = NULL;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c9_u)), FALSE);
  return c9_y;
}

static const mxArray *c9_b_emlrt_marshallOut(uint32_T c9_u)
{
  const mxArray *c9_y = NULL;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", &c9_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c9_y;
}

static void c9_b_info_helper(const mxArray **c9_info)
{
  const mxArray *c9_rhs64 = NULL;
  const mxArray *c9_lhs64 = NULL;
  const mxArray *c9_rhs65 = NULL;
  const mxArray *c9_lhs65 = NULL;
  const mxArray *c9_rhs66 = NULL;
  const mxArray *c9_lhs66 = NULL;
  const mxArray *c9_rhs67 = NULL;
  const mxArray *c9_lhs67 = NULL;
  const mxArray *c9_rhs68 = NULL;
  const mxArray *c9_lhs68 = NULL;
  const mxArray *c9_rhs69 = NULL;
  const mxArray *c9_lhs69 = NULL;
  const mxArray *c9_rhs70 = NULL;
  const mxArray *c9_lhs70 = NULL;
  const mxArray *c9_rhs71 = NULL;
  const mxArray *c9_lhs71 = NULL;
  const mxArray *c9_rhs72 = NULL;
  const mxArray *c9_lhs72 = NULL;
  const mxArray *c9_rhs73 = NULL;
  const mxArray *c9_lhs73 = NULL;
  const mxArray *c9_rhs74 = NULL;
  const mxArray *c9_lhs74 = NULL;
  const mxArray *c9_rhs75 = NULL;
  const mxArray *c9_lhs75 = NULL;
  const mxArray *c9_rhs76 = NULL;
  const mxArray *c9_lhs76 = NULL;
  const mxArray *c9_rhs77 = NULL;
  const mxArray *c9_lhs77 = NULL;
  const mxArray *c9_rhs78 = NULL;
  const mxArray *c9_lhs78 = NULL;
  const mxArray *c9_rhs79 = NULL;
  const mxArray *c9_lhs79 = NULL;
  const mxArray *c9_rhs80 = NULL;
  const mxArray *c9_lhs80 = NULL;
  const mxArray *c9_rhs81 = NULL;
  const mxArray *c9_lhs81 = NULL;
  const mxArray *c9_rhs82 = NULL;
  const mxArray *c9_lhs82 = NULL;
  const mxArray *c9_rhs83 = NULL;
  const mxArray *c9_lhs83 = NULL;
  const mxArray *c9_rhs84 = NULL;
  const mxArray *c9_lhs84 = NULL;
  const mxArray *c9_rhs85 = NULL;
  const mxArray *c9_lhs85 = NULL;
  const mxArray *c9_rhs86 = NULL;
  const mxArray *c9_lhs86 = NULL;
  const mxArray *c9_rhs87 = NULL;
  const mxArray *c9_lhs87 = NULL;
  const mxArray *c9_rhs88 = NULL;
  const mxArray *c9_lhs88 = NULL;
  const mxArray *c9_rhs89 = NULL;
  const mxArray *c9_lhs89 = NULL;
  const mxArray *c9_rhs90 = NULL;
  const mxArray *c9_lhs90 = NULL;
  const mxArray *c9_rhs91 = NULL;
  const mxArray *c9_lhs91 = NULL;
  const mxArray *c9_rhs92 = NULL;
  const mxArray *c9_lhs92 = NULL;
  const mxArray *c9_rhs93 = NULL;
  const mxArray *c9_lhs93 = NULL;
  const mxArray *c9_rhs94 = NULL;
  const mxArray *c9_lhs94 = NULL;
  const mxArray *c9_rhs95 = NULL;
  const mxArray *c9_lhs95 = NULL;
  const mxArray *c9_rhs96 = NULL;
  const mxArray *c9_lhs96 = NULL;
  const mxArray *c9_rhs97 = NULL;
  const mxArray *c9_lhs97 = NULL;
  const mxArray *c9_rhs98 = NULL;
  const mxArray *c9_lhs98 = NULL;
  const mxArray *c9_rhs99 = NULL;
  const mxArray *c9_lhs99 = NULL;
  const mxArray *c9_rhs100 = NULL;
  const mxArray *c9_lhs100 = NULL;
  const mxArray *c9_rhs101 = NULL;
  const mxArray *c9_lhs101 = NULL;
  const mxArray *c9_rhs102 = NULL;
  const mxArray *c9_lhs102 = NULL;
  const mxArray *c9_rhs103 = NULL;
  const mxArray *c9_lhs103 = NULL;
  const mxArray *c9_rhs104 = NULL;
  const mxArray *c9_lhs104 = NULL;
  const mxArray *c9_rhs105 = NULL;
  const mxArray *c9_lhs105 = NULL;
  const mxArray *c9_rhs106 = NULL;
  const mxArray *c9_lhs106 = NULL;
  const mxArray *c9_rhs107 = NULL;
  const mxArray *c9_lhs107 = NULL;
  const mxArray *c9_rhs108 = NULL;
  const mxArray *c9_lhs108 = NULL;
  const mxArray *c9_rhs109 = NULL;
  const mxArray *c9_lhs109 = NULL;
  const mxArray *c9_rhs110 = NULL;
  const mxArray *c9_lhs110 = NULL;
  const mxArray *c9_rhs111 = NULL;
  const mxArray *c9_lhs111 = NULL;
  const mxArray *c9_rhs112 = NULL;
  const mxArray *c9_lhs112 = NULL;
  const mxArray *c9_rhs113 = NULL;
  const mxArray *c9_lhs113 = NULL;
  const mxArray *c9_rhs114 = NULL;
  const mxArray *c9_lhs114 = NULL;
  const mxArray *c9_rhs115 = NULL;
  const mxArray *c9_lhs115 = NULL;
  const mxArray *c9_rhs116 = NULL;
  const mxArray *c9_lhs116 = NULL;
  const mxArray *c9_rhs117 = NULL;
  const mxArray *c9_lhs117 = NULL;
  const mxArray *c9_rhs118 = NULL;
  const mxArray *c9_lhs118 = NULL;
  const mxArray *c9_rhs119 = NULL;
  const mxArray *c9_lhs119 = NULL;
  const mxArray *c9_rhs120 = NULL;
  const mxArray *c9_lhs120 = NULL;
  const mxArray *c9_rhs121 = NULL;
  const mxArray *c9_lhs121 = NULL;
  const mxArray *c9_rhs122 = NULL;
  const mxArray *c9_lhs122 = NULL;
  const mxArray *c9_rhs123 = NULL;
  const mxArray *c9_lhs123 = NULL;
  const mxArray *c9_rhs124 = NULL;
  const mxArray *c9_lhs124 = NULL;
  const mxArray *c9_rhs125 = NULL;
  const mxArray *c9_lhs125 = NULL;
  const mxArray *c9_rhs126 = NULL;
  const mxArray *c9_lhs126 = NULL;
  const mxArray *c9_rhs127 = NULL;
  const mxArray *c9_lhs127 = NULL;
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 64);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 64);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c9_rhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 65);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 65);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c9_rhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon"),
                  "context", "context", 66);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 66);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c9_rhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 67);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 67);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c9_rhs67, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs67, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 68);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 68);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c9_rhs68, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs68, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 69);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 69);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c9_rhs69, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs69, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 70);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 70);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c9_rhs70, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs70, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 71);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 71);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c9_rhs71, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs71, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 72);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 72);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 72);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c9_rhs72, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs72, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 73);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 73);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 73);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c9_rhs73, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs73, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 74);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 74);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c9_rhs74, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs74, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 75);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 75);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 75);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c9_rhs75, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs75, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 76);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_ixamax"), "name", "name",
                  76);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c9_rhs76, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs76, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 77);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 77);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 77);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c9_rhs77, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs77, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 78);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 78);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c9_rhs78, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs78, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m"),
                  "context", "context", 79);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_ixamax"), "name",
                  "name", 79);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080370U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c9_rhs79, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs79, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "context", "context", 80);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 80);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c9_rhs80, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs80, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "context", "context", 81);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xcabs1"), "name", "name",
                  81);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "resolved", "resolved", 81);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822306U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c9_rhs81, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs81, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "context", "context", 82);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("abs"), "name", "name", 82);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 82);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c9_rhs82, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs82, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 83);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 83);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 83);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c9_rhs83, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs83, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 84);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 84);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 84);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822312U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c9_rhs84, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs84, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "context", "context", 85);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 85);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c9_rhs85, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs85, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "context", "context", 86);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 86);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 86);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 86);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c9_rhs86, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs86, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 87);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xswap"), "name", "name",
                  87);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717474U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c9_rhs87, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs87, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 88);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 88);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c9_rhs88, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs88, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m"),
                  "context", "context", 89);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 89);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c9_rhs89, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs89, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m"),
                  "context", "context", 90);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_xswap"), "name",
                  "name", 90);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080386U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c9_rhs90, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs90, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 91);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 91);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c9_rhs91, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs91, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 92);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("abs"), "name", "name", 92);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 92);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 92);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c9_rhs92, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs92, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 93);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 93);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 93);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 93);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c9_rhs93, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs93, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 94);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 94);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 94);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 94);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822312U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c9_rhs94, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs94, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 95);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 95);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 95);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c9_rhs95, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs95, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 96);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 96);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 96);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c9_rhs96, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs96, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 97);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 97);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 97);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 97);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c9_rhs97, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs97, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 98);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_div"), "name", "name", 98);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 98);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c9_rhs98, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs98, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 99);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xgeru"), "name", "name",
                  99);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"),
                  "resolved", "resolved", 99);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717472U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c9_rhs99, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs99, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 100);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 100);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c9_rhs100, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs100, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs100), "rhs", "rhs",
                  100);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs100), "lhs", "lhs",
                  100);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 101);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xger"), "name", "name",
                  101);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 101);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m"), "resolved",
                  "resolved", 101);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c9_rhs101, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs101, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs101), "rhs", "rhs",
                  101);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs101), "lhs", "lhs",
                  101);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m"), "context",
                  "context", 102);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 102);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 102);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c9_rhs102, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs102, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs102), "rhs", "rhs",
                  102);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs102), "lhs", "lhs",
                  102);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold"),
                  "context", "context", 103);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmax"), "name", "name", 103);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 103);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 103);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c9_rhs103, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs103, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs103), "rhs", "rhs",
                  103);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs103), "lhs", "lhs",
                  103);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold"),
                  "context", "context", 104);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("min"), "name", "name", 104);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 104);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1311258918U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c9_rhs104, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs104, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs104), "rhs", "rhs",
                  104);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs104), "lhs", "lhs",
                  104);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 105);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 105);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 105);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c9_rhs105, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs105, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs105), "rhs", "rhs",
                  105);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs105), "lhs", "lhs",
                  105);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 106);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 106);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 106);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 106);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1358189740U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c9_rhs106, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs106, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs106), "rhs", "rhs",
                  106);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs106), "lhs", "lhs",
                  106);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 107);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 107);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 107);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c9_rhs107, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs107, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs107), "rhs", "rhs",
                  107);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs107), "lhs", "lhs",
                  107);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 108);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 108);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c9_rhs108, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs108, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs108), "rhs", "rhs",
                  108);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs108), "lhs", "lhs",
                  108);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold"),
                  "context", "context", 109);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mtimes"), "name", "name", 109);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 109);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 109);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c9_rhs109, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs109, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs109), "rhs", "rhs",
                  109);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs109), "lhs", "lhs",
                  109);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m"),
                  "context", "context", 110);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 110);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 110);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c9_rhs110, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs110, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs110), "rhs", "rhs",
                  110);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs110), "lhs", "lhs",
                  110);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m"),
                  "context", "context", 111);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_xger"), "name",
                  "name", 111);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 111);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xger.m"),
                  "resolved", "resolved", 111);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080376U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c9_rhs111, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs111, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs111), "rhs", "rhs",
                  111);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs111), "lhs", "lhs",
                  111);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xger.m"),
                  "context", "context", 112);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_xgerx"), "name",
                  "name", 112);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "resolved", "resolved", 112);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c9_rhs112, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs112, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs112), "rhs", "rhs",
                  112);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs112), "lhs", "lhs",
                  112);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 113);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 113);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c9_rhs113, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs113, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs113), "rhs", "rhs",
                  113);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs113), "lhs", "lhs",
                  113);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 114);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("abs"), "name", "name", 114);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 114);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 114);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c9_rhs114, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs114, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs114), "rhs", "rhs",
                  114);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs114), "lhs", "lhs",
                  114);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 115);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 115);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c9_rhs115, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs115, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs115), "rhs", "rhs",
                  115);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs115), "lhs", "lhs",
                  115);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 116);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 116);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 116);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c9_rhs116, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs116, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs116), "rhs", "rhs",
                  116);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs116), "lhs", "lhs",
                  116);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 117);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 117);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 117);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c9_rhs117, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs117, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs117), "rhs", "rhs",
                  117);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs117), "lhs", "lhs",
                  117);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 118);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 118);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 118);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 118);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c9_rhs118, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs118, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs118), "rhs", "rhs",
                  118);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs118), "lhs", "lhs",
                  118);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 119);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 119);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 119);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 119);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c9_rhs119, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs119, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs119), "rhs", "rhs",
                  119);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs119), "lhs", "lhs",
                  119);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular"),
                  "context", "context", 120);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_warning"), "name", "name",
                  120);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 120);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822402U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c9_rhs120, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs120, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs120), "rhs", "rhs",
                  120);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs120), "lhs", "lhs",
                  120);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 121);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 121);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 121);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c9_rhs121, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs121, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs121), "rhs", "rhs",
                  121);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs121), "lhs", "lhs",
                  121);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 122);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 122);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 122);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c9_rhs122, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs122, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs122), "rhs", "rhs",
                  122);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs122), "lhs", "lhs",
                  122);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 123);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xtrsm"), "name", "name",
                  123);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"),
                  "resolved", "resolved", 123);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080378U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c9_rhs123, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs123, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs123), "rhs", "rhs",
                  123);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs123), "lhs", "lhs",
                  123);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 124);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 124);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 124);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c9_rhs124, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs124, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs124), "rhs", "rhs",
                  124);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs124), "lhs", "lhs",
                  124);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"),
                  "context", "context", 125);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 125);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c9_rhs125, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs125, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs125), "rhs", "rhs",
                  125);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs125), "lhs", "lhs",
                  125);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"),
                  "context", "context", 126);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 126);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 126);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c9_rhs126, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs126, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs126), "rhs", "rhs",
                  126);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs126), "lhs", "lhs",
                  126);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"),
                  "context", "context", 127);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_xtrsm"), "name",
                  "name", 127);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "resolved", "resolved", 127);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c9_rhs127, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs127, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs127), "rhs", "rhs",
                  127);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs127), "lhs", "lhs",
                  127);
  sf_mex_destroy(&c9_rhs64);
  sf_mex_destroy(&c9_lhs64);
  sf_mex_destroy(&c9_rhs65);
  sf_mex_destroy(&c9_lhs65);
  sf_mex_destroy(&c9_rhs66);
  sf_mex_destroy(&c9_lhs66);
  sf_mex_destroy(&c9_rhs67);
  sf_mex_destroy(&c9_lhs67);
  sf_mex_destroy(&c9_rhs68);
  sf_mex_destroy(&c9_lhs68);
  sf_mex_destroy(&c9_rhs69);
  sf_mex_destroy(&c9_lhs69);
  sf_mex_destroy(&c9_rhs70);
  sf_mex_destroy(&c9_lhs70);
  sf_mex_destroy(&c9_rhs71);
  sf_mex_destroy(&c9_lhs71);
  sf_mex_destroy(&c9_rhs72);
  sf_mex_destroy(&c9_lhs72);
  sf_mex_destroy(&c9_rhs73);
  sf_mex_destroy(&c9_lhs73);
  sf_mex_destroy(&c9_rhs74);
  sf_mex_destroy(&c9_lhs74);
  sf_mex_destroy(&c9_rhs75);
  sf_mex_destroy(&c9_lhs75);
  sf_mex_destroy(&c9_rhs76);
  sf_mex_destroy(&c9_lhs76);
  sf_mex_destroy(&c9_rhs77);
  sf_mex_destroy(&c9_lhs77);
  sf_mex_destroy(&c9_rhs78);
  sf_mex_destroy(&c9_lhs78);
  sf_mex_destroy(&c9_rhs79);
  sf_mex_destroy(&c9_lhs79);
  sf_mex_destroy(&c9_rhs80);
  sf_mex_destroy(&c9_lhs80);
  sf_mex_destroy(&c9_rhs81);
  sf_mex_destroy(&c9_lhs81);
  sf_mex_destroy(&c9_rhs82);
  sf_mex_destroy(&c9_lhs82);
  sf_mex_destroy(&c9_rhs83);
  sf_mex_destroy(&c9_lhs83);
  sf_mex_destroy(&c9_rhs84);
  sf_mex_destroy(&c9_lhs84);
  sf_mex_destroy(&c9_rhs85);
  sf_mex_destroy(&c9_lhs85);
  sf_mex_destroy(&c9_rhs86);
  sf_mex_destroy(&c9_lhs86);
  sf_mex_destroy(&c9_rhs87);
  sf_mex_destroy(&c9_lhs87);
  sf_mex_destroy(&c9_rhs88);
  sf_mex_destroy(&c9_lhs88);
  sf_mex_destroy(&c9_rhs89);
  sf_mex_destroy(&c9_lhs89);
  sf_mex_destroy(&c9_rhs90);
  sf_mex_destroy(&c9_lhs90);
  sf_mex_destroy(&c9_rhs91);
  sf_mex_destroy(&c9_lhs91);
  sf_mex_destroy(&c9_rhs92);
  sf_mex_destroy(&c9_lhs92);
  sf_mex_destroy(&c9_rhs93);
  sf_mex_destroy(&c9_lhs93);
  sf_mex_destroy(&c9_rhs94);
  sf_mex_destroy(&c9_lhs94);
  sf_mex_destroy(&c9_rhs95);
  sf_mex_destroy(&c9_lhs95);
  sf_mex_destroy(&c9_rhs96);
  sf_mex_destroy(&c9_lhs96);
  sf_mex_destroy(&c9_rhs97);
  sf_mex_destroy(&c9_lhs97);
  sf_mex_destroy(&c9_rhs98);
  sf_mex_destroy(&c9_lhs98);
  sf_mex_destroy(&c9_rhs99);
  sf_mex_destroy(&c9_lhs99);
  sf_mex_destroy(&c9_rhs100);
  sf_mex_destroy(&c9_lhs100);
  sf_mex_destroy(&c9_rhs101);
  sf_mex_destroy(&c9_lhs101);
  sf_mex_destroy(&c9_rhs102);
  sf_mex_destroy(&c9_lhs102);
  sf_mex_destroy(&c9_rhs103);
  sf_mex_destroy(&c9_lhs103);
  sf_mex_destroy(&c9_rhs104);
  sf_mex_destroy(&c9_lhs104);
  sf_mex_destroy(&c9_rhs105);
  sf_mex_destroy(&c9_lhs105);
  sf_mex_destroy(&c9_rhs106);
  sf_mex_destroy(&c9_lhs106);
  sf_mex_destroy(&c9_rhs107);
  sf_mex_destroy(&c9_lhs107);
  sf_mex_destroy(&c9_rhs108);
  sf_mex_destroy(&c9_lhs108);
  sf_mex_destroy(&c9_rhs109);
  sf_mex_destroy(&c9_lhs109);
  sf_mex_destroy(&c9_rhs110);
  sf_mex_destroy(&c9_lhs110);
  sf_mex_destroy(&c9_rhs111);
  sf_mex_destroy(&c9_lhs111);
  sf_mex_destroy(&c9_rhs112);
  sf_mex_destroy(&c9_lhs112);
  sf_mex_destroy(&c9_rhs113);
  sf_mex_destroy(&c9_lhs113);
  sf_mex_destroy(&c9_rhs114);
  sf_mex_destroy(&c9_lhs114);
  sf_mex_destroy(&c9_rhs115);
  sf_mex_destroy(&c9_lhs115);
  sf_mex_destroy(&c9_rhs116);
  sf_mex_destroy(&c9_lhs116);
  sf_mex_destroy(&c9_rhs117);
  sf_mex_destroy(&c9_lhs117);
  sf_mex_destroy(&c9_rhs118);
  sf_mex_destroy(&c9_lhs118);
  sf_mex_destroy(&c9_rhs119);
  sf_mex_destroy(&c9_lhs119);
  sf_mex_destroy(&c9_rhs120);
  sf_mex_destroy(&c9_lhs120);
  sf_mex_destroy(&c9_rhs121);
  sf_mex_destroy(&c9_lhs121);
  sf_mex_destroy(&c9_rhs122);
  sf_mex_destroy(&c9_lhs122);
  sf_mex_destroy(&c9_rhs123);
  sf_mex_destroy(&c9_lhs123);
  sf_mex_destroy(&c9_rhs124);
  sf_mex_destroy(&c9_lhs124);
  sf_mex_destroy(&c9_rhs125);
  sf_mex_destroy(&c9_lhs125);
  sf_mex_destroy(&c9_rhs126);
  sf_mex_destroy(&c9_lhs126);
  sf_mex_destroy(&c9_rhs127);
  sf_mex_destroy(&c9_lhs127);
}

static void c9_c_info_helper(const mxArray **c9_info)
{
  const mxArray *c9_rhs128 = NULL;
  const mxArray *c9_lhs128 = NULL;
  const mxArray *c9_rhs129 = NULL;
  const mxArray *c9_lhs129 = NULL;
  const mxArray *c9_rhs130 = NULL;
  const mxArray *c9_lhs130 = NULL;
  const mxArray *c9_rhs131 = NULL;
  const mxArray *c9_lhs131 = NULL;
  const mxArray *c9_rhs132 = NULL;
  const mxArray *c9_lhs132 = NULL;
  const mxArray *c9_rhs133 = NULL;
  const mxArray *c9_lhs133 = NULL;
  const mxArray *c9_rhs134 = NULL;
  const mxArray *c9_lhs134 = NULL;
  const mxArray *c9_rhs135 = NULL;
  const mxArray *c9_lhs135 = NULL;
  const mxArray *c9_rhs136 = NULL;
  const mxArray *c9_lhs136 = NULL;
  const mxArray *c9_rhs137 = NULL;
  const mxArray *c9_lhs137 = NULL;
  const mxArray *c9_rhs138 = NULL;
  const mxArray *c9_lhs138 = NULL;
  const mxArray *c9_rhs139 = NULL;
  const mxArray *c9_lhs139 = NULL;
  const mxArray *c9_rhs140 = NULL;
  const mxArray *c9_lhs140 = NULL;
  const mxArray *c9_rhs141 = NULL;
  const mxArray *c9_lhs141 = NULL;
  const mxArray *c9_rhs142 = NULL;
  const mxArray *c9_lhs142 = NULL;
  const mxArray *c9_rhs143 = NULL;
  const mxArray *c9_lhs143 = NULL;
  const mxArray *c9_rhs144 = NULL;
  const mxArray *c9_lhs144 = NULL;
  const mxArray *c9_rhs145 = NULL;
  const mxArray *c9_lhs145 = NULL;
  const mxArray *c9_rhs146 = NULL;
  const mxArray *c9_lhs146 = NULL;
  const mxArray *c9_rhs147 = NULL;
  const mxArray *c9_lhs147 = NULL;
  const mxArray *c9_rhs148 = NULL;
  const mxArray *c9_lhs148 = NULL;
  const mxArray *c9_rhs149 = NULL;
  const mxArray *c9_lhs149 = NULL;
  const mxArray *c9_rhs150 = NULL;
  const mxArray *c9_lhs150 = NULL;
  const mxArray *c9_rhs151 = NULL;
  const mxArray *c9_lhs151 = NULL;
  const mxArray *c9_rhs152 = NULL;
  const mxArray *c9_lhs152 = NULL;
  const mxArray *c9_rhs153 = NULL;
  const mxArray *c9_lhs153 = NULL;
  const mxArray *c9_rhs154 = NULL;
  const mxArray *c9_lhs154 = NULL;
  const mxArray *c9_rhs155 = NULL;
  const mxArray *c9_lhs155 = NULL;
  const mxArray *c9_rhs156 = NULL;
  const mxArray *c9_lhs156 = NULL;
  const mxArray *c9_rhs157 = NULL;
  const mxArray *c9_lhs157 = NULL;
  const mxArray *c9_rhs158 = NULL;
  const mxArray *c9_lhs158 = NULL;
  const mxArray *c9_rhs159 = NULL;
  const mxArray *c9_lhs159 = NULL;
  const mxArray *c9_rhs160 = NULL;
  const mxArray *c9_lhs160 = NULL;
  const mxArray *c9_rhs161 = NULL;
  const mxArray *c9_lhs161 = NULL;
  const mxArray *c9_rhs162 = NULL;
  const mxArray *c9_lhs162 = NULL;
  const mxArray *c9_rhs163 = NULL;
  const mxArray *c9_lhs163 = NULL;
  const mxArray *c9_rhs164 = NULL;
  const mxArray *c9_lhs164 = NULL;
  const mxArray *c9_rhs165 = NULL;
  const mxArray *c9_lhs165 = NULL;
  const mxArray *c9_rhs166 = NULL;
  const mxArray *c9_lhs166 = NULL;
  const mxArray *c9_rhs167 = NULL;
  const mxArray *c9_lhs167 = NULL;
  const mxArray *c9_rhs168 = NULL;
  const mxArray *c9_lhs168 = NULL;
  const mxArray *c9_rhs169 = NULL;
  const mxArray *c9_lhs169 = NULL;
  const mxArray *c9_rhs170 = NULL;
  const mxArray *c9_lhs170 = NULL;
  const mxArray *c9_rhs171 = NULL;
  const mxArray *c9_lhs171 = NULL;
  const mxArray *c9_rhs172 = NULL;
  const mxArray *c9_lhs172 = NULL;
  const mxArray *c9_rhs173 = NULL;
  const mxArray *c9_lhs173 = NULL;
  const mxArray *c9_rhs174 = NULL;
  const mxArray *c9_lhs174 = NULL;
  const mxArray *c9_rhs175 = NULL;
  const mxArray *c9_lhs175 = NULL;
  const mxArray *c9_rhs176 = NULL;
  const mxArray *c9_lhs176 = NULL;
  const mxArray *c9_rhs177 = NULL;
  const mxArray *c9_lhs177 = NULL;
  const mxArray *c9_rhs178 = NULL;
  const mxArray *c9_lhs178 = NULL;
  const mxArray *c9_rhs179 = NULL;
  const mxArray *c9_lhs179 = NULL;
  const mxArray *c9_rhs180 = NULL;
  const mxArray *c9_lhs180 = NULL;
  const mxArray *c9_rhs181 = NULL;
  const mxArray *c9_lhs181 = NULL;
  const mxArray *c9_rhs182 = NULL;
  const mxArray *c9_lhs182 = NULL;
  const mxArray *c9_rhs183 = NULL;
  const mxArray *c9_lhs183 = NULL;
  const mxArray *c9_rhs184 = NULL;
  const mxArray *c9_lhs184 = NULL;
  const mxArray *c9_rhs185 = NULL;
  const mxArray *c9_lhs185 = NULL;
  const mxArray *c9_rhs186 = NULL;
  const mxArray *c9_lhs186 = NULL;
  const mxArray *c9_rhs187 = NULL;
  const mxArray *c9_lhs187 = NULL;
  const mxArray *c9_rhs188 = NULL;
  const mxArray *c9_lhs188 = NULL;
  const mxArray *c9_rhs189 = NULL;
  const mxArray *c9_lhs189 = NULL;
  const mxArray *c9_rhs190 = NULL;
  const mxArray *c9_lhs190 = NULL;
  const mxArray *c9_rhs191 = NULL;
  const mxArray *c9_lhs191 = NULL;
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 128);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 128);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 128);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 128);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c9_rhs128, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs128, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs128), "rhs", "rhs",
                  128);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs128), "lhs", "lhs",
                  128);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 129);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 129);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 129);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 129);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c9_rhs129, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs129, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs129), "rhs", "rhs",
                  129);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs129), "lhs", "lhs",
                  129);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 130);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 130);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 130);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c9_rhs130, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs130, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs130), "rhs", "rhs",
                  130);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs130), "lhs", "lhs",
                  130);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 131);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 131);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 131);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c9_rhs131, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs131, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs131), "rhs", "rhs",
                  131);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs131), "lhs", "lhs",
                  131);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 132);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 132);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 132);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c9_rhs132, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs132, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs132), "rhs", "rhs",
                  132);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs132), "lhs", "lhs",
                  132);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 133);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 133);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 133);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 133);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c9_rhs133, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs133, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs133), "rhs", "rhs",
                  133);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs133), "lhs", "lhs",
                  133);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 134);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 134);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 134);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c9_rhs134, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs134, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs134), "rhs", "rhs",
                  134);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs134), "lhs", "lhs",
                  134);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 135);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("intmin"), "name", "name", 135);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 135);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c9_rhs135, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs135, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs135), "rhs", "rhs",
                  135);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs135), "lhs", "lhs",
                  135);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 136);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_div"), "name", "name", 136);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 136);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c9_rhs136, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs136, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs136), "rhs", "rhs",
                  136);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs136), "lhs", "lhs",
                  136);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "context",
                  "context", 137);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_qrsolve"), "name", "name",
                  137);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "resolved",
                  "resolved", 137);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c9_rhs137, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs137, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs137), "rhs", "rhs",
                  137);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs137), "lhs", "lhs",
                  137);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 138);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("min"), "name", "name", 138);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 138);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1311258918U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c9_rhs138, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs138, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs138), "rhs", "rhs",
                  138);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs138), "lhs", "lhs",
                  138);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 139);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xgeqp3"), "name", "name",
                  139);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeqp3.m"),
                  "resolved", "resolved", 139);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822404U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c9_rhs139, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs139, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs139), "rhs", "rhs",
                  139);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs139), "lhs", "lhs",
                  139);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgeqp3.m"),
                  "context", "context", 140);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_lapack_xgeqp3"), "name",
                  "name", 140);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeqp3.m"),
                  "resolved", "resolved", 140);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822408U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c9_rhs140, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs140, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs140), "rhs", "rhs",
                  140);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs140), "lhs", "lhs",
                  140);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeqp3.m"),
                  "context", "context", 141);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_matlab_zgeqp3"), "name",
                  "name", 141);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "resolved", "resolved", 141);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1290002966U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c9_rhs141, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs141, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs141), "rhs", "rhs",
                  141);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs141), "lhs", "lhs",
                  141);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 142);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 142);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 142);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c9_rhs142, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs142, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs142), "rhs", "rhs",
                  142);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs142), "lhs", "lhs",
                  142);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 143);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("min"), "name", "name", 143);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 143);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 143);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1311258918U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c9_rhs143, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs143, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs143), "rhs", "rhs",
                  143);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs143), "lhs", "lhs",
                  143);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 144);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 144);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 144);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c9_rhs144, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs144, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs144), "rhs", "rhs",
                  144);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs144), "lhs", "lhs",
                  144);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 145);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("colon"), "name", "name", 145);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 145);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 145);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1366165842U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c9_rhs145, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs145, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs145), "rhs", "rhs",
                  145);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs145), "lhs", "lhs",
                  145);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 146);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eps"), "name", "name", 146);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 146);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c9_rhs146, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs146, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs146), "rhs", "rhs",
                  146);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs146), "lhs", "lhs",
                  146);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 147);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("sqrt"), "name", "name", 147);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 147);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c9_rhs147, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs147, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs147), "rhs", "rhs",
                  147);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs147), "lhs", "lhs",
                  147);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 148);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_error"), "name", "name",
                  148);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 148);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1343833958U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c9_rhs148, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs148, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs148), "rhs", "rhs",
                  148);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs148), "lhs", "lhs",
                  148);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 149);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 149);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 149);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822338U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c9_rhs149, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs149, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs149), "rhs", "rhs",
                  149);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs149), "lhs", "lhs",
                  149);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 150);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 150);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 150);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c9_rhs150, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs150, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs150), "rhs", "rhs",
                  150);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs150), "lhs", "lhs",
                  150);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 151);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  151);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 151);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717472U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c9_rhs151, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs151, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs151), "rhs", "rhs",
                  151);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs151), "lhs", "lhs",
                  151);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 152);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 152);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 152);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c9_rhs152, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs152, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs152), "rhs", "rhs",
                  152);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs152), "lhs", "lhs",
                  152);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m"),
                  "context", "context", 153);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 153);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 153);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c9_rhs153, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs153, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs153), "rhs", "rhs",
                  153);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs153), "lhs", "lhs",
                  153);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m"),
                  "context", "context", 154);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_xnrm2"), "name",
                  "name", 154);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m"),
                  "resolved", "resolved", 154);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080384U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c9_rhs154, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs154, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs154), "rhs", "rhs",
                  154);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs154), "lhs", "lhs",
                  154);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m"),
                  "context", "context", 155);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("abs"), "name", "name", 155);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 155);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c9_rhs155, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs155, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs155), "rhs", "rhs",
                  155);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs155), "lhs", "lhs",
                  155);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m"),
                  "context", "context", 156);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("realmin"), "name", "name", 156);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 156);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1307654842U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c9_rhs156, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs156, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs156), "rhs", "rhs",
                  156);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs156), "lhs", "lhs",
                  156);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m"),
                  "context", "context", 157);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 157);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 157);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c9_rhs157, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs157, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs157), "rhs", "rhs",
                  157);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs157), "lhs", "lhs",
                  157);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m"),
                  "context", "context", 158);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 158);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c9_rhs158, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs158, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs158), "rhs", "rhs",
                  158);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs158), "lhs", "lhs",
                  158);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m"),
                  "context", "context", 159);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 159);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 159);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 159);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c9_rhs159, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs159, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs159), "rhs", "rhs",
                  159);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs159), "lhs", "lhs",
                  159);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m"),
                  "context", "context", 160);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 160);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 160);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 160);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c9_rhs160, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs160, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs160), "rhs", "rhs",
                  160);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs160), "lhs", "lhs",
                  160);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xnrm2.m"),
                  "context", "context", 161);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 161);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 161);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 161);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c9_rhs161, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs161, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs161), "rhs", "rhs",
                  161);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs161), "lhs", "lhs",
                  161);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 162);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 162);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 162);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 162);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c9_rhs162, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs162, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs162), "rhs", "rhs",
                  162);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs162), "lhs", "lhs",
                  162);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 163);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 163);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c9_rhs163, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs163, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs163), "rhs", "rhs",
                  163);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs163), "lhs", "lhs",
                  163);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 164);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 164);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 164);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 164);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c9_rhs164, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs164, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs164), "rhs", "rhs",
                  164);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs164), "lhs", "lhs",
                  164);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 165);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 165);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 165);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 165);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c9_rhs165, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs165, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs165), "rhs", "rhs",
                  165);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs165), "lhs", "lhs",
                  165);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 166);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 166);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 166);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 166);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c9_rhs166, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs166, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs166), "rhs", "rhs",
                  166);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs166), "lhs", "lhs",
                  166);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 167);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_ixamax"), "name", "name",
                  167);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "resolved", "resolved", 167);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c9_rhs167, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs167, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs167), "rhs", "rhs",
                  167);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs167), "lhs", "lhs",
                  167);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 168);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xswap"), "name", "name",
                  168);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"),
                  "resolved", "resolved", 168);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717474U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c9_rhs168, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs168, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs168), "rhs", "rhs",
                  168);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs168), "lhs", "lhs",
                  168);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 169);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_matlab_zlarfg"), "name",
                  "name", 169);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "resolved", "resolved", 169);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822422U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c9_rhs169, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs169, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs169), "rhs", "rhs",
                  169);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs169), "lhs", "lhs",
                  169);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 170);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 170);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 170);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c9_rhs170, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs170, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs170), "rhs", "rhs",
                  170);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs170), "lhs", "lhs",
                  170);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 171);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  171);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 171);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717472U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c9_rhs171, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs171, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs171), "rhs", "rhs",
                  171);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs171), "lhs", "lhs",
                  171);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 172);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_dlapy2"), "name", "name",
                  172);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_dlapy2.m"), "resolved",
                  "resolved", 172);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1350414254U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c9_rhs172, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs172, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs172), "rhs", "rhs",
                  172);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs172), "lhs", "lhs",
                  172);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 173);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("realmin"), "name", "name", 173);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 173);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 173);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1307654842U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c9_rhs173, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs173, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs173), "rhs", "rhs",
                  173);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs173), "lhs", "lhs",
                  173);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 174);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eps"), "name", "name", 174);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 174);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 174);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c9_rhs174, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs174, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs174), "rhs", "rhs",
                  174);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs174), "lhs", "lhs",
                  174);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 175);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("abs"), "name", "name", 175);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 175);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 175);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c9_rhs175, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs175, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs175), "rhs", "rhs",
                  175);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs175), "lhs", "lhs",
                  175);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 176);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 176);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 176);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c9_rhs176, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs176, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs176), "rhs", "rhs",
                  176);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs176), "lhs", "lhs",
                  176);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 177);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 177);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 177);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c9_rhs177, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs177, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs177), "rhs", "rhs",
                  177);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs177), "lhs", "lhs",
                  177);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 178);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xscal"), "name", "name",
                  178);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 178);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"),
                  "resolved", "resolved", 178);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717472U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c9_rhs178, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs178, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs178), "rhs", "rhs",
                  178);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs178), "lhs", "lhs",
                  178);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 179);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 179);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 179);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 179);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c9_rhs179, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs179, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs179), "rhs", "rhs",
                  179);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs179), "lhs", "lhs",
                  179);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 180);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 180);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 180);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 180);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c9_rhs180, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs180, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs180), "rhs", "rhs",
                  180);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs180), "lhs", "lhs",
                  180);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 181);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 181);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 181);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 181);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c9_rhs181, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs181, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs181), "rhs", "rhs",
                  181);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs181), "lhs", "lhs",
                  181);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 182);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_xscal"), "name",
                  "name", 182);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 182);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "resolved", "resolved", 182);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080384U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c9_rhs182, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs182, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs182), "rhs", "rhs",
                  182);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs182), "lhs", "lhs",
                  182);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 183);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 183);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 183);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 183);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c9_rhs183, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs183, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs183), "rhs", "rhs",
                  183);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs183), "lhs", "lhs",
                  183);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 184);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 184);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 184);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 184);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 184);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 184);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 184);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 184);
  sf_mex_assign(&c9_rhs184, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs184, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs184), "rhs", "rhs",
                  184);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs184), "lhs", "lhs",
                  184);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 185);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 185);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 185);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 185);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 185);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 185);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 185);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 185);
  sf_mex_assign(&c9_rhs185, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs185, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs185), "rhs", "rhs",
                  185);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs185), "lhs", "lhs",
                  185);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 186);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 186);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 186);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 186);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 186);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 186);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 186);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 186);
  sf_mex_assign(&c9_rhs186, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs186, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs186), "rhs", "rhs",
                  186);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs186), "lhs", "lhs",
                  186);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 187);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 187);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 187);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 187);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 187);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 187);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 187);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 187);
  sf_mex_assign(&c9_rhs187, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs187, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs187), "rhs", "rhs",
                  187);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs187), "lhs", "lhs",
                  187);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 188);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mtimes"), "name", "name", 188);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 188);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 188);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 188);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 188);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 188);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 188);
  sf_mex_assign(&c9_rhs188, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs188, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs188), "rhs", "rhs",
                  188);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs188), "lhs", "lhs",
                  188);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 189);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_div"), "name", "name", 189);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 189);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 189);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 189);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 189);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 189);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 189);
  sf_mex_assign(&c9_rhs189, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs189, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs189), "rhs", "rhs",
                  189);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs189), "lhs", "lhs",
                  189);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"),
                  "context", "context", 190);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 190);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 190);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 190);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 190);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 190);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 190);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 190);
  sf_mex_assign(&c9_rhs190, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs190, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs190), "rhs", "rhs",
                  190);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs190), "lhs", "lhs",
                  190);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m!below_threshold"),
                  "context", "context", 191);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("length"), "name", "name", 191);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 191);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 191);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1303149806U), "fileTimeLo",
                  "fileTimeLo", 191);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 191);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 191);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 191);
  sf_mex_assign(&c9_rhs191, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs191, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs191), "rhs", "rhs",
                  191);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs191), "lhs", "lhs",
                  191);
  sf_mex_destroy(&c9_rhs128);
  sf_mex_destroy(&c9_lhs128);
  sf_mex_destroy(&c9_rhs129);
  sf_mex_destroy(&c9_lhs129);
  sf_mex_destroy(&c9_rhs130);
  sf_mex_destroy(&c9_lhs130);
  sf_mex_destroy(&c9_rhs131);
  sf_mex_destroy(&c9_lhs131);
  sf_mex_destroy(&c9_rhs132);
  sf_mex_destroy(&c9_lhs132);
  sf_mex_destroy(&c9_rhs133);
  sf_mex_destroy(&c9_lhs133);
  sf_mex_destroy(&c9_rhs134);
  sf_mex_destroy(&c9_lhs134);
  sf_mex_destroy(&c9_rhs135);
  sf_mex_destroy(&c9_lhs135);
  sf_mex_destroy(&c9_rhs136);
  sf_mex_destroy(&c9_lhs136);
  sf_mex_destroy(&c9_rhs137);
  sf_mex_destroy(&c9_lhs137);
  sf_mex_destroy(&c9_rhs138);
  sf_mex_destroy(&c9_lhs138);
  sf_mex_destroy(&c9_rhs139);
  sf_mex_destroy(&c9_lhs139);
  sf_mex_destroy(&c9_rhs140);
  sf_mex_destroy(&c9_lhs140);
  sf_mex_destroy(&c9_rhs141);
  sf_mex_destroy(&c9_lhs141);
  sf_mex_destroy(&c9_rhs142);
  sf_mex_destroy(&c9_lhs142);
  sf_mex_destroy(&c9_rhs143);
  sf_mex_destroy(&c9_lhs143);
  sf_mex_destroy(&c9_rhs144);
  sf_mex_destroy(&c9_lhs144);
  sf_mex_destroy(&c9_rhs145);
  sf_mex_destroy(&c9_lhs145);
  sf_mex_destroy(&c9_rhs146);
  sf_mex_destroy(&c9_lhs146);
  sf_mex_destroy(&c9_rhs147);
  sf_mex_destroy(&c9_lhs147);
  sf_mex_destroy(&c9_rhs148);
  sf_mex_destroy(&c9_lhs148);
  sf_mex_destroy(&c9_rhs149);
  sf_mex_destroy(&c9_lhs149);
  sf_mex_destroy(&c9_rhs150);
  sf_mex_destroy(&c9_lhs150);
  sf_mex_destroy(&c9_rhs151);
  sf_mex_destroy(&c9_lhs151);
  sf_mex_destroy(&c9_rhs152);
  sf_mex_destroy(&c9_lhs152);
  sf_mex_destroy(&c9_rhs153);
  sf_mex_destroy(&c9_lhs153);
  sf_mex_destroy(&c9_rhs154);
  sf_mex_destroy(&c9_lhs154);
  sf_mex_destroy(&c9_rhs155);
  sf_mex_destroy(&c9_lhs155);
  sf_mex_destroy(&c9_rhs156);
  sf_mex_destroy(&c9_lhs156);
  sf_mex_destroy(&c9_rhs157);
  sf_mex_destroy(&c9_lhs157);
  sf_mex_destroy(&c9_rhs158);
  sf_mex_destroy(&c9_lhs158);
  sf_mex_destroy(&c9_rhs159);
  sf_mex_destroy(&c9_lhs159);
  sf_mex_destroy(&c9_rhs160);
  sf_mex_destroy(&c9_lhs160);
  sf_mex_destroy(&c9_rhs161);
  sf_mex_destroy(&c9_lhs161);
  sf_mex_destroy(&c9_rhs162);
  sf_mex_destroy(&c9_lhs162);
  sf_mex_destroy(&c9_rhs163);
  sf_mex_destroy(&c9_lhs163);
  sf_mex_destroy(&c9_rhs164);
  sf_mex_destroy(&c9_lhs164);
  sf_mex_destroy(&c9_rhs165);
  sf_mex_destroy(&c9_lhs165);
  sf_mex_destroy(&c9_rhs166);
  sf_mex_destroy(&c9_lhs166);
  sf_mex_destroy(&c9_rhs167);
  sf_mex_destroy(&c9_lhs167);
  sf_mex_destroy(&c9_rhs168);
  sf_mex_destroy(&c9_lhs168);
  sf_mex_destroy(&c9_rhs169);
  sf_mex_destroy(&c9_lhs169);
  sf_mex_destroy(&c9_rhs170);
  sf_mex_destroy(&c9_lhs170);
  sf_mex_destroy(&c9_rhs171);
  sf_mex_destroy(&c9_lhs171);
  sf_mex_destroy(&c9_rhs172);
  sf_mex_destroy(&c9_lhs172);
  sf_mex_destroy(&c9_rhs173);
  sf_mex_destroy(&c9_lhs173);
  sf_mex_destroy(&c9_rhs174);
  sf_mex_destroy(&c9_lhs174);
  sf_mex_destroy(&c9_rhs175);
  sf_mex_destroy(&c9_lhs175);
  sf_mex_destroy(&c9_rhs176);
  sf_mex_destroy(&c9_lhs176);
  sf_mex_destroy(&c9_rhs177);
  sf_mex_destroy(&c9_lhs177);
  sf_mex_destroy(&c9_rhs178);
  sf_mex_destroy(&c9_lhs178);
  sf_mex_destroy(&c9_rhs179);
  sf_mex_destroy(&c9_lhs179);
  sf_mex_destroy(&c9_rhs180);
  sf_mex_destroy(&c9_lhs180);
  sf_mex_destroy(&c9_rhs181);
  sf_mex_destroy(&c9_lhs181);
  sf_mex_destroy(&c9_rhs182);
  sf_mex_destroy(&c9_lhs182);
  sf_mex_destroy(&c9_rhs183);
  sf_mex_destroy(&c9_lhs183);
  sf_mex_destroy(&c9_rhs184);
  sf_mex_destroy(&c9_lhs184);
  sf_mex_destroy(&c9_rhs185);
  sf_mex_destroy(&c9_lhs185);
  sf_mex_destroy(&c9_rhs186);
  sf_mex_destroy(&c9_lhs186);
  sf_mex_destroy(&c9_rhs187);
  sf_mex_destroy(&c9_lhs187);
  sf_mex_destroy(&c9_rhs188);
  sf_mex_destroy(&c9_lhs188);
  sf_mex_destroy(&c9_rhs189);
  sf_mex_destroy(&c9_lhs189);
  sf_mex_destroy(&c9_rhs190);
  sf_mex_destroy(&c9_lhs190);
  sf_mex_destroy(&c9_rhs191);
  sf_mex_destroy(&c9_lhs191);
}

static void c9_d_info_helper(const mxArray **c9_info)
{
  const mxArray *c9_rhs192 = NULL;
  const mxArray *c9_lhs192 = NULL;
  const mxArray *c9_rhs193 = NULL;
  const mxArray *c9_lhs193 = NULL;
  const mxArray *c9_rhs194 = NULL;
  const mxArray *c9_lhs194 = NULL;
  const mxArray *c9_rhs195 = NULL;
  const mxArray *c9_lhs195 = NULL;
  const mxArray *c9_rhs196 = NULL;
  const mxArray *c9_lhs196 = NULL;
  const mxArray *c9_rhs197 = NULL;
  const mxArray *c9_lhs197 = NULL;
  const mxArray *c9_rhs198 = NULL;
  const mxArray *c9_lhs198 = NULL;
  const mxArray *c9_rhs199 = NULL;
  const mxArray *c9_lhs199 = NULL;
  const mxArray *c9_rhs200 = NULL;
  const mxArray *c9_lhs200 = NULL;
  const mxArray *c9_rhs201 = NULL;
  const mxArray *c9_lhs201 = NULL;
  const mxArray *c9_rhs202 = NULL;
  const mxArray *c9_lhs202 = NULL;
  const mxArray *c9_rhs203 = NULL;
  const mxArray *c9_lhs203 = NULL;
  const mxArray *c9_rhs204 = NULL;
  const mxArray *c9_lhs204 = NULL;
  const mxArray *c9_rhs205 = NULL;
  const mxArray *c9_lhs205 = NULL;
  const mxArray *c9_rhs206 = NULL;
  const mxArray *c9_lhs206 = NULL;
  const mxArray *c9_rhs207 = NULL;
  const mxArray *c9_lhs207 = NULL;
  const mxArray *c9_rhs208 = NULL;
  const mxArray *c9_lhs208 = NULL;
  const mxArray *c9_rhs209 = NULL;
  const mxArray *c9_lhs209 = NULL;
  const mxArray *c9_rhs210 = NULL;
  const mxArray *c9_lhs210 = NULL;
  const mxArray *c9_rhs211 = NULL;
  const mxArray *c9_lhs211 = NULL;
  const mxArray *c9_rhs212 = NULL;
  const mxArray *c9_lhs212 = NULL;
  const mxArray *c9_rhs213 = NULL;
  const mxArray *c9_lhs213 = NULL;
  const mxArray *c9_rhs214 = NULL;
  const mxArray *c9_lhs214 = NULL;
  const mxArray *c9_rhs215 = NULL;
  const mxArray *c9_lhs215 = NULL;
  const mxArray *c9_rhs216 = NULL;
  const mxArray *c9_lhs216 = NULL;
  const mxArray *c9_rhs217 = NULL;
  const mxArray *c9_lhs217 = NULL;
  const mxArray *c9_rhs218 = NULL;
  const mxArray *c9_lhs218 = NULL;
  const mxArray *c9_rhs219 = NULL;
  const mxArray *c9_lhs219 = NULL;
  const mxArray *c9_rhs220 = NULL;
  const mxArray *c9_lhs220 = NULL;
  const mxArray *c9_rhs221 = NULL;
  const mxArray *c9_lhs221 = NULL;
  const mxArray *c9_rhs222 = NULL;
  const mxArray *c9_lhs222 = NULL;
  const mxArray *c9_rhs223 = NULL;
  const mxArray *c9_lhs223 = NULL;
  const mxArray *c9_rhs224 = NULL;
  const mxArray *c9_lhs224 = NULL;
  const mxArray *c9_rhs225 = NULL;
  const mxArray *c9_lhs225 = NULL;
  const mxArray *c9_rhs226 = NULL;
  const mxArray *c9_lhs226 = NULL;
  const mxArray *c9_rhs227 = NULL;
  const mxArray *c9_lhs227 = NULL;
  const mxArray *c9_rhs228 = NULL;
  const mxArray *c9_lhs228 = NULL;
  const mxArray *c9_rhs229 = NULL;
  const mxArray *c9_lhs229 = NULL;
  const mxArray *c9_rhs230 = NULL;
  const mxArray *c9_lhs230 = NULL;
  const mxArray *c9_rhs231 = NULL;
  const mxArray *c9_lhs231 = NULL;
  const mxArray *c9_rhs232 = NULL;
  const mxArray *c9_lhs232 = NULL;
  const mxArray *c9_rhs233 = NULL;
  const mxArray *c9_lhs233 = NULL;
  const mxArray *c9_rhs234 = NULL;
  const mxArray *c9_lhs234 = NULL;
  const mxArray *c9_rhs235 = NULL;
  const mxArray *c9_lhs235 = NULL;
  const mxArray *c9_rhs236 = NULL;
  const mxArray *c9_lhs236 = NULL;
  const mxArray *c9_rhs237 = NULL;
  const mxArray *c9_lhs237 = NULL;
  const mxArray *c9_rhs238 = NULL;
  const mxArray *c9_lhs238 = NULL;
  const mxArray *c9_rhs239 = NULL;
  const mxArray *c9_lhs239 = NULL;
  const mxArray *c9_rhs240 = NULL;
  const mxArray *c9_lhs240 = NULL;
  const mxArray *c9_rhs241 = NULL;
  const mxArray *c9_lhs241 = NULL;
  const mxArray *c9_rhs242 = NULL;
  const mxArray *c9_lhs242 = NULL;
  const mxArray *c9_rhs243 = NULL;
  const mxArray *c9_lhs243 = NULL;
  const mxArray *c9_rhs244 = NULL;
  const mxArray *c9_lhs244 = NULL;
  const mxArray *c9_rhs245 = NULL;
  const mxArray *c9_lhs245 = NULL;
  const mxArray *c9_rhs246 = NULL;
  const mxArray *c9_lhs246 = NULL;
  const mxArray *c9_rhs247 = NULL;
  const mxArray *c9_lhs247 = NULL;
  const mxArray *c9_rhs248 = NULL;
  const mxArray *c9_lhs248 = NULL;
  const mxArray *c9_rhs249 = NULL;
  const mxArray *c9_lhs249 = NULL;
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 192);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_matlab_zlarf"), "name",
                  "name", 192);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 192);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "resolved", "resolved", 192);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822422U), "fileTimeLo",
                  "fileTimeLo", 192);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 192);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 192);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 192);
  sf_mex_assign(&c9_rhs192, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs192, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs192), "rhs", "rhs",
                  192);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs192), "lhs", "lhs",
                  192);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 193);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 193);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 193);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 193);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 193);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 193);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 193);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 193);
  sf_mex_assign(&c9_rhs193, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs193, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs193), "rhs", "rhs",
                  193);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs193), "lhs", "lhs",
                  193);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 194);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 194);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 194);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 194);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 194);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 194);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 194);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 194);
  sf_mex_assign(&c9_rhs194, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs194, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs194), "rhs", "rhs",
                  194);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs194), "lhs", "lhs",
                  194);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 195);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("isequal"), "name", "name", 195);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 195);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 195);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822358U), "fileTimeLo",
                  "fileTimeLo", 195);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 195);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 195);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 195);
  sf_mex_assign(&c9_rhs195, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs195, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs195), "rhs", "rhs",
                  195);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs195), "lhs", "lhs",
                  195);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 196);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 196);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 196);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 196);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822386U), "fileTimeLo",
                  "fileTimeLo", 196);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 196);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 196);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 196);
  sf_mex_assign(&c9_rhs196, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs196, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs196), "rhs", "rhs",
                  196);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs196), "lhs", "lhs",
                  196);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 197);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 197);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 197);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 197);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 197);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 197);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 197);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 197);
  sf_mex_assign(&c9_rhs197, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs197, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs197), "rhs", "rhs",
                  197);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs197), "lhs", "lhs",
                  197);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 198);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 198);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 198);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 198);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 198);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 198);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 198);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 198);
  sf_mex_assign(&c9_rhs198, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs198, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs198), "rhs", "rhs",
                  198);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs198), "lhs", "lhs",
                  198);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 199);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 199);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 199);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 199);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 199);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 199);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 199);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 199);
  sf_mex_assign(&c9_rhs199, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs199, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs199), "rhs", "rhs",
                  199);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs199), "lhs", "lhs",
                  199);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 200);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 200);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 200);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 200);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 200);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 200);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 200);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 200);
  sf_mex_assign(&c9_rhs200, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs200, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs200), "rhs", "rhs",
                  200);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs200), "lhs", "lhs",
                  200);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m!ilazlc"),
                  "context", "context", 201);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 201);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 201);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 201);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 201);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 201);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 201);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 201);
  sf_mex_assign(&c9_rhs201, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs201, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs201), "rhs", "rhs",
                  201);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs201), "lhs", "lhs",
                  201);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m!ilazlc"),
                  "context", "context", 202);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 202);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 202);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 202);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 202);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 202);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 202);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 202);
  sf_mex_assign(&c9_rhs202, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs202, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs202), "rhs", "rhs",
                  202);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs202), "lhs", "lhs",
                  202);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m!ilazlc"),
                  "context", "context", 203);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 203);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 203);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 203);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 203);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 203);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 203);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 203);
  sf_mex_assign(&c9_rhs203, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs203, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs203), "rhs", "rhs",
                  203);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs203), "lhs", "lhs",
                  203);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m!ilazlc"),
                  "context", "context", 204);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 204);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 204);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 204);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 204);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 204);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 204);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 204);
  sf_mex_assign(&c9_rhs204, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs204, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs204), "rhs", "rhs",
                  204);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs204), "lhs", "lhs",
                  204);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m!ilazlc"),
                  "context", "context", 205);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 205);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 205);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 205);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 205);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 205);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 205);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 205);
  sf_mex_assign(&c9_rhs205, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs205, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs205), "rhs", "rhs",
                  205);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs205), "lhs", "lhs",
                  205);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 206);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xgemv"), "name", "name",
                  206);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 206);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"),
                  "resolved", "resolved", 206);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 206);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 206);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 206);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 206);
  sf_mex_assign(&c9_rhs206, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs206, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs206), "rhs", "rhs",
                  206);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs206), "lhs", "lhs",
                  206);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"), "context",
                  "context", 207);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 207);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 207);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 207);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 207);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 207);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 207);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 207);
  sf_mex_assign(&c9_rhs207, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs207, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs207), "rhs", "rhs",
                  207);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs207), "lhs", "lhs",
                  207);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 208);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 208);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 208);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 208);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 208);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 208);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 208);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 208);
  sf_mex_assign(&c9_rhs208, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs208, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs208), "rhs", "rhs",
                  208);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs208), "lhs", "lhs",
                  208);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 209);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 209);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 209);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 209);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 209);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 209);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 209);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 209);
  sf_mex_assign(&c9_rhs209, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs209, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs209), "rhs", "rhs",
                  209);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs209), "lhs", "lhs",
                  209);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 210);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_xgemv"), "name",
                  "name", 210);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 210);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "resolved", "resolved", 210);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 210);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 210);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 210);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 210);
  sf_mex_assign(&c9_rhs210, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs210, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs210), "rhs", "rhs",
                  210);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs210), "lhs", "lhs",
                  210);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 211);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 211);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 211);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 211);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 211);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 211);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 211);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 211);
  sf_mex_assign(&c9_rhs211, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs211, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs211), "rhs", "rhs",
                  211);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs211), "lhs", "lhs",
                  211);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 212);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 212);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 212);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 212);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 212);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 212);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 212);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 212);
  sf_mex_assign(&c9_rhs212, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs212, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs212), "rhs", "rhs",
                  212);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs212), "lhs", "lhs",
                  212);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 213);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 213);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 213);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 213);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 213);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 213);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 213);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 213);
  sf_mex_assign(&c9_rhs213, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs213, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs213), "rhs", "rhs",
                  213);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs213), "lhs", "lhs",
                  213);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 214);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 214);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 214);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 214);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 214);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 214);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 214);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 214);
  sf_mex_assign(&c9_rhs214, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs214, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs214), "rhs", "rhs",
                  214);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs214), "lhs", "lhs",
                  214);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 215);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 215);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 215);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 215);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 215);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 215);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 215);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 215);
  sf_mex_assign(&c9_rhs215, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs215, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs215), "rhs", "rhs",
                  215);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs215), "lhs", "lhs",
                  215);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 216);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 216);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 216);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 216);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 216);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 216);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 216);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 216);
  sf_mex_assign(&c9_rhs216, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs216, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs216), "rhs", "rhs",
                  216);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs216), "lhs", "lhs",
                  216);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 217);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.conjtimes"),
                  "name", "name", 217);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 217);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/conjtimes.m"),
                  "resolved", "resolved", 217);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1360286186U), "fileTimeLo",
                  "fileTimeLo", 217);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 217);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 217);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 217);
  sf_mex_assign(&c9_rhs217, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs217, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs217), "rhs", "rhs",
                  217);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs217), "lhs", "lhs",
                  217);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"),
                  "context", "context", 218);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xgerc"), "name", "name",
                  218);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 218);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgerc.m"),
                  "resolved", "resolved", 218);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717472U), "fileTimeLo",
                  "fileTimeLo", 218);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 218);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 218);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 218);
  sf_mex_assign(&c9_rhs218, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs218, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs218), "rhs", "rhs",
                  218);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs218), "lhs", "lhs",
                  218);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgerc.m"), "context",
                  "context", 219);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 219);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 219);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 219);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 219);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 219);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 219);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 219);
  sf_mex_assign(&c9_rhs219, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs219, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs219), "rhs", "rhs",
                  219);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs219), "lhs", "lhs",
                  219);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgerc.m"), "context",
                  "context", 220);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xger"), "name", "name",
                  220);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 220);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m"), "resolved",
                  "resolved", 220);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 220);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 220);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 220);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 220);
  sf_mex_assign(&c9_rhs220, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs220, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs220), "rhs", "rhs",
                  220);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs220), "lhs", "lhs",
                  220);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 221);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("abs"), "name", "name", 221);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 221);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 221);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 221);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 221);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 221);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 221);
  sf_mex_assign(&c9_rhs221, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs221, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs221), "rhs", "rhs",
                  221);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs221), "lhs", "lhs",
                  221);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"),
                  "context", "context", 222);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mtimes"), "name", "name", 222);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 222);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 222);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 222);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 222);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 222);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 222);
  sf_mex_assign(&c9_rhs222, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs222, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs222), "rhs", "rhs",
                  222);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs222), "lhs", "lhs",
                  222);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 223);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("max"), "name", "name", 223);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 223);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "resolved",
                  "resolved", 223);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1311258916U), "fileTimeLo",
                  "fileTimeLo", 223);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 223);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 223);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 223);
  sf_mex_assign(&c9_rhs223, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs223, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs223), "rhs", "rhs",
                  223);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs223), "lhs", "lhs",
                  223);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "context",
                  "context", 224);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 224);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 224);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 224);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 224);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 224);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 224);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 224);
  sf_mex_assign(&c9_rhs224, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs224, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs224), "rhs", "rhs",
                  224);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs224), "lhs", "lhs",
                  224);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 225);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xcabs1"), "name", "name",
                  225);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 225);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "resolved", "resolved", 225);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822306U), "fileTimeLo",
                  "fileTimeLo", 225);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 225);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 225);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 225);
  sf_mex_assign(&c9_rhs225, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs225, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs225), "rhs", "rhs",
                  225);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs225), "lhs", "lhs",
                  225);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 226);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("mtimes"), "name", "name", 226);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 226);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 226);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 226);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 226);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 226);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 226);
  sf_mex_assign(&c9_rhs226, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs226, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs226), "rhs", "rhs",
                  226);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs226), "lhs", "lhs",
                  226);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 227);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eps"), "name", "name", 227);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 227);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 227);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 227);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 227);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 227);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 227);
  sf_mex_assign(&c9_rhs227, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs227, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs227), "rhs", "rhs",
                  227);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs227), "lhs", "lhs",
                  227);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 228);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_flt2str"), "name", "name",
                  228);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 228);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "resolved",
                  "resolved", 228);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 228);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 228);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 228);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 228);
  sf_mex_assign(&c9_rhs228, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs228, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs228), "rhs", "rhs",
                  228);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs228), "lhs", "lhs",
                  228);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "context",
                  "context", 229);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "name", "name", 229);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 229);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m"), "resolved",
                  "resolved", 229);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1319733568U), "fileTimeLo",
                  "fileTimeLo", 229);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 229);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 229);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 229);
  sf_mex_assign(&c9_rhs229, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs229, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs229), "rhs", "rhs",
                  229);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs229), "lhs", "lhs",
                  229);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 230);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_warning"), "name", "name",
                  230);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 230);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 230);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822402U), "fileTimeLo",
                  "fileTimeLo", 230);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 230);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 230);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 230);
  sf_mex_assign(&c9_rhs230, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs230, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs230), "rhs", "rhs",
                  230);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs230), "lhs", "lhs",
                  230);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 231);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 231);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 231);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 231);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 231);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 231);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 231);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 231);
  sf_mex_assign(&c9_rhs231, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs231, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs231), "rhs", "rhs",
                  231);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs231), "lhs", "lhs",
                  231);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 232);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.conjtimes"),
                  "name", "name", 232);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 232);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/conjtimes.m"),
                  "resolved", "resolved", 232);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1360286186U), "fileTimeLo",
                  "fileTimeLo", 232);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 232);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 232);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 232);
  sf_mex_assign(&c9_rhs232, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs232, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs232), "rhs", "rhs",
                  232);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs232), "lhs", "lhs",
                  232);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_qrsolve.m"), "context",
                  "context", 233);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_div"), "name", "name", 233);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 233);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 233);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 233);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 233);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 233);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 233);
  sf_mex_assign(&c9_rhs233, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs233, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs233), "rhs", "rhs",
                  233);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs233), "lhs", "lhs",
                  233);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 234);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 234);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 234);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 234);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 234);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 234);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 234);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 234);
  sf_mex_assign(&c9_rhs234, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs234, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs234), "rhs", "rhs",
                  234);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs234), "lhs", "lhs",
                  234);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 235);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 235);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 235);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 235);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 235);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 235);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 235);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 235);
  sf_mex_assign(&c9_rhs235, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs235, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs235), "rhs", "rhs",
                  235);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs235), "lhs", "lhs",
                  235);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 236);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  236);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 236);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 236);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 236);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 236);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 236);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 236);
  sf_mex_assign(&c9_rhs236, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs236, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs236), "rhs", "rhs",
                  236);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs236), "lhs", "lhs",
                  236);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 237);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 237);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 237);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 237);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 237);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 237);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 237);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 237);
  sf_mex_assign(&c9_rhs237, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs237, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs237), "rhs", "rhs",
                  237);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs237), "lhs", "lhs",
                  237);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 238);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 238);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 238);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 238);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 238);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 238);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 238);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 238);
  sf_mex_assign(&c9_rhs238, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs238, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs238), "rhs", "rhs",
                  238);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs238), "lhs", "lhs",
                  238);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 239);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 239);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 239);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 239);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 239);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 239);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 239);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 239);
  sf_mex_assign(&c9_rhs239, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs239, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs239), "rhs", "rhs",
                  239);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs239), "lhs", "lhs",
                  239);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 240);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 240);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 240);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 240);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 240);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 240);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 240);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 240);
  sf_mex_assign(&c9_rhs240, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs240, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs240), "rhs", "rhs",
                  240);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs240), "lhs", "lhs",
                  240);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 241);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 241);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 241);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 241);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 241);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 241);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 241);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 241);
  sf_mex_assign(&c9_rhs241, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs241, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs241), "rhs", "rhs",
                  241);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs241), "lhs", "lhs",
                  241);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 242);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 242);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 242);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 242);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 242);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 242);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 242);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 242);
  sf_mex_assign(&c9_rhs242, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs242, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs242), "rhs", "rhs",
                  242);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs242), "lhs", "lhs",
                  242);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 243);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 243);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 243);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 243);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 243);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 243);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 243);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 243);
  sf_mex_assign(&c9_rhs243, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs243, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs243), "rhs", "rhs",
                  243);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs243), "lhs", "lhs",
                  243);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 244);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 244);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 244);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 244);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 244);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 244);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 244);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 244);
  sf_mex_assign(&c9_rhs244, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs244, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs244), "rhs", "rhs",
                  244);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs244), "lhs", "lhs",
                  244);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 245);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 245);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 245);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 245);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 245);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 245);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 245);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 245);
  sf_mex_assign(&c9_rhs245, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs245, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs245), "rhs", "rhs",
                  245);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs245), "lhs", "lhs",
                  245);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 246);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 246);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 246);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 246);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 246);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 246);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 246);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 246);
  sf_mex_assign(&c9_rhs246, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs246, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs246), "rhs", "rhs",
                  246);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs246), "lhs", "lhs",
                  246);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 247);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 247);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 247);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 247);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 247);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 247);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 247);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 247);
  sf_mex_assign(&c9_rhs247, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs247, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs247), "rhs", "rhs",
                  247);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs247), "lhs", "lhs",
                  247);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(""), "context", "context", 248);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("tanh"), "name", "name", 248);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 248);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tanh.m"), "resolved",
                  "resolved", 248);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1343833988U), "fileTimeLo",
                  "fileTimeLo", 248);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 248);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 248);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 248);
  sf_mex_assign(&c9_rhs248, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs248, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs248), "rhs", "rhs",
                  248);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs248), "lhs", "lhs",
                  248);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tanh.m"), "context",
                  "context", 249);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("eml_scalar_tanh"), "name",
                  "name", 249);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 249);
  sf_mex_addfield(*c9_info, c9_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_tanh.m"),
                  "resolved", "resolved", 249);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(1286822340U), "fileTimeLo",
                  "fileTimeLo", 249);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 249);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 249);
  sf_mex_addfield(*c9_info, c9_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 249);
  sf_mex_assign(&c9_rhs249, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c9_lhs249, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_rhs249), "rhs", "rhs",
                  249);
  sf_mex_addfield(*c9_info, sf_mex_duplicatearraysafe(&c9_lhs249), "lhs", "lhs",
                  249);
  sf_mex_destroy(&c9_rhs192);
  sf_mex_destroy(&c9_lhs192);
  sf_mex_destroy(&c9_rhs193);
  sf_mex_destroy(&c9_lhs193);
  sf_mex_destroy(&c9_rhs194);
  sf_mex_destroy(&c9_lhs194);
  sf_mex_destroy(&c9_rhs195);
  sf_mex_destroy(&c9_lhs195);
  sf_mex_destroy(&c9_rhs196);
  sf_mex_destroy(&c9_lhs196);
  sf_mex_destroy(&c9_rhs197);
  sf_mex_destroy(&c9_lhs197);
  sf_mex_destroy(&c9_rhs198);
  sf_mex_destroy(&c9_lhs198);
  sf_mex_destroy(&c9_rhs199);
  sf_mex_destroy(&c9_lhs199);
  sf_mex_destroy(&c9_rhs200);
  sf_mex_destroy(&c9_lhs200);
  sf_mex_destroy(&c9_rhs201);
  sf_mex_destroy(&c9_lhs201);
  sf_mex_destroy(&c9_rhs202);
  sf_mex_destroy(&c9_lhs202);
  sf_mex_destroy(&c9_rhs203);
  sf_mex_destroy(&c9_lhs203);
  sf_mex_destroy(&c9_rhs204);
  sf_mex_destroy(&c9_lhs204);
  sf_mex_destroy(&c9_rhs205);
  sf_mex_destroy(&c9_lhs205);
  sf_mex_destroy(&c9_rhs206);
  sf_mex_destroy(&c9_lhs206);
  sf_mex_destroy(&c9_rhs207);
  sf_mex_destroy(&c9_lhs207);
  sf_mex_destroy(&c9_rhs208);
  sf_mex_destroy(&c9_lhs208);
  sf_mex_destroy(&c9_rhs209);
  sf_mex_destroy(&c9_lhs209);
  sf_mex_destroy(&c9_rhs210);
  sf_mex_destroy(&c9_lhs210);
  sf_mex_destroy(&c9_rhs211);
  sf_mex_destroy(&c9_lhs211);
  sf_mex_destroy(&c9_rhs212);
  sf_mex_destroy(&c9_lhs212);
  sf_mex_destroy(&c9_rhs213);
  sf_mex_destroy(&c9_lhs213);
  sf_mex_destroy(&c9_rhs214);
  sf_mex_destroy(&c9_lhs214);
  sf_mex_destroy(&c9_rhs215);
  sf_mex_destroy(&c9_lhs215);
  sf_mex_destroy(&c9_rhs216);
  sf_mex_destroy(&c9_lhs216);
  sf_mex_destroy(&c9_rhs217);
  sf_mex_destroy(&c9_lhs217);
  sf_mex_destroy(&c9_rhs218);
  sf_mex_destroy(&c9_lhs218);
  sf_mex_destroy(&c9_rhs219);
  sf_mex_destroy(&c9_lhs219);
  sf_mex_destroy(&c9_rhs220);
  sf_mex_destroy(&c9_lhs220);
  sf_mex_destroy(&c9_rhs221);
  sf_mex_destroy(&c9_lhs221);
  sf_mex_destroy(&c9_rhs222);
  sf_mex_destroy(&c9_lhs222);
  sf_mex_destroy(&c9_rhs223);
  sf_mex_destroy(&c9_lhs223);
  sf_mex_destroy(&c9_rhs224);
  sf_mex_destroy(&c9_lhs224);
  sf_mex_destroy(&c9_rhs225);
  sf_mex_destroy(&c9_lhs225);
  sf_mex_destroy(&c9_rhs226);
  sf_mex_destroy(&c9_lhs226);
  sf_mex_destroy(&c9_rhs227);
  sf_mex_destroy(&c9_lhs227);
  sf_mex_destroy(&c9_rhs228);
  sf_mex_destroy(&c9_lhs228);
  sf_mex_destroy(&c9_rhs229);
  sf_mex_destroy(&c9_lhs229);
  sf_mex_destroy(&c9_rhs230);
  sf_mex_destroy(&c9_lhs230);
  sf_mex_destroy(&c9_rhs231);
  sf_mex_destroy(&c9_lhs231);
  sf_mex_destroy(&c9_rhs232);
  sf_mex_destroy(&c9_lhs232);
  sf_mex_destroy(&c9_rhs233);
  sf_mex_destroy(&c9_lhs233);
  sf_mex_destroy(&c9_rhs234);
  sf_mex_destroy(&c9_lhs234);
  sf_mex_destroy(&c9_rhs235);
  sf_mex_destroy(&c9_lhs235);
  sf_mex_destroy(&c9_rhs236);
  sf_mex_destroy(&c9_lhs236);
  sf_mex_destroy(&c9_rhs237);
  sf_mex_destroy(&c9_lhs237);
  sf_mex_destroy(&c9_rhs238);
  sf_mex_destroy(&c9_lhs238);
  sf_mex_destroy(&c9_rhs239);
  sf_mex_destroy(&c9_lhs239);
  sf_mex_destroy(&c9_rhs240);
  sf_mex_destroy(&c9_lhs240);
  sf_mex_destroy(&c9_rhs241);
  sf_mex_destroy(&c9_lhs241);
  sf_mex_destroy(&c9_rhs242);
  sf_mex_destroy(&c9_lhs242);
  sf_mex_destroy(&c9_rhs243);
  sf_mex_destroy(&c9_lhs243);
  sf_mex_destroy(&c9_rhs244);
  sf_mex_destroy(&c9_lhs244);
  sf_mex_destroy(&c9_rhs245);
  sf_mex_destroy(&c9_lhs245);
  sf_mex_destroy(&c9_rhs246);
  sf_mex_destroy(&c9_lhs246);
  sf_mex_destroy(&c9_rhs247);
  sf_mex_destroy(&c9_lhs247);
  sf_mex_destroy(&c9_rhs248);
  sf_mex_destroy(&c9_lhs248);
  sf_mex_destroy(&c9_rhs249);
  sf_mex_destroy(&c9_lhs249);
}

static void c9_diag(SFc9_simulationInstanceStruct *chartInstance, real_T
                    c9_v_data[20], int32_T c9_v_sizes[2], real_T c9_d_data[400],
                    int32_T c9_d_sizes[2])
{
  int32_T c9_nv;
  int32_T c9_m;
  int32_T c9_iv4[2];
  int32_T c9_iv5[2];
  int32_T c9_iv6[2];
  int32_T c9_d;
  int32_T c9_b_d;
  int32_T c9_loop_ub;
  int32_T c9_i42;
  int32_T c9_b_nv;
  int32_T c9_b;
  int32_T c9_b_b;
  boolean_T c9_overflow;
  int32_T c9_j;
  int32_T c9_b_j;
  c9_nv = (int32_T)(real_T)c9_v_sizes[1];
  c9_m = c9_nv;
  c9_iv4[0] = c9_m;
  c9_iv4[1] = c9_m;
  c9_iv5[0] = c9_iv4[0];
  c9_iv5[1] = c9_iv4[1];
  c9_d_sizes[0] = c9_iv5[0];
  c9_iv6[0] = c9_iv4[0];
  c9_iv6[1] = c9_iv4[1];
  c9_d_sizes[1] = c9_iv6[1];
  c9_d = c9_d_sizes[0];
  c9_b_d = c9_d_sizes[1];
  c9_loop_ub = c9_iv4[0] * c9_iv4[1] - 1;
  for (c9_i42 = 0; c9_i42 <= c9_loop_ub; c9_i42++) {
    c9_d_data[c9_i42] = 0.0;
  }

  c9_b_nv = c9_nv;
  c9_b = c9_b_nv;
  c9_b_b = c9_b;
  if (1 > c9_b_b) {
    c9_overflow = FALSE;
  } else {
    c9_overflow = (c9_b_b > 2147483646);
  }

  if (c9_overflow) {
    c9_check_forloop_overflow_error(chartInstance, TRUE);
  }

  for (c9_j = 1; c9_j <= c9_b_nv; c9_j++) {
    c9_b_j = c9_j;
    c9_d_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1, c9_d_sizes[0], 1, 0) +
               c9_d_sizes[0] * (_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
      c9_d_sizes[1], 2, 0) - 1)) - 1] = c9_v_data[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      c9_b_j, 1, c9_v_sizes[1], 1, 0) - 1];
  }
}

static void c9_check_forloop_overflow_error(SFc9_simulationInstanceStruct
  *chartInstance, boolean_T c9_overflow)
{
  int32_T c9_i43;
  static char_T c9_cv2[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c9_u[34];
  const mxArray *c9_y = NULL;
  int32_T c9_i44;
  static char_T c9_cv3[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c9_b_u[23];
  const mxArray *c9_b_y = NULL;
  for (c9_i43 = 0; c9_i43 < 34; c9_i43++) {
    c9_u[c9_i43] = c9_cv2[c9_i43];
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 10, 0U, 1U, 0U, 2, 1, 34), FALSE);
  for (c9_i44 = 0; c9_i44 < 23; c9_i44++) {
    c9_b_u[c9_i44] = c9_cv3[c9_i44];
  }

  c9_b_y = NULL;
  sf_mex_assign(&c9_b_y, sf_mex_create("y", c9_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U, 14,
    c9_y, 14, c9_b_y));
}

static void c9_eye(SFc9_simulationInstanceStruct *chartInstance, real_T
                   c9_varargin_1, real_T c9_I_data[400], int32_T c9_I_sizes[2])
{
  real_T c9_arg1;
  real_T c9_m;
  real_T c9_b_n;
  int32_T c9_d;
  int32_T c9_b_m[2];
  int32_T c9_c_m[2];
  int32_T c9_I;
  int32_T c9_b_I;
  int32_T c9_loop_ub;
  int32_T c9_i45;
  int32_T c9_b_d;
  int32_T c9_b;
  int32_T c9_b_b;
  boolean_T c9_overflow;
  int32_T c9_k;
  int32_T c9_b_k;
  c9_arg1 = c9_varargin_1;
  c9_m = c9_arg1;
  c9_b_n = c9_arg1;
  c9_eml_assert_valid_size_arg(chartInstance, c9_m);
  c9_d = (int32_T)c9_m;
  c9_b_m[0] = (int32_T)c9_m;
  c9_b_m[1] = (int32_T)c9_b_n;
  c9_I_sizes[0] = c9_b_m[0];
  c9_c_m[0] = (int32_T)c9_m;
  c9_c_m[1] = (int32_T)c9_b_n;
  c9_I_sizes[1] = c9_c_m[1];
  c9_I = c9_I_sizes[0];
  c9_b_I = c9_I_sizes[1];
  c9_loop_ub = (int32_T)c9_m * (int32_T)c9_b_n - 1;
  for (c9_i45 = 0; c9_i45 <= c9_loop_ub; c9_i45++) {
    c9_I_data[c9_i45] = 0.0;
  }

  if (c9_d > 0) {
    c9_b_d = c9_d;
    c9_b = c9_b_d;
    c9_b_b = c9_b;
    if (1 > c9_b_b) {
      c9_overflow = FALSE;
    } else {
      c9_overflow = (c9_b_b > 2147483646);
    }

    if (c9_overflow) {
      c9_check_forloop_overflow_error(chartInstance, TRUE);
    }

    for (c9_k = 1; c9_k <= c9_b_d; c9_k++) {
      c9_b_k = c9_k;
      c9_I_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_k, 1, c9_I_sizes[0], 1, 0)
                 + c9_I_sizes[0] * (_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_k, 1,
        c9_I_sizes[1], 2, 0) - 1)) - 1] = 1.0;
    }
  }
}

static void c9_eml_assert_valid_size_arg(SFc9_simulationInstanceStruct
  *chartInstance, real_T c9_varargin_1)
{
  real_T c9_arg;
  boolean_T c9_p;
  boolean_T c9_b2;
  int32_T c9_i46;
  static char_T c9_cv4[28] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'N', 'o', 'n', 'I', 'n', 't', 'e', 'g', 'e', 'r', 'I', 'n',
    'p', 'u', 't' };

  char_T c9_u[28];
  const mxArray *c9_y = NULL;
  int32_T c9_b_u;
  const mxArray *c9_b_y = NULL;
  int32_T c9_c_u;
  const mxArray *c9_c_y = NULL;
  c9_arg = c9_varargin_1;
  if (c9_arg != c9_arg) {
    c9_p = FALSE;
  } else {
    c9_p = TRUE;
  }

  c9_b2 = c9_p;
  if (c9_b2) {
  } else {
    for (c9_i46 = 0; c9_i46 < 28; c9_i46++) {
      c9_u[c9_i46] = c9_cv4[c9_i46];
    }

    c9_y = NULL;
    sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 10, 0U, 1U, 0U, 2, 1, 28),
                  FALSE);
    c9_b_u = MIN_int32_T;
    c9_b_y = NULL;
    sf_mex_assign(&c9_b_y, sf_mex_create("y", &c9_b_u, 6, 0U, 0U, 0U, 0), FALSE);
    c9_c_u = MAX_int32_T;
    c9_c_y = NULL;
    sf_mex_assign(&c9_c_y, sf_mex_create("y", &c9_c_u, 6, 0U, 0U, 0U, 0), FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 3U,
      14, c9_y, 14, c9_b_y, 14, c9_c_y));
  }
}

static void c9_mldivide(SFc9_simulationInstanceStruct *chartInstance, real_T
  c9_A_data[400], int32_T c9_A_sizes[2], real_T c9_B_data[400], int32_T
  c9_B_sizes[2], real_T c9_Y_data[400], int32_T c9_Y_sizes[2])
{
  int32_T c9_i47;
  static char_T c9_cv5[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'd', 'i', 'm', 'a', 'g', 'r', 'e', 'e' };

  char_T c9_u[21];
  const mxArray *c9_y = NULL;
  boolean_T c9_b3;
  boolean_T c9_b4;
  boolean_T c9_b5;
  boolean_T c9_b6;
  real_T c9_dv1[2];
  int32_T c9_iv7[2];
  int32_T c9_iv8[2];
  int32_T c9_Y;
  int32_T c9_b_Y;
  int32_T c9_loop_ub;
  int32_T c9_i48;
  int32_T c9_b_A_sizes[2];
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_b_loop_ub;
  int32_T c9_i49;
  real_T c9_b_A_data[400];
  int32_T c9_b_B_sizes[2];
  int32_T c9_B;
  int32_T c9_b_B;
  int32_T c9_c_loop_ub;
  int32_T c9_i50;
  real_T c9_b_B_data[400];
  int32_T c9_c_A_sizes[2];
  int32_T c9_c_A;
  int32_T c9_d_A;
  int32_T c9_d_loop_ub;
  int32_T c9_i51;
  real_T c9_c_A_data[400];
  int32_T c9_c_B_sizes[2];
  int32_T c9_c_B;
  int32_T c9_d_B;
  int32_T c9_e_loop_ub;
  int32_T c9_i52;
  real_T c9_c_B_data[400];
  boolean_T guard1 = FALSE;
  if ((real_T)c9_B_sizes[0] == (real_T)c9_A_sizes[0]) {
  } else {
    for (c9_i47 = 0; c9_i47 < 21; c9_i47++) {
      c9_u[c9_i47] = c9_cv5[c9_i47];
    }

    c9_y = NULL;
    sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 10, 0U, 1U, 0U, 2, 1, 21),
                  FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
      14, c9_y));
  }

  c9_b3 = (c9_A_sizes[0] == 0);
  c9_b4 = (c9_A_sizes[1] == 0);
  guard1 = FALSE;
  if (c9_b3 || c9_b4) {
    guard1 = TRUE;
  } else {
    c9_b5 = (c9_B_sizes[0] == 0);
    c9_b6 = (c9_B_sizes[1] == 0);
    if (c9_b5 || c9_b6) {
      guard1 = TRUE;
    } else if ((real_T)c9_A_sizes[0] == (real_T)c9_A_sizes[1]) {
      c9_b_A_sizes[0] = c9_A_sizes[0];
      c9_b_A_sizes[1] = c9_A_sizes[1];
      c9_A = c9_b_A_sizes[0];
      c9_b_A = c9_b_A_sizes[1];
      c9_b_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
      for (c9_i49 = 0; c9_i49 <= c9_b_loop_ub; c9_i49++) {
        c9_b_A_data[c9_i49] = c9_A_data[c9_i49];
      }

      c9_b_B_sizes[0] = c9_B_sizes[0];
      c9_b_B_sizes[1] = c9_B_sizes[1];
      c9_B = c9_b_B_sizes[0];
      c9_b_B = c9_b_B_sizes[1];
      c9_c_loop_ub = c9_B_sizes[0] * c9_B_sizes[1] - 1;
      for (c9_i50 = 0; c9_i50 <= c9_c_loop_ub; c9_i50++) {
        c9_b_B_data[c9_i50] = c9_B_data[c9_i50];
      }

      c9_eml_lusolve(chartInstance, c9_b_A_data, c9_b_A_sizes, c9_b_B_data,
                     c9_b_B_sizes, c9_Y_data, c9_Y_sizes);
    } else {
      c9_c_A_sizes[0] = c9_A_sizes[0];
      c9_c_A_sizes[1] = c9_A_sizes[1];
      c9_c_A = c9_c_A_sizes[0];
      c9_d_A = c9_c_A_sizes[1];
      c9_d_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
      for (c9_i51 = 0; c9_i51 <= c9_d_loop_ub; c9_i51++) {
        c9_c_A_data[c9_i51] = c9_A_data[c9_i51];
      }

      c9_c_B_sizes[0] = c9_B_sizes[0];
      c9_c_B_sizes[1] = c9_B_sizes[1];
      c9_c_B = c9_c_B_sizes[0];
      c9_d_B = c9_c_B_sizes[1];
      c9_e_loop_ub = c9_B_sizes[0] * c9_B_sizes[1] - 1;
      for (c9_i52 = 0; c9_i52 <= c9_e_loop_ub; c9_i52++) {
        c9_c_B_data[c9_i52] = c9_B_data[c9_i52];
      }

      c9_eml_qrsolve(chartInstance, c9_c_A_data, c9_c_A_sizes, c9_c_B_data,
                     c9_c_B_sizes, c9_Y_data, c9_Y_sizes);
    }
  }

  if (guard1 == TRUE) {
    c9_eml_scalar_eg(chartInstance);
    c9_dv1[0] = (real_T)c9_A_sizes[1];
    c9_dv1[1] = (real_T)c9_B_sizes[1];
    c9_iv7[0] = (int32_T)c9_dv1[0];
    c9_iv7[1] = (int32_T)c9_dv1[1];
    c9_Y_sizes[0] = c9_iv7[0];
    c9_iv8[0] = (int32_T)c9_dv1[0];
    c9_iv8[1] = (int32_T)c9_dv1[1];
    c9_Y_sizes[1] = c9_iv8[1];
    c9_Y = c9_Y_sizes[0];
    c9_b_Y = c9_Y_sizes[1];
    c9_loop_ub = (int32_T)c9_dv1[0] * (int32_T)c9_dv1[1] - 1;
    for (c9_i48 = 0; c9_i48 <= c9_loop_ub; c9_i48++) {
      c9_Y_data[c9_i48] = 0.0;
    }
  }
}

static void c9_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_eml_lusolve(SFc9_simulationInstanceStruct *chartInstance, real_T
  c9_A_data[400], int32_T c9_A_sizes[2], real_T c9_B_data[400], int32_T
  c9_B_sizes[2], real_T c9_X_data[400], int32_T c9_X_sizes[2])
{
  int32_T c9_b_A_sizes[2];
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_loop_ub;
  int32_T c9_i53;
  real_T c9_b_A_data[400];
  int32_T c9_b_n;
  int32_T c9_m;
  int32_T c9_c_n;
  int32_T c9_lda;
  int32_T c9_b_m;
  int32_T c9_d_n;
  int32_T c9_b_lda;
  int32_T c9_c_m;
  int32_T c9_e_n;
  int32_T c9_c_lda;
  int32_T c9_varargin_1;
  int32_T c9_varargin_2;
  int32_T c9_b_varargin_2;
  int32_T c9_varargin_3;
  int32_T c9_x;
  int32_T c9_y;
  int32_T c9_b_x;
  int32_T c9_b_y;
  int32_T c9_xk;
  int32_T c9_yk;
  int32_T c9_c_x;
  int32_T c9_c_y;
  int32_T c9_minval;
  int32_T c9_d;
  int32_T c9_b;
  int32_T c9_b_b;
  int32_T c9_a;
  int32_T c9_bi;
  int32_T c9_ipiv_sizes[2];
  int32_T c9_ipiv_data[20];
  int32_T c9_info;
  int32_T c9_b_a;
  int32_T c9_ldap1;
  int32_T c9_b_varargin_1;
  int32_T c9_c_varargin_2;
  int32_T c9_d_varargin_2;
  int32_T c9_b_varargin_3;
  int32_T c9_d_x;
  int32_T c9_d_y;
  int32_T c9_e_x;
  int32_T c9_e_y;
  int32_T c9_b_xk;
  int32_T c9_b_yk;
  int32_T c9_f_x;
  int32_T c9_f_y;
  int32_T c9_i54;
  int32_T c9_c_b;
  int32_T c9_d_b;
  boolean_T c9_overflow;
  int32_T c9_j;
  int32_T c9_b_j;
  int32_T c9_c_a;
  int32_T c9_jm1;
  int32_T c9_d_a;
  int32_T c9_e_b;
  int32_T c9_mmj;
  int32_T c9_e_a;
  int32_T c9_f_b;
  int32_T c9_c;
  int32_T c9_g_b;
  int32_T c9_jj;
  int32_T c9_f_a;
  int32_T c9_jp1j;
  int32_T c9_g_a;
  int32_T c9_b_c;
  int32_T c9_f_n;
  int32_T c9_ix0;
  int32_T c9_g_n;
  int32_T c9_b_ix0;
  int32_T c9_idxmax;
  int32_T c9_h_n;
  int32_T c9_c_ix0;
  int32_T c9_var;
  ptrdiff_t c9_n_t;
  ptrdiff_t c9_incx_t;
  double * c9_xix0_t;
  ptrdiff_t c9_idxmax_t;
  int32_T c9_h_a;
  int32_T c9_jpiv_offset;
  int32_T c9_i_a;
  int32_T c9_h_b;
  int32_T c9_jpiv;
  int32_T c9_j_a;
  int32_T c9_i_b;
  int32_T c9_c_c;
  int32_T c9_j_b;
  int32_T c9_jrow;
  int32_T c9_k_a;
  int32_T c9_k_b;
  int32_T c9_jprow;
  int32_T c9_b_jp1j;
  int32_T c9_l_a;
  int32_T c9_d_c;
  int32_T c9_m_a;
  int32_T c9_l_b;
  int32_T c9_i55;
  int32_T c9_n_a;
  int32_T c9_m_b;
  int32_T c9_o_a;
  int32_T c9_n_b;
  boolean_T c9_b_overflow;
  int32_T c9_b_i;
  int32_T c9_c_i;
  real_T c9_g_x;
  real_T c9_g_y;
  real_T c9_z;
  int32_T c9_p_a;
  int32_T c9_o_b;
  int32_T c9_e_c;
  int32_T c9_q_a;
  int32_T c9_p_b;
  int32_T c9_f_c;
  int32_T c9_r_a;
  int32_T c9_q_b;
  int32_T c9_g_c;
  int32_T c9_d_m;
  int32_T c9_i_n;
  int32_T c9_d_ix0;
  int32_T c9_iy0;
  int32_T c9_incy;
  int32_T c9_ia0;
  int32_T c9_d_lda;
  real_T c9_d4;
  int32_T c9_b_info;
  int32_T c9_c_info;
  int32_T c9_d_info;
  int32_T c9_X;
  int32_T c9_b_X;
  int32_T c9_b_loop_ub;
  int32_T c9_i56;
  int32_T c9_nb;
  int32_T c9_j_n;
  int32_T c9_r_b;
  int32_T c9_s_b;
  boolean_T c9_c_overflow;
  int32_T c9_d_i;
  int32_T c9_e_i;
  int32_T c9_ip;
  int32_T c9_b_nb;
  int32_T c9_t_b;
  int32_T c9_u_b;
  boolean_T c9_d_overflow;
  int32_T c9_c_j;
  int32_T c9_d_j;
  real_T c9_temp;
  int32_T c9_c_A_sizes[2];
  int32_T c9_c_A;
  int32_T c9_d_A;
  int32_T c9_c_loop_ub;
  int32_T c9_i57;
  real_T c9_c_A_data[400];
  int32_T c9_d_A_sizes[2];
  int32_T c9_e_A;
  int32_T c9_f_A;
  int32_T c9_d_loop_ub;
  int32_T c9_i58;
  real_T c9_d_A_data[400];
  c9_b_A_sizes[0] = c9_A_sizes[0];
  c9_b_A_sizes[1] = c9_A_sizes[1];
  c9_A = c9_b_A_sizes[0];
  c9_b_A = c9_b_A_sizes[1];
  c9_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
  for (c9_i53 = 0; c9_i53 <= c9_loop_ub; c9_i53++) {
    c9_b_A_data[c9_i53] = c9_A_data[c9_i53];
  }

  c9_b_n = (int32_T)(real_T)c9_b_A_sizes[1];
  c9_m = c9_b_n;
  c9_c_n = c9_b_n;
  c9_lda = c9_b_n;
  c9_b_m = c9_m;
  c9_d_n = c9_c_n;
  c9_b_lda = c9_lda;
  c9_c_m = c9_b_m;
  c9_e_n = c9_d_n;
  c9_c_lda = c9_b_lda;
  c9_realmin(chartInstance);
  c9_eps(chartInstance);
  c9_varargin_1 = c9_c_m;
  c9_varargin_2 = c9_e_n;
  c9_b_varargin_2 = c9_varargin_1;
  c9_varargin_3 = c9_varargin_2;
  c9_x = c9_b_varargin_2;
  c9_y = c9_varargin_3;
  c9_b_x = c9_x;
  c9_b_y = c9_y;
  c9_b_eml_scalar_eg(chartInstance);
  c9_xk = c9_b_x;
  c9_yk = c9_b_y;
  c9_c_x = c9_xk;
  c9_c_y = c9_yk;
  c9_b_eml_scalar_eg(chartInstance);
  c9_minval = muIntScalarMin_sint32(c9_c_x, c9_c_y);
  c9_d = c9_minval;
  c9_b = c9_d;
  c9_b_b = c9_b;
  c9_a = c9_b_b;
  c9_bi = c9_a;
  c9_eml_signed_integer_colon(chartInstance, c9_bi, c9_ipiv_data, c9_ipiv_sizes);
  c9_info = 0;
  if (c9_c_m < 1) {
  } else if (c9_e_n < 1) {
  } else {
    c9_b_a = c9_c_lda + 1;
    c9_ldap1 = c9_b_a;
    c9_b_varargin_1 = c9_c_m;
    c9_c_varargin_2 = c9_e_n;
    c9_d_varargin_2 = c9_b_varargin_1 - 1;
    c9_b_varargin_3 = c9_c_varargin_2;
    c9_d_x = c9_d_varargin_2;
    c9_d_y = c9_b_varargin_3;
    c9_e_x = c9_d_x;
    c9_e_y = c9_d_y;
    c9_b_eml_scalar_eg(chartInstance);
    c9_b_xk = c9_e_x;
    c9_b_yk = c9_e_y;
    c9_f_x = c9_b_xk;
    c9_f_y = c9_b_yk;
    c9_b_eml_scalar_eg(chartInstance);
    c9_i54 = muIntScalarMin_sint32(c9_f_x, c9_f_y);
    c9_c_b = c9_i54;
    c9_d_b = c9_c_b;
    if (1 > c9_d_b) {
      c9_overflow = FALSE;
    } else {
      c9_overflow = (c9_d_b > 2147483646);
    }

    if (c9_overflow) {
      c9_check_forloop_overflow_error(chartInstance, TRUE);
    }

    for (c9_j = 1; c9_j <= c9_i54; c9_j++) {
      c9_b_j = c9_j;
      c9_c_a = c9_b_j - 1;
      c9_jm1 = c9_c_a;
      c9_d_a = c9_c_m;
      c9_e_b = c9_b_j;
      c9_mmj = c9_d_a - c9_e_b;
      c9_e_a = c9_jm1;
      c9_f_b = c9_ldap1;
      c9_c = c9_e_a * c9_f_b;
      c9_g_b = c9_c + 1;
      c9_jj = c9_g_b;
      c9_f_a = c9_jj + 1;
      c9_jp1j = c9_f_a;
      c9_g_a = c9_mmj;
      c9_b_c = c9_g_a;
      c9_f_n = c9_b_c + 1;
      c9_ix0 = c9_jj;
      c9_g_n = c9_f_n;
      c9_b_ix0 = c9_ix0;
      if (c9_g_n < 1) {
        c9_idxmax = 0;
      } else {
        c9_h_n = c9_g_n;
        c9_c_ix0 = c9_b_ix0;
        c9_var = c9_h_n;
        c9_n_t = (ptrdiff_t)(c9_var);
        c9_incx_t = (ptrdiff_t)(1);
        c9_xix0_t = (double *)(&c9_b_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("",
          c9_c_ix0, 1, c9_b_A_sizes[0] * c9_b_A_sizes[1], 1, 0) - 1]);
        c9_idxmax_t = idamax(&c9_n_t, c9_xix0_t, &c9_incx_t);
        c9_idxmax = (int32_T)(c9_idxmax_t);
      }

      c9_h_a = c9_idxmax - 1;
      c9_jpiv_offset = c9_h_a;
      c9_i_a = c9_jj;
      c9_h_b = c9_jpiv_offset;
      c9_jpiv = c9_i_a + c9_h_b;
      if (c9_b_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_jpiv, 1, c9_b_A_sizes[0]
           * c9_b_A_sizes[1], 1, 0) - 1] != 0.0) {
        if (c9_jpiv_offset != 0) {
          c9_j_a = c9_b_j;
          c9_i_b = c9_jpiv_offset;
          c9_c_c = c9_j_a + c9_i_b;
          c9_ipiv_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1, c9_ipiv_sizes
            [1], 1, 0) - 1] = c9_c_c;
          c9_j_b = c9_jm1 + 1;
          c9_jrow = c9_j_b;
          c9_k_a = c9_jrow;
          c9_k_b = c9_jpiv_offset;
          c9_jprow = c9_k_a + c9_k_b;
          c9_c_eml_xswap(chartInstance, c9_e_n, c9_b_A_data, c9_b_A_sizes,
                         c9_jrow, c9_c_lda, c9_jprow, c9_c_lda);
        }

        c9_b_jp1j = c9_jp1j;
        c9_l_a = c9_mmj;
        c9_d_c = c9_l_a;
        c9_m_a = c9_jp1j;
        c9_l_b = c9_d_c - 1;
        c9_i55 = c9_m_a + c9_l_b;
        c9_n_a = c9_b_jp1j;
        c9_m_b = c9_i55;
        c9_o_a = c9_n_a;
        c9_n_b = c9_m_b;
        if (c9_o_a > c9_n_b) {
          c9_b_overflow = FALSE;
        } else {
          c9_b_overflow = (c9_n_b > 2147483646);
        }

        if (c9_b_overflow) {
          c9_check_forloop_overflow_error(chartInstance, TRUE);
        }

        for (c9_b_i = c9_b_jp1j; c9_b_i <= c9_i55; c9_b_i++) {
          c9_c_i = c9_b_i;
          c9_g_x = c9_b_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1,
            c9_b_A_sizes[0] * c9_b_A_sizes[1], 1, 0) - 1];
          c9_g_y = c9_b_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_jj, 1,
            c9_b_A_sizes[0] * c9_b_A_sizes[1], 1, 0) - 1];
          c9_z = c9_g_x / c9_g_y;
          c9_b_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1, c9_b_A_sizes[0]
            * c9_b_A_sizes[1], 1, 0) - 1] = c9_z;
        }
      } else {
        c9_info = c9_b_j;
      }

      c9_p_a = c9_e_n;
      c9_o_b = c9_b_j;
      c9_e_c = c9_p_a - c9_o_b;
      c9_q_a = c9_jj;
      c9_p_b = c9_c_lda;
      c9_f_c = c9_q_a + c9_p_b;
      c9_r_a = c9_jj;
      c9_q_b = c9_ldap1;
      c9_g_c = c9_r_a + c9_q_b;
      c9_d_m = c9_mmj;
      c9_i_n = c9_e_c;
      c9_d_ix0 = c9_jp1j;
      c9_iy0 = c9_f_c;
      c9_incy = c9_c_lda;
      c9_ia0 = c9_g_c;
      c9_d_lda = c9_c_lda;
      c9_d4 = -1.0;
      c9_c_eml_xger(chartInstance, c9_d_m, c9_i_n, c9_d4, c9_d_ix0, c9_iy0,
                    c9_incy, c9_b_A_data, c9_b_A_sizes, c9_ia0, c9_d_lda);
    }

    if (c9_info == 0) {
      if (c9_c_m <= c9_e_n) {
        if (!(c9_b_A_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_m, 1,
               c9_b_A_sizes[0], 1, 0) + c9_b_A_sizes[0] *
                           (_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_m, 1,
                c9_b_A_sizes[1], 2, 0) - 1)) - 1] != 0.0)) {
          c9_info = c9_c_m;
        }
      }
    }
  }

  c9_b_info = c9_info;
  c9_c_info = c9_b_info;
  c9_d_info = c9_c_info;
  if (c9_d_info > 0) {
    c9_eml_warning(chartInstance);
  }

  c9_eml_scalar_eg(chartInstance);
  c9_X_sizes[0] = c9_B_sizes[0];
  c9_X_sizes[1] = c9_B_sizes[1];
  c9_X = c9_X_sizes[0];
  c9_b_X = c9_X_sizes[1];
  c9_b_loop_ub = c9_B_sizes[0] * c9_B_sizes[1] - 1;
  for (c9_i56 = 0; c9_i56 <= c9_b_loop_ub; c9_i56++) {
    c9_X_data[c9_i56] = c9_B_data[c9_i56];
  }

  c9_nb = (int32_T)(real_T)c9_B_sizes[1];
  c9_j_n = c9_b_n;
  c9_r_b = c9_j_n;
  c9_s_b = c9_r_b;
  if (1 > c9_s_b) {
    c9_c_overflow = FALSE;
  } else {
    c9_c_overflow = (c9_s_b > 2147483646);
  }

  if (c9_c_overflow) {
    c9_check_forloop_overflow_error(chartInstance, TRUE);
  }

  for (c9_d_i = 1; c9_d_i <= c9_j_n; c9_d_i++) {
    c9_e_i = c9_d_i;
    if (c9_ipiv_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_e_i, 1, c9_ipiv_sizes[1],
         1, 0) - 1] != c9_e_i) {
      c9_ip = c9_ipiv_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_e_i, 1,
        c9_ipiv_sizes[1], 1, 0) - 1];
      c9_b_nb = c9_nb;
      c9_t_b = c9_b_nb;
      c9_u_b = c9_t_b;
      if (1 > c9_u_b) {
        c9_d_overflow = FALSE;
      } else {
        c9_d_overflow = (c9_u_b > 2147483646);
      }

      if (c9_d_overflow) {
        c9_check_forloop_overflow_error(chartInstance, TRUE);
      }

      for (c9_c_j = 1; c9_c_j <= c9_b_nb; c9_c_j++) {
        c9_d_j = c9_c_j;
        c9_temp = c9_X_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_e_i, 1,
          c9_X_sizes[0], 1, 0) + c9_X_sizes[0] * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
          c9_d_j, 1, c9_X_sizes[1], 2, 0) - 1)) - 1];
        c9_X_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_e_i, 1, c9_X_sizes[0], 1,
                    0) + c9_X_sizes[0] * (_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_d_j,
          1, c9_X_sizes[1], 2, 0) - 1)) - 1] = c9_X_data
          [(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_ip, 1, c9_X_sizes[0], 1, 0) +
            c9_X_sizes[0] * (_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_d_j, 1,
              c9_X_sizes[1], 2, 0) - 1)) - 1];
        c9_X_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_ip, 1, c9_X_sizes[0], 1, 0)
                   + c9_X_sizes[0] * (_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_d_j, 1,
          c9_X_sizes[1], 2, 0) - 1)) - 1] = c9_temp;
      }
    }
  }

  c9_c_A_sizes[0] = c9_b_A_sizes[0];
  c9_c_A_sizes[1] = c9_b_A_sizes[1];
  c9_c_A = c9_c_A_sizes[0];
  c9_d_A = c9_c_A_sizes[1];
  c9_c_loop_ub = c9_b_A_sizes[0] * c9_b_A_sizes[1] - 1;
  for (c9_i57 = 0; c9_i57 <= c9_c_loop_ub; c9_i57++) {
    c9_c_A_data[c9_i57] = c9_b_A_data[c9_i57];
  }

  c9_c_eml_xtrsm(chartInstance, c9_b_n, c9_nb, c9_c_A_data, c9_c_A_sizes, c9_b_n,
                 c9_X_data, c9_X_sizes, c9_b_n);
  c9_d_A_sizes[0] = c9_b_A_sizes[0];
  c9_d_A_sizes[1] = c9_b_A_sizes[1];
  c9_e_A = c9_d_A_sizes[0];
  c9_f_A = c9_d_A_sizes[1];
  c9_d_loop_ub = c9_b_A_sizes[0] * c9_b_A_sizes[1] - 1;
  for (c9_i58 = 0; c9_i58 <= c9_d_loop_ub; c9_i58++) {
    c9_d_A_data[c9_i58] = c9_b_A_data[c9_i58];
  }

  c9_d_eml_xtrsm(chartInstance, c9_b_n, c9_nb, c9_d_A_data, c9_d_A_sizes, c9_b_n,
                 c9_X_data, c9_X_sizes, c9_b_n);
}

static void c9_realmin(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_eps(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_b_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_eml_signed_integer_colon(SFc9_simulationInstanceStruct
  *chartInstance, int32_T c9_b, int32_T c9_y_data[20], int32_T c9_y_sizes[2])
{
  int32_T c9_b_b;
  int32_T c9_c_b;
  int32_T c9_b_n;
  int32_T c9_a;
  int32_T c9_c;
  int32_T c9_b_a;
  int32_T c9_span;
  int32_T c9_c_a;
  int32_T c9_nm1;
  int32_T c9_d_a;
  int32_T c9_d_b;
  int32_T c9_e_a;
  int32_T c9_tmp_sizes[2];
  int32_T c9_iv9[2];
  int32_T c9_i59;
  int32_T c9_i60;
  int32_T c9_loop_ub;
  int32_T c9_i61;
  int32_T c9_tmp_data[20];
  int32_T c9_i62;
  int32_T c9_yk;
  int32_T c9_c_n;
  int32_T c9_k;
  int32_T c9_b_k;
  int32_T c9_f_a;
  c9_b_b = c9_b;
  c9_c_b = c9_b_b;
  if (c9_c_b < 1) {
    c9_b_n = 0;
  } else {
    c9_a = c9_c_b;
    c9_c = c9_a;
    c9_b_a = c9_c - 1;
    c9_span = c9_b_a;
    c9_c_a = c9_span;
    c9_nm1 = c9_c_a;
    c9_d_a = c9_nm1;
    c9_d_b = c9_d_a;
    c9_e_a = c9_d_b + 1;
    c9_b_n = c9_e_a;
  }

  c9_tmp_sizes[0] = 1;
  c9_iv9[0] = 1;
  c9_iv9[1] = c9_b_n;
  c9_tmp_sizes[1] = c9_iv9[1];
  c9_i59 = c9_tmp_sizes[0];
  c9_i60 = c9_tmp_sizes[1];
  c9_loop_ub = c9_b_n - 1;
  for (c9_i61 = 0; c9_i61 <= c9_loop_ub; c9_i61++) {
    c9_tmp_data[c9_i61] = 0;
  }

  for (c9_i62 = 0; c9_i62 < 2; c9_i62++) {
    c9_y_sizes[c9_i62] = c9_tmp_sizes[c9_i62];
  }

  if (c9_b_n > 0) {
    (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_y_sizes[1], 1, 0);
    c9_y_data[0] = 1;
    c9_yk = 1;
    c9_c_n = c9_b_n;
    for (c9_k = 2; c9_k <= c9_c_n; c9_k++) {
      c9_b_k = c9_k;
      c9_f_a = c9_yk + 1;
      c9_yk = c9_f_a;
      c9_y_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_k, 1, c9_y_sizes[1], 1, 0)
        - 1] = c9_yk;
    }
  }
}

static boolean_T c9_eml_use_refblas(SFc9_simulationInstanceStruct *chartInstance)
{
  return FALSE;
}

static real_T c9_abs(SFc9_simulationInstanceStruct *chartInstance, real_T c9_x)
{
  real_T c9_b_x;
  c9_b_x = c9_x;
  return muDoubleScalarAbs(c9_b_x);
}

static void c9_eml_xswap(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0, int32_T
  c9_incx, int32_T c9_iy0, int32_T c9_incy, real_T c9_b_x_data[400], int32_T
  c9_b_x_sizes[2])
{
  int32_T c9_x;
  int32_T c9_b_x;
  int32_T c9_loop_ub;
  int32_T c9_i63;
  c9_b_x_sizes[0] = c9_x_sizes[0];
  c9_b_x_sizes[1] = c9_x_sizes[1];
  c9_x = c9_b_x_sizes[0];
  c9_b_x = c9_b_x_sizes[1];
  c9_loop_ub = c9_x_sizes[0] * c9_x_sizes[1] - 1;
  for (c9_i63 = 0; c9_i63 <= c9_loop_ub; c9_i63++) {
    c9_b_x_data[c9_i63] = c9_x_data[c9_i63];
  }

  c9_c_eml_xswap(chartInstance, c9_b_n, c9_b_x_data, c9_b_x_sizes, c9_ix0,
                 c9_incx, c9_iy0, c9_incy);
}

static void c9_eml_xger(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_alpha1, int32_T c9_ix0, int32_T c9_iy0,
  int32_T c9_incy, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T c9_ia0,
  int32_T c9_lda, real_T c9_b_A_data[400], int32_T c9_b_A_sizes[2])
{
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_loop_ub;
  int32_T c9_i64;
  c9_b_A_sizes[0] = c9_A_sizes[0];
  c9_b_A_sizes[1] = c9_A_sizes[1];
  c9_A = c9_b_A_sizes[0];
  c9_b_A = c9_b_A_sizes[1];
  c9_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
  for (c9_i64 = 0; c9_i64 <= c9_loop_ub; c9_i64++) {
    c9_b_A_data[c9_i64] = c9_A_data[c9_i64];
  }

  c9_c_eml_xger(chartInstance, c9_m, c9_b_n, c9_alpha1, c9_ix0, c9_iy0, c9_incy,
                c9_b_A_data, c9_b_A_sizes, c9_ia0, c9_lda);
}

static void c9_c_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_eml_warning(SFc9_simulationInstanceStruct *chartInstance)
{
  int32_T c9_i65;
  static char_T c9_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c9_u[27];
  const mxArray *c9_y = NULL;
  for (c9_i65 = 0; c9_i65 < 27; c9_i65++) {
    c9_u[c9_i65] = c9_varargin_1[c9_i65];
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 10, 0U, 1U, 0U, 2, 1, 27), FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
    14, c9_y));
}

static void c9_eml_xtrsm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B_data[400], int32_T c9_B_sizes[2], int32_T c9_ldb, real_T
  c9_b_B_data[400], int32_T c9_b_B_sizes[2])
{
  int32_T c9_B;
  int32_T c9_b_B;
  int32_T c9_loop_ub;
  int32_T c9_i66;
  int32_T c9_b_A_sizes[2];
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_b_loop_ub;
  int32_T c9_i67;
  real_T c9_b_A_data[400];
  c9_b_B_sizes[0] = c9_B_sizes[0];
  c9_b_B_sizes[1] = c9_B_sizes[1];
  c9_B = c9_b_B_sizes[0];
  c9_b_B = c9_b_B_sizes[1];
  c9_loop_ub = c9_B_sizes[0] * c9_B_sizes[1] - 1;
  for (c9_i66 = 0; c9_i66 <= c9_loop_ub; c9_i66++) {
    c9_b_B_data[c9_i66] = c9_B_data[c9_i66];
  }

  c9_b_A_sizes[0] = c9_A_sizes[0];
  c9_b_A_sizes[1] = c9_A_sizes[1];
  c9_A = c9_b_A_sizes[0];
  c9_b_A = c9_b_A_sizes[1];
  c9_b_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
  for (c9_i67 = 0; c9_i67 <= c9_b_loop_ub; c9_i67++) {
    c9_b_A_data[c9_i67] = c9_A_data[c9_i67];
  }

  c9_c_eml_xtrsm(chartInstance, c9_m, c9_b_n, c9_b_A_data, c9_b_A_sizes, c9_lda,
                 c9_b_B_data, c9_b_B_sizes, c9_ldb);
}

static void c9_d_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_b_eml_xtrsm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B_data[400], int32_T c9_B_sizes[2], int32_T c9_ldb, real_T
  c9_b_B_data[400], int32_T c9_b_B_sizes[2])
{
  int32_T c9_B;
  int32_T c9_b_B;
  int32_T c9_loop_ub;
  int32_T c9_i68;
  int32_T c9_b_A_sizes[2];
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_b_loop_ub;
  int32_T c9_i69;
  real_T c9_b_A_data[400];
  c9_b_B_sizes[0] = c9_B_sizes[0];
  c9_b_B_sizes[1] = c9_B_sizes[1];
  c9_B = c9_b_B_sizes[0];
  c9_b_B = c9_b_B_sizes[1];
  c9_loop_ub = c9_B_sizes[0] * c9_B_sizes[1] - 1;
  for (c9_i68 = 0; c9_i68 <= c9_loop_ub; c9_i68++) {
    c9_b_B_data[c9_i68] = c9_B_data[c9_i68];
  }

  c9_b_A_sizes[0] = c9_A_sizes[0];
  c9_b_A_sizes[1] = c9_A_sizes[1];
  c9_A = c9_b_A_sizes[0];
  c9_b_A = c9_b_A_sizes[1];
  c9_b_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
  for (c9_i69 = 0; c9_i69 <= c9_b_loop_ub; c9_i69++) {
    c9_b_A_data[c9_i69] = c9_A_data[c9_i69];
  }

  c9_d_eml_xtrsm(chartInstance, c9_m, c9_b_n, c9_b_A_data, c9_b_A_sizes, c9_lda,
                 c9_b_B_data, c9_b_B_sizes, c9_ldb);
}

static void c9_eml_qrsolve(SFc9_simulationInstanceStruct *chartInstance, real_T
  c9_A_data[400], int32_T c9_A_sizes[2], real_T c9_B_data[400], int32_T
  c9_B_sizes[2], real_T c9_Y_data[400], int32_T c9_Y_sizes[2])
{
  real_T c9_m;
  real_T c9_b_n;
  real_T c9_nb;
  real_T c9_varargin_1;
  real_T c9_varargin_2;
  real_T c9_b_varargin_2;
  real_T c9_varargin_3;
  real_T c9_x;
  real_T c9_y;
  real_T c9_b_x;
  real_T c9_b_y;
  real_T c9_xk;
  real_T c9_yk;
  real_T c9_c_x;
  real_T c9_c_y;
  real_T c9_mn;
  int32_T c9_b_m;
  int32_T c9_c_n;
  int32_T c9_b_varargin_1;
  int32_T c9_c_varargin_2;
  int32_T c9_d_varargin_2;
  int32_T c9_b_varargin_3;
  int32_T c9_d_x;
  int32_T c9_d_y;
  int32_T c9_e_x;
  int32_T c9_e_y;
  int32_T c9_b_xk;
  int32_T c9_b_yk;
  int32_T c9_f_x;
  int32_T c9_f_y;
  int32_T c9_b_mn;
  int32_T c9_iv10[2];
  int32_T c9_vn2_sizes;
  int32_T c9_loop_ub;
  int32_T c9_i70;
  real_T c9_vn2_data[20];
  int32_T c9_tau_sizes;
  int32_T c9_d;
  int32_T c9_b;
  int32_T c9_b_b;
  int32_T c9_a;
  int32_T c9_bi;
  int32_T c9_jpvt_sizes[2];
  int32_T c9_jpvt_data[20];
  boolean_T c9_b7;
  boolean_T c9_b8;
  int32_T c9_work_sizes;
  int32_T c9_b_loop_ub;
  int32_T c9_i71;
  real_T c9_work_data[20];
  int32_T c9_c_loop_ub;
  int32_T c9_i72;
  int32_T c9_vn1_sizes;
  int32_T c9_k;
  int32_T c9_d_n;
  int32_T c9_c_b;
  int32_T c9_d_b;
  boolean_T c9_overflow;
  int32_T c9_j;
  int32_T c9_b_j;
  int32_T c9_b_A_sizes[2];
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_d_loop_ub;
  int32_T c9_i73;
  real_T c9_b_A_data[400];
  real_T c9_vn1_data[20];
  int32_T c9_b_a;
  int32_T c9_e_b;
  int32_T c9_c_mn;
  int32_T c9_f_b;
  int32_T c9_g_b;
  boolean_T c9_b_overflow;
  int32_T c9_b_i;
  int32_T c9_c_i;
  int32_T c9_c_a;
  int32_T c9_im1;
  int32_T c9_d_a;
  int32_T c9_ip1;
  int32_T c9_e_a;
  int32_T c9_h_b;
  int32_T c9_c;
  int32_T c9_f_a;
  int32_T c9_i_b;
  int32_T c9_i_i;
  int32_T c9_g_a;
  int32_T c9_j_b;
  int32_T c9_nmi;
  int32_T c9_h_a;
  int32_T c9_k_b;
  int32_T c9_mmi;
  int32_T c9_l_b;
  int32_T c9_mmip1;
  int32_T c9_m_b;
  int32_T c9_nmip1;
  int32_T c9_i_a;
  int32_T c9_b_vn1_sizes;
  int32_T c9_e_loop_ub;
  int32_T c9_i74;
  real_T c9_b_vn1_data[20];
  int32_T c9_n_b;
  int32_T c9_pvt;
  int32_T c9_j_a;
  int32_T c9_b_c;
  int32_T c9_k_a;
  int32_T c9_o_b;
  int32_T c9_c_c;
  int32_T c9_p_b;
  int32_T c9_pvtcol;
  int32_T c9_l_a;
  int32_T c9_q_b;
  int32_T c9_d_c;
  int32_T c9_r_b;
  int32_T c9_mcol;
  int32_T c9_itemp;
  real_T c9_atmp;
  int32_T c9_m_a;
  int32_T c9_e_c;
  real_T c9_b_atmp;
  real_T c9_d5;
  real_T c9_tau_data[20];
  real_T c9_c_atmp;
  real_T c9_d6;
  real_T c9_d7;
  int32_T c9_n_a;
  int32_T c9_s_b;
  int32_T c9_f_c;
  int32_T c9_o_a;
  int32_T c9_t_b;
  int32_T c9_i_ip1;
  int32_T c9_b_ip1;
  int32_T c9_e_n;
  int32_T c9_p_a;
  int32_T c9_u_b;
  int32_T c9_q_a;
  int32_T c9_v_b;
  boolean_T c9_c_overflow;
  int32_T c9_c_j;
  int32_T c9_r_a;
  int32_T c9_g_c;
  int32_T c9_s_a;
  int32_T c9_w_b;
  int32_T c9_h_c;
  int32_T c9_t_a;
  int32_T c9_x_b;
  int32_T c9_i_j;
  real_T c9_temp1;
  real_T c9_u_a;
  real_T c9_y_b;
  real_T c9_g_y;
  real_T c9_temp2;
  real_T c9_v_a;
  real_T c9_ab_b;
  real_T c9_h_y;
  real_T c9_w_a;
  real_T c9_bb_b;
  int32_T c9_x_a;
  int32_T c9_i_c;
  int32_T c9_c_A_sizes[2];
  int32_T c9_c_A;
  int32_T c9_d_A;
  int32_T c9_f_loop_ub;
  int32_T c9_i75;
  real_T c9_c_A_data[400];
  real_T c9_y_a;
  real_T c9_cb_b;
  real_T c9_i_y;
  real_T c9_rankR;
  real_T c9_c_varargin_1;
  real_T c9_e_varargin_2;
  real_T c9_f_varargin_2;
  real_T c9_c_varargin_3;
  real_T c9_g_x;
  real_T c9_j_y;
  real_T c9_h_x;
  real_T c9_k_y;
  real_T c9_c_xk;
  real_T c9_c_yk;
  real_T c9_i_x;
  real_T c9_l_y;
  real_T c9_maxval;
  real_T c9_j_x;
  real_T c9_k_x;
  real_T c9_l_x;
  real_T c9_m_y;
  real_T c9_m_x;
  real_T c9_n_x;
  real_T c9_n_y;
  real_T c9_b_d;
  real_T c9_ab_a;
  real_T c9_db_b;
  real_T c9_o_y;
  real_T c9_bb_a;
  real_T c9_tol;
  real_T c9_d_mn;
  int32_T c9_i76;
  int32_T c9_b_k;
  real_T c9_c_k;
  real_T c9_o_x;
  real_T c9_p_x;
  real_T c9_q_x;
  real_T c9_p_y;
  real_T c9_r_x;
  real_T c9_s_x;
  real_T c9_q_y;
  real_T c9_c_d;
  char_T c9_cv6[14];
  int32_T c9_i77;
  char_T c9_cv7[14];
  real_T c9_dv2[2];
  int32_T c9_iv11[2];
  int32_T c9_iv12[2];
  int32_T c9_Y;
  int32_T c9_b_Y;
  int32_T c9_g_loop_ub;
  int32_T c9_i78;
  real_T c9_e_mn;
  int32_T c9_i79;
  int32_T c9_d_j;
  real_T c9_e_j;
  real_T c9_tauj;
  real_T c9_b_nb;
  int32_T c9_i80;
  int32_T c9_d_k;
  real_T c9_wj;
  real_T c9_d8;
  real_T c9_c_m;
  int32_T c9_i81;
  int32_T c9_d_i;
  real_T c9_e_i;
  real_T c9_cb_a;
  real_T c9_eb_b;
  real_T c9_z;
  real_T c9_db_a;
  real_T c9_fb_b;
  real_T c9_d9;
  real_T c9_d_m;
  int32_T c9_i82;
  int32_T c9_f_i;
  real_T c9_eb_a;
  real_T c9_gb_b;
  real_T c9_r_y;
  real_T c9_rr;
  real_T c9_c_nb;
  int32_T c9_i83;
  int32_T c9_e_k;
  real_T c9_b_rr;
  int32_T c9_i84;
  int32_T c9_g_i;
  real_T c9_c_rr;
  int32_T c9_i85;
  int32_T c9_f_j;
  int32_T c9_pj;
  real_T c9_t_x;
  real_T c9_s_y;
  real_T c9_b_z;
  real_T c9_d10;
  int32_T c9_i86;
  int32_T c9_h_i;
  real_T c9_fb_a;
  real_T c9_hb_b;
  real_T c9_t_y;
  boolean_T exitg1;
  c9_m = (real_T)c9_A_sizes[0];
  c9_b_n = (real_T)c9_A_sizes[1];
  c9_nb = (real_T)c9_B_sizes[1];
  c9_varargin_1 = c9_m;
  c9_varargin_2 = c9_b_n;
  c9_b_varargin_2 = c9_varargin_1;
  c9_varargin_3 = c9_varargin_2;
  c9_x = c9_b_varargin_2;
  c9_y = c9_varargin_3;
  c9_b_x = c9_x;
  c9_b_y = c9_y;
  c9_c_eml_scalar_eg(chartInstance);
  c9_xk = c9_b_x;
  c9_yk = c9_b_y;
  c9_c_x = c9_xk;
  c9_c_y = c9_yk;
  c9_c_eml_scalar_eg(chartInstance);
  c9_mn = muDoubleScalarMin(c9_c_x, c9_c_y);
  c9_b_m = (int32_T)(real_T)c9_A_sizes[0];
  c9_c_n = (int32_T)(real_T)c9_A_sizes[1];
  c9_b_varargin_1 = c9_b_m;
  c9_c_varargin_2 = c9_c_n;
  c9_d_varargin_2 = c9_b_varargin_1;
  c9_b_varargin_3 = c9_c_varargin_2;
  c9_d_x = c9_d_varargin_2;
  c9_d_y = c9_b_varargin_3;
  c9_e_x = c9_d_x;
  c9_e_y = c9_d_y;
  c9_b_eml_scalar_eg(chartInstance);
  c9_b_xk = c9_e_x;
  c9_b_yk = c9_e_y;
  c9_f_x = c9_b_xk;
  c9_f_y = c9_b_yk;
  c9_b_eml_scalar_eg(chartInstance);
  c9_b_mn = muIntScalarMin_sint32(c9_f_x, c9_f_y);
  c9_d_eml_scalar_eg(chartInstance);
  c9_iv10[0] = c9_b_mn;
  c9_iv10[1] = 1;
  c9_vn2_sizes = c9_iv10[0];
  c9_loop_ub = c9_iv10[0] - 1;
  for (c9_i70 = 0; c9_i70 <= c9_loop_ub; c9_i70++) {
    c9_vn2_data[c9_i70] = 0.0;
  }

  c9_tau_sizes = c9_vn2_sizes;
  c9_d = c9_c_n;
  c9_b = c9_d;
  c9_b_b = c9_b;
  c9_a = c9_b_b;
  c9_bi = c9_a;
  c9_eml_signed_integer_colon(chartInstance, c9_bi, c9_jpvt_data, c9_jpvt_sizes);
  c9_b7 = (c9_A_sizes[0] == 0);
  c9_b8 = (c9_A_sizes[1] == 0);
  if (c9_b7 || c9_b8) {
  } else {
    c9_d_eml_scalar_eg(chartInstance);
    c9_iv10[0] = c9_c_n;
    c9_iv10[1] = 1;
    c9_work_sizes = c9_iv10[0];
    c9_b_loop_ub = c9_iv10[0] - 1;
    for (c9_i71 = 0; c9_i71 <= c9_b_loop_ub; c9_i71++) {
      c9_work_data[c9_i71] = 0.0;
    }

    c9_eps(chartInstance);
    c9_sqrt(chartInstance, 2.2204460492503131E-16);
    c9_d_eml_scalar_eg(chartInstance);
    c9_iv10[0] = c9_c_n;
    c9_iv10[1] = 1;
    c9_vn2_sizes = c9_iv10[0];
    c9_c_loop_ub = c9_iv10[0] - 1;
    for (c9_i72 = 0; c9_i72 <= c9_c_loop_ub; c9_i72++) {
      c9_vn2_data[c9_i72] = 0.0;
    }

    c9_vn1_sizes = c9_vn2_sizes;
    c9_vn2_sizes = c9_vn1_sizes;
    c9_k = 1;
    c9_d_n = c9_c_n;
    c9_c_b = c9_d_n;
    c9_d_b = c9_c_b;
    if (1 > c9_d_b) {
      c9_overflow = FALSE;
    } else {
      c9_overflow = (c9_d_b > 2147483646);
    }

    if (c9_overflow) {
      c9_check_forloop_overflow_error(chartInstance, TRUE);
    }

    for (c9_j = 1; c9_j <= c9_d_n; c9_j++) {
      c9_b_j = c9_j;
      c9_b_A_sizes[0] = c9_A_sizes[0];
      c9_b_A_sizes[1] = c9_A_sizes[1];
      c9_A = c9_b_A_sizes[0];
      c9_b_A = c9_b_A_sizes[1];
      c9_d_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
      for (c9_i73 = 0; c9_i73 <= c9_d_loop_ub; c9_i73++) {
        c9_b_A_data[c9_i73] = c9_A_data[c9_i73];
      }

      c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1, c9_vn1_sizes, 1, 0)
        - 1] = c9_eml_xnrm2(chartInstance, c9_b_m, c9_b_A_data, c9_b_A_sizes,
                            c9_k);
      c9_vn2_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1, c9_vn2_sizes, 1, 0)
        - 1] = c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
        c9_vn1_sizes, 1, 0) - 1];
      c9_b_a = c9_k;
      c9_e_b = c9_b_m;
      c9_k = c9_b_a + c9_e_b;
    }

    c9_c_mn = c9_b_mn;
    c9_f_b = c9_c_mn;
    c9_g_b = c9_f_b;
    if (1 > c9_g_b) {
      c9_b_overflow = FALSE;
    } else {
      c9_b_overflow = (c9_g_b > 2147483646);
    }

    if (c9_b_overflow) {
      c9_check_forloop_overflow_error(chartInstance, TRUE);
    }

    for (c9_b_i = 1; c9_b_i <= c9_c_mn; c9_b_i++) {
      c9_c_i = c9_b_i;
      c9_c_a = c9_c_i - 1;
      c9_im1 = c9_c_a;
      c9_d_a = c9_c_i + 1;
      c9_ip1 = c9_d_a;
      c9_e_a = c9_im1;
      c9_h_b = c9_b_m;
      c9_c = c9_e_a * c9_h_b;
      c9_f_a = c9_c_i;
      c9_i_b = c9_c;
      c9_i_i = c9_f_a + c9_i_b;
      c9_g_a = c9_c_n;
      c9_j_b = c9_c_i;
      c9_nmi = c9_g_a - c9_j_b;
      c9_h_a = c9_b_m;
      c9_k_b = c9_c_i;
      c9_mmi = c9_h_a - c9_k_b;
      c9_l_b = c9_mmi + 1;
      c9_mmip1 = c9_l_b;
      c9_m_b = c9_nmi + 1;
      c9_nmip1 = c9_m_b;
      c9_i_a = c9_im1;
      c9_b_vn1_sizes = c9_vn1_sizes;
      c9_e_loop_ub = c9_vn1_sizes - 1;
      for (c9_i74 = 0; c9_i74 <= c9_e_loop_ub; c9_i74++) {
        c9_b_vn1_data[c9_i74] = c9_vn1_data[c9_i74];
      }

      c9_n_b = c9_eml_ixamax(chartInstance, c9_nmip1, c9_b_vn1_data, *(int32_T (*)
        [1])&c9_b_vn1_sizes, c9_c_i);
      c9_pvt = c9_i_a + c9_n_b;
      if (c9_pvt != c9_c_i) {
        c9_j_a = c9_pvt - 1;
        c9_b_c = c9_j_a;
        c9_k_a = c9_b_m;
        c9_o_b = c9_b_c;
        c9_c_c = c9_k_a * c9_o_b;
        c9_p_b = c9_c_c + 1;
        c9_pvtcol = c9_p_b;
        c9_l_a = c9_b_m;
        c9_q_b = c9_im1;
        c9_d_c = c9_l_a * c9_q_b;
        c9_r_b = c9_d_c + 1;
        c9_mcol = c9_r_b;
        c9_d_eml_xswap(chartInstance, c9_b_m, c9_A_data, c9_A_sizes, c9_pvtcol,
                       c9_mcol);
        c9_itemp = c9_jpvt_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_pvt, 1,
          c9_jpvt_sizes[1], 1, 0) - 1];
        c9_jpvt_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_pvt, 1, c9_jpvt_sizes[1],
          1, 0) - 1] = c9_jpvt_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1,
          c9_jpvt_sizes[1], 1, 0) - 1];
        c9_jpvt_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1, c9_jpvt_sizes[1],
          1, 0) - 1] = c9_itemp;
        c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_pvt, 1, c9_vn1_sizes, 1,
          0) - 1] = c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1,
          c9_vn1_sizes, 1, 0) - 1];
        c9_vn2_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_pvt, 1, c9_vn2_sizes, 1,
          0) - 1] = c9_vn2_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1,
          c9_vn2_sizes, 1, 0) - 1];
      }

      c9_atmp = c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_i_i, 1, c9_A_sizes
        [0] * c9_A_sizes[1], 1, 0) - 1];
      if (c9_c_i < c9_b_m) {
        c9_m_a = c9_i_i + 1;
        c9_e_c = c9_m_a;
        c9_b_atmp = c9_atmp;
        c9_d5 = c9_c_eml_matlab_zlarfg(chartInstance, c9_mmip1, &c9_b_atmp,
          c9_A_data, c9_A_sizes, c9_e_c);
        c9_atmp = c9_b_atmp;
        c9_tau_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1, c9_tau_sizes, 1,
          0) - 1] = c9_d5;
      } else {
        c9_c_atmp = c9_atmp;
        c9_d6 = c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_i_i, 1, c9_A_sizes
          [0] * c9_A_sizes[1], 1, 0) - 1];
        c9_d7 = c9_d_eml_matlab_zlarfg(chartInstance, &c9_c_atmp, &c9_d6);
        c9_atmp = c9_c_atmp;
        c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_i_i, 1, c9_A_sizes[0] *
          c9_A_sizes[1], 1, 0) - 1] = c9_d6;
        c9_tau_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1, c9_tau_sizes, 1,
          0) - 1] = c9_d7;
      }

      c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_i_i, 1, c9_A_sizes[0] *
        c9_A_sizes[1], 1, 0) - 1] = c9_atmp;
      if (c9_c_i < c9_c_n) {
        c9_atmp = c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_i_i, 1,
          c9_A_sizes[0] * c9_A_sizes[1], 1, 0) - 1];
        c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_i_i, 1, c9_A_sizes[0] *
          c9_A_sizes[1], 1, 0) - 1] = 1.0;
        c9_n_a = c9_c_i;
        c9_s_b = c9_b_m;
        c9_f_c = c9_n_a * c9_s_b;
        c9_o_a = c9_c_i;
        c9_t_b = c9_f_c;
        c9_i_ip1 = c9_o_a + c9_t_b;
        c9_b_eml_matlab_zlarf(chartInstance, c9_mmip1, c9_nmi, c9_i_i,
                              c9_tau_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i,
          1, c9_tau_sizes, 1, 0) - 1], c9_A_data, c9_A_sizes, c9_i_ip1, c9_b_m,
                              c9_work_data, &c9_work_sizes);
        c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_i_i, 1, c9_A_sizes[0] *
          c9_A_sizes[1], 1, 0) - 1] = c9_atmp;
      }

      c9_b_ip1 = c9_ip1;
      c9_e_n = c9_c_n;
      c9_p_a = c9_b_ip1;
      c9_u_b = c9_e_n;
      c9_q_a = c9_p_a;
      c9_v_b = c9_u_b;
      if (c9_q_a > c9_v_b) {
        c9_c_overflow = FALSE;
      } else {
        c9_c_overflow = (c9_v_b > 2147483646);
      }

      if (c9_c_overflow) {
        c9_check_forloop_overflow_error(chartInstance, TRUE);
      }

      for (c9_c_j = c9_b_ip1; c9_c_j <= c9_e_n; c9_c_j++) {
        c9_b_j = c9_c_j;
        c9_r_a = c9_b_j - 1;
        c9_g_c = c9_r_a;
        c9_s_a = c9_b_m;
        c9_w_b = c9_g_c;
        c9_h_c = c9_s_a * c9_w_b;
        c9_t_a = c9_c_i;
        c9_x_b = c9_h_c;
        c9_i_j = c9_t_a + c9_x_b;
        if (c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1, c9_vn1_sizes,
             1, 0) - 1] != 0.0) {
          c9_temp1 = c9_abs(chartInstance, c9_A_data
                            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_i, 1,
            c9_A_sizes[0], 1, 0) + c9_A_sizes[0] * (_SFD_EML_ARRAY_BOUNDS_CHECK(
            "", c9_b_j, 1, c9_A_sizes[1], 2, 0) - 1)) - 1]) /
            c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1, c9_vn1_sizes,
            1, 0) - 1];
          c9_u_a = c9_temp1;
          c9_y_b = c9_temp1;
          c9_g_y = c9_u_a * c9_y_b;
          c9_temp1 = 1.0 - c9_g_y;
          if (c9_temp1 < 0.0) {
            c9_temp1 = 0.0;
          }

          c9_temp2 = c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
            c9_vn1_sizes, 1, 0) - 1] / c9_vn2_data[_SFD_EML_ARRAY_BOUNDS_CHECK(
            "", c9_b_j, 1, c9_vn2_sizes, 1, 0) - 1];
          c9_v_a = c9_temp2;
          c9_ab_b = c9_temp2;
          c9_h_y = c9_v_a * c9_ab_b;
          c9_w_a = c9_temp1;
          c9_bb_b = c9_h_y;
          c9_temp2 = c9_w_a * c9_bb_b;
          if (c9_temp2 <= 1.4901161193847656E-8) {
            if (c9_c_i < c9_b_m) {
              c9_x_a = c9_i_j + 1;
              c9_i_c = c9_x_a;
              c9_c_A_sizes[0] = c9_A_sizes[0];
              c9_c_A_sizes[1] = c9_A_sizes[1];
              c9_c_A = c9_c_A_sizes[0];
              c9_d_A = c9_c_A_sizes[1];
              c9_f_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
              for (c9_i75 = 0; c9_i75 <= c9_f_loop_ub; c9_i75++) {
                c9_c_A_data[c9_i75] = c9_A_data[c9_i75];
              }

              c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
                c9_vn1_sizes, 1, 0) - 1] = c9_c_eml_xnrm2(chartInstance, c9_mmi,
                c9_c_A_data, c9_c_A_sizes, c9_i_c);
              c9_vn2_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
                c9_vn2_sizes, 1, 0) - 1] =
                c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
                c9_vn1_sizes, 1, 0) - 1];
            } else {
              c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
                c9_vn1_sizes, 1, 0) - 1] = 0.0;
              c9_vn2_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
                c9_vn2_sizes, 1, 0) - 1] = 0.0;
            }
          } else {
            c9_y_a = c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1,
              c9_vn1_sizes, 1, 0) - 1];
            c9_cb_b = c9_temp1;
            c9_b_sqrt(chartInstance, &c9_cb_b);
            c9_i_y = c9_y_a * c9_cb_b;
            c9_vn1_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_j, 1, c9_vn1_sizes,
              1, 0) - 1] = c9_i_y;
          }
        }
      }
    }
  }

  c9_rankR = 0.0;
  if (c9_mn > 0.0) {
    c9_eps(chartInstance);
    (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_A_sizes[0], 1, 0);
    (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_A_sizes[1], 2, 0);
    c9_c_varargin_1 = c9_m;
    c9_e_varargin_2 = c9_b_n;
    c9_f_varargin_2 = c9_c_varargin_1;
    c9_c_varargin_3 = c9_e_varargin_2;
    c9_g_x = c9_f_varargin_2;
    c9_j_y = c9_c_varargin_3;
    c9_h_x = c9_g_x;
    c9_k_y = c9_j_y;
    c9_c_eml_scalar_eg(chartInstance);
    c9_c_xk = c9_h_x;
    c9_c_yk = c9_k_y;
    c9_i_x = c9_c_xk;
    c9_l_y = c9_c_yk;
    c9_c_eml_scalar_eg(chartInstance);
    c9_maxval = muDoubleScalarMax(c9_i_x, c9_l_y);
    c9_j_x = c9_A_data[0];
    c9_k_x = c9_j_x;
    c9_l_x = c9_k_x;
    c9_m_y = muDoubleScalarAbs(c9_l_x);
    c9_m_x = 0.0;
    c9_n_x = c9_m_x;
    c9_n_y = muDoubleScalarAbs(c9_n_x);
    c9_b_d = c9_m_y + c9_n_y;
    c9_ab_a = c9_maxval;
    c9_db_b = c9_b_d;
    c9_o_y = c9_ab_a * c9_db_b;
    c9_bb_a = c9_o_y;
    c9_tol = c9_bb_a * 2.2204460492503131E-16;
    c9_d_mn = c9_mn;
    c9_i76 = (int32_T)c9_d_mn - 1;
    c9_b_k = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (c9_b_k <= c9_i76)) {
      c9_c_k = 1.0 + (real_T)c9_b_k;
      c9_o_x = c9_A_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)c9_c_k, 1, c9_A_sizes[0], 1, 0) + c9_A_sizes[0] * ((int32_T)
        (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_c_k, 1, c9_A_sizes[1],
        2, 0) - 1)) - 1];
      c9_p_x = c9_o_x;
      c9_q_x = c9_p_x;
      c9_p_y = muDoubleScalarAbs(c9_q_x);
      c9_r_x = 0.0;
      c9_s_x = c9_r_x;
      c9_q_y = muDoubleScalarAbs(c9_s_x);
      c9_c_d = c9_p_y + c9_q_y;
      if (c9_c_d <= c9_tol) {
        c9_eml_flt2str(chartInstance, c9_tol, c9_cv6);
        for (c9_i77 = 0; c9_i77 < 14; c9_i77++) {
          c9_cv7[c9_i77] = c9_cv6[c9_i77];
        }

        c9_b_eml_warning(chartInstance, c9_rankR, c9_cv7);
        exitg1 = TRUE;
      } else {
        c9_rankR++;
        c9_b_k++;
      }
    }
  }

  c9_eml_scalar_eg(chartInstance);
  c9_dv2[0] = c9_b_n;
  c9_dv2[1] = c9_nb;
  c9_iv11[0] = (int32_T)c9_dv2[0];
  c9_iv11[1] = (int32_T)c9_dv2[1];
  c9_Y_sizes[0] = c9_iv11[0];
  c9_iv12[0] = (int32_T)c9_dv2[0];
  c9_iv12[1] = (int32_T)c9_dv2[1];
  c9_Y_sizes[1] = c9_iv12[1];
  c9_Y = c9_Y_sizes[0];
  c9_b_Y = c9_Y_sizes[1];
  c9_g_loop_ub = (int32_T)c9_dv2[0] * (int32_T)c9_dv2[1] - 1;
  for (c9_i78 = 0; c9_i78 <= c9_g_loop_ub; c9_i78++) {
    c9_Y_data[c9_i78] = 0.0;
  }

  c9_e_mn = c9_mn;
  c9_i79 = (int32_T)c9_e_mn - 1;
  for (c9_d_j = 0; c9_d_j <= c9_i79; c9_d_j++) {
    c9_e_j = 1.0 + (real_T)c9_d_j;
    c9_tauj = c9_tau_data[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)c9_e_j, 1, c9_tau_sizes, 1, 0) - 1];
    if (c9_tauj != 0.0) {
      c9_b_nb = c9_nb;
      c9_i80 = (int32_T)c9_b_nb - 1;
      for (c9_d_k = 0; c9_d_k <= c9_i80; c9_d_k++) {
        c9_c_k = 1.0 + (real_T)c9_d_k;
        c9_wj = c9_B_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)c9_e_j, 1, c9_B_sizes[0], 1, 0) + c9_B_sizes[0] * ((int32_T)
          (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_c_k, 1,
          c9_B_sizes[1], 2, 0) - 1)) - 1];
        c9_d8 = c9_e_j + 1.0;
        c9_c_m = c9_m;
        c9_i81 = (int32_T)(c9_c_m + (1.0 - c9_d8)) - 1;
        for (c9_d_i = 0; c9_d_i <= c9_i81; c9_d_i++) {
          c9_e_i = c9_d8 + (real_T)c9_d_i;
          c9_cb_a = c9_A_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
                                (int32_T)c9_e_i, 1, c9_A_sizes[0], 1, 0) +
                               c9_A_sizes[0] * ((int32_T)(real_T)
            _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_j, 1, c9_A_sizes[1], 2,
            0) - 1)) - 1];
          c9_eb_b = c9_B_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
                                (int32_T)c9_e_i, 1, c9_B_sizes[0], 1, 0) +
                               c9_B_sizes[0] * ((int32_T)(real_T)
            _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_c_k, 1, c9_B_sizes[1], 2,
            0) - 1)) - 1];
          c9_z = c9_cb_a * c9_eb_b;
          c9_wj += c9_z;
        }

        c9_db_a = c9_tauj;
        c9_fb_b = c9_wj;
        c9_wj = c9_db_a * c9_fb_b;
        if (c9_wj != 0.0) {
          c9_B_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                      c9_e_j, 1, c9_B_sizes[0], 1, 0) + c9_B_sizes[0] *
                     ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            c9_c_k, 1, c9_B_sizes[1], 2, 0) - 1)) - 1] = c9_B_data[((int32_T)
            (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_j, 1,
            c9_B_sizes[0], 1, 0) + c9_B_sizes[0] * ((int32_T)(real_T)
            _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_c_k, 1, c9_B_sizes[1], 2,
            0) - 1)) - 1] - c9_wj;
          c9_d9 = c9_e_j + 1.0;
          c9_d_m = c9_m;
          c9_i82 = (int32_T)(c9_d_m + (1.0 - c9_d9)) - 1;
          for (c9_f_i = 0; c9_f_i <= c9_i82; c9_f_i++) {
            c9_e_i = c9_d9 + (real_T)c9_f_i;
            c9_eb_a = c9_A_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
                                  (int32_T)c9_e_i, 1, c9_A_sizes[0], 1, 0) +
                                 c9_A_sizes[0] * ((int32_T)(real_T)
              _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_j, 1, c9_A_sizes[1],
              2, 0) - 1)) - 1];
            c9_gb_b = c9_wj;
            c9_r_y = c9_eb_a * c9_gb_b;
            c9_B_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                        c9_e_i, 1, c9_B_sizes[0], 1, 0) + c9_B_sizes[0] *
                       ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)c9_c_k, 1, c9_B_sizes[1], 2, 0) - 1)) - 1] = c9_B_data
              [((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_i,
                 1, c9_B_sizes[0], 1, 0) + c9_B_sizes[0] * ((int32_T)(real_T)
                 _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_c_k, 1, c9_B_sizes
                  [1], 2, 0) - 1)) - 1] - c9_r_y;
          }
        }
      }
    }
  }

  c9_rr = c9_rankR;
  c9_c_nb = c9_nb;
  c9_i83 = (int32_T)c9_c_nb - 1;
  for (c9_e_k = 0; c9_e_k <= c9_i83; c9_e_k++) {
    c9_c_k = 1.0 + (real_T)c9_e_k;
    c9_b_rr = c9_rr;
    c9_i84 = (int32_T)c9_b_rr;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c9_b_rr, mxDOUBLE_CLASS, c9_i84);
    for (c9_g_i = 0; c9_g_i < c9_i84; c9_g_i++) {
      c9_e_i = 1.0 + (real_T)c9_g_i;
      c9_Y_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_jpvt_data[(int32_T)(real_T)
                  _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_i, 1,
        c9_jpvt_sizes[1], 1, 0) - 1], 1, c9_Y_sizes[0], 1, 0) + c9_Y_sizes[0] *
                 ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        c9_c_k, 1, c9_Y_sizes[1], 2, 0) - 1)) - 1] = c9_B_data[((int32_T)(real_T)
        _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_i, 1, c9_B_sizes[0], 1, 0)
        + c9_B_sizes[0] * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)c9_c_k, 1, c9_B_sizes[1], 2, 0) - 1)) - 1];
    }

    c9_c_rr = c9_rr;
    c9_i85 = (int32_T)-(1.0 + (-1.0 - c9_c_rr));
    _SFD_FOR_LOOP_VECTOR_CHECK(c9_c_rr, -1.0, 1.0, mxDOUBLE_CLASS, c9_i85);
    for (c9_f_j = 0; c9_f_j < c9_i85; c9_f_j++) {
      c9_e_j = c9_c_rr + -(real_T)c9_f_j;
      c9_pj = c9_jpvt_data[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)c9_e_j, 1, c9_jpvt_sizes[1], 1, 0) - 1];
      c9_t_x = c9_Y_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_pj, 1, c9_Y_sizes[0],
        1, 0) + c9_Y_sizes[0] * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
                            (int32_T)c9_c_k, 1, c9_Y_sizes[1], 2, 0) - 1)) - 1];
      c9_s_y = c9_A_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)c9_e_j, 1, c9_A_sizes[0], 1, 0) + c9_A_sizes[0] * ((int32_T)
        (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_j, 1, c9_A_sizes[1],
        2, 0) - 1)) - 1];
      c9_b_z = c9_t_x / c9_s_y;
      c9_Y_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_pj, 1, c9_Y_sizes[0], 1, 0)
                 + c9_Y_sizes[0] * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK
                  ("", (int32_T)c9_c_k, 1, c9_Y_sizes[1], 2, 0) - 1)) - 1] =
        c9_b_z;
      c9_d10 = c9_e_j - 1.0;
      c9_i86 = (int32_T)c9_d10 - 1;
      for (c9_h_i = 0; c9_h_i <= c9_i86; c9_h_i++) {
        c9_e_i = 1.0 + (real_T)c9_h_i;
        c9_fb_a = c9_Y_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_pj, 1,
          c9_Y_sizes[0], 1, 0) + c9_Y_sizes[0] * ((int32_T)(real_T)
          _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_c_k, 1, c9_Y_sizes[1], 2,
          0) - 1)) - 1];
        c9_hb_b = c9_A_data[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)c9_e_i, 1, c9_A_sizes[0], 1, 0) + c9_A_sizes[0] * ((int32_T)
          (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_j, 1,
          c9_A_sizes[1], 2, 0) - 1)) - 1];
        c9_t_y = c9_fb_a * c9_hb_b;
        c9_Y_data[(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_jpvt_data[(int32_T)(real_T)
                    _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_i, 1,
          c9_jpvt_sizes[1], 1, 0) - 1], 1, c9_Y_sizes[0], 1, 0) + c9_Y_sizes[0] *
                   ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          c9_c_k, 1, c9_Y_sizes[1], 2, 0) - 1)) - 1] = c9_Y_data
          [(_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_jpvt_data[(int32_T)(real_T)
             _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_e_i, 1, c9_jpvt_sizes[1],
              1, 0) - 1], 1, c9_Y_sizes[0], 1, 0) + c9_Y_sizes[0] * ((int32_T)
             (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_c_k, 1,
              c9_Y_sizes[1], 2, 0) - 1)) - 1] - c9_t_y;
      }
    }
  }
}

static real_T c9_sqrt(SFc9_simulationInstanceStruct *chartInstance, real_T c9_x)
{
  real_T c9_b_x;
  c9_b_x = c9_x;
  c9_b_sqrt(chartInstance, &c9_b_x);
  return c9_b_x;
}

static real_T c9_eml_xnrm2(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0)
{
  real_T c9_y;
  int32_T c9_c_n;
  int32_T c9_b_ix0;
  int32_T c9_d_n;
  int32_T c9_c_ix0;
  int32_T c9_var;
  ptrdiff_t c9_n_t;
  ptrdiff_t c9_incx_t;
  double * c9_xix0_t;
  c9_c_n = c9_b_n;
  c9_b_ix0 = c9_ix0;
  if (c9_c_n < 1) {
    c9_y = 0.0;
  } else {
    c9_d_n = c9_c_n;
    c9_c_ix0 = c9_b_ix0;
    c9_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_var);
    c9_incx_t = (ptrdiff_t)(1);
    c9_xix0_t = (double *)(&c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ix0,
      1, c9_x_sizes[0] * c9_x_sizes[1], 1, 0) - 1]);
    c9_y = dnrm2(&c9_n_t, c9_xix0_t, &c9_incx_t);
  }

  return c9_y;
}

static int32_T c9_eml_ixamax(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_b_n, real_T c9_x_data[20], int32_T c9_x_sizes[1], int32_T c9_ix0)
{
  int32_T c9_idxmax;
  int32_T c9_c_n;
  int32_T c9_b_ix0;
  int32_T c9_d_n;
  int32_T c9_c_ix0;
  int32_T c9_var;
  ptrdiff_t c9_n_t;
  ptrdiff_t c9_incx_t;
  double * c9_xix0_t;
  ptrdiff_t c9_idxmax_t;
  c9_c_n = c9_b_n;
  c9_b_ix0 = c9_ix0;
  if (c9_c_n < 1) {
    c9_idxmax = 0;
  } else {
    c9_d_n = c9_c_n;
    c9_c_ix0 = c9_b_ix0;
    c9_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_var);
    c9_incx_t = (ptrdiff_t)(1);
    c9_xix0_t = (double *)(&c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ix0,
      1, c9_x_sizes[0], 1, 0) - 1]);
    c9_idxmax_t = idamax(&c9_n_t, c9_xix0_t, &c9_incx_t);
    c9_idxmax = (int32_T)(c9_idxmax_t);
  }

  return c9_idxmax;
}

static void c9_b_eml_xswap(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0, int32_T
  c9_iy0, real_T c9_b_x_data[400], int32_T c9_b_x_sizes[2])
{
  int32_T c9_x;
  int32_T c9_b_x;
  int32_T c9_loop_ub;
  int32_T c9_i87;
  c9_b_x_sizes[0] = c9_x_sizes[0];
  c9_b_x_sizes[1] = c9_x_sizes[1];
  c9_x = c9_b_x_sizes[0];
  c9_b_x = c9_b_x_sizes[1];
  c9_loop_ub = c9_x_sizes[0] * c9_x_sizes[1] - 1;
  for (c9_i87 = 0; c9_i87 <= c9_loop_ub; c9_i87++) {
    c9_b_x_data[c9_i87] = c9_x_data[c9_i87];
  }

  c9_d_eml_xswap(chartInstance, c9_b_n, c9_b_x_data, c9_b_x_sizes, c9_ix0,
                 c9_iy0);
}

static void c9_eml_matlab_zlarfg(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_b_n, real_T c9_alpha1, real_T c9_x_data[400], int32_T c9_x_sizes[2],
  int32_T c9_ix0, real_T *c9_b_alpha1, real_T c9_b_x_data[400], int32_T
  c9_b_x_sizes[2], real_T *c9_tau)
{
  int32_T c9_x;
  int32_T c9_b_x;
  int32_T c9_loop_ub;
  int32_T c9_i88;
  *c9_b_alpha1 = c9_alpha1;
  c9_b_x_sizes[0] = c9_x_sizes[0];
  c9_b_x_sizes[1] = c9_x_sizes[1];
  c9_x = c9_b_x_sizes[0];
  c9_b_x = c9_b_x_sizes[1];
  c9_loop_ub = c9_x_sizes[0] * c9_x_sizes[1] - 1;
  for (c9_i88 = 0; c9_i88 <= c9_loop_ub; c9_i88++) {
    c9_b_x_data[c9_i88] = c9_x_data[c9_i88];
  }

  *c9_tau = c9_c_eml_matlab_zlarfg(chartInstance, c9_b_n, c9_b_alpha1,
    c9_b_x_data, c9_b_x_sizes, c9_ix0);
}

static void c9_eml_xscal(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_a, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T
  c9_ix0, real_T c9_b_x_data[400], int32_T c9_b_x_sizes[2])
{
  int32_T c9_x;
  int32_T c9_b_x;
  int32_T c9_loop_ub;
  int32_T c9_i89;
  c9_b_x_sizes[0] = c9_x_sizes[0];
  c9_b_x_sizes[1] = c9_x_sizes[1];
  c9_x = c9_b_x_sizes[0];
  c9_b_x = c9_b_x_sizes[1];
  c9_loop_ub = c9_x_sizes[0] * c9_x_sizes[1] - 1;
  for (c9_i89 = 0; c9_i89 <= c9_loop_ub; c9_i89++) {
    c9_b_x_data[c9_i89] = c9_x_data[c9_i89];
  }

  c9_c_eml_xscal(chartInstance, c9_b_n, c9_a, c9_b_x_data, c9_b_x_sizes, c9_ix0);
}

static void c9_b_eml_matlab_zlarfg(SFc9_simulationInstanceStruct *chartInstance,
  real_T c9_alpha1, real_T c9_x, real_T *c9_b_alpha1, real_T *c9_b_x, real_T
  *c9_tau)
{
  *c9_b_alpha1 = c9_alpha1;
  *c9_b_x = c9_x;
  *c9_tau = c9_d_eml_matlab_zlarfg(chartInstance, c9_b_alpha1, c9_b_x);
}

static void c9_b_eml_xnrm2(SFc9_simulationInstanceStruct *chartInstance)
{
}

static real_T c9_b_eml_xscal(SFc9_simulationInstanceStruct *chartInstance,
  real_T c9_x)
{
  real_T c9_b_x;
  c9_b_x = c9_x;
  c9_d_eml_xscal(chartInstance, &c9_b_x);
  return c9_b_x;
}

static void c9_e_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_eml_matlab_zlarf(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_m, int32_T c9_b_n, int32_T c9_iv0, real_T c9_tau, real_T c9_C_data
  [400], int32_T c9_C_sizes[2], int32_T c9_ic0, int32_T c9_ldc, real_T
  c9_work_data[20], int32_T c9_work_sizes[1], real_T c9_b_C_data[400], int32_T
  c9_b_C_sizes[2], real_T c9_b_work_data[20], int32_T c9_b_work_sizes[1])
{
  int32_T c9_C;
  int32_T c9_b_C;
  int32_T c9_loop_ub;
  int32_T c9_i90;
  int32_T c9_b_loop_ub;
  int32_T c9_i91;
  c9_b_C_sizes[0] = c9_C_sizes[0];
  c9_b_C_sizes[1] = c9_C_sizes[1];
  c9_C = c9_b_C_sizes[0];
  c9_b_C = c9_b_C_sizes[1];
  c9_loop_ub = c9_C_sizes[0] * c9_C_sizes[1] - 1;
  for (c9_i90 = 0; c9_i90 <= c9_loop_ub; c9_i90++) {
    c9_b_C_data[c9_i90] = c9_C_data[c9_i90];
  }

  c9_b_work_sizes[0] = c9_work_sizes[0];
  c9_b_loop_ub = c9_work_sizes[0] - 1;
  for (c9_i91 = 0; c9_i91 <= c9_b_loop_ub; c9_i91++) {
    c9_b_work_data[c9_i91] = c9_work_data[c9_i91];
  }

  c9_b_eml_matlab_zlarf(chartInstance, c9_m, c9_b_n, c9_iv0, c9_tau, c9_b_C_data,
                        c9_b_C_sizes, c9_ic0, c9_ldc, c9_b_work_data,
                        c9_b_work_sizes);
}

static void c9_eml_xgemv(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_ia0, int32_T c9_lda, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T
  c9_ix0, real_T c9_y_data[20], int32_T c9_y_sizes[1], real_T c9_b_y_data[20],
  int32_T c9_b_y_sizes[1])
{
  int32_T c9_loop_ub;
  int32_T c9_i92;
  int32_T c9_b_A_sizes[2];
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_b_loop_ub;
  int32_T c9_i93;
  real_T c9_b_A_data[400];
  int32_T c9_b_x_sizes[2];
  int32_T c9_x;
  int32_T c9_b_x;
  int32_T c9_c_loop_ub;
  int32_T c9_i94;
  real_T c9_b_x_data[400];
  c9_b_y_sizes[0] = c9_y_sizes[0];
  c9_loop_ub = c9_y_sizes[0] - 1;
  for (c9_i92 = 0; c9_i92 <= c9_loop_ub; c9_i92++) {
    c9_b_y_data[c9_i92] = c9_y_data[c9_i92];
  }

  c9_b_A_sizes[0] = c9_A_sizes[0];
  c9_b_A_sizes[1] = c9_A_sizes[1];
  c9_A = c9_b_A_sizes[0];
  c9_b_A = c9_b_A_sizes[1];
  c9_b_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
  for (c9_i93 = 0; c9_i93 <= c9_b_loop_ub; c9_i93++) {
    c9_b_A_data[c9_i93] = c9_A_data[c9_i93];
  }

  c9_b_x_sizes[0] = c9_x_sizes[0];
  c9_b_x_sizes[1] = c9_x_sizes[1];
  c9_x = c9_b_x_sizes[0];
  c9_b_x = c9_b_x_sizes[1];
  c9_c_loop_ub = c9_x_sizes[0] * c9_x_sizes[1] - 1;
  for (c9_i94 = 0; c9_i94 <= c9_c_loop_ub; c9_i94++) {
    c9_b_x_data[c9_i94] = c9_x_data[c9_i94];
  }

  c9_b_eml_xgemv(chartInstance, c9_m, c9_b_n, c9_b_A_data, c9_b_A_sizes, c9_ia0,
                 c9_lda, c9_b_x_data, c9_b_x_sizes, c9_ix0, c9_b_y_data,
                 c9_b_y_sizes);
}

static void c9_b_eml_xger(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_alpha1, int32_T c9_ix0, real_T c9_y_data[20],
  int32_T c9_y_sizes[1], real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_ia0, int32_T c9_lda, real_T c9_b_A_data[400], int32_T c9_b_A_sizes[2])
{
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_loop_ub;
  int32_T c9_i95;
  int32_T c9_b_y_sizes;
  int32_T c9_b_loop_ub;
  int32_T c9_i96;
  real_T c9_b_y_data[20];
  c9_b_A_sizes[0] = c9_A_sizes[0];
  c9_b_A_sizes[1] = c9_A_sizes[1];
  c9_A = c9_b_A_sizes[0];
  c9_b_A = c9_b_A_sizes[1];
  c9_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
  for (c9_i95 = 0; c9_i95 <= c9_loop_ub; c9_i95++) {
    c9_b_A_data[c9_i95] = c9_A_data[c9_i95];
  }

  c9_b_y_sizes = c9_y_sizes[0];
  c9_b_loop_ub = c9_y_sizes[0] - 1;
  for (c9_i96 = 0; c9_i96 <= c9_b_loop_ub; c9_i96++) {
    c9_b_y_data[c9_i96] = c9_y_data[c9_i96];
  }

  c9_d_eml_xger(chartInstance, c9_m, c9_b_n, c9_alpha1, c9_ix0, c9_b_y_data,
                *(int32_T (*)[1])&c9_b_y_sizes, c9_b_A_data, c9_b_A_sizes,
                c9_ia0, c9_lda);
}

static real_T c9_c_eml_xnrm2(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0)
{
  real_T c9_y;
  int32_T c9_c_n;
  int32_T c9_b_ix0;
  int32_T c9_d_n;
  int32_T c9_c_ix0;
  int32_T c9_var;
  ptrdiff_t c9_n_t;
  ptrdiff_t c9_incx_t;
  double * c9_xix0_t;
  c9_c_n = c9_b_n;
  c9_b_ix0 = c9_ix0;
  if (c9_c_n < 1) {
    c9_y = 0.0;
  } else {
    c9_d_n = c9_c_n;
    c9_c_ix0 = c9_b_ix0;
    c9_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_var);
    c9_incx_t = (ptrdiff_t)(1);
    c9_xix0_t = (double *)(&c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ix0,
      1, c9_x_sizes[0] * c9_x_sizes[1], 1, 0) - 1]);
    c9_y = dnrm2(&c9_n_t, c9_xix0_t, &c9_incx_t);
  }

  return c9_y;
}

static void c9_eml_flt2str(SFc9_simulationInstanceStruct *chartInstance, real_T
  c9_x, char_T c9_str[14])
{
  int32_T c9_i97;
  static char_T c9_cv8[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c9_u[8];
  const mxArray *c9_y = NULL;
  real_T c9_b_u;
  const mxArray *c9_b_y = NULL;
  real_T c9_c_u;
  const mxArray *c9_c_y = NULL;
  real_T c9_d_u;
  const mxArray *c9_d_y = NULL;
  for (c9_i97 = 0; c9_i97 < 8; c9_i97++) {
    c9_u[c9_i97] = c9_cv8[c9_i97];
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 10, 0U, 1U, 0U, 2, 1, 8), FALSE);
  c9_b_u = 14.0;
  c9_b_y = NULL;
  sf_mex_assign(&c9_b_y, sf_mex_create("y", &c9_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  c9_c_u = 6.0;
  c9_c_y = NULL;
  sf_mex_assign(&c9_c_y, sf_mex_create("y", &c9_c_u, 0, 0U, 0U, 0U, 0), FALSE);
  c9_d_u = c9_x;
  c9_d_y = NULL;
  sf_mex_assign(&c9_d_y, sf_mex_create("y", &c9_d_u, 0, 0U, 0U, 0U, 0), FALSE);
  c9_e_emlrt_marshallIn(chartInstance, sf_mex_call_debug("sprintf", 1U, 2U, 14,
    sf_mex_call_debug("sprintf", 1U, 3U, 14, c9_y, 14, c9_b_y, 14, c9_c_y), 14,
    c9_d_y), "sprintf", c9_str);
}

static void c9_b_eml_warning(SFc9_simulationInstanceStruct *chartInstance,
  real_T c9_varargin_2, char_T c9_varargin_3[14])
{
  int32_T c9_i98;
  static char_T c9_varargin_1[32] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 'r', 'a', 'n', 'k', 'D', 'e', 'f', 'i', 'c', 'i',
    'e', 'n', 't', 'M', 'a', 't', 'r', 'i', 'x' };

  char_T c9_u[32];
  const mxArray *c9_y = NULL;
  real_T c9_b_u;
  const mxArray *c9_b_y = NULL;
  int32_T c9_i99;
  char_T c9_c_u[14];
  const mxArray *c9_c_y = NULL;
  for (c9_i98 = 0; c9_i98 < 32; c9_i98++) {
    c9_u[c9_i98] = c9_varargin_1[c9_i98];
  }

  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", c9_u, 10, 0U, 1U, 0U, 2, 1, 32), FALSE);
  c9_b_u = c9_varargin_2;
  c9_b_y = NULL;
  sf_mex_assign(&c9_b_y, sf_mex_create("y", &c9_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  for (c9_i99 = 0; c9_i99 < 14; c9_i99++) {
    c9_c_u[c9_i99] = c9_varargin_3[c9_i99];
  }

  c9_c_y = NULL;
  sf_mex_assign(&c9_c_y, sf_mex_create("y", c9_c_u, 10, 0U, 1U, 0U, 2, 1, 14),
                FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 3U,
    14, c9_y, 14, c9_b_y, 14, c9_c_y));
}

static void c9_f_eml_scalar_eg(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_eml_xgemm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_k, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B[4], int32_T c9_ldb, real_T c9_C_data[20], int32_T
  c9_C_sizes[1], int32_T c9_ldc, real_T c9_b_C_data[20], int32_T c9_b_C_sizes[1])
{
  int32_T c9_loop_ub;
  int32_T c9_i100;
  int32_T c9_b_A_sizes[2];
  int32_T c9_A;
  int32_T c9_b_A;
  int32_T c9_b_loop_ub;
  int32_T c9_i101;
  real_T c9_b_A_data[400];
  int32_T c9_i102;
  real_T c9_b_B[4];
  c9_b_C_sizes[0] = c9_C_sizes[0];
  c9_loop_ub = c9_C_sizes[0] - 1;
  for (c9_i100 = 0; c9_i100 <= c9_loop_ub; c9_i100++) {
    c9_b_C_data[c9_i100] = c9_C_data[c9_i100];
  }

  c9_b_A_sizes[0] = c9_A_sizes[0];
  c9_b_A_sizes[1] = c9_A_sizes[1];
  c9_A = c9_b_A_sizes[0];
  c9_b_A = c9_b_A_sizes[1];
  c9_b_loop_ub = c9_A_sizes[0] * c9_A_sizes[1] - 1;
  for (c9_i101 = 0; c9_i101 <= c9_b_loop_ub; c9_i101++) {
    c9_b_A_data[c9_i101] = c9_A_data[c9_i101];
  }

  for (c9_i102 = 0; c9_i102 < 4; c9_i102++) {
    c9_b_B[c9_i102] = c9_B[c9_i102];
  }

  c9_b_eml_xgemm(chartInstance, c9_m, c9_k, c9_b_A_data, c9_b_A_sizes, c9_lda,
                 c9_b_B, c9_ldb, c9_b_C_data, c9_b_C_sizes, c9_ldc);
}

static void c9_below_threshold(SFc9_simulationInstanceStruct *chartInstance)
{
}

static void c9_tanh(SFc9_simulationInstanceStruct *chartInstance, real_T
                    c9_x_data[20], int32_T c9_x_sizes[1], real_T c9_b_x_data[20],
                    int32_T c9_b_x_sizes[1])
{
  int32_T c9_loop_ub;
  int32_T c9_i103;
  c9_b_x_sizes[0] = c9_x_sizes[0];
  c9_loop_ub = c9_x_sizes[0] - 1;
  for (c9_i103 = 0; c9_i103 <= c9_loop_ub; c9_i103++) {
    c9_b_x_data[c9_i103] = c9_x_data[c9_i103];
  }

  c9_b_tanh(chartInstance, c9_b_x_data, c9_b_x_sizes);
}

static void c9_e_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_sprintf, const char_T *c9_identifier, char_T c9_y[14])
{
  emlrtMsgIdentifier c9_thisId;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_sprintf), &c9_thisId, c9_y);
  sf_mex_destroy(&c9_sprintf);
}

static void c9_f_emlrt_marshallIn(SFc9_simulationInstanceStruct *chartInstance,
  const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId, char_T c9_y[14])
{
  char_T c9_cv9[14];
  int32_T c9_i104;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), c9_cv9, 1, 10, 0U, 1, 0U, 2, 1,
                14);
  for (c9_i104 = 0; c9_i104 < 14; c9_i104++) {
    c9_y[c9_i104] = c9_cv9[c9_i104];
  }

  sf_mex_destroy(&c9_u);
}

static const mxArray *c9_e_sf_marshallOut(void *chartInstanceVoid, void
  *c9_inData)
{
  const mxArray *c9_mxArrayOutData = NULL;
  int32_T c9_u;
  const mxArray *c9_y = NULL;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_mxArrayOutData = NULL;
  c9_u = *(int32_T *)c9_inData;
  c9_y = NULL;
  sf_mex_assign(&c9_y, sf_mex_create("y", &c9_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c9_mxArrayOutData, c9_y, FALSE);
  return c9_mxArrayOutData;
}

static int32_T c9_g_emlrt_marshallIn(SFc9_simulationInstanceStruct
  *chartInstance, const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId)
{
  int32_T c9_y;
  int32_T c9_i105;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), &c9_i105, 1, 6, 0U, 0, 0U, 0);
  c9_y = c9_i105;
  sf_mex_destroy(&c9_u);
  return c9_y;
}

static void c9_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c9_mxArrayInData, const char_T *c9_varName, void *c9_outData)
{
  const mxArray *c9_b_sfEvent;
  const char_T *c9_identifier;
  emlrtMsgIdentifier c9_thisId;
  int32_T c9_y;
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)chartInstanceVoid;
  c9_b_sfEvent = sf_mex_dup(c9_mxArrayInData);
  c9_identifier = c9_varName;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_y = c9_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c9_b_sfEvent),
    &c9_thisId);
  sf_mex_destroy(&c9_b_sfEvent);
  *(int32_T *)c9_outData = c9_y;
  sf_mex_destroy(&c9_mxArrayInData);
}

static uint8_T c9_h_emlrt_marshallIn(SFc9_simulationInstanceStruct
  *chartInstance, const mxArray *c9_b_is_active_c9_simulation, const char_T
  *c9_identifier)
{
  uint8_T c9_y;
  emlrtMsgIdentifier c9_thisId;
  c9_thisId.fIdentifier = c9_identifier;
  c9_thisId.fParent = NULL;
  c9_y = c9_i_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c9_b_is_active_c9_simulation), &c9_thisId);
  sf_mex_destroy(&c9_b_is_active_c9_simulation);
  return c9_y;
}

static uint8_T c9_i_emlrt_marshallIn(SFc9_simulationInstanceStruct
  *chartInstance, const mxArray *c9_u, const emlrtMsgIdentifier *c9_parentId)
{
  uint8_T c9_y;
  uint8_T c9_u0;
  sf_mex_import(c9_parentId, sf_mex_dup(c9_u), &c9_u0, 1, 3, 0U, 0, 0U, 0);
  c9_y = c9_u0;
  sf_mex_destroy(&c9_u);
  return c9_y;
}

static void c9_c_eml_xswap(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0, int32_T
  c9_incx, int32_T c9_iy0, int32_T c9_incy)
{
  int32_T c9_c_n;
  int32_T c9_b_ix0;
  int32_T c9_b_incx;
  int32_T c9_b_iy0;
  int32_T c9_b_incy;
  int32_T c9_d_n;
  int32_T c9_c_ix0;
  int32_T c9_c_incx;
  int32_T c9_c_iy0;
  int32_T c9_c_incy;
  int32_T c9_ix;
  int32_T c9_iy;
  int32_T c9_x;
  int32_T c9_b_x;
  int32_T c9_ixinc;
  int32_T c9_c_x;
  int32_T c9_d_x;
  int32_T c9_iyinc;
  int32_T c9_e_n;
  int32_T c9_b;
  int32_T c9_b_b;
  boolean_T c9_overflow;
  int32_T c9_k;
  int32_T c9_e_x[1];
  real_T c9_temp;
  int32_T c9_f_x[1];
  int32_T c9_g_x[1];
  int32_T c9_h_x[1];
  int32_T c9_a;
  int32_T c9_c_b;
  int32_T c9_b_a;
  int32_T c9_d_b;
  c9_c_n = c9_b_n;
  c9_b_ix0 = c9_ix0;
  c9_b_incx = c9_incx;
  c9_b_iy0 = c9_iy0;
  c9_b_incy = c9_incy;
  c9_d_n = c9_c_n;
  c9_c_ix0 = c9_b_ix0;
  c9_c_incx = c9_b_incx;
  c9_c_iy0 = c9_b_iy0;
  c9_c_incy = c9_b_incy;
  c9_ix = c9_c_ix0;
  c9_iy = c9_c_iy0;
  c9_x = c9_c_incx;
  c9_b_x = c9_x;
  c9_ixinc = c9_b_x;
  c9_c_x = c9_c_incy;
  c9_d_x = c9_c_x;
  c9_iyinc = c9_d_x;
  c9_e_n = c9_d_n;
  c9_b = c9_e_n;
  c9_b_b = c9_b;
  if (1 > c9_b_b) {
    c9_overflow = FALSE;
  } else {
    c9_overflow = (c9_b_b > 2147483646);
  }

  if (c9_overflow) {
    c9_check_forloop_overflow_error(chartInstance, TRUE);
  }

  for (c9_k = 1; c9_k <= c9_e_n; c9_k++) {
    c9_e_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_temp = c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_ix, 1, c9_x_sizes[0] *
      c9_x_sizes[1], 1, 0) - 1];
    c9_f_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_g_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_ix, 1, c9_x_sizes[0] *
      c9_x_sizes[1], 1, 0) - 1] = c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      c9_iy, 1, c9_x_sizes[0] * c9_x_sizes[1], 1, 0) - 1];
    c9_h_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_iy, 1, c9_x_sizes[0] *
      c9_x_sizes[1], 1, 0) - 1] = c9_temp;
    c9_a = c9_ix;
    c9_c_b = c9_ixinc;
    c9_ix = c9_a + c9_c_b;
    c9_b_a = c9_iy;
    c9_d_b = c9_iyinc;
    c9_iy = c9_b_a + c9_d_b;
  }
}

static void c9_c_eml_xger(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_alpha1, int32_T c9_ix0, int32_T c9_iy0,
  int32_T c9_incy, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T c9_ia0,
  int32_T c9_lda)
{
  int32_T c9_b_m;
  int32_T c9_c_n;
  int32_T c9_b_ix0;
  int32_T c9_b_iy0;
  int32_T c9_b_incy;
  int32_T c9_b_ia0;
  int32_T c9_b_lda;
  int32_T c9_c_m;
  int32_T c9_d_n;
  real_T c9_b_alpha1;
  int32_T c9_c_ix0;
  int32_T c9_c_iy0;
  int32_T c9_c_incy;
  int32_T c9_c_ia0;
  int32_T c9_c_lda;
  int32_T c9_var;
  ptrdiff_t c9_m_t;
  int32_T c9_b_var;
  ptrdiff_t c9_n_t;
  ptrdiff_t c9_incx_t;
  int32_T c9_c_var;
  ptrdiff_t c9_incy_t;
  int32_T c9_d_var;
  ptrdiff_t c9_lda_t;
  double * c9_alpha1_t;
  int32_T c9_A[1];
  double * c9_Aia0_t;
  int32_T c9_b_A[1];
  double * c9_Aix0_t;
  int32_T c9_c_A[1];
  double * c9_Aiy0_t;
  c9_b_m = c9_m;
  c9_c_n = c9_b_n;
  c9_b_ix0 = c9_ix0;
  c9_b_iy0 = c9_iy0;
  c9_b_incy = c9_incy;
  c9_b_ia0 = c9_ia0;
  c9_b_lda = c9_lda;
  if (c9_b_m < 1) {
  } else if (c9_c_n < 1) {
  } else {
    c9_c_m = c9_b_m;
    c9_d_n = c9_c_n;
    c9_b_alpha1 = -1.0;
    c9_c_ix0 = c9_b_ix0;
    c9_c_iy0 = c9_b_iy0;
    c9_c_incy = c9_b_incy;
    c9_c_ia0 = c9_b_ia0;
    c9_c_lda = c9_b_lda;
    c9_var = c9_c_m;
    c9_m_t = (ptrdiff_t)(c9_var);
    c9_b_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_b_var);
    c9_incx_t = (ptrdiff_t)(1);
    c9_c_var = c9_c_incy;
    c9_incy_t = (ptrdiff_t)(c9_c_var);
    c9_d_var = c9_c_lda;
    c9_lda_t = (ptrdiff_t)(c9_d_var);
    c9_alpha1_t = (double *)(&c9_b_alpha1);
    c9_A[0] = c9_A_sizes[0] * c9_A_sizes[1];
    c9_Aia0_t = (double *)(&c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ia0,
      1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0) - 1]);
    c9_b_A[0] = c9_A_sizes[0] * c9_A_sizes[1];
    c9_Aix0_t = (double *)(&c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ix0,
      1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0) - 1]);
    c9_c_A[0] = c9_A_sizes[0] * c9_A_sizes[1];
    c9_Aiy0_t = (double *)(&c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_iy0,
      1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0) - 1]);
    dger(&c9_m_t, &c9_n_t, c9_alpha1_t, c9_Aix0_t, &c9_incx_t, c9_Aiy0_t,
         &c9_incy_t, c9_Aia0_t, &c9_lda_t);
  }
}

static void c9_c_eml_xtrsm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B_data[400], int32_T c9_B_sizes[2], int32_T c9_ldb)
{
  int32_T c9_b_m;
  int32_T c9_c_n;
  int32_T c9_b_lda;
  int32_T c9_b_ldb;
  int32_T c9_c_m;
  int32_T c9_d_n;
  real_T c9_alpha1;
  int32_T c9_c_lda;
  int32_T c9_c_ldb;
  char_T c9_DIAGA;
  char_T c9_TRANSA;
  char_T c9_UPLO;
  char_T c9_SIDE;
  int32_T c9_var;
  ptrdiff_t c9_m_t;
  int32_T c9_b_var;
  ptrdiff_t c9_n_t;
  int32_T c9_c_var;
  ptrdiff_t c9_lda_t;
  int32_T c9_d_var;
  ptrdiff_t c9_ldb_t;
  double * c9_Aia0_t;
  int32_T c9_B[1];
  double * c9_Bib0_t;
  double * c9_alpha1_t;
  c9_b_m = c9_m;
  c9_c_n = c9_b_n;
  c9_b_lda = c9_lda;
  c9_b_ldb = c9_ldb;
  if (c9_b_m < 1) {
  } else if (c9_c_n < 1) {
  } else {
    c9_c_m = c9_b_m;
    c9_d_n = c9_c_n;
    c9_alpha1 = 1.0;
    c9_c_lda = c9_b_lda;
    c9_c_ldb = c9_b_ldb;
    c9_DIAGA = 'U';
    c9_TRANSA = 'N';
    c9_UPLO = 'L';
    c9_SIDE = 'L';
    c9_var = c9_c_m;
    c9_m_t = (ptrdiff_t)(c9_var);
    c9_b_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_b_var);
    c9_c_var = c9_c_lda;
    c9_lda_t = (ptrdiff_t)(c9_c_var);
    c9_d_var = c9_c_ldb;
    c9_ldb_t = (ptrdiff_t)(c9_d_var);
    _SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0);
    c9_Aia0_t = (double *)(&c9_A_data[0]);
    _SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_B_sizes[0] * c9_B_sizes[1], 1, 0);
    c9_B[0] = c9_B_sizes[0] * c9_B_sizes[1];
    c9_Bib0_t = (double *)(&c9_B_data[0]);
    c9_alpha1_t = (double *)(&c9_alpha1);
    dtrsm(&c9_SIDE, &c9_UPLO, &c9_TRANSA, &c9_DIAGA, &c9_m_t, &c9_n_t,
          c9_alpha1_t, c9_Aia0_t, &c9_lda_t, c9_Bib0_t, &c9_ldb_t);
  }
}

static void c9_d_eml_xtrsm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B_data[400], int32_T c9_B_sizes[2], int32_T c9_ldb)
{
  int32_T c9_b_m;
  int32_T c9_c_n;
  int32_T c9_b_lda;
  int32_T c9_b_ldb;
  int32_T c9_c_m;
  int32_T c9_d_n;
  real_T c9_alpha1;
  int32_T c9_c_lda;
  int32_T c9_c_ldb;
  char_T c9_DIAGA;
  char_T c9_TRANSA;
  char_T c9_UPLO;
  char_T c9_SIDE;
  int32_T c9_var;
  ptrdiff_t c9_m_t;
  int32_T c9_b_var;
  ptrdiff_t c9_n_t;
  int32_T c9_c_var;
  ptrdiff_t c9_lda_t;
  int32_T c9_d_var;
  ptrdiff_t c9_ldb_t;
  double * c9_Aia0_t;
  int32_T c9_B[1];
  double * c9_Bib0_t;
  double * c9_alpha1_t;
  c9_b_m = c9_m;
  c9_c_n = c9_b_n;
  c9_b_lda = c9_lda;
  c9_b_ldb = c9_ldb;
  if (c9_b_m < 1) {
  } else if (c9_c_n < 1) {
  } else {
    c9_c_m = c9_b_m;
    c9_d_n = c9_c_n;
    c9_alpha1 = 1.0;
    c9_c_lda = c9_b_lda;
    c9_c_ldb = c9_b_ldb;
    c9_DIAGA = 'N';
    c9_TRANSA = 'N';
    c9_UPLO = 'U';
    c9_SIDE = 'L';
    c9_var = c9_c_m;
    c9_m_t = (ptrdiff_t)(c9_var);
    c9_b_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_b_var);
    c9_c_var = c9_c_lda;
    c9_lda_t = (ptrdiff_t)(c9_c_var);
    c9_d_var = c9_c_ldb;
    c9_ldb_t = (ptrdiff_t)(c9_d_var);
    _SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0);
    c9_Aia0_t = (double *)(&c9_A_data[0]);
    _SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_B_sizes[0] * c9_B_sizes[1], 1, 0);
    c9_B[0] = c9_B_sizes[0] * c9_B_sizes[1];
    c9_Bib0_t = (double *)(&c9_B_data[0]);
    c9_alpha1_t = (double *)(&c9_alpha1);
    dtrsm(&c9_SIDE, &c9_UPLO, &c9_TRANSA, &c9_DIAGA, &c9_m_t, &c9_n_t,
          c9_alpha1_t, c9_Aia0_t, &c9_lda_t, c9_Bib0_t, &c9_ldb_t);
  }
}

static void c9_b_sqrt(SFc9_simulationInstanceStruct *chartInstance, real_T *c9_x)
{
  *c9_x = muDoubleScalarSqrt(*c9_x);
}

static void c9_d_eml_xswap(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T c9_ix0, int32_T
  c9_iy0)
{
  int32_T c9_c_n;
  int32_T c9_b_ix0;
  int32_T c9_b_iy0;
  int32_T c9_d_n;
  int32_T c9_c_ix0;
  int32_T c9_c_iy0;
  int32_T c9_ix;
  int32_T c9_iy;
  int32_T c9_e_n;
  int32_T c9_b;
  int32_T c9_b_b;
  boolean_T c9_overflow;
  int32_T c9_k;
  int32_T c9_x[1];
  real_T c9_temp;
  int32_T c9_b_x[1];
  int32_T c9_c_x[1];
  int32_T c9_d_x[1];
  int32_T c9_a;
  int32_T c9_b_a;
  c9_c_n = c9_b_n;
  c9_b_ix0 = c9_ix0;
  c9_b_iy0 = c9_iy0;
  c9_d_n = c9_c_n;
  c9_c_ix0 = c9_b_ix0;
  c9_c_iy0 = c9_b_iy0;
  c9_ix = c9_c_ix0;
  c9_iy = c9_c_iy0;
  c9_e_n = c9_d_n;
  c9_b = c9_e_n;
  c9_b_b = c9_b;
  if (1 > c9_b_b) {
    c9_overflow = FALSE;
  } else {
    c9_overflow = (c9_b_b > 2147483646);
  }

  if (c9_overflow) {
    c9_check_forloop_overflow_error(chartInstance, TRUE);
  }

  for (c9_k = 1; c9_k <= c9_e_n; c9_k++) {
    c9_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_temp = c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_ix, 1, c9_x_sizes[0] *
      c9_x_sizes[1], 1, 0) - 1];
    c9_b_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_c_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_ix, 1, c9_x_sizes[0] *
      c9_x_sizes[1], 1, 0) - 1] = c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      c9_iy, 1, c9_x_sizes[0] * c9_x_sizes[1], 1, 0) - 1];
    c9_d_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_iy, 1, c9_x_sizes[0] *
      c9_x_sizes[1], 1, 0) - 1] = c9_temp;
    c9_a = c9_ix + 1;
    c9_ix = c9_a;
    c9_b_a = c9_iy + 1;
    c9_iy = c9_b_a;
  }
}

static real_T c9_c_eml_matlab_zlarfg(SFc9_simulationInstanceStruct
  *chartInstance, int32_T c9_b_n, real_T *c9_alpha1, real_T c9_x_data[400],
  int32_T c9_x_sizes[2], int32_T c9_ix0)
{
  real_T c9_tau;
  int32_T c9_nm1;
  int32_T c9_b_x_sizes[2];
  int32_T c9_loop_ub;
  int32_T c9_i106;
  int32_T c9_b_loop_ub;
  int32_T c9_i107;
  real_T c9_b_x_data[400];
  real_T c9_xnorm;
  real_T c9_x1;
  real_T c9_x2;
  real_T c9_a;
  real_T c9_b;
  real_T c9_beta1;
  real_T c9_x;
  real_T c9_b_x;
  real_T c9_y;
  int32_T c9_knt;
  int32_T c9_b_a;
  real_T c9_d11;
  real_T c9_c_a;
  real_T c9_d_a;
  real_T c9_c_x;
  real_T c9_d_x;
  real_T c9_b_y;
  int32_T c9_c_x_sizes[2];
  int32_T c9_c_loop_ub;
  int32_T c9_i108;
  int32_T c9_d_loop_ub;
  int32_T c9_i109;
  real_T c9_c_x_data[400];
  real_T c9_b_x1;
  real_T c9_b_x2;
  real_T c9_e_a;
  real_T c9_b_b;
  real_T c9_e_x;
  real_T c9_c_y;
  real_T c9_d_y;
  int32_T c9_b_knt;
  int32_T c9_c_b;
  int32_T c9_d_b;
  boolean_T c9_overflow;
  int32_T c9_k;
  real_T c9_f_a;
  real_T c9_f_x;
  real_T c9_e_y;
  real_T c9_f_y;
  c9_tau = 0.0;
  if (c9_b_n <= 0) {
  } else {
    c9_nm1 = c9_b_n - 1;
    c9_b_x_sizes[0] = c9_x_sizes[0];
    c9_b_x_sizes[1] = c9_x_sizes[1];
    c9_loop_ub = c9_x_sizes[1] - 1;
    for (c9_i106 = 0; c9_i106 <= c9_loop_ub; c9_i106++) {
      c9_b_loop_ub = c9_x_sizes[0] - 1;
      for (c9_i107 = 0; c9_i107 <= c9_b_loop_ub; c9_i107++) {
        c9_b_x_data[c9_i107 + c9_b_x_sizes[0] * c9_i106] = c9_x_data[c9_i107 +
          c9_x_sizes[0] * c9_i106];
      }
    }

    c9_xnorm = c9_eml_xnrm2(chartInstance, c9_nm1, c9_b_x_data, c9_b_x_sizes,
      c9_ix0);
    if (c9_xnorm != 0.0) {
      c9_x1 = *c9_alpha1;
      c9_x2 = c9_xnorm;
      c9_a = c9_x1;
      c9_b = c9_x2;
      c9_beta1 = muDoubleScalarHypot(c9_a, c9_b);
      if (*c9_alpha1 >= 0.0) {
        c9_beta1 = -c9_beta1;
      }

      c9_realmin(chartInstance);
      c9_eps(chartInstance);
      c9_x = c9_beta1;
      c9_b_x = c9_x;
      c9_y = muDoubleScalarAbs(c9_b_x);
      if (c9_y < 1.0020841800044864E-292) {
        c9_knt = 0;
        do {
          c9_b_a = c9_knt + 1;
          c9_knt = c9_b_a;
          c9_d11 = 9.9792015476736E+291;
          c9_c_eml_xscal(chartInstance, c9_nm1, c9_d11, c9_x_data, c9_x_sizes,
                         c9_ix0);
          c9_c_a = c9_beta1;
          c9_beta1 = c9_c_a * 9.9792015476736E+291;
          c9_d_a = *c9_alpha1;
          *c9_alpha1 = c9_d_a * 9.9792015476736E+291;
          c9_c_x = c9_beta1;
          c9_d_x = c9_c_x;
          c9_b_y = muDoubleScalarAbs(c9_d_x);
        } while (!(c9_b_y >= 1.0020841800044864E-292));

        c9_c_x_sizes[0] = c9_x_sizes[0];
        c9_c_x_sizes[1] = c9_x_sizes[1];
        c9_c_loop_ub = c9_x_sizes[1] - 1;
        for (c9_i108 = 0; c9_i108 <= c9_c_loop_ub; c9_i108++) {
          c9_d_loop_ub = c9_x_sizes[0] - 1;
          for (c9_i109 = 0; c9_i109 <= c9_d_loop_ub; c9_i109++) {
            c9_c_x_data[c9_i109 + c9_c_x_sizes[0] * c9_i108] = c9_x_data[c9_i109
              + c9_x_sizes[0] * c9_i108];
          }
        }

        c9_xnorm = c9_eml_xnrm2(chartInstance, c9_nm1, c9_c_x_data, c9_c_x_sizes,
          c9_ix0);
        c9_b_x1 = *c9_alpha1;
        c9_b_x2 = c9_xnorm;
        c9_e_a = c9_b_x1;
        c9_b_b = c9_b_x2;
        c9_beta1 = muDoubleScalarHypot(c9_e_a, c9_b_b);
        if (*c9_alpha1 >= 0.0) {
          c9_beta1 = -c9_beta1;
        }

        c9_e_x = c9_beta1 - *c9_alpha1;
        c9_c_y = c9_beta1;
        c9_tau = c9_e_x / c9_c_y;
        c9_d_y = *c9_alpha1 - c9_beta1;
        *c9_alpha1 = 1.0 / c9_d_y;
        c9_c_eml_xscal(chartInstance, c9_nm1, *c9_alpha1, c9_x_data, c9_x_sizes,
                       c9_ix0);
        c9_b_knt = c9_knt;
        c9_c_b = c9_b_knt;
        c9_d_b = c9_c_b;
        if (1 > c9_d_b) {
          c9_overflow = FALSE;
        } else {
          c9_overflow = (c9_d_b > 2147483646);
        }

        if (c9_overflow) {
          c9_check_forloop_overflow_error(chartInstance, TRUE);
        }

        for (c9_k = 1; c9_k <= c9_b_knt; c9_k++) {
          c9_f_a = c9_beta1;
          c9_beta1 = c9_f_a * 1.0020841800044864E-292;
        }

        *c9_alpha1 = c9_beta1;
      } else {
        c9_f_x = c9_beta1 - *c9_alpha1;
        c9_e_y = c9_beta1;
        c9_tau = c9_f_x / c9_e_y;
        c9_f_y = *c9_alpha1 - c9_beta1;
        *c9_alpha1 = 1.0 / c9_f_y;
        c9_c_eml_xscal(chartInstance, c9_nm1, *c9_alpha1, c9_x_data, c9_x_sizes,
                       c9_ix0);
        *c9_alpha1 = c9_beta1;
      }
    }
  }

  return c9_tau;
}

static void c9_c_eml_xscal(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_b_n, real_T c9_a, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T
  c9_ix0)
{
  int32_T c9_c_n;
  real_T c9_b_a;
  int32_T c9_b_ix0;
  int32_T c9_d_n;
  real_T c9_c_a;
  int32_T c9_c_ix0;
  int32_T c9_var;
  ptrdiff_t c9_n_t;
  ptrdiff_t c9_incx_t;
  int32_T c9_x[1];
  double * c9_xix0_t;
  double * c9_a_t;
  c9_c_n = c9_b_n;
  c9_b_a = c9_a;
  c9_b_ix0 = c9_ix0;
  if (c9_c_n < 1) {
  } else {
    c9_d_n = c9_c_n;
    c9_c_a = c9_b_a;
    c9_c_ix0 = c9_b_ix0;
    c9_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_var);
    c9_incx_t = (ptrdiff_t)(1);
    c9_x[0] = c9_x_sizes[0] * c9_x_sizes[1];
    c9_xix0_t = (double *)(&c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ix0,
      1, c9_x_sizes[0] * c9_x_sizes[1], 1, 0) - 1]);
    c9_a_t = (double *)(&c9_c_a);
    dscal(&c9_n_t, c9_a_t, c9_xix0_t, &c9_incx_t);
  }
}

static real_T c9_d_eml_matlab_zlarfg(SFc9_simulationInstanceStruct
  *chartInstance, real_T *c9_alpha1, real_T *c9_x)
{
  real_T c9_tau;
  c9_tau = 0.0;
  c9_b_eml_xnrm2(chartInstance);
  return c9_tau;
}

static void c9_d_eml_xscal(SFc9_simulationInstanceStruct *chartInstance, real_T *
  c9_x)
{
}

static void c9_b_eml_matlab_zlarf(SFc9_simulationInstanceStruct *chartInstance,
  int32_T c9_m, int32_T c9_b_n, int32_T c9_iv0, real_T c9_tau, real_T c9_C_data
  [400], int32_T c9_C_sizes[2], int32_T c9_ic0, int32_T c9_ldc, real_T
  c9_work_data[20], int32_T c9_work_sizes[1])
{
  int32_T c9_lastv;
  int32_T c9_a;
  int32_T c9_c;
  int32_T c9_b;
  int32_T c9_b_c;
  int32_T c9_b_a;
  int32_T c9_b_b;
  int32_T c9_b_i;
  int32_T c9_C[1];
  int32_T c9_c_a;
  int32_T c9_d_a;
  int32_T c9_b_m;
  int32_T c9_c_n;
  int32_T c9_ia0;
  int32_T c9_lda;
  int32_T c9_lastc;
  int32_T c9_e_a;
  int32_T c9_c_c;
  int32_T c9_f_a;
  int32_T c9_c_b;
  int32_T c9_d_c;
  int32_T c9_g_a;
  int32_T c9_d_b;
  int32_T c9_coltop;
  int32_T c9_h_a;
  int32_T c9_e_c;
  int32_T c9_i_a;
  int32_T c9_e_b;
  int32_T c9_colbottom;
  int32_T c9_b_coltop;
  int32_T c9_b_colbottom;
  int32_T c9_j_a;
  int32_T c9_f_b;
  int32_T c9_k_a;
  int32_T c9_g_b;
  boolean_T c9_overflow;
  int32_T c9_ia;
  int32_T c9_b_ia;
  int32_T c9_b_C[1];
  int32_T c9_l_a;
  int32_T c9_b_C_sizes[2];
  int32_T c9_loop_ub;
  int32_T c9_i110;
  int32_T c9_b_loop_ub;
  int32_T c9_i111;
  real_T c9_b_C_data[400];
  int32_T c9_c_C_sizes[2];
  int32_T c9_c_loop_ub;
  int32_T c9_i112;
  int32_T c9_d_loop_ub;
  int32_T c9_i113;
  real_T c9_c_C_data[400];
  int32_T c9_c_m;
  int32_T c9_d_n;
  real_T c9_alpha1;
  int32_T c9_ix0;
  int32_T c9_b_ia0;
  int32_T c9_b_lda;
  int32_T c9_b_work_sizes;
  int32_T c9_e_loop_ub;
  int32_T c9_i114;
  real_T c9_b_work_data[20];
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  if (c9_tau != 0.0) {
    c9_lastv = c9_m;
    c9_a = c9_lastv;
    c9_c = c9_a;
    c9_b = c9_c - 1;
    c9_b_c = c9_b;
    c9_b_a = c9_iv0;
    c9_b_b = c9_b_c;
    c9_b_i = c9_b_a + c9_b_b;
    exitg3 = FALSE;
    while ((exitg3 == FALSE) && (c9_lastv > 0)) {
      c9_C[0] = c9_C_sizes[0] * c9_C_sizes[1];
      if (c9_C_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_i, 1, c9_C_sizes[0] *
           c9_C_sizes[1], 1, 0) - 1] == 0.0) {
        c9_c_a = c9_lastv - 1;
        c9_lastv = c9_c_a;
        c9_d_a = c9_b_i - 1;
        c9_b_i = c9_d_a;
      } else {
        exitg3 = TRUE;
      }
    }

    c9_b_m = c9_lastv;
    c9_c_n = c9_b_n;
    c9_ia0 = c9_ic0;
    c9_lda = c9_ldc;
    c9_lastc = c9_c_n;
    exitg2 = FALSE;
    while ((exitg2 == FALSE) && (c9_lastc > 0)) {
      c9_e_a = c9_lastc;
      c9_c_c = c9_e_a;
      c9_f_a = c9_c_c - 1;
      c9_c_b = c9_lda;
      c9_d_c = c9_f_a * c9_c_b;
      c9_g_a = c9_ia0;
      c9_d_b = c9_d_c;
      c9_coltop = c9_g_a + c9_d_b;
      c9_h_a = c9_b_m;
      c9_e_c = c9_h_a;
      c9_i_a = c9_coltop;
      c9_e_b = c9_e_c - 1;
      c9_colbottom = c9_i_a + c9_e_b;
      c9_b_coltop = c9_coltop;
      c9_b_colbottom = c9_colbottom;
      c9_j_a = c9_b_coltop;
      c9_f_b = c9_b_colbottom;
      c9_k_a = c9_j_a;
      c9_g_b = c9_f_b;
      if (c9_k_a > c9_g_b) {
        c9_overflow = FALSE;
      } else {
        c9_overflow = (c9_g_b > 2147483646);
      }

      if (c9_overflow) {
        c9_check_forloop_overflow_error(chartInstance, TRUE);
      }

      c9_ia = c9_b_coltop;
      do {
        exitg1 = 0;
        if (c9_ia <= c9_b_colbottom) {
          c9_b_ia = c9_ia;
          c9_b_C[0] = c9_C_sizes[0] * c9_C_sizes[1];
          if (c9_C_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_ia, 1, c9_C_sizes[0]
               * c9_C_sizes[1], 1, 0) - 1] != 0.0) {
            exitg1 = 1;
          } else {
            c9_ia++;
          }
        } else {
          c9_l_a = c9_lastc - 1;
          c9_lastc = c9_l_a;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = TRUE;
      }
    }
  } else {
    c9_lastv = 0;
    c9_lastc = 0;
  }

  if (c9_lastv > 0) {
    c9_b_C_sizes[0] = c9_C_sizes[0];
    c9_b_C_sizes[1] = c9_C_sizes[1];
    c9_loop_ub = c9_C_sizes[1] - 1;
    for (c9_i110 = 0; c9_i110 <= c9_loop_ub; c9_i110++) {
      c9_b_loop_ub = c9_C_sizes[0] - 1;
      for (c9_i111 = 0; c9_i111 <= c9_b_loop_ub; c9_i111++) {
        c9_b_C_data[c9_i111 + c9_b_C_sizes[0] * c9_i110] = c9_C_data[c9_i111 +
          c9_C_sizes[0] * c9_i110];
      }
    }

    c9_c_C_sizes[0] = c9_C_sizes[0];
    c9_c_C_sizes[1] = c9_C_sizes[1];
    c9_c_loop_ub = c9_C_sizes[1] - 1;
    for (c9_i112 = 0; c9_i112 <= c9_c_loop_ub; c9_i112++) {
      c9_d_loop_ub = c9_C_sizes[0] - 1;
      for (c9_i113 = 0; c9_i113 <= c9_d_loop_ub; c9_i113++) {
        c9_c_C_data[c9_i113 + c9_c_C_sizes[0] * c9_i112] = c9_C_data[c9_i113 +
          c9_C_sizes[0] * c9_i112];
      }
    }

    c9_b_eml_xgemv(chartInstance, c9_lastv, c9_lastc, c9_b_C_data, c9_b_C_sizes,
                   c9_ic0, c9_ldc, c9_c_C_data, c9_c_C_sizes, c9_iv0,
                   c9_work_data, c9_work_sizes);
    c9_c_m = c9_lastv;
    c9_d_n = c9_lastc;
    c9_alpha1 = -c9_tau;
    c9_ix0 = c9_iv0;
    c9_b_ia0 = c9_ic0;
    c9_b_lda = c9_ldc;
    c9_b_work_sizes = c9_work_sizes[0];
    c9_e_loop_ub = c9_work_sizes[0] - 1;
    for (c9_i114 = 0; c9_i114 <= c9_e_loop_ub; c9_i114++) {
      c9_b_work_data[c9_i114] = c9_work_data[c9_i114];
    }

    c9_d_eml_xger(chartInstance, c9_c_m, c9_d_n, c9_alpha1, c9_ix0,
                  c9_b_work_data, *(int32_T (*)[1])&c9_b_work_sizes, c9_C_data,
                  c9_C_sizes, c9_b_ia0, c9_b_lda);
  }
}

static void c9_b_eml_xgemv(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_ia0, int32_T c9_lda, real_T c9_x_data[400], int32_T c9_x_sizes[2], int32_T
  c9_ix0, real_T c9_y_data[20], int32_T c9_y_sizes[1])
{
  int32_T c9_b_m;
  int32_T c9_c_n;
  int32_T c9_b_ia0;
  int32_T c9_b_lda;
  int32_T c9_b_ix0;
  int32_T c9_c_m;
  int32_T c9_d_n;
  real_T c9_alpha1;
  int32_T c9_c_ia0;
  int32_T c9_c_lda;
  int32_T c9_c_ix0;
  real_T c9_beta1;
  char_T c9_TRANSA;
  int32_T c9_var;
  ptrdiff_t c9_m_t;
  int32_T c9_b_var;
  ptrdiff_t c9_n_t;
  int32_T c9_c_var;
  ptrdiff_t c9_lda_t;
  ptrdiff_t c9_incx_t;
  ptrdiff_t c9_incy_t;
  double * c9_alpha1_t;
  double * c9_beta1_t;
  double * c9_yiy0_t;
  double * c9_Aia0_t;
  double * c9_xix0_t;
  c9_b_m = c9_m;
  c9_c_n = c9_b_n;
  c9_b_ia0 = c9_ia0;
  c9_b_lda = c9_lda;
  c9_b_ix0 = c9_ix0;
  if (c9_c_n < 1) {
  } else {
    c9_c_m = c9_b_m;
    c9_d_n = c9_c_n;
    c9_alpha1 = 1.0;
    c9_c_ia0 = c9_b_ia0;
    c9_c_lda = c9_b_lda;
    c9_c_ix0 = c9_b_ix0;
    c9_beta1 = 0.0;
    c9_TRANSA = 'C';
    c9_var = c9_c_m;
    c9_m_t = (ptrdiff_t)(c9_var);
    c9_b_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_b_var);
    c9_c_var = c9_c_lda;
    c9_lda_t = (ptrdiff_t)(c9_c_var);
    c9_incx_t = (ptrdiff_t)(1);
    c9_incy_t = (ptrdiff_t)(1);
    c9_alpha1_t = (double *)(&c9_alpha1);
    c9_beta1_t = (double *)(&c9_beta1);
    _SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_y_sizes[0], 1, 0);
    c9_yiy0_t = (double *)(&c9_y_data[0]);
    c9_Aia0_t = (double *)(&c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ia0,
      1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0) - 1]);
    c9_xix0_t = (double *)(&c9_x_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ix0,
      1, c9_x_sizes[0] * c9_x_sizes[1], 1, 0) - 1]);
    dgemv(&c9_TRANSA, &c9_m_t, &c9_n_t, c9_alpha1_t, c9_Aia0_t, &c9_lda_t,
          c9_xix0_t, &c9_incx_t, c9_beta1_t, c9_yiy0_t, &c9_incy_t);
  }
}

static void c9_d_eml_xger(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_b_n, real_T c9_alpha1, int32_T c9_ix0, real_T c9_y_data[20],
  int32_T c9_y_sizes[1], real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_ia0, int32_T c9_lda)
{
  int32_T c9_b_m;
  int32_T c9_c_n;
  real_T c9_b_alpha1;
  int32_T c9_b_ix0;
  int32_T c9_b_ia0;
  int32_T c9_b_lda;
  int32_T c9_c_m;
  int32_T c9_d_n;
  real_T c9_c_alpha1;
  int32_T c9_c_ix0;
  int32_T c9_c_ia0;
  int32_T c9_c_lda;
  int32_T c9_var;
  ptrdiff_t c9_m_t;
  int32_T c9_b_var;
  ptrdiff_t c9_n_t;
  ptrdiff_t c9_incx_t;
  ptrdiff_t c9_incy_t;
  int32_T c9_c_var;
  ptrdiff_t c9_lda_t;
  double * c9_alpha1_t;
  int32_T c9_A[1];
  double * c9_Aia0_t;
  int32_T c9_b_A[1];
  double * c9_Aix0_t;
  double * c9_yiy0_t;
  c9_b_m = c9_m;
  c9_c_n = c9_b_n;
  c9_b_alpha1 = c9_alpha1;
  c9_b_ix0 = c9_ix0;
  c9_b_ia0 = c9_ia0;
  c9_b_lda = c9_lda;
  if (c9_c_n < 1) {
  } else {
    c9_c_m = c9_b_m;
    c9_d_n = c9_c_n;
    c9_c_alpha1 = c9_b_alpha1;
    c9_c_ix0 = c9_b_ix0;
    c9_c_ia0 = c9_b_ia0;
    c9_c_lda = c9_b_lda;
    c9_var = c9_c_m;
    c9_m_t = (ptrdiff_t)(c9_var);
    c9_b_var = c9_d_n;
    c9_n_t = (ptrdiff_t)(c9_b_var);
    c9_incx_t = (ptrdiff_t)(1);
    c9_incy_t = (ptrdiff_t)(1);
    c9_c_var = c9_c_lda;
    c9_lda_t = (ptrdiff_t)(c9_c_var);
    c9_alpha1_t = (double *)(&c9_c_alpha1);
    c9_A[0] = c9_A_sizes[0] * c9_A_sizes[1];
    c9_Aia0_t = (double *)(&c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ia0,
      1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0) - 1]);
    c9_b_A[0] = c9_A_sizes[0] * c9_A_sizes[1];
    c9_Aix0_t = (double *)(&c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_c_ix0,
      1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0) - 1]);
    _SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_y_sizes[0], 1, 0);
    c9_yiy0_t = (double *)(&c9_y_data[0]);
    dger(&c9_m_t, &c9_n_t, c9_alpha1_t, c9_Aix0_t, &c9_incx_t, c9_yiy0_t,
         &c9_incy_t, c9_Aia0_t, &c9_lda_t);
  }
}

static void c9_b_eml_xgemm(SFc9_simulationInstanceStruct *chartInstance, int32_T
  c9_m, int32_T c9_k, real_T c9_A_data[400], int32_T c9_A_sizes[2], int32_T
  c9_lda, real_T c9_B[4], int32_T c9_ldb, real_T c9_C_data[20], int32_T
  c9_C_sizes[1], int32_T c9_ldc)
{
  int32_T c9_b_m;
  int32_T c9_b_k;
  int32_T c9_b_lda;
  int32_T c9_b_ldb;
  int32_T c9_b_ldc;
  int32_T c9_c_m;
  int32_T c9_c_k;
  int32_T c9_c_lda;
  int32_T c9_c_ldb;
  int32_T c9_c_ldc;
  int32_T c9_lda1;
  int32_T c9_ldb1;
  int32_T c9_ldc1;
  int32_T c9_c;
  int32_T c9_b;
  int32_T c9_lastColC;
  int32_T c9_b_ldc1;
  int32_T c9_b_lastColC;
  int32_T c9_d;
  int32_T c9_b_b;
  int32_T c9_b_d;
  int32_T c9_c_b;
  boolean_T c9_overflow;
  int32_T c9_cr;
  int32_T c9_b_cr;
  int32_T c9_a;
  int32_T c9_i115;
  int32_T c9_b_a;
  int32_T c9_d_b;
  int32_T c9_i116;
  int32_T c9_c_a;
  int32_T c9_e_b;
  int32_T c9_d_a;
  int32_T c9_f_b;
  boolean_T c9_b_overflow;
  int32_T c9_ic;
  int32_T c9_b_ic;
  int32_T c9_br;
  int32_T c9_c_ldc1;
  int32_T c9_c_lastColC;
  int32_T c9_c_d;
  int32_T c9_g_b;
  int32_T c9_d_d;
  int32_T c9_h_b;
  boolean_T c9_c_overflow;
  int32_T c9_c_cr;
  int32_T c9_ar;
  int32_T c9_e_a;
  int32_T c9_i117;
  int32_T c9_f_a;
  int32_T c9_i_b;
  int32_T c9_i118;
  int32_T c9_g_a;
  int32_T c9_j_b;
  int32_T c9_h_a;
  int32_T c9_k_b;
  boolean_T c9_d_overflow;
  int32_T c9_ib;
  int32_T c9_b_ib;
  real_T c9_temp;
  int32_T c9_ia;
  int32_T c9_i_a;
  int32_T c9_i119;
  int32_T c9_j_a;
  int32_T c9_l_b;
  int32_T c9_i120;
  int32_T c9_k_a;
  int32_T c9_m_b;
  int32_T c9_l_a;
  int32_T c9_n_b;
  boolean_T c9_e_overflow;
  int32_T c9_c_ic;
  int32_T c9_m_a;
  int32_T c9_n_a;
  int32_T c9_o_b;
  int32_T c9_o_a;
  int32_T c9_p_b;
  int32_T c9_d_m;
  int32_T c9_d_k;
  real_T c9_alpha1;
  int32_T c9_d_lda;
  int32_T c9_d_ldb;
  real_T c9_beta1;
  int32_T c9_d_ldc;
  char_T c9_TRANSB;
  char_T c9_TRANSA;
  int32_T c9_var;
  ptrdiff_t c9_m_t;
  ptrdiff_t c9_n_t;
  int32_T c9_b_var;
  ptrdiff_t c9_k_t;
  int32_T c9_c_var;
  ptrdiff_t c9_lda_t;
  int32_T c9_d_var;
  ptrdiff_t c9_ldb_t;
  int32_T c9_e_var;
  ptrdiff_t c9_ldc_t;
  double * c9_alpha1_t;
  double * c9_Aia0_t;
  double * c9_Bib0_t;
  double * c9_beta1_t;
  double * c9_Cic0_t;
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  c9_b_m = c9_m;
  c9_b_k = c9_k;
  c9_b_lda = c9_lda;
  c9_b_ldb = c9_ldb;
  c9_b_ldc = c9_ldc;
  if (c9_eml_use_refblas(chartInstance)) {
    c9_e_eml_scalar_eg(chartInstance);
    c9_e_eml_scalar_eg(chartInstance);
    c9_c_m = c9_b_m;
    c9_c_k = c9_b_k;
    c9_c_lda = c9_b_lda;
    c9_c_ldb = c9_b_ldb;
    c9_c_ldc = c9_b_ldc;
    if (c9_c_m == 0) {
    } else {
      c9_lda1 = c9_c_lda;
      c9_ldb1 = c9_c_ldb;
      c9_ldc1 = c9_c_ldc;
      c9_c = 0;
      c9_b = c9_c;
      c9_lastColC = c9_b;
      c9_b_ldc1 = c9_ldc1;
      c9_b_lastColC = c9_lastColC;
      c9_d = c9_b_ldc1;
      c9_b_b = c9_b_lastColC;
      c9_b_d = c9_d;
      c9_c_b = c9_b_b;
      guard2 = FALSE;
      if (c9_b_d == 0) {
        guard2 = TRUE;
      } else if (0 > c9_c_b) {
        guard2 = TRUE;
      } else {
        c9_overflow = (c9_c_b > MAX_int32_T - c9_b_d);
      }

      if (guard2 == TRUE) {
        c9_overflow = FALSE;
      }

      if (c9_overflow) {
        c9_check_forloop_overflow_error(chartInstance, TRUE);
      }

      c9_cr = 0;
      while ((c9_b_ldc1 > 0) && (c9_cr <= c9_b_lastColC)) {
        c9_b_cr = c9_cr;
        c9_a = c9_b_cr + 1;
        c9_i115 = c9_a;
        c9_b_a = c9_b_cr;
        c9_d_b = c9_c_m;
        c9_i116 = c9_b_a + c9_d_b;
        c9_c_a = c9_i115;
        c9_e_b = c9_i116;
        c9_d_a = c9_c_a;
        c9_f_b = c9_e_b;
        if (c9_d_a > c9_f_b) {
          c9_b_overflow = FALSE;
        } else {
          c9_b_overflow = (c9_f_b > 2147483646);
        }

        if (c9_b_overflow) {
          c9_check_forloop_overflow_error(chartInstance, TRUE);
        }

        for (c9_ic = c9_i115; c9_ic <= c9_i116; c9_ic++) {
          c9_b_ic = c9_ic;
          c9_C_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_ic, 1, c9_C_sizes[0], 1,
            0) - 1] = 0.0;
        }

        c9_cr += c9_b_ldc1;
      }

      c9_br = 0;
      c9_c_ldc1 = c9_ldc1;
      c9_c_lastColC = c9_lastColC;
      c9_c_d = c9_c_ldc1;
      c9_g_b = c9_c_lastColC;
      c9_d_d = c9_c_d;
      c9_h_b = c9_g_b;
      guard1 = FALSE;
      if (c9_d_d == 0) {
        guard1 = TRUE;
      } else if (0 > c9_h_b) {
        guard1 = TRUE;
      } else {
        c9_c_overflow = (c9_h_b > MAX_int32_T - c9_d_d);
      }

      if (guard1 == TRUE) {
        c9_c_overflow = FALSE;
      }

      if (c9_c_overflow) {
        c9_check_forloop_overflow_error(chartInstance, TRUE);
      }

      c9_c_cr = 0;
      while ((c9_c_ldc1 > 0) && (c9_c_cr <= c9_c_lastColC)) {
        c9_b_cr = c9_c_cr;
        c9_ar = 0;
        c9_e_a = c9_br + 1;
        c9_i117 = c9_e_a;
        c9_f_a = c9_br;
        c9_i_b = c9_c_k;
        c9_i118 = c9_f_a + c9_i_b;
        c9_g_a = c9_i117;
        c9_j_b = c9_i118;
        c9_h_a = c9_g_a;
        c9_k_b = c9_j_b;
        if (c9_h_a > c9_k_b) {
          c9_d_overflow = FALSE;
        } else {
          c9_d_overflow = (c9_k_b > 2147483646);
        }

        if (c9_d_overflow) {
          c9_check_forloop_overflow_error(chartInstance, TRUE);
        }

        for (c9_ib = c9_i117; c9_ib <= c9_i118; c9_ib++) {
          c9_b_ib = c9_ib;
          if (c9_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_ib, 1, 4, 1, 0) - 1] !=
              0.0) {
            c9_temp = c9_B[c9_b_ib - 1];
            c9_ia = c9_ar;
            c9_i_a = c9_b_cr + 1;
            c9_i119 = c9_i_a;
            c9_j_a = c9_b_cr;
            c9_l_b = c9_c_m;
            c9_i120 = c9_j_a + c9_l_b;
            c9_k_a = c9_i119;
            c9_m_b = c9_i120;
            c9_l_a = c9_k_a;
            c9_n_b = c9_m_b;
            if (c9_l_a > c9_n_b) {
              c9_e_overflow = FALSE;
            } else {
              c9_e_overflow = (c9_n_b > 2147483646);
            }

            if (c9_e_overflow) {
              c9_check_forloop_overflow_error(chartInstance, TRUE);
            }

            for (c9_c_ic = c9_i119; c9_c_ic <= c9_i120; c9_c_ic++) {
              c9_b_ic = c9_c_ic;
              c9_m_a = c9_ia + 1;
              c9_ia = c9_m_a;
              c9_C_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_ic, 1, c9_C_sizes[0],
                1, 0) - 1] = c9_C_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_b_ic,
                1, c9_C_sizes[0], 1, 0) - 1] + c9_temp *
                c9_A_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c9_ia, 1, c9_A_sizes[0]
                * c9_A_sizes[1], 1, 0) - 1];
            }
          }

          c9_n_a = c9_ar;
          c9_o_b = c9_lda1;
          c9_ar = c9_n_a + c9_o_b;
        }

        c9_o_a = c9_br;
        c9_p_b = c9_ldb1;
        c9_br = c9_o_a + c9_p_b;
        c9_c_cr += c9_c_ldc1;
      }
    }
  } else {
    c9_below_threshold(chartInstance);
    if (c9_b_m < 1) {
    } else if (c9_b_k < 1) {
    } else {
      c9_d_m = c9_b_m;
      c9_d_k = c9_b_k;
      c9_alpha1 = 1.0;
      c9_d_lda = c9_b_lda;
      c9_d_ldb = c9_b_ldb;
      c9_beta1 = 0.0;
      c9_d_ldc = c9_b_ldc;
      c9_TRANSB = 'N';
      c9_TRANSA = 'N';
      c9_var = c9_d_m;
      c9_m_t = (ptrdiff_t)(c9_var);
      c9_n_t = (ptrdiff_t)(1);
      c9_b_var = c9_d_k;
      c9_k_t = (ptrdiff_t)(c9_b_var);
      c9_c_var = c9_d_lda;
      c9_lda_t = (ptrdiff_t)(c9_c_var);
      c9_d_var = c9_d_ldb;
      c9_ldb_t = (ptrdiff_t)(c9_d_var);
      c9_e_var = c9_d_ldc;
      c9_ldc_t = (ptrdiff_t)(c9_e_var);
      c9_alpha1_t = (double *)(&c9_alpha1);
      _SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_A_sizes[0] * c9_A_sizes[1], 1, 0);
      c9_Aia0_t = (double *)(&c9_A_data[0]);
      c9_Bib0_t = (double *)(&c9_B[0]);
      c9_beta1_t = (double *)(&c9_beta1);
      _SFD_EML_ARRAY_BOUNDS_CHECK("", 1, 1, c9_C_sizes[0], 1, 0);
      c9_Cic0_t = (double *)(&c9_C_data[0]);
      dgemm(&c9_TRANSA, &c9_TRANSB, &c9_m_t, &c9_n_t, &c9_k_t, c9_alpha1_t,
            c9_Aia0_t, &c9_lda_t, c9_Bib0_t, &c9_ldb_t, c9_beta1_t, c9_Cic0_t,
            &c9_ldc_t);
    }
  }
}

static void c9_b_tanh(SFc9_simulationInstanceStruct *chartInstance, real_T
                      c9_x_data[20], int32_T c9_x_sizes[1])
{
  real_T c9_d12;
  int32_T c9_i121;
  int32_T c9_k;
  real_T c9_b_k;
  real_T c9_x;
  real_T c9_b_x;
  c9_d12 = (real_T)c9_x_sizes[0];
  c9_i121 = (int32_T)c9_d12 - 1;
  for (c9_k = 0; c9_k <= c9_i121; c9_k++) {
    c9_b_k = 1.0 + (real_T)c9_k;
    c9_x = c9_x_data[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      c9_b_k, 1, c9_x_sizes[0], 1, 0) - 1];
    c9_b_x = c9_x;
    c9_b_x = muDoubleScalarTanh(c9_b_x);
    c9_x_data[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)c9_b_k,
      1, c9_x_sizes[0], 1, 0) - 1] = c9_b_x;
  }
}

static void init_dsm_address_info(SFc9_simulationInstanceStruct *chartInstance)
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

void sf_c9_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1987238611U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(147245680U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2264701406U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2655748693U);
}

mxArray *sf_c9_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("Lu03muUfFEvEfeDapiDzm");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

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
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c9_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c9_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c9_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"v_tilde\",},{M[8],M[0],T\"is_active_c9_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c9_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc9_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc9_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           9,
           1,
           1,
           5,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"Gam_est");
          _SFD_SET_DATA_PROPS(1,2,0,1,"v_tilde");
          _SFD_SET_DATA_PROPS(2,10,0,0,"n");
          _SFD_SET_DATA_PROPS(3,10,0,0,"k_xi");
          _SFD_SET_DATA_PROPS(4,10,0,0,"i");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,178);

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c9_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c9_sf_marshallOut,(MexInFcnForType)c9_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c9_sf_marshallOut,(MexInFcnForType)c9_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c9_sf_marshallOut,(MexInFcnForType)c9_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c9_sf_marshallOut,(MexInFcnForType)c9_sf_marshallIn);

        {
          real_T *c9_v_tilde;
          real_T (*c9_Gam_est)[4];
          c9_v_tilde = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          c9_Gam_est = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c9_Gam_est);
          _SFD_SET_DATA_VALUE_PTR(1U, c9_v_tilde);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c9_n);
          _SFD_SET_DATA_VALUE_PTR(3U, &chartInstance->c9_k_xi);
          _SFD_SET_DATA_VALUE_PTR(4U, &chartInstance->c9_i);
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
  return "MqMLod0FmHsYrBUgfZ425F";
}

static void sf_opaque_initialize_c9_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc9_simulationInstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c9_simulation((SFc9_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c9_simulation((SFc9_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c9_simulation(void *chartInstanceVar)
{
  enable_c9_simulation((SFc9_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c9_simulation(void *chartInstanceVar)
{
  disable_c9_simulation((SFc9_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c9_simulation(void *chartInstanceVar)
{
  sf_c9_simulation((SFc9_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c9_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c9_simulation
    ((SFc9_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c9_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c9_simulation(SimStruct* S, const mxArray *
  st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c9_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c9_simulation((SFc9_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c9_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c9_simulation(S);
}

static void sf_opaque_set_sim_state_c9_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c9_simulation(S, st);
}

static void sf_opaque_terminate_c9_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc9_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c9_simulation((SFc9_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc9_simulation((SFc9_simulationInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c9_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c9_simulation((SFc9_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c9_simulation(SimStruct *S)
{
  /* Actual parameters from chart:
     i k_xi n
   */
  const char_T *rtParamNames[] = { "i", "k_xi", "n" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for i*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for k_xi*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);

  /* registration for n*/
  ssRegDlgParamAsRunTimeParam(S, 2, 2, rtParamNames[2], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      9);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,9,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,9,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,9);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,9,1);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,9,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 1; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,9);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2479061222U));
  ssSetChecksum1(S,(170206450U));
  ssSetChecksum2(S,(1066122733U));
  ssSetChecksum3(S,(14096751U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c9_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c9_simulation(SimStruct *S)
{
  SFc9_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc9_simulationInstanceStruct *)utMalloc(sizeof
    (SFc9_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc9_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c9_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c9_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c9_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c9_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c9_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c9_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c9_simulation;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c9_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c9_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c9_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c9_simulation;
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

void c9_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c9_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c9_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c9_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c9_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
