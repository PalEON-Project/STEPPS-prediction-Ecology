// manual gradient for the prediction model

#define EIGEN_DONT_PARALLELIZE

#include <string>
#include <vector>

#include <stan/callbacks/logger.hpp>
#include <stan/io/cmd_line.hpp>
#include <stan/io/dump.hpp>
#include <stan/io/writer.hpp>
#include <stan/io/reader.hpp>
#include <stan/callbacks/logger.hpp>
#include <stan/model/prob_grad.hpp>
#include <stan/services/util/mcmc_writer.hpp>

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/CXX11/Tensor>

#include <boost/math/special_functions/digamma.hpp>

#include "timer.hpp"

using namespace std;

namespace pred_model_namespace {

  using namespace stan::math;
  using namespace stan::io;
  using namespace Eigen;

  using stan::math::digamma;
  using stan::math::lgamma;
  using stan::math::digamma;

  typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
  typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
  typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;


    static double const ksi_lb = 0;
    static double const ksi_ub = 5; // used to be 20!
    static double const omega_lb = 0;
    static double const omega_ub = 1;
    //static double const lambda_lb = 0;
    //static double const lambda_ub = 1;

    static double const mu_std = 2; //from 5 // used to be 20!
  

  class pred_model : public stan::model::prob_grad {
private:
    int K;
    int N;
    int T;
    int N_knots;
    int N_cores;
    double res;
    vector<vector<int> > y;
    vector_d rho;
    vector_d eta;
    vector_d gamma;
    vector_d sum_w_pot;
    vector_d phi;
    vector<int> idx_cores;
    matrix_d d_knots;
    matrix_d d_inter;
    vector<matrix_d> w;
    matrix_d lag;
    vector_d eta2;
    vector_d zeros;
    vector_d ones;
    matrix_d Eye_knots;
    vector<matrix_d> C_s;
    vector<Eigen::LLT<matrix_d> > llt_s;
    vector<matrix_d> C_s_L;
    vector<matrix_d> C_s_inv;
    vector<matrix_d> c_s;
    vector<matrix_d> H_s;
    //matrix_d M;
    vector<matrix_d> M_H_s;

    vector_d lambda;

    Eigen::HouseholderQR<vector_d> qr;
    matrix_d fatQ;
    vector_d Q;
    row_vector_d QT;

public:
    pred_model(stan::io::var_context& context__,
	       unsigned int random_seed,	       
	       std::ostream* pstream__ = 0)
        : prob_grad::prob_grad(0) {
        static const char* function__ = "pred_model_namespace::pred_model(%1%)";
        (void) function__; // dummy call to supress warning
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        context__.validate_dims("data initialization", "K", "int", context__.to_vec());
        K = int(0);
        vals_i__ = context__.vals_i("K");
        pos__ = 0;
        K = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N", "int", context__.to_vec());
        N = int(0);
        vals_i__ = context__.vals_i("N");
        pos__ = 0;
        N = vals_i__[pos__++];
        context__.validate_dims("data initialization", "T", "int", context__.to_vec());
        T = int(0);
        vals_i__ = context__.vals_i("T");
        pos__ = 0;
        T = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N_knots", "int", context__.to_vec());
        N_knots = int(0);
        vals_i__ = context__.vals_i("N_knots");
        pos__ = 0;
        N_knots = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N_cores", "int", context__.to_vec());
        N_cores = int(0);
        vals_i__ = context__.vals_i("N_cores");
        pos__ = 0;
        N_cores = vals_i__[pos__++];
        context__.validate_dims("data initialization", "y", "int", context__.to_vec((N_cores * T),K));
        stan::math::validate_non_negative_index("y", "(N_cores * T)", (N_cores * T));
        stan::math::validate_non_negative_index("y", "K", K);
        y = std::vector<std::vector<int> >((N_cores * T),std::vector<int>(K,int(0)));
        vals_i__ = context__.vals_i("y");
        pos__ = 0;
        size_t y_limit_1__ = K;
        for (size_t i_1__ = 0; i_1__ < y_limit_1__; ++i_1__) {
            size_t y_limit_0__ = (N_cores * T);
            for (size_t i_0__ = 0; i_0__ < y_limit_0__; ++i_0__) {
                y[i_0__][i_1__] = vals_i__[pos__++];
            }
        }
        res = double(0.0);
        vals_r__ = context__.vals_r("res");
        pos__ = 0;
        res = vals_r__[pos__++];
        stan::math::validate_non_negative_index("rho", "K", K);
        rho = vector_d(K);
        context__.validate_dims("data initialization", "rho", "vector_d", context__.to_vec(K));
        vals_r__ = context__.vals_r("rho");
        pos__ = 0;
        size_t rho_i_vec_lim__ = K;
        for (size_t i_vec__ = 0; i_vec__ < rho_i_vec_lim__; ++i_vec__) {
            rho[i_vec__] = vals_r__[pos__++];
        }
        stan::math::validate_non_negative_index("eta", "K", K);
        eta = vector_d(K);
        context__.validate_dims("data initialization", "eta", "vector_d", context__.to_vec(K));
        vals_r__ = context__.vals_r("eta");
        pos__ = 0;
        size_t eta_i_vec_lim__ = K;
        for (size_t i_vec__ = 0; i_vec__ < eta_i_vec_lim__; ++i_vec__) {
            eta[i_vec__] = vals_r__[pos__++];
        }
	lambda = rho;
        context__.validate_dims("data initialization", "gamma", "double", context__.to_vec(K));
        gamma = vector_d(K);
        vals_r__ = context__.vals_r("gamma");
        pos__ = 0;
	for (int k = 0; k < K; ++k) {
	  gamma[k] = vals_r__[pos__++];
	}
        // context__.validate_dims("data initialization", "psi", "double", context__.to_vec());
        // psi = double(0);
        // vals_r__ = context__.vals_r("psi");
        // pos__ = 0;
        // psi = vals_r__[pos__++];
        stan::math::validate_non_negative_index("sum_w_pot", "K", K);
        sum_w_pot = vector_d(K);
        context__.validate_dims("data initialization", "sum_w_pot", "vector_d", context__.to_vec(K));
        vals_r__ = context__.vals_r("sum_w_pot");
        pos__ = 0;
        size_t sum_w_pot_i_vec_lim__ = K;
        for (size_t i_vec__ = 0; i_vec__ < sum_w_pot_i_vec_lim__; ++i_vec__) {
            sum_w_pot[i_vec__] = vals_r__[pos__++];
        }
        stan::math::validate_non_negative_index("phi", "K", K);
        phi = vector_d(K);
        context__.validate_dims("data initialization", "phi", "vector_d", context__.to_vec(K));
        vals_r__ = context__.vals_r("phi");
        pos__ = 0;
        size_t phi_i_vec_lim__ = K;
        for (size_t i_vec__ = 0; i_vec__ < phi_i_vec_lim__; ++i_vec__) {
            phi[i_vec__] = vals_r__[pos__++];
        }
        context__.validate_dims("data initialization", "idx_cores", "int", context__.to_vec(N_cores));
        stan::math::validate_non_negative_index("idx_cores", "N_cores", N_cores);
        idx_cores = std::vector<int>(N_cores,int(0));
        vals_i__ = context__.vals_i("idx_cores");
        pos__ = 0;
        size_t idx_cores_limit_0__ = N_cores;
        for (size_t i_0__ = 0; i_0__ < idx_cores_limit_0__; ++i_0__) {
            idx_cores[i_0__] = vals_i__[pos__++];
        }
        // context__.validate_dims("data initialization", "d", "matrix_d", context__.to_vec(N,N));
        // stan::math::validate_non_negative_index("d", "N", N);
        // stan::math::validate_non_negative_index("d", "N", N);
        // d = matrix_d(N,N);
        // vals_r__ = context__.vals_r("d");
        // pos__ = 0;
        // size_t d_m_mat_lim__ = N;
        // size_t d_n_mat_lim__ = N;
        // for (size_t n_mat__ = 0; n_mat__ < d_n_mat_lim__; ++n_mat__) {
        //     for (size_t m_mat__ = 0; m_mat__ < d_m_mat_lim__; ++m_mat__) {
        //         d(m_mat__,n_mat__) = vals_r__[pos__++];
        //     }
        // }
        context__.validate_dims("data initialization", "d_knots", "matrix_d", context__.to_vec(N_knots,N_knots));
        stan::math::validate_non_negative_index("d_knots", "N_knots", N_knots);
        stan::math::validate_non_negative_index("d_knots", "N_knots", N_knots);
        d_knots = matrix_d(N_knots,N_knots);
        vals_r__ = context__.vals_r("d_knots");
        pos__ = 0;
        size_t d_knots_m_mat_lim__ = N_knots;
        size_t d_knots_n_mat_lim__ = N_knots;
        for (size_t n_mat__ = 0; n_mat__ < d_knots_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < d_knots_m_mat_lim__; ++m_mat__) {
                d_knots(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "d_inter", "matrix_d", context__.to_vec(N,N_knots));
        stan::math::validate_non_negative_index("d_inter", "N", N);
        stan::math::validate_non_negative_index("d_inter", "N_knots", N_knots);
        d_inter = matrix_d(N,N_knots);
        vals_r__ = context__.vals_r("d_inter");
        pos__ = 0;
        size_t d_inter_m_mat_lim__ = N;
        size_t d_inter_n_mat_lim__ = N_knots;
        for (size_t n_mat__ = 0; n_mat__ < d_inter_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < d_inter_m_mat_lim__; ++m_mat__) {
                d_inter(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "w", "matrix_d", context__.to_vec(K,N_cores,N));
        stan::math::validate_non_negative_index("w", "K", K);
        stan::math::validate_non_negative_index("w", "N_cores", N_cores);
        stan::math::validate_non_negative_index("w", "N", N);

        w = vector<matrix_d>(K,matrix_d(N_cores,N));
        vals_r__ = context__.vals_r("w");
        pos__ = 0;
        size_t w_m_mat_lim__ = N_cores;
        size_t w_n_mat_lim__ = N;
	for (size_t k = 0; k < K; ++k) {
	  for (size_t n_mat__ = 0; n_mat__ < w_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < w_m_mat_lim__; ++m_mat__) {
	      w[k](m_mat__,n_mat__) = vals_r__[pos__++];
            }
	  }
	}

        context__.validate_dims("data initialization", "lag", "matrix_d", context__.to_vec(T,T));
        stan::math::validate_non_negative_index("lag", "T", T);
        stan::math::validate_non_negative_index("lag", "T", T);
        lag = matrix_d(T,T);
        vals_r__ = context__.vals_r("lag");
        pos__ = 0;
        size_t lag_m_mat_lim__ = T;
        size_t lag_n_mat_lim__ = T;
        for (size_t n_mat__ = 0; n_mat__ < lag_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < lag_m_mat_lim__; ++m_mat__) {
                lag(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        // context__.validate_dims("data initialization", "N_p", "int", context__.to_vec());
        // N_p = int(0);
        // vals_i__ = context__.vals_i("N_p");
        // pos__ = 0;
        // N_p = vals_i__[pos__++];
        // context__.validate_dims("data initialization", "P", "double", context__.to_vec());
        // P = double(0);
        // vals_r__ = context__.vals_r("P");
        // pos__ = 0;
        // P = vals_r__[pos__++];

        // validate data
        try {
	  check_greater_or_equal(function__,"K",K,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of K: ") + std::string(e.what()));
        };
        try {
	  check_greater_or_equal(function__,"N",N,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of N: ") + std::string(e.what()));
        };
        try {
	  check_greater_or_equal(function__,"T",T,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of T: ") + std::string(e.what()));
        };
        try {
	  check_greater_or_equal(function__,"N_knots",N_knots,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of N_knots: ") + std::string(e.what()));
        };
        try {
	  check_greater_or_equal(function__,"N_cores",N_cores,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of N_cores: ") + std::string(e.what()));
        };
        // try {
	//   check_greater_or_equal(function__,"res",res,0);
        // } catch (std::domain_error& e) {
        //     throw std::domain_error(std::string("Invalid value of res: ") + std::string(e.what()));
        // };
        try {
	  check_greater_or_equal(function__,"rho",rho,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of rho: ") + std::string(e.what()));
        };
        try {
	  check_greater_or_equal(function__,"eta",eta,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of eta: ") + std::string(e.what()));
        };
        try {
	  check_greater_or_equal(function__,"gamma",gamma,0);
	  check_less_or_equal(function__,"gamma",gamma,1);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of gamma: ") + std::string(e.what()));
        };
        try {
	  check_greater_or_equal(function__,"sum_w_pot",sum_w_pot,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of sum_w_pot: ") + std::string(e.what()));
        };
        try {
	  check_greater_or_equal(function__,"phi",phi,0);
        } catch (std::domain_error& e) {
            throw std::domain_error(std::string("Invalid value of phi: ") + std::string(e.what()));
        };
        // try {
        //     check_greater_or_equal(function__,N_p,0,"N_p");
        // } catch (std::domain_error& e) {
        //     throw std::domain_error(std::string("Invalid value of N_p: ") + std::string(e.what()));
        // };
        stan::math::validate_non_negative_index("eta2", "K", K);
        eta2 = vector_d(K);
        stan::math::validate_non_negative_index("zeros", "N_knots", N_knots);
        zeros = vector_d(N_knots);
	stan::math::validate_non_negative_index("ones", "N", N);
	ones = vector_d(N);
        stan::math::validate_non_negative_index("Eye_knots", "N_knots", N_knots);
        stan::math::validate_non_negative_index("Eye_knots", "N_knots", N_knots);
        Eye_knots = matrix_d(N_knots,N_knots);
        stan::math::validate_non_negative_index("C_s", "K", K);
        stan::math::validate_non_negative_index("C_s", "N_knots", N_knots);
        stan::math::validate_non_negative_index("C_s", "N_knots", N_knots);
        C_s = std::vector<matrix_d>(K,matrix_d(N_knots,N_knots));
        stan::math::validate_non_negative_index("C_s_L", "K", K);
        stan::math::validate_non_negative_index("C_s_L", "N_knots", N_knots);
        stan::math::validate_non_negative_index("C_s_L", "N_knots", N_knots);
        C_s_L = std::vector<matrix_d>(K,matrix_d(N_knots,N_knots));
        stan::math::validate_non_negative_index("C_s_inv", "K", K);
        stan::math::validate_non_negative_index("C_s_inv", "N_knots", N_knots);
        stan::math::validate_non_negative_index("C_s_inv", "N_knots", N_knots);
        C_s_inv = std::vector<matrix_d>(K,matrix_d(N_knots,N_knots));
        stan::math::validate_non_negative_index("c_s", "K", K);
        stan::math::validate_non_negative_index("c_s", "N", N);
        stan::math::validate_non_negative_index("c_s", "N_knots", N_knots);
        c_s = std::vector<matrix_d>(K,matrix_d(N,N_knots));
        stan::math::validate_non_negative_index("H_s", "K", K);
        stan::math::validate_non_negative_index("H_s", "N", N);
        stan::math::validate_non_negative_index("H_s", "N_knots", N_knots);
        H_s = std::vector<matrix_d>(K,matrix_d(N,N_knots));
        // stan::math::validate_non_negative_index("M", "(N * T)", (N * T));
        // stan::math::validate_non_negative_index("M", "(N * T)", (N * T));
        // M = matrix_d((N * T),(N * T));
	stan::math::validate_non_negative_index("M_H_s", "N", N);
        stan::math::validate_non_negative_index("M_H_s", "N_knots", N_knots);
        M_H_s = std::vector<matrix_d>(K,matrix_d(N,N_knots));

        for (int k = 1; k <= K; ++k) {
          stan::math::assign(get_base1_lhs(eta2,k,"eta2",1), (get_base1(eta,k,"eta",1) * get_base1(eta,k,"eta",1)));
        }
        for (int i = 1; i <= N_knots; ++i) {
            stan::math::assign(get_base1_lhs(zeros,i,"zeros",1), 0.0);
        }
	for (int i = 1; i <= N; ++i) {
	  stan::math::assign(get_base1_lhs(ones,i,"ones",1), 1.0);
	}
        for (int i = 1; i <= N_knots; ++i) {
            for (int j = 1; j <= N_knots; ++j) {
                if (as_bool(logical_eq(i,j))) {
                    stan::math::assign(get_base1_lhs(Eye_knots,i,j,"Eye_knots",1), 1.0);
                } else {
                    stan::math::assign(get_base1_lhs(Eye_knots,i,j,"Eye_knots",1), 0);
                }
            }
        }

      qr   = ones.householderQr();
      fatQ = qr.householderQ();
      Q    = fatQ.col(0); // thin Q, as returned in R
      QT = Q.transpose();

      llt_s.resize(K);
      c_s.resize(K);
      C_s_L.resize(K);
      C_s_inv.resize(K);
      H_s.resize(K);
      M_H_s.resize(K);

      for (int k = 0; k < K; ++k) {
	C_s[k] = exp(-d_knots.array() / rho[k]).matrix();

	llt_s[k]    = C_s[k].llt();
	C_s_L[k]    = llt_s[k].matrixL();
	C_s_inv[k]  = llt_s[k].solve(Eigen::MatrixXd::Identity(N_knots,N_knots));
	c_s[k]      = exp(- 1.0/rho[k] * d_inter.array()).matrix();
	H_s[k]      = c_s[k] * C_s_inv[k];
	M_H_s[k] = H_s[k] - Q * ( QT * H_s[k]);
	//w = res * res * w;
      }

        // validate transformed data

        // set parameter ranges
        num_params_r__ = 0U;
        param_ranges_i__.clear();
        ++num_params_r__;
        num_params_r__ += K;
        num_params_r__ += K;
	num_params_r__ += (T-1) * K;
        num_params_r__ += N_knots * K;
        num_params_r__ += N_knots * (K * (T-1));
        num_params_r__ += (N * T) * K;
    }

    ~pred_model() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
			 std::stringstream* msgs) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;


        if (!(context__.contains_r("ksi")))
            throw std::runtime_error("variable ksi missing");
        vals_r__ = context__.vals_r("ksi");
        pos__ = 0U;
        context__.validate_dims("initialization", "ksi", "double", context__.to_vec());
        double ksi(0);
        ksi = vals_r__[pos__++];
        try { writer__.scalar_lub_unconstrain(ksi_lb,ksi_ub,ksi); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable ksi: ") + e.what()); }

        if (!(context__.contains_r("omega")))
            throw std::runtime_error("variable omega missing");
        vals_r__ = context__.vals_r("omega");
        pos__ = 0U;
        context__.validate_dims("initialization", "omega", "vector_d", context__.to_vec(K));
        vector_d omega(K);
        for (int j1__ = 0U; j1__ < K; ++j1__)
            omega(j1__) = vals_r__[pos__++];
        try { writer__.vector_lub_unconstrain(omega_lb,omega_ub,omega); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable omega: ") + e.what()); }

        if (!(context__.contains_r("mu")))
            throw std::runtime_error("variable mu missing");
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        context__.validate_dims("initialization", "mu", "vector_d", context__.to_vec(K));
        vector_d mu(K);
        for (int j1__ = 0U; j1__ < K; ++j1__)
            mu(j1__) = vals_r__[pos__++];
        try { writer__.vector_unconstrain(mu); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable mu: ") + e.what()); }

        if (!(context__.contains_r("mu_t")))
            throw std::runtime_error("variable mu_t missing");
        vals_r__ = context__.vals_r("mu_t");
        pos__ = 0U;
        context__.validate_dims("initialization", "mu_t", "vector_d", context__.to_vec(K,T-1));
        std::vector<vector_d> mu_t(K, vector_d(T-1));
        for (int j1__ = 0U; j1__ < T-1; ++j1__)
            for (int i0__ = 0U; i0__ < K; ++i0__)
                mu_t[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < K; ++i0__)
            try { writer__.vector_unconstrain(mu_t[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable mu_t: ") + e.what()); }

        if (!(context__.contains_r("alpha_s")))
            throw std::runtime_error("variable alpha_s missing");
        vals_r__ = context__.vals_r("alpha_s");
        pos__ = 0U;
        context__.validate_dims("initialization", "alpha_s", "vector_d", context__.to_vec(K,N_knots));
        std::vector<vector_d> alpha_s(K,vector_d(N_knots));
        for (int j1__ = 0U; j1__ < N_knots; ++j1__)
            for (int i0__ = 0U; i0__ < K; ++i0__)
                alpha_s[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < K; ++i0__)
            try { writer__.vector_unconstrain(alpha_s[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable alpha_s: ") + e.what()); }

        if (!(context__.contains_r("alpha_t")))
            throw std::runtime_error("variable alpha_t missing");
        vals_r__ = context__.vals_r("alpha_t");
        pos__ = 0U;
        context__.validate_dims("initialization", "alpha_t", "vector_d", context__.to_vec((K * (T-1)),N_knots));
        std::vector<vector_d> alpha_t((K * (T-1)),vector_d(N_knots));
        for (int j1__ = 0U; j1__ < N_knots; ++j1__)
	  for (int i0__ = 0U; i0__ < (K * (T-1)); ++i0__)
                alpha_t[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < (K * (T-1)); ++i0__)
            try { writer__.vector_unconstrain(alpha_t[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable alpha_t: ") + e.what()); }

        if (!(context__.contains_r("g")))
            throw std::runtime_error("variable g missing");
        vals_r__ = context__.vals_r("g");
        pos__ = 0U;
        context__.validate_dims("initialization", "g", "vector_d", context__.to_vec(K,(N * T)));
        std::vector<vector_d> g(K,vector_d((N * T)));
        for (int j1__ = 0U; j1__ < (N * T); ++j1__)
            for (int i0__ = 0U; i0__ < K; ++i0__)
                g[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < K; ++i0__)
            try { writer__.vector_unconstrain(g[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable g: ") + e.what()); }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

   //  void transform_inits(const stan::io::var_context& context__,
   //                       std::vector<int>& params_i__,
   //                       std::vector<double>& params_r__,
   // 			 std::stringstream& msgs) const {
   // void transform_inits(const stan::io::var_context& context,
   //                       Eigen::VectorXd& params_r,
   //                       ostream* output) const {
   //    std::vector<double> params_r_vec;
   //    std::vector<int> params_i_vec;
   //    std::stringstream msg;
   //    transform_inits(context, params_i_vec, params_r_vec, &msg);
   //    params_r.resize(params_r_vec.size());
   //    for (int i = 0; i < params_r.size(); ++i)
   //      params_r(i) = params_r_vec[i];
   //  }


    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////


    inline double lub_transform(const double x, const double lb, const double ub,
				double &lja, double &ja, double &dj) const
    {
      double inv_logit_x;
      if (x > 0) {
        double exp_minus_x = exp(-x);
	double exp_minus_x_p1 = exp_minus_x + 1.0;
        inv_logit_x = 1.0 / (1.0 + exp_minus_x);
        lja = log(ub - lb) - x - 2 * stan::math::log1p(exp_minus_x);
	ja  = (ub - lb) * exp_minus_x / (exp_minus_x_p1 * exp_minus_x_p1);
	dj  = -1.0 + 2.0 * exp_minus_x / (1 + exp_minus_x);
        if ((x < std::numeric_limits<double>::infinity())
            && (inv_logit_x == 1))
	  inv_logit_x = 1 - 1e-15;
      } else {
        double exp_x = exp(x);
	double exp_x_p1 = exp_x + 1.0;
        inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
        lja = log(ub - lb) + x - 2 * stan::math::log1p(exp_x);
	ja  = (ub - lb) * exp_x / (exp_x_p1 * exp_x_p1);
	dj  = -1.0 + 2.0 / (exp_x + 1.0);
        if ((x > -std::numeric_limits<double>::infinity())
            && (inv_logit_x== 0))
	  inv_logit_x = 1e-15;
      }
      return lb + (ub - lb) * inv_logit_x;
    }

    double normal_log_double(const vector_d& y, const double mu, const double sigma) const {
      double lp = 0.0;
      double inv_sigma = 1.0/sigma;

      int size_y = y.size();

      for (int n = 0; n < size_y; n++) {
	const double y_minus_mu_over_sigma = (y[n] - mu) * inv_sigma;
	const double y_minus_mu_over_sigma_squared = y_minus_mu_over_sigma * y_minus_mu_over_sigma;
	lp -= 0.5 * y_minus_mu_over_sigma_squared;
      }

      return lp;
    }

    double normal_log_double(const double y, const double mu, const double sigma, const int sigma_fixed) const {
      double lp = 0.0;
      double inv_sigma = 1.0/sigma;
      double log_sigma = log(sigma);

      const double y_minus_mu_over_sigma = (y - mu) * inv_sigma;
      const double y_minus_mu_over_sigma_squared = y_minus_mu_over_sigma * y_minus_mu_over_sigma;

      lp -= 0.5 * y_minus_mu_over_sigma_squared;

      if (sigma_fixed != 1){
	lp -= log_sigma;
      }

      return lp;
    }

    double multi_normal_cholesky_log_double(const vector_d& y,
					    const vector_d& mu,
					    const matrix_d& L, bool const_matrix=false) const {

      double lp = 0.0;

      if (!const_matrix) {
        vector_d L_log_diag = L.diagonal().array().log().matrix();
        lp -= sum(L_log_diag);
      }

      vector_d y_minus_mu = y.array() - mu.array();
      vector_d half = mdivide_left_tri_low(L, y_minus_mu);
      lp -= 0.5 * dot_self(half);

      return lp;
    }

    double multi_normal_prec_log_double(const vector_d& y,
                                        const vector_d& mu,
                                        const matrix_d& Sigma) const {

      double lp = 0.0;
      if (y.rows() == 0)
        return lp;

      LDLT_factor<double,Eigen::Dynamic,Eigen::Dynamic> ldlt_Sigma(Sigma);

      lp += 0.5*log_determinant_ldlt(ldlt_Sigma);

      vector_d y_minus_mu(y.size());
      for (int i = 0; i < y.size(); i++)
        y_minus_mu(i) = y(i)-mu(i);

      lp -= 0.5 * (y_minus_mu.transpose() * Sigma * y_minus_mu).trace();

      return lp;
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////


    template <bool propto, bool jacobian, bool lp_only>
    double log_prob_grad(vector<double>& params_r,
                         vector<double>& gradient,
                         ostream* pstream = 0) const {

      double lp = 0.0;

      //
      // unpack model parameters and constrain
      //

      vector<int> params_i;
      stan::io::reader<double> in(params_r, params_i);

      double ksi, ksi_lja, ksi_ja, ksi_dj;
      ksi = lub_transform(in.scalar(), ksi_lb, ksi_ub, ksi_lja, ksi_ja, ksi_dj);

      if (jacobian)
      	lp += ksi_lja;

      vector_d omega(K), omega_lja(K), omega_ja(K), omega_dj(K);
      for (int i=0; i<K; ++i)
      	omega[i] = lub_transform(in.scalar(), omega_lb, omega_ub, omega_lja[i], omega_ja[i], omega_dj[i]);

      if (jacobian)
      	for (int i=0; i<K; ++i)
      	  lp += omega_lja[i];

      vector_d  mu = in.vector_constrain(K,lp);

      vector<vector_d> mu_t(K);
      for (int k = 0; k < K; ++k) {
      	if (jacobian)
      	  mu_t[k] = in.vector_constrain(T-1,lp);
      	else
      	  mu_t[k] = in.vector_constrain(T-1);
      }

      vector<vector_d> alpha_s(K);
      for (int k = 0; k<K; ++k) {
      	if (jacobian)
      	  alpha_s[k] = in.vector_constrain(N_knots,lp);
      	else
      	  alpha_s[k] = in.vector_constrain(N_knots);
      }

      vector<vector_d> alpha_t(K*(T-1));
      for (int k = 0; k<(K*(T-1)); ++k) {
      	if (jacobian)
      	  alpha_t[k] = in.vector_constrain(N_knots,lp);
      	else
      	  alpha_t[k] = in.vector_constrain(N_knots);
      }

      vector<vector_d> g(K);
      for (int k = 0; k < K; ++k) {
      	if (jacobian)
      	  g[k] = in.vector_constrain((N * T),lp);
      	else
      	  g[k] = in.vector_constrain((N * T));
      }

      //
      // compute log probability
      //

      // priors
      lp += normal_log_double(mu, 0, mu_std);

      vector<double> omega2(K);
      for (int k=0; k<K; ++k) {
      	omega2[k] = omega[k] * omega[k];
      }

      vector<vector_d> qvar(K);
      vector<vector_d> var_g(K);
      vector<vector_d> mu_g(K);
      matrix_d         exp_g(N*T, K);

      for (int k=0; k < K; ++k) {
	var_g[k].resize(N*T);
        mu_g[k].resize(N*T);
	qvar[k].resize(N*T);
      }

      double lp_thrd[K];
      for (int k = 0; k<K; ++k){
	lp_thrd[k] = 0;
      }

      #pragma omp parallel for
      for (int k = 0; k<K; ++k){
      	lp_thrd[k] = multi_normal_cholesky_log_double(alpha_s[k], zeros, eta[k] * C_s_L[k], true);

	vector_d M_Halpha_s = M_H_s[k] * alpha_s[k];

      	for (int i = 0; i<N; ++i){
      	  for (int t = 0; t<T; ++t){
	    mu_g[k][i * T + t] = M_Halpha_s[i];
      	  }
      	}

	//first nonzero alpha_t
	lp_thrd[k] += multi_normal_cholesky_log_double(alpha_t[k*(T-1)], zeros, eta[k] * C_s_L[k], false);

	for (int t=1; t<(T-1); ++t){
	  lp_thrd[k] += multi_normal_cholesky_log_double(alpha_t[k*(T-1) + t], omega[k] * alpha_t[k*(T-1) + t-1], sqrt((1 - omega2[k])) * eta[k] * C_s_L[k]);
      	}

	vector<vector_d> qQinv_alpha(T-1);
	for (int t=0; t<(T-1); ++t){
	  qQinv_alpha[t] = H_s[k] * alpha_t[k*(T-1) + t];
	}

      	// time-varying mean
	lp_thrd[k] += normal_log_double(mu_t[k][0], 0, ksi, 0);

      	for (int t=1; t<T-1; ++t){
      	  lp_thrd[k] += normal_log_double(mu_t[k][t], mu_t[k][t-1], ksi,  0);
      	}

	row_vector_d c_i;
	row_vector_d q_i;

      	for (int i = 0; i<N; ++i){
      	  c_i  = c_s[k].row(i);
      	  double cvar = eta2[k] * c_i * C_s_inv[k] * c_i.transpose();

      	  //q_i  = q_s[k].row(i);
      	  //qvar[k][i] = c_i * C_s_inv[k] * c_i.transpose();//q_i * Q_s_inv[k] * q_i.transpose();
	  qvar[k][i] = cvar;

      	  mu_g[k][i * T] = mu[k] + mu_g[k][i * T];
      	  var_g[k][i * T] = eta2[k] - cvar;

      	  for (int t=1; t<T; ++t){
      	    mu_g[k][i * T + t] = mu[k] + mu_t[k][t-1] + mu_g[k][i * T + t] +  qQinv_alpha[t-1][i];
      	    var_g[k][i * T + t] = eta2[k] - cvar +  eta2[k] * (1-omega2[k]) -  (1-omega2[k]) * qvar[k][i];
      	  }
      	}

      	for (int i = 0; i<N*T; ++i){
	  if (var_g[k][i] <= 1.e-8) {
	    var_g[k][i] = 0.0;
	  }
	  var_g[k][i] += 1.e-5;//1;
	  lp_thrd[k] += normal_log_double(g[k][i], mu_g[k][i], sqrt(var_g[k][i]), 0);
      	}

      	// exp_g.col(k) = exp(g[k]);
        for (int i=0; i<N*T; ++i)
          exp_g(i,k) = exp(g[k][i]);
      }

      for (int k=0; k<K; ++k)
      	lp += lp_thrd[k];

      //std::cout<<"LP after g : " << lp << std::endl;

      vector_d sum_exp_g = exp_g.rowwise().sum();
      matrix_d r(N*T,K);

      #pragma omp parallel for collapse(2)
      for (int k = 0; k < K; ++k)
      	for (int i = 0; i < N*T; ++i)
      	  r(i,k) = exp_g(i,k) / sum_exp_g(i);

      // #pragma omp parallel for
      // for (int i = 0; i < N*T; ++i)
      // 	r(i,K-1) = 1. / (1. + sum_exp_g(i));

      // #pragma omp parallel for
      // for (int k = 0; k < K; ++k)
      // 	for (int i = 0; i < N*T; ++i)
      // 	  r(i,k) = exp_g(i,k) / (1. + sum_exp_g(i));

      // #pragma omp parallel for
      // for (int i = 0; i < N*T; ++i)
      // 	r(i,K-1) = 1. / (1. + sum_exp_g(i));

      matrix_d r_new(N_cores * T, K);
      matrix_d out_sum(N_cores*T,K);
      vector<vector_d> sum_w(K);

      for (int k = 0; k < K; ++k) {
	sum_w[k] = vector_d(N*T);
      }

      out_sum.fill(0);

      #pragma omp parallel for collapse(2)
      for (int k = 0; k < K; ++k) {
	for (int i = 0; i < N_cores; ++i) {
	  for (int t = 0; t < T; ++t) {
	    int idx_core = (idx_cores[i] - 1) * T + t;
	    r_new(i * T + t, k) = gamma[k] * r(idx_core, k);
	    for (int j = 0; j < N; ++j) {
	      if (j != idx_cores[i] -1) {
		out_sum(i * T + t, k) += res * res * w[k](i,j) * r(j * T + t, k);
	      }
	    }

	    sum_w[k][i*T+t] = sum_w_pot[k];
	    r_new(i * T + t, k) += out_sum(i *T + t, k) * (1 - gamma[k]) / sum_w[k][i*T+t];
	  }
	}
      }

      vector<int> N_grains(N_cores*T);
      vector_d    A(N_cores*T);
      matrix_d    kappa(N_cores*T,K);
      vector<int> y_row_sum(N_cores*T);

      // lp! so don't use omp here
      for (int i = 0; i < N_cores * T; ++i) {
      	y_row_sum[i] = 0.0;
      	for (int k = 0; k < K; ++k)
      	  y_row_sum[i] += y[i][k];

      	if (y_row_sum[i] > 0) {
      	  for (int k = 0; k < K; ++k){
      	    kappa(i,k) = phi(k) * r_new(i,k);
      	  }

      	  A[i] = kappa.row(i).sum();
      	  N_grains[i] = y_row_sum[i];

      	  lp += lgamma(N_grains[i] + 1.0) + lgamma(A[i]) - lgamma(N_grains[i] + A[i]);
      	  for (int k = 0; k < K; ++k) {
      	    lp += -lgamma(y[i][k] + 1) + lgamma(y[i][k] + kappa(i,k)) - lgamma(kappa(i,k));
	  }
      	}
      }

      if (lp_only) {
        return lp;
      }

      fill(gradient.begin(), gradient.end(), 0.0);

      double inv_ksi1 = 1.0 / ksi;
      double inv_ksi2 = 1.0 / (ksi * ksi);
      double inv_ksi3 = 1.0 / (ksi * ksi * ksi);

      for (int k=0; k<K; ++k){
	gradient[0] += -inv_ksi1 + inv_ksi3 * mu_t[k][0] * mu_t[k][0];
	for (int t=1; t<T-1; ++t)
	  gradient[0] += -inv_ksi1 + inv_ksi3 * (mu_t[k][t] - mu_t[k][t-1]) * (mu_t[k][t] - mu_t[k][t-1]);
      }

      // ksi jacobian adjustment
      gradient[0] = gradient[0] * ksi_ja + ksi_dj;

      // partial of mu_t normal

      #pragma omp parallel for
      for (int k=0; k<K; ++k) {
	for (int t=1; t<T-1; ++t){
	  int idx_mut = 1 + 2*K + k*(T-1) + t;
	  gradient[idx_mut] -= inv_ksi2 * (mu_t[k][t] - mu_t[k][t-1]);

	  if (t<T-2){
	    gradient[idx_mut] +=  inv_ksi2 * (mu_t[k][t+1] - mu_t[k][t]);
	  }
	}

	// mu_t[k][0]
	gradient[1 + 2*K + k*(T-1)] -= inv_ksi2 * mu_t[k][0];
	gradient[1 + 2*K + k*(T-1)] += inv_ksi2 * (mu_t[k][1] - mu_t[k][0]);

	// partial of g wrt mu_t
      	for (int n=0; n<N; ++n){
      	  for (int t=0; t<T-1; ++t){
      	    double idx = n*T + t+1;
      	    gradient[1 + 2*K + k*(T-1) + t] += (g[k](idx) - mu_g[k][idx]) / var_g[k][idx];
      	  }
      	}

	// partial of MVN for alpha_s wrt to alpha_s
	vector_d tmp = 1 / eta2[k] * C_s_inv[k] * alpha_s[k];
	for (int v=0; v<N_knots; ++v) {
	  int idx_alpha_s = 1 + K*2 + K*(T-1) + k*N_knots + v;
	  gradient[idx_alpha_s] -= tmp[v];
	}

	double alpha_Qinv_alpha_first;
	alpha_Qinv_alpha_first = (alpha_t[k*(T-1)]).transpose() * C_s_inv[k] * alpha_t[k*(T-1)];

	for (int t=1; t<(T-1); ++t){

	  double alpha_Qinv_alpha_full;
	  alpha_Qinv_alpha_full = (alpha_t[k*(T-1)+t] - omega[k] * alpha_t[k*(T-1)+t-1]).transpose() * C_s_inv[k] * (alpha_t[k*(T-1)+t] - omega[k] * alpha_t[k*(T-1)+t-1]);

	  double alpha_Qinv_alpha_d1;
	  alpha_Qinv_alpha_d1 =  ( (alpha_t[k*(T-1)+t] - omega[k] * alpha_t[k*(T-1)+t-1]).transpose() * C_s_inv[k] *  alpha_t[k*(T-1)+t-1]);

	  gradient[1 + k] +=  N_knots * omega[k] / ( 1 - omega2[k]);
	  gradient[1 + k] +=  - omega[k] / eta2[k] / (1 - omega2[k]) / (1 - omega2[k]) * alpha_Qinv_alpha_full;
	  gradient[1 + k] += 1 / eta2[k] / (1 - omega2[k]) *  alpha_Qinv_alpha_d1;
	}

	// partial of mu prior
	gradient[1 + K + k] -=  mu[k] / (mu_std * mu_std);

	// partials of g normal
	for (int t=0; t<T; ++t){
	  for (int n=0; n<N; ++n){

	    double A      = g[k][n*T+t] - mu_g[k][n*T+t];
	    double B      = var_g[k][n*T+t];
	    double AoverB = A/B;

	    gradient[1 + K + k] += AoverB;

	    int idx_g        = 1 + K*2 + K*(T-1) + K*N_knots + K*(T-1)*N_knots + k*N*T + n*T + t;

	    // g wrt g
	    gradient[idx_g] -= AoverB;

	    if (t > 0){
	      // wrt sigma
	      gradient[1 + k] += (1 / B - AoverB * AoverB ) * omega[k] * (eta2[k] - qvar[k][n]);
	    } 

	    for (int v=0; v<N_knots; ++v){
	      int idx_alpha_s  = 1 + K*2 + K*(T-1) +  k*N_knots + v;
	      gradient[idx_alpha_s] += AoverB * M_H_s[k](n,v);
	    }

	    if (t > 0){
	      for (int v=0; v<N_knots; ++v){
	      	int idx_alpha_t = 1 + K*2 + K*(T-1) + K*N_knots + (k*(T-1) + t-1)*N_knots + v;
	      	gradient[idx_alpha_t] += AoverB * H_s[k](n,v);
	      }
	    }


      	  } // n
      	} // t

	// sigma jacobian adjustments
	gradient[1 + k] = gradient[1 + k] * omega_ja[k] + omega_dj[k];

	// partial of MVN for alpha_t (wrt alpha_t)
	vector_d Qinv_alphat_next;
	vector_d Qinv_alphat_now;

	for (int t=0; t<(T-1); ++t) {

	  if (t>0) {
	    Qinv_alphat_now = - 1 / eta2[k] / (1 - omega2[k]) *  C_s_inv[k] * (alpha_t[k*(T-1)+t] - omega[k] * alpha_t[k*(T-1)+t-1]);
	  } else {
	    Qinv_alphat_now = - 1 / eta2[k] * C_s_inv[k] * alpha_t[k*(T-1)+t];
	  }

	  if (t < (T-2)) {
	    Qinv_alphat_next =  omega[k] * 1/ eta2[k] / (1-omega2[k]) * C_s_inv[k] * (alpha_t[k*(T-1)+t+1] - omega[k] * alpha_t[k*(T-1)+t]);
	  }

	  for (int v=0; v<N_knots; ++v){
	    int idx_alpha_t = 1 + K*2 + K*(T-1) + K*N_knots + (k*(T-1) + t)*N_knots + v;
	    gradient[idx_alpha_t] += Qinv_alphat_now[v];
	    if (t < T-2) {
	      gradient[idx_alpha_t] += Qinv_alphat_next[v];
	    }
	  }
	}
      } // k

      // partial of dirmult
      Eigen::Tensor<double, 3> fac3(N,K,N_cores);

      #pragma omp parallel for collapse(2)      
      for (int i=0; i<N_cores; ++i) {
	for (int m=0; m<K; ++m) {
	  for (int c=0; c<N; ++c) {
	    const double drnew_case2 = (1-gamma[m]) * res * res / sum_w_pot[m];
	    fac3(c,m,i) = (idx_cores[i]-1 == c) ? gamma[m] : w[m](i,c) * drnew_case2;
	  } // c
	} // m
      } // i
      
      #pragma omp parallel for collapse(2)
      for (int k=0; k<K; ++k) {
	for (int t=0; t<T; ++t) {
	  const int idx1 = 1 + K*2 + K*(T-1) + K*N_knots + K*(T-1)*N_knots + k*N*T;
	  for (int i=0; i<N_cores; ++i) {
	    const int idx_core = idx_cores[i] - 1;
	    const int si = i*T + t;
	    if (y_row_sum[si] > 0.0) {
	      double dirmultp1 = -digamma(A[si] + N_grains[si]) + digamma(A[si]);
	      for (int m=0; m<K; ++m) {
		const double drnew_case2 = (1-gamma[m]) * res * res / sum_w_pot[m];
		const double dirmultp2 = digamma(y[si][m] + kappa(si,m)) - digamma(kappa(si,m));
		const double fac1 = (dirmultp1 + dirmultp2) * phi[m];

		if (m == k) {
		  for (int c=0; c<N; ++c) {
		    int C = c*T + t;
		    const double sumgp1 = sum_exp_g[C];
		    const double sumgp1inv2 = 1.0 / (sumgp1*sumgp1);
		    const double dr = exp_g(C,m) * (sumgp1 - exp_g(C,m)) * sumgp1inv2;
		    gradient[idx1 + C] += fac1 * fac3(c,m,i) * dr;
		  } // c
		} else {
		  for (int c=0; c<N; ++c) {
		    int C = c*T + t;
		    const double sumgp1 = sum_exp_g[C];
		    const double sumgp1inv2 = 1.0 / (sumgp1*sumgp1);
		    const double dr = -exp_g(C,m) * exp_g(C,k) * sumgp1inv2;
		    gradient[idx1 + C] += fac1 * fac3(c,m,i) * dr;
		  } // c
		} // if
	      } // m
	    } // if
	  } // i
	} // t
      } // k

      return lp;
    }

    template <bool propto, bool jacobian>
    double log_prob(vector<double>& params_r,
                    vector<int>& params_i,
                    std::ostream* msgs = 0) const {
      vector<double> gradient;
      return log_prob_grad<propto, jacobian, true>(params_r, gradient, msgs);
    }

    template <bool propto, bool jacobian, typename T__>
    double log_prob(vector<double>& params_r,
                    vector<int>& params_i,
                    std::ostream* msgs = 0) const {

      vector<double> gradient;
      return log_prob_grad<propto, jacobian, true>(params_r, gradient, msgs);
    }

    template <bool propto, bool jacobian, typename T__>
    double log_prob(Eigen::Matrix<T__,Eigen::Dynamic,1>& params_r,
                    std::ostream* msgs = 0) const {
      vector<T__> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      vector<double> gradient;
      return log_prob_grad<propto, jacobian, true>(vec_params_r, gradient, msgs);
    }

    // double log_prob(std::vector<stan::math::var>&, std::vector<double>&, std::ostream*&) const {

    // }

    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("ksi");
        names__.push_back("omega");
        // names__.push_back("lambda");
        names__.push_back("mu");
        names__.push_back("mu_t");
        names__.push_back("alpha_s");
        names__.push_back("alpha_t");
        names__.push_back("g");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dims__.push_back(T);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dims__.push_back(N_knots);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((K * T));
        dims__.push_back(N_knots);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dims__.push_back((N * T));
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "pred_model_namespace::write_array(%1%)";
        (void) function__; // dummy call to supress warning
        // read-transform, write parameters
        double ksi = in__.scalar_lub_constrain(ksi_lb,ksi_ub);
        vector_d omega = in__.vector_lub_constrain(omega_lb,omega_ub,K);
        // vector_d lambda = in__.vector_lub_constrain(lambda_lb,lambda_ub,K);
        vector_d mu = in__.vector_constrain(K);
        vector<vector_d> mu_t;
        size_t dim_mu_t_0__ = K;
        for (size_t k_0__ = 0; k_0__ < dim_mu_t_0__; ++k_0__) {
            mu_t.push_back(in__.vector_constrain(T-1));
        }
        vector<vector_d> alpha_s;
        size_t dim_alpha_s_0__ = K;
        for (size_t k_0__ = 0; k_0__ < dim_alpha_s_0__; ++k_0__) {
            alpha_s.push_back(in__.vector_constrain(N_knots));
        }
        vector<vector_d> alpha_t;
        size_t dim_alpha_t_0__ = (K * (T-1));
        for (size_t k_0__ = 0; k_0__ < dim_alpha_t_0__; ++k_0__) {
            alpha_t.push_back(in__.vector_constrain(N_knots));
        }
        vector<vector_d> g;
        size_t dim_g_0__ = K;
        for (size_t k_0__ = 0; k_0__ < dim_g_0__; ++k_0__) {
            g.push_back(in__.vector_constrain((N * T)));
        }
        vars__.push_back(ksi);
        for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(omega[k_0__]);
        }
        // for (int k_0__ = 0; k_0__ < K; ++k_0__) {
        //     vars__.push_back(lambda[k_0__]);
        // }
        for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(mu[k_0__]);
        }
        for (int k_1__ = 0; k_1__ < T-1; ++k_1__) {
            for (int k_0__ = 0; k_0__ < K; ++k_0__) {
                vars__.push_back(mu_t[k_0__][k_1__]);
            }
        }
        for (int k_1__ = 0; k_1__ < N_knots; ++k_1__) {
            for (int k_0__ = 0; k_0__ < K; ++k_0__) {
                vars__.push_back(alpha_s[k_0__][k_1__]);
            }
        }
        for (int k_1__ = 0; k_1__ < N_knots; ++k_1__) {
	  for (int k_0__ = 0; k_0__ < (K * (T-1)); ++k_0__) {
                vars__.push_back(alpha_t[k_0__][k_1__]);
            }
        }
        for (int k_1__ = 0; k_1__ < (N * T); ++k_1__) {
            for (int k_0__ = 0; k_0__ < K; ++k_0__) {
                vars__.push_back(g[k_0__][k_1__]);
            }
        }

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__; // dummy call to supress warning
        stan::math::accumulator<double> lp_accum__;



        // validate transformed parameters

        // write transformed parameters

        if (!include_gqs__) return;
        // declare and define generated quantities


        // validate generated quantities

        // write generated quantities
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "pred_model";
    }

    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "ksi";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "omega" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        // for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
        //     param_name_stream__.str(std::string());
        //     param_name_stream__ << "lambda" << '.' << k_0__;
        //     param_names__.push_back(param_name_stream__.str());
        // }
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_1__ = 1; k_1__ <= T-1; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu_t" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "alpha_s" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
	  for (int k_0__ = 1; k_0__ <= (K * (T-1)); ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "alpha_t" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= (N * T); ++k_1__) {
            for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "g" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "ksi";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "omega" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        // for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
        //     param_name_stream__.str(std::string());
        //     param_name_stream__ << "lambda" << '.' << k_0__;
        //     param_names__.push_back(param_name_stream__.str());
        // }
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_1__ = 1; k_1__ <= T-1; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu_t" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "alpha_s" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
	  for (int k_0__ = 1; k_0__ <= (K * (T-1)); ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "alpha_t" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_1__ = 1; k_1__ <= (N * T); ++k_1__) {
            for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "g" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }

}; // model

} // namespace

typedef pred_model_namespace::pred_model stan_model;

namespace stan {
  namespace model {

    using namespace std;

    template <bool propto, bool jacobian_adjust_transform>
    double log_prob_grad(const stan_model& model,
                         std::vector<double>& params_r,
                         std::vector<int>& params_i,
                         std::vector<double>& gradient,
                         std::ostream* msgs = 0) {
      gradient.resize(params_r.size());
      return model.log_prob_grad<propto, jacobian_adjust_transform, false>(params_r, gradient, msgs);
    }

    void gradient(const stan_model& model,
                  const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
                  double& f,
                  Eigen::Matrix<double, Eigen::Dynamic, 1>& grad_f,
		  std::ostream* msgs = 0) {
      vector<double> params_r(x.rows());
      vector<double> grad(x.rows());
      for (size_t i = 0; i < params_r.size(); i++)
        params_r[i] = x[i];
      paleon::Timer timer_lpg(1);
      timer_lpg.tic(0);
      f = model.log_prob_grad<true, true, false>(params_r, grad, msgs);
      timer_lpg.toc(0);
      timer_lpg.echo("=> gradient ");
      for (size_t i = 0; i < params_r.size(); i++)
        grad_f[i] = grad[i];
    }

    void gradient(const stan_model& model,
                  const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
                  double& f,
                  Eigen::Matrix<double, Eigen::Dynamic, 1>& grad_f,
		  callbacks::logger& logger) {
      gradient(model, x, f, grad_f);
    }
  }
}
