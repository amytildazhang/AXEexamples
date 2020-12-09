// Generated by rstantools.  Do not edit by hand.

/*
    AXEexamples is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AXEexamples is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AXEexamples.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_carslp_cv_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_carslp_cv");
    reader.add_event(50, 48, "end", "model_carslp_cv");
    return reader;
}
#include <stan_meta_header.hpp>
class model_carslp_cv
  : public stan::model::model_base_crtp<model_carslp_cv> {
private:
        int n;
        int p;
        matrix_d X;
        std::vector<int> y_train;
        vector_d log_offset;
        matrix_d W;
        int test_idx;
        vector_d zeros;
        matrix_d D;
public:
    model_carslp_cv(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_carslp_cv(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_carslp_cv_namespace::model_carslp_cv";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            check_greater_or_equal(function__, "n", n, 1);
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "p", "int", context__.to_vec());
            p = int(0);
            vals_i__ = context__.vals_i("p");
            pos__ = 0;
            p = vals_i__[pos__++];
            check_greater_or_equal(function__, "p", p, 1);
            current_statement_begin__ = 6;
            validate_non_negative_index("X", "n", n);
            validate_non_negative_index("X", "p", p);
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(n,p));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, p);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = p;
            size_t X_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("y_train", "(n - 1)", (n - 1));
            context__.validate_dims("data initialization", "y_train", "int", context__.to_vec((n - 1)));
            y_train = std::vector<int>((n - 1), int(0));
            vals_i__ = context__.vals_i("y_train");
            pos__ = 0;
            size_t y_train_k_0_max__ = (n - 1);
            for (size_t k_0__ = 0; k_0__ < y_train_k_0_max__; ++k_0__) {
                y_train[k_0__] = vals_i__[pos__++];
            }
            size_t y_train_i_0_max__ = (n - 1);
            for (size_t i_0__ = 0; i_0__ < y_train_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "y_train[i_0__]", y_train[i_0__], 0);
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("log_offset", "n", n);
            context__.validate_dims("data initialization", "log_offset", "vector_d", context__.to_vec(n));
            log_offset = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("log_offset");
            pos__ = 0;
            size_t log_offset_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < log_offset_j_1_max__; ++j_1__) {
                log_offset(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("W", "n", n);
            validate_non_negative_index("W", "n", n);
            context__.validate_dims("data initialization", "W", "matrix_d", context__.to_vec(n,n));
            W = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, n);
            vals_r__ = context__.vals_r("W");
            pos__ = 0;
            size_t W_j_2_max__ = n;
            size_t W_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < W_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < W_j_1_max__; ++j_1__) {
                    W(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            check_greater_or_equal(function__, "W", W, 0);
            check_less_or_equal(function__, "W", W, 1);
            current_statement_begin__ = 11;
            context__.validate_dims("data initialization", "test_idx", "int", context__.to_vec());
            test_idx = int(0);
            vals_i__ = context__.vals_i("test_idx");
            pos__ = 0;
            test_idx = vals_i__[pos__++];
            check_greater_or_equal(function__, "test_idx", test_idx, 1);
            check_less_or_equal(function__, "test_idx", test_idx, n);
            // initialize transformed data variables
            current_statement_begin__ = 14;
            validate_non_negative_index("zeros", "n", n);
            zeros = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            stan::math::fill(zeros, DUMMY_VAR__);
            current_statement_begin__ = 15;
            validate_non_negative_index("D", "n", n);
            validate_non_negative_index("D", "n", n);
            D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, n);
            stan::math::fill(D, DUMMY_VAR__);
            // execute transformed data statements
            {
            current_statement_begin__ = 17;
            validate_non_negative_index("W_rowsums", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> W_rowsums(n);
            stan::math::initialize(W_rowsums, DUMMY_VAR__);
            stan::math::fill(W_rowsums, DUMMY_VAR__);
            current_statement_begin__ = 18;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 19;
                stan::model::assign(W_rowsums, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            sum(stan::model::rvalue(W, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "W")), 
                            "assigning variable W_rowsums");
            }
            current_statement_begin__ = 21;
            stan::math::assign(D, diag_matrix(W_rowsums));
            }
            current_statement_begin__ = 23;
            stan::math::assign(zeros, rep_vector(0, n));
            // validate transformed data
            current_statement_begin__ = 15;
            check_greater_or_equal(function__, "D", D, 0);
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 26;
            validate_non_negative_index("beta", "p", p);
            num_params_r__ += p;
            current_statement_begin__ = 27;
            validate_non_negative_index("phi", "n", n);
            num_params_r__ += n;
            current_statement_begin__ = 28;
            num_params_r__ += 1;
            current_statement_begin__ = 29;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_carslp_cv() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 26;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "p", p);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(p);
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 27;
        if (!(context__.contains_r("phi")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable phi missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("phi");
        pos__ = 0U;
        validate_non_negative_index("phi", "n", n);
        context__.validate_dims("parameter initialization", "phi", "vector_d", context__.to_vec(n));
        Eigen::Matrix<double, Eigen::Dynamic, 1> phi(n);
        size_t phi_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < phi_j_1_max__; ++j_1__) {
            phi(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(phi);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable phi: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 28;
        if (!(context__.contains_r("tau")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable tau missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("tau");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "tau", "double", context__.to_vec());
        double tau(0);
        tau = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, tau);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable tau: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 29;
        if (!(context__.contains_r("alpha")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable alpha missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("alpha");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "alpha", "double", context__.to_vec());
        double alpha(0);
        alpha = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0, 1, alpha);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable alpha: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 26;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(p, lp__);
            else
                beta = in__.vector_constrain(p);
            current_statement_begin__ = 27;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> phi;
            (void) phi;  // dummy to suppress unused var warning
            if (jacobian__)
                phi = in__.vector_constrain(n, lp__);
            else
                phi = in__.vector_constrain(n);
            current_statement_begin__ = 28;
            local_scalar_t__ tau;
            (void) tau;  // dummy to suppress unused var warning
            if (jacobian__)
                tau = in__.scalar_lb_constrain(0, lp__);
            else
                tau = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 29;
            local_scalar_t__ alpha;
            (void) alpha;  // dummy to suppress unused var warning
            if (jacobian__)
                alpha = in__.scalar_lub_constrain(0, 1, lp__);
            else
                alpha = in__.scalar_lub_constrain(0, 1);
            // transformed parameters
            current_statement_begin__ = 33;
            validate_non_negative_index("mu", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu(n);
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 34;
            validate_non_negative_index("mu_train", "(n - 1)", (n - 1));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu_train((n - 1));
            stan::math::initialize(mu_train, DUMMY_VAR__);
            stan::math::fill(mu_train, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 36;
            stan::math::assign(mu, add(add(multiply(X, beta), phi), log_offset));
            current_statement_begin__ = 37;
            stan::math::assign(mu_train, stan::model::rvalue(mu, stan::model::cons_list(stan::model::index_min_max(1, (n - 1)), stan::model::nil_index_list()), "mu"));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 33;
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(mu(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: mu" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mu: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 34;
            size_t mu_train_j_1_max__ = (n - 1);
            for (size_t j_1__ = 0; j_1__ < mu_train_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(mu_train(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: mu_train" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mu_train: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 40;
            lp_accum__.add(multi_normal_prec_log<propto__>(phi, zeros, multiply(tau, subtract(D, multiply(alpha, W)))));
            current_statement_begin__ = 41;
            lp_accum__.add(normal_log<propto__>(beta, 0, 1));
            current_statement_begin__ = 42;
            lp_accum__.add(gamma_log<propto__>(tau, 2, 2));
            current_statement_begin__ = 43;
            lp_accum__.add(poisson_log_log<propto__>(y_train, mu_train));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("beta");
        names__.push_back("phi");
        names__.push_back("tau");
        names__.push_back("alpha");
        names__.push_back("mu");
        names__.push_back("mu_train");
        names__.push_back("y_test");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((n - 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
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
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_carslp_cv_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(p);
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> phi = in__.vector_constrain(n);
        size_t phi_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < phi_j_1_max__; ++j_1__) {
            vars__.push_back(phi(j_1__));
        }
        double tau = in__.scalar_lb_constrain(0);
        vars__.push_back(tau);
        double alpha = in__.scalar_lub_constrain(0, 1);
        vars__.push_back(alpha);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 33;
            validate_non_negative_index("mu", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> mu(n);
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 34;
            validate_non_negative_index("mu_train", "(n - 1)", (n - 1));
            Eigen::Matrix<double, Eigen::Dynamic, 1> mu_train((n - 1));
            stan::math::initialize(mu_train, DUMMY_VAR__);
            stan::math::fill(mu_train, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 36;
            stan::math::assign(mu, add(add(multiply(X, beta), phi), log_offset));
            current_statement_begin__ = 37;
            stan::math::assign(mu_train, stan::model::rvalue(mu, stan::model::cons_list(stan::model::index_min_max(1, (n - 1)), stan::model::nil_index_list()), "mu"));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t mu_j_1_max__ = n;
                for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                    vars__.push_back(mu(j_1__));
                }
                size_t mu_train_j_1_max__ = (n - 1);
                for (size_t j_1__ = 0; j_1__ < mu_train_j_1_max__; ++j_1__) {
                    vars__.push_back(mu_train(j_1__));
                }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 47;
            int y_test;
            (void) y_test;  // dummy to suppress unused var warning
            stan::math::fill(y_test, std::numeric_limits<int>::min());
            stan::math::assign(y_test,poisson_log_rng(get_base1(mu, n, "mu", 1), base_rng__));
            // validate, write generated quantities
            current_statement_begin__ = 47;
            vars__.push_back(y_test);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
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
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_carslp_cv";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t phi_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < phi_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "phi" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t mu_train_j_1_max__ = (n - 1);
            for (size_t j_1__ = 0; j_1__ < mu_train_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu_train" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "y_test";
        param_names__.push_back(param_name_stream__.str());
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t phi_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < phi_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "phi" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t mu_train_j_1_max__ = (n - 1);
            for (size_t j_1__ = 0; j_1__ < mu_train_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu_train" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "y_test";
        param_names__.push_back(param_name_stream__.str());
    }
}; // model
}  // namespace
typedef model_carslp_cv_namespace::model_carslp_cv stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
