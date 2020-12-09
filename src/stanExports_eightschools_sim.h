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
namespace model_eightschools_sim_namespace {
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
    reader.add_event(0, 0, "start", "model_eightschools_sim");
    reader.add_event(46, 44, "end", "model_eightschools_sim");
    return reader;
}
#include <stan_meta_header.hpp>
class model_eightschools_sim
  : public stan::model::model_base_crtp<model_eightschools_sim> {
private:
        int J;
        int N;
        int H;
        std::vector<double> y;
        std::vector<double> sigma;
        std::vector<int> sigma_idx;
        std::vector<int> train_idx;
        std::vector<int> test_idx;
        std::vector<int> school_train;
        std::vector<int> school_test;
public:
    model_eightschools_sim(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_eightschools_sim(stan::io::var_context& context__,
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
        static const char* function__ = "model_eightschools_sim_namespace::model_eightschools_sim";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "J", "int", context__.to_vec());
            J = int(0);
            vals_i__ = context__.vals_i("J");
            pos__ = 0;
            J = vals_i__[pos__++];
            check_greater_or_equal(function__, "J", J, 1);
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 1);
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "H", "int", context__.to_vec());
            H = int(0);
            vals_i__ = context__.vals_i("H");
            pos__ = 0;
            H = vals_i__[pos__++];
            check_greater_or_equal(function__, "H", H, 0);
            check_less_or_equal(function__, "H", H, N);
            current_statement_begin__ = 8;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "double", context__.to_vec(N));
            y = std::vector<double>(N, double(0));
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                y[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("sigma", "J", J);
            context__.validate_dims("data initialization", "sigma", "double", context__.to_vec(J));
            sigma = std::vector<double>(J, double(0));
            vals_r__ = context__.vals_r("sigma");
            pos__ = 0;
            size_t sigma_k_0_max__ = J;
            for (size_t k_0__ = 0; k_0__ < sigma_k_0_max__; ++k_0__) {
                sigma[k_0__] = vals_r__[pos__++];
            }
            size_t sigma_i_0_max__ = J;
            for (size_t i_0__ = 0; i_0__ < sigma_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "sigma[i_0__]", sigma[i_0__], 0);
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("sigma_idx", "N", N);
            context__.validate_dims("data initialization", "sigma_idx", "int", context__.to_vec(N));
            sigma_idx = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("sigma_idx");
            pos__ = 0;
            size_t sigma_idx_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < sigma_idx_k_0_max__; ++k_0__) {
                sigma_idx[k_0__] = vals_i__[pos__++];
            }
            size_t sigma_idx_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < sigma_idx_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "sigma_idx[i_0__]", sigma_idx[i_0__], 0);
                check_less_or_equal(function__, "sigma_idx[i_0__]", sigma_idx[i_0__], J);
            }
            current_statement_begin__ = 13;
            validate_non_negative_index("train_idx", "(N - H)", (N - H));
            context__.validate_dims("data initialization", "train_idx", "int", context__.to_vec((N - H)));
            train_idx = std::vector<int>((N - H), int(0));
            vals_i__ = context__.vals_i("train_idx");
            pos__ = 0;
            size_t train_idx_k_0_max__ = (N - H);
            for (size_t k_0__ = 0; k_0__ < train_idx_k_0_max__; ++k_0__) {
                train_idx[k_0__] = vals_i__[pos__++];
            }
            size_t train_idx_i_0_max__ = (N - H);
            for (size_t i_0__ = 0; i_0__ < train_idx_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "train_idx[i_0__]", train_idx[i_0__], 1);
                check_less_or_equal(function__, "train_idx[i_0__]", train_idx[i_0__], N);
            }
            current_statement_begin__ = 14;
            validate_non_negative_index("test_idx", "H", H);
            context__.validate_dims("data initialization", "test_idx", "int", context__.to_vec(H));
            test_idx = std::vector<int>(H, int(0));
            vals_i__ = context__.vals_i("test_idx");
            pos__ = 0;
            size_t test_idx_k_0_max__ = H;
            for (size_t k_0__ = 0; k_0__ < test_idx_k_0_max__; ++k_0__) {
                test_idx[k_0__] = vals_i__[pos__++];
            }
            size_t test_idx_i_0_max__ = H;
            for (size_t i_0__ = 0; i_0__ < test_idx_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "test_idx[i_0__]", test_idx[i_0__], 1);
                check_less_or_equal(function__, "test_idx[i_0__]", test_idx[i_0__], N);
            }
            // initialize transformed data variables
            current_statement_begin__ = 18;
            validate_non_negative_index("school_train", "(N - H)", (N - H));
            school_train = std::vector<int>((N - H), int(0));
            stan::math::fill(school_train, std::numeric_limits<int>::min());
            stan::math::assign(school_train,stan::model::rvalue(sigma_idx, stan::model::cons_list(stan::model::index_multi(train_idx), stan::model::nil_index_list()), "sigma_idx"));
            current_statement_begin__ = 19;
            validate_non_negative_index("school_test", "H", H);
            school_test = std::vector<int>(H, int(0));
            stan::math::fill(school_test, std::numeric_limits<int>::min());
            stan::math::assign(school_test,stan::model::rvalue(sigma_idx, stan::model::cons_list(stan::model::index_multi(test_idx), stan::model::nil_index_list()), "sigma_idx"));
            // execute transformed data statements
            // validate transformed data
            current_statement_begin__ = 18;
            size_t school_train_i_0_max__ = (N - H);
            for (size_t i_0__ = 0; i_0__ < school_train_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "school_train[i_0__]", school_train[i_0__], 1);
                check_less_or_equal(function__, "school_train[i_0__]", school_train[i_0__], J);
            }
            current_statement_begin__ = 19;
            size_t school_test_i_0_max__ = H;
            for (size_t i_0__ = 0; i_0__ < school_test_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "school_test[i_0__]", school_test[i_0__], 1);
                check_less_or_equal(function__, "school_test[i_0__]", school_test[i_0__], J);
            }
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 24;
            num_params_r__ += 1;
            current_statement_begin__ = 25;
            num_params_r__ += 1;
            current_statement_begin__ = 26;
            validate_non_negative_index("eta", "J", J);
            num_params_r__ += J;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_eightschools_sim() { }
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
        current_statement_begin__ = 24;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu", "double", context__.to_vec());
        double mu(0);
        mu = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 25;
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
        current_statement_begin__ = 26;
        if (!(context__.contains_r("eta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable eta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("eta");
        pos__ = 0U;
        validate_non_negative_index("eta", "J", J);
        context__.validate_dims("parameter initialization", "eta", "vector_d", context__.to_vec(J));
        Eigen::Matrix<double, Eigen::Dynamic, 1> eta(J);
        size_t eta_j_1_max__ = J;
        for (size_t j_1__ = 0; j_1__ < eta_j_1_max__; ++j_1__) {
            eta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(eta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable eta: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 24;
            local_scalar_t__ mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.scalar_constrain(lp__);
            else
                mu = in__.scalar_constrain();
            current_statement_begin__ = 25;
            local_scalar_t__ tau;
            (void) tau;  // dummy to suppress unused var warning
            if (jacobian__)
                tau = in__.scalar_lb_constrain(0, lp__);
            else
                tau = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 26;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta;
            (void) eta;  // dummy to suppress unused var warning
            if (jacobian__)
                eta = in__.vector_constrain(J, lp__);
            else
                eta = in__.vector_constrain(J);
            // transformed parameters
            current_statement_begin__ = 29;
            validate_non_negative_index("theta", "J", J);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> theta(J);
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            stan::math::assign(theta,multiply(tau, eta));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 29;
            size_t theta_j_1_max__ = J;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(theta(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: theta" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable theta: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 33;
            lp_accum__.add(normal_log<propto__>(eta, 0, 1));
            current_statement_begin__ = 34;
            lp_accum__.add(normal_log<propto__>(stan::model::rvalue(y, stan::model::cons_list(stan::model::index_multi(train_idx), stan::model::nil_index_list()), "y"), add(mu, stan::model::rvalue(theta, stan::model::cons_list(stan::model::index_multi(school_train), stan::model::nil_index_list()), "theta")), stan::model::rvalue(sigma, stan::model::cons_list(stan::model::index_multi(school_train), stan::model::nil_index_list()), "sigma")));
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
        names__.push_back("mu");
        names__.push_back("tau");
        names__.push_back("eta");
        names__.push_back("theta");
        names__.push_back("log_lik");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(J);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(J);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
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
        static const char* function__ = "model_eightschools_sim_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double mu = in__.scalar_constrain();
        vars__.push_back(mu);
        double tau = in__.scalar_lb_constrain(0);
        vars__.push_back(tau);
        Eigen::Matrix<double, Eigen::Dynamic, 1> eta = in__.vector_constrain(J);
        size_t eta_j_1_max__ = J;
        for (size_t j_1__ = 0; j_1__ < eta_j_1_max__; ++j_1__) {
            vars__.push_back(eta(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 29;
            validate_non_negative_index("theta", "J", J);
            Eigen::Matrix<double, Eigen::Dynamic, 1> theta(J);
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            stan::math::assign(theta,multiply(tau, eta));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t theta_j_1_max__ = J;
                for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                    vars__.push_back(theta(j_1__));
                }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 38;
            validate_non_negative_index("log_lik", "N", N);
            Eigen::Matrix<double, Eigen::Dynamic, 1> log_lik(N);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 40;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 41;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            normal_log(get_base1(y, n, "y", 1), (mu + get_base1(theta, get_base1(sigma_idx, n, "sigma_idx", 1), "theta", 1)), get_base1(sigma, get_base1(sigma_idx, n, "sigma_idx", 1), "sigma", 1)), 
                            "assigning variable log_lik");
            }
            // validate, write generated quantities
            current_statement_begin__ = 38;
            size_t log_lik_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                vars__.push_back(log_lik(j_1__));
            }
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
        return "model_eightschools_sim";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        size_t eta_j_1_max__ = J;
        for (size_t j_1__ = 0; j_1__ < eta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "eta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t theta_j_1_max__ = J;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "theta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        size_t log_lik_j_1_max__ = N;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        size_t eta_j_1_max__ = J;
        for (size_t j_1__ = 0; j_1__ < eta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "eta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t theta_j_1_max__ = J;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "theta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        size_t log_lik_j_1_max__ = N;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_eightschools_sim_namespace::model_eightschools_sim stan_model;
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