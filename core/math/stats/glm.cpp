/* Copyright (c) 2008-2023 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Covered Software is provided under this License on an "as is"
 * basis, without warranty of any kind, either expressed, implied, or
 * statutory, including, without limitation, warranties that the
 * Covered Software is free of defects, merchantable, fit for a
 * particular purpose or non-infringing.
 * See the Mozilla Public License v. 2.0 for more details.
 *
 * For more details, see http://www.mrtrix.org/.
 */

#include "math/stats/glm.h"

#include "debug.h"
#include "thread_queue.h"
#include "file/matrix.h"
#include "math/betainc.h"
#include "math/erfinv.h"
#include "math/welch_satterthwaite.h"
#include "misc/bitset.h"

//#define GLM_ALL_STATS_DEBUG
//#define GLM_TEST_DEBUG

namespace MR
{
  namespace Math
  {
    namespace Stats
    {
      namespace GLM
      {



        const char* const column_ones_description =
            "In some software packages, a column of ones is automatically added to the "
            "GLM design matrix; the purpose of this column is to estimate the \"global "
            "intercept\", which is the predicted value of the observed variable if all "
            "explanatory variables were to be zero. However there are rare situations "
            "where including such a column would not be appropriate for a particular "
            "experimental design. Hence, in MRtrix3 statistical inference commands, "
            "it is up to the user to determine whether or not this column of ones should "
            "be included in their design matrix, and add it explicitly if necessary. "
            "The contrast matrix must also reflect the presence of this additional column.";


        App::OptionGroup glm_options (const std::string& element_name)
        {
          using namespace App;
          OptionGroup result = OptionGroup ("Options related to the General Linear Model (GLM)")

            + Option ("variance", "define variance groups for the G-statistic; "
                                  "measurements for which the expected variance is equivalent should contain the same index")
              + Argument ("file").type_file_in()

            + Option ("ftests", "perform F-tests; input text file should contain, for each F-test, a row containing "
                                "ones and zeros, where ones indicate the rows of the contrast matrix to be included "
                                "in the F-test.")
              + Argument ("path").type_file_in()

            + Option ("fonly", "only assess F-tests; do not perform statistical inference on entries in the contrast matrix")

            + Option ("column", "add a column to the design matrix corresponding to subject " + element_name + "-wise values "
                                "(note that the contrast matrix must include an additional column for each use of this option); "
                                "the text file provided via this option should contain a file name for each subject").allow_multiple()
              + Argument ("path").type_file_in();

          return result;
        }




        void check_design (const matrix_type& design, const bool extra_factors)
        {
          Eigen::ColPivHouseholderQR<matrix_type> decomp;
          decomp.setThreshold (1e-5);
          decomp = decomp.compute (design);
          if (decomp.rank() < design.cols()) {
            if (extra_factors) {
              CONSOLE ("Design matrix is rank-deficient before addition of element-wise columns");
            } else {
              WARN ("Design matrix is rank-deficient; processing may proceed, but manually checking your matrix is advised");
            }
          } else {
            const default_type cond = Math::condition_number (design);
            if (cond > 100.0) {
              if (extra_factors) {
                CONSOLE ("Design matrix conditioning is poor (condition number: " + str(cond, 6) + ") before the addition of element-wise columns");
              } else {
                WARN ("Design matrix conditioning is poor (condition number: " + str(cond, 6) + "); model fitting may be highly influenced by noise");
              }
            } else {
              CONSOLE (std::string ("Design matrix condition number") + (extra_factors ? " (without element-wise columns)" : "") + ": " + str(cond, 6));
            }
          }
        }



        index_array_type load_variance_groups (const index_type num_inputs)
        {
          auto opt = App::get_options ("variance");
          if (!opt.size())
            return index_array_type();
          try {
            auto data = File::Matrix::load_vector<index_type> (opt[0][0]);
            if (index_type(data.size()) != num_inputs)
              throw Exception ("Number of entries in variance group file \"" + std::string(opt[0][0]) + "\" (" + str(data.size()) + ") does not match number of inputs (" + str(num_inputs) + ")");
            const index_type min_coeff = data.minCoeff();
            const index_type max_coeff = data.maxCoeff();
            if (min_coeff > 1)
              throw Exception ("Minimum coefficient needs to be either zero or one");
            if (max_coeff == min_coeff) {
              WARN ("Only a single variance group is defined in file \"" + opt[0][0] + "\"; variance groups will not be used");
              return index_array_type();
            }
            vector<index_type> count_per_group (max_coeff + 1, 0);
            for (index_type i = 0; i != index_type(data.size()); ++i)
              count_per_group[data[i]]++;
            for (index_type vg_index = min_coeff; vg_index <= max_coeff; ++vg_index) {
              if (!count_per_group[vg_index])
                throw Exception ("No entries found for variance group " + str(vg_index));
            }
            if (min_coeff)
              data.array() -= 1;
            return data.array();
          } catch (Exception& e) {
            throw Exception (e, "unable to read file \"" + opt[0][0] + "\" as variance group data");
          }
        }




        vector<Hypothesis> load_hypotheses (const std::string& file_path)
        {
          vector<Hypothesis> hypotheses;
          const matrix_type contrast_matrix = File::Matrix::load_matrix (file_path);
          for (index_type row = 0; row != index_type(contrast_matrix.rows()); ++row)
            hypotheses.emplace_back (Hypothesis (contrast_matrix.row (row), row));
          auto opt = App::get_options ("ftests");
          if (opt.size()) {
            const matrix_type ftest_matrix = File::Matrix::load_matrix (opt[0][0]);
            if (ftest_matrix.cols() != contrast_matrix.rows())
              throw Exception ("Number of columns in F-test matrix (" + str(ftest_matrix.cols()) + ") does not match number of rows in contrast matrix (" + str(contrast_matrix.rows()) + ")");
            if (!((ftest_matrix.array() == 0.0) + (ftest_matrix.array() == 1.0)).all())
              throw Exception ("F-test array must contain ones and zeros only");
            for (index_type ftest_index = 0; ftest_index != index_type(ftest_matrix.rows()); ++ftest_index) {
              if (!ftest_matrix.row (ftest_index).count())
                throw Exception ("Row " + str(ftest_index+1) + " of F-test matrix does not contain any ones");
              matrix_type this_f_matrix (ftest_matrix.row (ftest_index).count(), contrast_matrix.cols());
              index_type ftest_row = 0;
              for (index_type contrast_row = 0; contrast_row != index_type(contrast_matrix.rows()); ++contrast_row) {
                if (ftest_matrix (ftest_index, contrast_row))
                  this_f_matrix.row (ftest_row++) = contrast_matrix.row (contrast_row);
              }
              hypotheses.emplace_back (Hypothesis (this_f_matrix, ftest_index));
            }
            if (App::get_options ("fonly").size()) {
              vector<Hypothesis> new_hypotheses;
              for (index_type index = contrast_matrix.rows(); index != hypotheses.size(); ++index)
                new_hypotheses.push_back (std::move (hypotheses[index]));
              std::swap (hypotheses, new_hypotheses);
            }
          } else if (App::get_options ("fonly").size()) {
            throw Exception ("Cannot perform F-tests exclusively (-fonly option): No F-test matrix was provided (-ftests option)");
          }
          return hypotheses;
        }






        matrix_type solve_betas (const matrix_type& measurements, const matrix_type& design)
        {
          return design.jacobiSvd (Eigen::ComputeThinU | Eigen::ComputeThinV).solve (measurements);
        }



        vector_type abs_effect_size (const matrix_type& measurements, const matrix_type& design, const Hypothesis& hypothesis)
        {
          if (hypothesis.is_F())
            return vector_type::Constant (measurements.rows(), std::numeric_limits<default_type>::quiet_NaN());
          else
            return hypothesis.matrix() * solve_betas (measurements, design);
        }

        matrix_type abs_effect_size (const matrix_type& measurements, const matrix_type& design, const vector<Hypothesis>& hypotheses)
        {
          matrix_type result (measurements.cols(), hypotheses.size());
          for (index_type ic = 0; ic != hypotheses.size(); ++ic)
            result.col (ic) = abs_effect_size (measurements, design, hypotheses[ic]);
          return result;
        }



        matrix_type stdev (const matrix_type& measurements, const matrix_type& design)
        {
          const matrix_type residuals = measurements - design * solve_betas (measurements, design);
          const matrix_type sse = residuals.colwise().squaredNorm();
          return (sse.array() / value_type(design.rows()-Math::rank (design))).sqrt();
        }


        matrix_type stdev (const matrix_type& measurements, const matrix_type& design, const index_array_type& variance_groups)
        {
          assert (measurements.rows() == design.rows());
          if (!variance_groups.size())
            return stdev (measurements, design);
          assert (measurements.rows() == variance_groups.rows());
          // Residual-forming matrix
          const matrix_type R (matrix_type::Identity (design.rows(), design.rows()) - (design * Math::pinv (design)));
          // Residuals
          const matrix_type e (R * measurements);
          // One standard deviation per element per variance group
          // Rows are variance groups, columns are elements
          const index_type num_vgs = variance_groups.array().maxCoeff() + 1;
          // Sum of residual-forming matrix diagonal elements within each variance group
          //   will be equivalent across elements
          vector_type Rnn_sums (vector_type::Zero (num_vgs));
          for (index_type i = 0; i != index_type(measurements.rows()); ++i)
            Rnn_sums[variance_groups[i]] += R.diagonal()[i];
          // For each variance group, get the sum of squared residuals within that group
          matrix_type result (num_vgs, measurements.cols());
          for (index_type ie = 0; ie != index_type(measurements.cols()); ++ie) {
            vector_type sse (vector_type::Zero (num_vgs));
            for (index_type i = 0; i != index_type(measurements.rows()); ++i)
              sse[variance_groups[i]] += Math::pow2 (e (i, ie));
            // (Rnn_sum / sse) is the inverse of the estimated variance
            result.col (ie) = (sse.array() / Rnn_sums.array()).sqrt();
          }
          return result;
        }



        vector_type std_effect_size (const matrix_type& measurements, const matrix_type& design, const Hypothesis& hypothesis)
        {
          if (hypothesis.is_F())
            return vector_type::Constant (measurements.cols(), std::numeric_limits<default_type>::quiet_NaN());
          return abs_effect_size (measurements, design, hypothesis).array() / stdev (measurements, design).array().col(0);
        }

        matrix_type std_effect_size (const matrix_type& measurements, const matrix_type& design, const vector<Hypothesis>& hypotheses)
        {
          const vector_type stdev_reciprocal = vector_type::Ones (measurements.cols()) / stdev (measurements, design).array().col(0);
          matrix_type result (measurements.cols(), hypotheses.size());
          for (index_type ic = 0; ic != hypotheses.size(); ++ic)
            result.col (ic) = abs_effect_size (measurements, design, hypotheses[ic]) * stdev_reciprocal;
          return result;
        }






        void all_stats (const matrix_type& measurements,
                        const matrix_type& design,
                        const vector<Hypothesis>& hypotheses,
                        const index_array_type& variance_groups,
                        matrix_type& betas,
                        matrix_type& abs_effect_size,
                        matrix_type& std_effect_size,
                        matrix_type& stdev)
        {
#ifndef GLM_ALL_STATS_DEBUG
          // If this function is being invoked from the other version of all_stats(),
          //   on an element-by-element basis, don't interfere with the progress bar
          //   that's being displayed by that outer looping function
          std::unique_ptr<ProgressBar> progress;
          if (measurements.cols() > 1)
            progress.reset (new ProgressBar ("Calculating basic properties of default permutation", 5));
#endif
          betas = solve_betas (measurements, design);
#ifdef GLM_ALL_STATS_DEBUG
          std::cerr << "Betas: " << betas.rows() << " x " << betas.cols() << ", max " << betas.array().maxCoeff() << "\n";
#else
          if (progress)
            ++*progress;
#endif
          abs_effect_size.resize (measurements.cols(), hypotheses.size());
          for (index_type ic = 0; ic != hypotheses.size(); ++ic) {
            if (hypotheses[ic].is_F()) {
              abs_effect_size.col (ic).fill (std::numeric_limits<default_type>::quiet_NaN());
            } else {
              abs_effect_size.col (ic) = (hypotheses[ic].matrix() * betas).row (0);
            }
          }
#ifdef GLM_ALL_STATS_DEBUG
          std::cerr << "abs_effect_size: " << abs_effect_size.rows() << " x " << abs_effect_size.cols() << ", max " << abs_effect_size.array().maxCoeff() << "\n";
#else
          if (progress)
            ++*progress;
#endif
          // Explicit calculation of residuals before SSE, rather than in a single
          //   step, appears to be necessary for compatibility with Eigen 3.2.0
          const matrix_type residuals = (measurements - design * betas);
#ifdef GLM_ALL_STATS_DEBUG
          std::cerr << "Residuals: " << residuals.rows() << " x " << residuals.cols() << ", max " << residuals.array().maxCoeff() << "\n";
#else
          if (progress)
            ++*progress;
#endif
          stdev = GLM::stdev (measurements, design, variance_groups);
#ifdef GLM_ALL_STATS_DEBUG
          std::cerr << "stdev: " << stdev.rows() << " x " << stdev.cols() << ", max " << stdev.maxCoeff() << "\n";
#else
          if (progress)
            ++*progress;
#endif
          if (variance_groups.size())
            std_effect_size = matrix_type::Constant (measurements.cols(), hypotheses.size(), std::numeric_limits<default_type>::quiet_NaN());
          else
            std_effect_size = abs_effect_size.array().colwise() / stdev.transpose().array().col(0);
#ifdef GLM_ALL_STATS_DEBUG
          std::cerr << "std_effect_size: " << std_effect_size.rows() << " x " << std_effect_size.cols() << ", max " << std_effect_size.array().maxCoeff() << "\n";
#endif
        }



        void all_stats (const matrix_type& measurements,
                        const matrix_type& fixed_design,
                        const vector<CohortDataImport>& extra_data,
                        const vector<Hypothesis>& hypotheses,
                        const index_array_type& variance_groups,
                        vector_type& cond,
                        matrix_type& betas,
                        matrix_type& abs_effect_size,
                        matrix_type& std_effect_size,
                        matrix_type& stdev)
        {
          if (extra_data.empty() && measurements.allFinite()) {
            all_stats (measurements, fixed_design, hypotheses, variance_groups, betas, abs_effect_size, std_effect_size, stdev);
            return;
          }

          class Source
          {
            public:
              Source (const index_type num_elements) :
                  num_elements (num_elements),
                  counter (0),
                  progress (new ProgressBar ("Calculating basic properties of default permutation", num_elements)) { }
              bool operator() (index_type& element_index)
              {
                element_index = counter++;
                if (element_index >= num_elements) {
                  progress.reset();
                  return false;
                }
                assert (progress);
                ++(*progress);
                return true;
              }
            private:
              const index_type num_elements;
              index_type counter;
              std::unique_ptr<ProgressBar> progress;
          };

          class Functor
          {
            public:
              Functor (const matrix_type& data, const matrix_type& design_fixed, const vector<CohortDataImport>& extra_data, const vector<Hypothesis>& hypotheses, const index_array_type& variance_groups,
                       vector_type& cond, matrix_type& betas, matrix_type& abs_effect_size, matrix_type& std_effect_size, matrix_type& stdev) :
                  data (data),
                  design_fixed (design_fixed),
                  extra_data (extra_data),
                  hypotheses (hypotheses),
                  variance_groups (variance_groups),
                  global_cond (cond),
                  global_betas (betas),
                  global_abs_effect_size (abs_effect_size),
                  global_std_effect_size (std_effect_size),
                  global_stdev (stdev),
                  num_vgs (variance_groups.size() ? variance_groups.maxCoeff()+1 : 1)
              {
                assert (index_type(design_fixed.cols()) + extra_data.size() == index_type(hypotheses[0].cols()));
              }
              bool operator() (const index_type& element_index)
              {
                const matrix_type element_data = data.col (element_index);
                matrix_type element_design (design_fixed.rows(), design_fixed.cols() + extra_data.size());
                element_design.leftCols (design_fixed.cols()) = design_fixed;
                // For each element-wise design matrix column,
                //   acquire the data for this particular element, without permutation
                for (index_type col = 0; col != extra_data.size(); ++col)
                  element_design.col (design_fixed.cols() + col) = (extra_data[col]) (element_index);
                // For each element-wise design matrix, remove any NaN values
                //   present in either the input data or imported from the element-wise design matrix column data
                index_type valid_rows = 0;
                for (index_type row = 0; row != index_type(data.rows()); ++row) {
                  if (std::isfinite (element_data(row)) && element_design.row (row).allFinite())
                    ++valid_rows;
                }
                default_type condition_number = 0.0;
                if (valid_rows == data.rows()) { // No NaNs present
                  condition_number = Math::condition_number (element_design);
                  if (!std::isfinite (condition_number) || condition_number > 1e5) {
                    zero();
                  } else {
                    Math::Stats::GLM::all_stats (element_data, element_design, hypotheses, variance_groups,
                                                 local_betas, local_abs_effect_size, local_std_effect_size, local_stdev);
                  }
                } else if (valid_rows >= element_design.cols()) {
                  // Need to reduce the data and design matrices to contain only finite data
                  matrix_type element_data_finite (valid_rows, 1);
                  matrix_type element_design_finite (valid_rows, element_design.cols());
                  index_array_type variance_groups_finite (variance_groups.size() ? valid_rows : 0);
                  index_type output_row = 0;
                  for (index_type row = 0; row != index_type(data.rows()); ++row) {
                    if (std::isfinite (element_data(row)) && element_design.row (row).allFinite()) {
                      element_data_finite(output_row, 0) = element_data(row);
                      element_design_finite.row (output_row) = element_design.row (row);
                      if (variance_groups.size())
                        variance_groups_finite[output_row] = variance_groups[row];
                      ++output_row;
                    }
                  }
                  assert (output_row == valid_rows);
                  assert (element_data_finite.allFinite());
                  assert (element_design_finite.allFinite());
                  condition_number = Math::condition_number (element_design_finite);
                  if (!std::isfinite (condition_number) || condition_number > 1e5) {
                    zero();
                  } else {
                    Math::Stats::GLM::all_stats (element_data_finite, element_design_finite, hypotheses, variance_groups_finite,
                                                 local_betas, local_abs_effect_size, local_std_effect_size, local_stdev);
                  }
                } else { // Insufficient data to fit model at all
                  zero();
                }
                global_cond[element_index] = condition_number;
                global_betas.col (element_index) = local_betas;
                global_abs_effect_size.row (element_index) = local_abs_effect_size.row (0);
                global_std_effect_size.row (element_index) = local_std_effect_size.row (0);
                global_stdev.col (element_index) = local_stdev;
                return true;
              }
            private:
              const matrix_type& data;
              const matrix_type& design_fixed;
              const vector<CohortDataImport>& extra_data;
              const vector<Hypothesis>& hypotheses;
              const index_array_type& variance_groups;
              vector_type& global_cond;
              matrix_type& global_betas;
              matrix_type& global_abs_effect_size;
              matrix_type& global_std_effect_size;
              matrix_type& global_stdev;
              matrix_type local_betas, local_abs_effect_size, local_std_effect_size, local_stdev;
              const index_type num_vgs;

              void zero () {
                local_betas = matrix_type::Zero (global_betas.rows(), 1);
                local_abs_effect_size = matrix_type::Zero (1, hypotheses.size());
                local_std_effect_size = matrix_type::Zero (1, hypotheses.size());
                local_stdev = matrix_type::Zero (num_vgs, 1);
                for (index_type ih = 0; ih != hypotheses.size(); ++ih) {
                  if (hypotheses[ih].is_F())
                    local_abs_effect_size (0, ih) = local_std_effect_size (0, ih) = std::numeric_limits<default_type>::quiet_NaN();
                }
              }
          };

          Source source (measurements.cols());
          Functor functor (measurements, fixed_design, extra_data, hypotheses, variance_groups,
                           cond, betas, abs_effect_size, std_effect_size, stdev);
          Thread::run_queue (source, Thread::batch (index_type()), Thread::multi (functor));
        }










        // Same model partitioning as is used in FSL randomise
        template <class MatrixType>
        Hypothesis::Partition Hypothesis::partition (const MatrixType& design) const
        {
          // eval() calls necessary for older versions of Eigen / compiler to work:
          //   can't seem to map Eigen template result to const matrix_type& as the Math::pinv() input
          const matrix_type D = (design.transpose() * design).inverse();
          // Note: Cu is transposed with respect to how contrast matrices are stored elsewhere
          const matrix_type Cu = Eigen::FullPivLU<matrix_type> (c).kernel();
          const matrix_type inv_cDc = (c * D * c.transpose()).inverse();
          // Note: Cv is transposed with respect to convention just as Cu is
          const matrix_type Cv = Cu - c.transpose() * inv_cDc * c * D * Cu;
          const matrix_type X = design * D * c.transpose() * inv_cDc;
          // .inverse() leads to NaNs with no nuisance regressors
          const matrix_type Z = Cv.isZero() ?
                                matrix_type::Zero (design.rows(), 1) :
                                (design * D * Cv * (Cv.transpose() * D * Cv).inverse()).eval();
          return Partition (X, Z);
        }
        template Hypothesis::Partition Hypothesis::partition (const matrix_type&) const;
        template Hypothesis::Partition Hypothesis::partition (const matrix_type::ConstRowsBlockXpr&) const;



        void Hypothesis::check_nonzero() const
        {
          if (c.isZero())
            throw Exception ("Cannot specify a contrast that consists entirely of zeroes");
        }



        matrix_type Hypothesis::check_rank (const matrix_type& in, const index_type index) const
        {
          // FullPivLU.image() provides column-space of matrix;
          //   here we want the row-space (since it's degeneracy in contrast matrix rows
          //   that has led to the rank-deficiency, whereas we can't exclude factor columns).
          //   Hence the transposing.
          Eigen::FullPivLU<matrix_type> decomp (in.transpose());
          if (decomp.rank() == in.rows())
            return in;
          WARN ("F-test " + str(index+1) + " is rank-deficient; row-space matrix decomposition will instead be used");
          INFO ("Original matrix: " + str(in));
          const matrix_type result = decomp.image (in.transpose()).transpose();
          INFO ("Decomposed matrix: " + str(result));
          return result;
        }







        SharedFixedBase::SharedFixedBase (const matrix_type& measurements,
                                          const matrix_type& design,
                                          const vector<Hypothesis>& hypotheses) :
            pinvM (Math::pinv (design)),
            Rm (matrix_type::Identity (measurements.rows(), measurements.rows()) - (design*pinvM))
        {
          for (const auto& h : hypotheses)
            partitions.emplace_back (h.partition (design));
        }



        SharedVariableBase::SharedVariableBase (const vector<CohortDataImport>& importers,
                                                const bool nans_in_data,
                                                const bool nans_in_columns) :
            importers (importers),
            nans_in_data (nans_in_data),
            nans_in_columns (nans_in_columns) { }




        SharedHeteroscedasticBase::SharedHeteroscedasticBase (const vector<Hypothesis>& hypotheses,
                                                              const index_array_type& variance_groups) :
            VG (variance_groups),
            num_vgs (variance_groups.maxCoeff() + 1),
            gamma_weights (vector_type::Zero (hypotheses.size()))
        {
          for (index_type ih = 0; ih != hypotheses.size(); ++ih) {
            const index_type s = hypotheses[ih].rank();
            gamma_weights[ih] = 2.0*(s-1) / default_type(s*(s+2));
          }
        }








        TestFixedHomoscedastic::Shared::Shared (const matrix_type& measurements, const matrix_type& design, const vector<Hypothesis>& hypotheses) :
            SharedFixedBase (measurements, design, hypotheses)
        {
          for (size_t ih = 0; ih != hypotheses.size(); ++ih) {
            XtX.emplace_back (partitions[ih].X.transpose() * partitions[ih].X);
            dof.push_back (measurements.rows() - partitions[ih].rank_x - partitions[ih].rank_z);
            one_over_dof.push_back (1.0 / default_type(dof[ih]));
#ifdef GLM_TEST_DEBUG
            VAR (ih);
            VAR (hypotheses[ih].matrix());
            VAR (partitions[ih].X);
            VAR (XtX[ih]);
            VAR (dof[ih]);
            VAR (one_over_dof[ih]);
#endif
          }
        }



        TestFixedHomoscedastic::TestFixedHomoscedastic (const matrix_type& measurements, const matrix_type& design, const vector<Hypothesis>& hypotheses) :
            TestBase (measurements, design, hypotheses),
#ifdef NDEBUG
            Sy (measurements.rows(), measurements.cols()),
            lambdas (design.cols(), measurements.cols()),
            residuals (measurements.rows(), measurements.cols()),
            sse (measurements.cols())
#else
            Sy (matrix_type::Constant (measurements.rows(), measurements.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            lambdas (matrix_type::Constant (design.cols(), measurements.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            residuals (matrix_type::Constant (measurements.rows(), measurements.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            sse (vector_type::Constant (measurements.cols(), std::numeric_limits<default_type>::signaling_NaN()))
#endif
        {
          shared.reset (new Shared (measurements, design, hypotheses));
          for (index_type ih = 0; ih != hypotheses.size(); ++ih)
#ifdef NDEBUG
            betas.emplace_back (matrix_type (hypotheses[ih].matrix().rows(), 1));
#else
            betas.emplace_back (matrix_type::Constant (hypotheses[ih].matrix().rows(), 1, std::numeric_limits<default_type>::signaling_NaN()));
#endif
        }



        TestFixedHomoscedastic::TestFixedHomoscedastic (const TestFixedHomoscedastic& that) :
            TestBase (that),
#ifdef NDEBUG
            Sy (num_inputs(), num_elements()),
            lambdas (num_factors(), num_elements()),
            residuals (num_inputs(), num_elements()),
            sse (num_elements())
#else
            Sy (matrix_type::Constant (num_inputs(), num_elements(), std::numeric_limits<default_type>::signaling_NaN())),
            lambdas (matrix_type::Constant (num_factors(), num_elements(), std::numeric_limits<default_type>::signaling_NaN())),
            residuals (matrix_type::Constant (num_inputs(), num_elements(), std::numeric_limits<default_type>::signaling_NaN())),
            sse (vector_type::Constant (num_elements(), std::numeric_limits<default_type>::signaling_NaN()))
#endif
        {
          for (index_type ih = 0; ih != num_hypotheses(); ++ih)
#ifdef NDEBUG
            betas.emplace_back (matrix_type (c[ih].matrix().rows(), 1));
#else
            betas.emplace_back (matrix_type::Constant (c[ih].matrix().rows(), 1, std::numeric_limits<default_type>::signaling_NaN()));
#endif
        }



        std::unique_ptr<TestBase> TestFixedHomoscedastic::__clone() const
        {
          std::unique_ptr<TestBase> result (new TestFixedHomoscedastic (*this));
          return result;
        }






        void TestFixedHomoscedastic::operator() (const shuffle_matrix_type& shuffling_matrix,
                                                 matrix_type& stats,
                                                 matrix_type& zstats)
        {
          assert (index_type(shuffling_matrix.rows()) == num_inputs());
          stats .resize (num_elements(), num_hypotheses());
          zstats.resize (num_elements(), num_hypotheses());

          // Freedman-Lane for fixed design matrix case
          // Each hypothesis needs to be handled explicitly on its own
          for (index_type ih = 0; ih != num_hypotheses(); ++ih) {

            // First, we perform permutation of the input data
            // In Freedman-Lane, the initial 'effective' regression against the nuisance
            //   variables, and permutation of the data, are done in a single step
#ifdef GLM_TEST_DEBUG
            VAR (shuffling_matrix.rows());
            VAR (shuffling_matrix.cols());
            VAR (S().partitions[ih].Rz.rows());
            VAR (S().partitions[ih].Rz.cols());
            VAR (y.rows());
            VAR (y.cols());
#endif
            Sy.noalias() = shuffling_matrix.cast<default_type>() * S().partitions[ih].Rz * y;
#ifdef GLM_TEST_DEBUG
            VAR (Sy.rows());
            VAR (Sy.cols());
            VAR (S().pinvM.rows());
            VAR (S().pinvM.cols());
#endif
            // Now, we regress this shuffled data against the full model
            lambdas.noalias() = S().pinvM * Sy;
#ifdef GLM_TEST_DEBUG
            VAR (lambdas.rows());
            VAR (lambdas.cols());
            VAR (S().Rm.rows());
            VAR (S().Rm.cols());
            VAR (S().XtX[ih].rows());
            VAR (S().XtX[ih].cols());
#endif
            sse = (S().Rm*Sy).colwise().squaredNorm();
#ifdef GLM_TEST_DEBUG
            VAR (sse.size());
#endif
            for (index_type ie = 0; ie != num_elements(); ++ie) {
              betas[ih].noalias() = c[ih].matrix() * lambdas.col (ie);
              const default_type F = ((betas[ih].transpose() * S().XtX[ih] * betas[ih]) (0,0) / c[ih].rank())
                                     / (S().one_over_dof[ih] * sse[ie]);
              if (!std::isfinite (F)) {
                stats  (ie, ih) = zstats (ie, ih) = value_type(0);
              } else if (c[ih].is_F()) {
                stats  (ie, ih) = F;
                zstats (ie, ih) =
#ifdef MRTRIX_USE_ZSTATISTIC_LOOKUP
                S().
#else
                Math::
#endif
                F2z (F, c[ih].rank(), S().dof[ih]);
              } else {
                assert (betas[ih].rows() == 1);
                stats  (ie, ih) = std::sqrt (F) * (betas[ih].sum() > 0.0 ? 1.0 : -1.0);
                zstats (ie, ih) =
#ifdef MRTRIX_USE_ZSTATISTIC_LOOKUP
                S().
#else
                Math::
#endif
                t2z (stats (ie, ih), S().dof[ih]);
              }
            }

          }
        }









        TestFixedHeteroscedastic::Shared::Shared (const matrix_type& measurements,
                                                  const matrix_type& design,
                                                  const vector<Hypothesis>& hypotheses,
                                                  const index_array_type& variance_groups) :
            SharedFixedBase (measurements, design, hypotheses),
            SharedHeteroscedasticBase (hypotheses, variance_groups),
            inputs_per_vg (num_vgs, 0),
            Rnn_sums (vector_type::Zero (num_vgs))
        {
          // Pre-calculate whatever can be pre-calculated for G-statistic
          for (index_type input = 0; input != measurements.rows(); ++input) {
            // Number of inputs belonging to each VG
            inputs_per_vg[variance_groups[input]]++;
            // Sum of diagonal entries of residual-forming matrix corresponding to each VG
            Rnn_sums[variance_groups[input]] += Rm.diagonal()[input];
          }
#ifdef GLM_TEST_DEBUG
          VAR (inputs_per_vg);
          VAR (Rnn_sums.transpose());
#endif
          // Reciprocals of the sums of diagonal entries of residual-forming matrix corresponding to each VG
          inv_Rnn_sums = Rnn_sums.inverse();
#ifdef GLM_TEST_DEBUG
          VAR (inv_Rnn_sums.transpose());
#endif
        }



        TestFixedHeteroscedastic::TestFixedHeteroscedastic (const matrix_type& measurements,
                                                            const matrix_type& design,
                                                            const vector<Hypothesis>& hypotheses,
                                                            const index_array_type& variance_groups) :
            TestBase (measurements, design, hypotheses),
#ifdef NDEBUG
            Sy (measurements.rows(), measurements.cols()),
            lambdas (design.cols(), measurements.cols()),
            sq_residuals (measurements.rows(), measurements.cols()),
            W (measurements.rows())

#else
            Sy (matrix_type::Constant (measurements.rows(), measurements.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            lambdas (matrix_type::Constant (design.cols(), measurements.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            sq_residuals (decltype(sq_residuals)::Constant (measurements.rows(), measurements.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            W (decltype(W)::Constant (measurements.rows(), std::numeric_limits<default_type>::signaling_NaN()))
#endif
        {
          shared.reset (new Shared (measurements, design, hypotheses, variance_groups));
          // Require shared to have been constructed first
#ifdef NDEBUG
          sse_per_vg.resize (num_variance_groups(), measurements.cols());
          Wterms.resize (num_variance_groups(), measurements.cols());
#else
          sse_per_vg = decltype(sse_per_vg)::Constant (num_variance_groups(), measurements.cols(), std::numeric_limits<default_type>::signaling_NaN());
          Wterms = decltype(Wterms)::Constant (num_variance_groups(), measurements.cols(), std::numeric_limits<default_type>::signaling_NaN());
#endif
        }



        TestFixedHeteroscedastic::TestFixedHeteroscedastic (const TestFixedHeteroscedastic& that) :
            TestBase (that),
#ifdef NDEBUG
            Sy (num_inputs(), num_elements()),
            lambdas (num_factors(), num_elements()),
            sq_residuals (num_inputs(), num_elements()),
            sse_per_vg (num_variance_groups(), num_elements()),
            Wterms (num_variance_groups(), num_elements()),
            W (num_inputs())
#else
            Sy (matrix_type::Constant (num_inputs(), num_elements(), std::numeric_limits<default_type>::signaling_NaN())),
            lambdas (matrix_type::Constant (num_factors(), num_elements(), std::numeric_limits<default_type>::signaling_NaN())),
            sq_residuals (decltype(sq_residuals)::Constant (num_inputs(), num_elements(), std::numeric_limits<default_type>::signaling_NaN())),
            sse_per_vg (decltype(sse_per_vg)::Constant (num_variance_groups(), num_elements(), std::numeric_limits<default_type>::signaling_NaN())),
            Wterms (decltype(Wterms)::Constant (num_variance_groups(), num_elements(), std::numeric_limits<default_type>::signaling_NaN())),
            W (decltype(W)::Constant (num_inputs(), std::numeric_limits<default_type>::signaling_NaN()))
#endif
        { }



        std::unique_ptr<TestBase> TestFixedHeteroscedastic::__clone() const
        {
          std::unique_ptr<TestBase> result (new TestFixedHeteroscedastic (*this));
          return result;
        }






        void TestFixedHeteroscedastic::operator() (const shuffle_matrix_type& shuffling_matrix, matrix_type& stats, matrix_type& zstats)
        {
          assert (index_type(shuffling_matrix.rows()) == num_inputs());
          stats.resize (num_elements(), num_hypotheses());
          zstats.resize (num_elements(), num_hypotheses());

#ifdef GLM_TEST_DEBUG
          //VAR (shuffling_matrix);
#endif

          for (index_type ih = 0; ih != num_hypotheses(); ++ih) {
            // First two steps are identical to the homoscedastic case
            Sy.noalias() = shuffling_matrix.cast<default_type>() * S().partitions[ih].Rz * y;
#ifdef GLM_TEST_DEBUG
            //VAR (Sy);
            VAR (Sy.rows());
            VAR (Sy.cols());
#endif
            lambdas.noalias() = S().pinvM * Sy;
#ifdef GLM_TEST_DEBUG
            //VAR (lambdas);
            VAR (lambdas.rows());
            VAR (lambdas.cols());
#endif
            // Compute sum of residuals per VG immediately
            // Variance groups appear across rows, and one column per element tested
            // Immediately calculate squared residuals; simplifies summation over variance groups
            sq_residuals = (S().Rm * Sy).array().square();
#ifdef GLM_TEST_DEBUG
            //VAR (sq_residuals);
#endif
            sse_per_vg.setZero();
            for (index_type input = 0; input != num_inputs(); ++input)
              sse_per_vg.row (S().VG[input]) += sq_residuals.row (input);
#ifdef GLM_TEST_DEBUG
            //VAR (sse_per_vg);
            VAR (sse_per_vg.rows());
            VAR (sse_per_vg.cols());
#endif
            // These terms are what appears in the weighting matrix based on the VG to which each input belongs;
            //   one row per variance group, one column per element to be tested
            Wterms = sse_per_vg.array().inverse().colwise() * S().Rnn_sums;
            for (index_type col = 0; col != num_elements(); ++col) {
              for (index_type row = 0; row != num_variance_groups(); ++row) {
                if (!std::isfinite (Wterms (row, col)))
                  Wterms (row, col) = 0.0;
              }
            }
#ifdef GLM_TEST_DEBUG
            //VAR (Wterms);
            VAR (Wterms.rows());
            VAR (Wterms.cols());
#endif
            for (index_type ie = 0; ie != num_elements(); ++ie) {
              // Need to construct the weights diagonal matrix; is unique for each element
              default_type W_trace (0.0);
              for (index_type input = 0; input != num_inputs(); ++input) {
                W[input] = Wterms (S().VG[input], ie);
                W_trace += W[input];
              }
#ifdef GLM_TEST_DEBUG
              VAR (W_trace);
#endif
              const default_type numerator = lambdas.col (ie).transpose() * c[ih].matrix().transpose() * (c[ih].matrix() * (M.transpose() * W.asDiagonal() * M).inverse() * c[ih].matrix().transpose()).inverse() * c[ih].matrix() * lambdas.col (ie);
#ifdef GLM_TEST_DEBUG
              VAR (numerator);
#endif
              default_type gamma (0.0);
              for (index_type vg_index = 0; vg_index != num_variance_groups(); ++vg_index)
                // Since Wnn is the same for every n in the variance group, can compute that summation as the product of:
                //   - the value inserted in W for that particular VG
                //   - the number of inputs that are a part of that VG
                gamma += S().inv_Rnn_sums[vg_index] * Math::pow2 (1.0 - ((Wterms(vg_index, ie) * S().inputs_per_vg[vg_index]) / W_trace));
              gamma = 1.0 + (S().gamma_weights[ih] * gamma);
#ifdef GLM_TEST_DEBUG
              VAR (gamma);
#endif
              const default_type denominator = gamma * c[ih].rank();
              const default_type G = numerator / denominator;
              if (!std::isfinite (G)) {
                stats  (ie, ih) = zstats (ie, ih) = value_type(0);
              } else {
                stats  (ie, ih) = c[ih].is_F() ?
                                  G :
                                  std::sqrt (G) * ((c[ih].matrix() * lambdas.col (ie)).sum() > 0.0 ? 1.0 : -1.0);
                if (c[ih].is_F() && c[ih].rank() > 1) {
                  const default_type dof = 2.0 * default_type(c[ih].rank() - 1) / (3.0 * (gamma - 1.0));
#ifdef GLM_TEST_DEBUG
                  VAR (dof);
#endif
                  zstats (ie, ih) = Math::G2z (G, c[ih].rank(), dof);
                } else {
                  const default_type dof = Math::welch_satterthwaite (Wterms.col (ie).inverse(), S().inputs_per_vg);
#ifdef GLM_TEST_DEBUG
                  VAR (dof);
#endif
                  zstats (ie, ih) = c[ih].is_F() ?
                                    Math::G2z (G, c[ih].rank(), dof) :
                                    Math::v2z (stats (ie, ih), dof);
                }
              }
            }

          }
        }










        TestVariableBase::TestVariableBase (const matrix_type& measurements,
                                            const matrix_type& design,
                                            const vector<Hypothesis>& hypotheses,
                                            const vector<CohortDataImport>& importers) :
            TestBase (measurements, design, hypotheses),
#ifdef NDEBUG
            dof (measurements.cols(), hypotheses.size()),
            extra_column_data (measurements.rows(), importers.size()),
            element_mask (measurements.rows()),
            permuted_mask (measurements.rows()),
            intermediate_shuffling_matrix (measurements.rows(), measurements.rows()),
            shuffling_matrix_masked (measurements.rows(), measurements.rows()),
            Mfull_masked (measurements.rows(), design.cols() + importers.size()),
            pinvMfull_masked (design.cols() + importers.size(), measurements.rows()),
            Rm (measurements.rows(), measurements.rows()),
            y_masked (measurements.rows()),
            Sy (measurements.rows()),
            lambda (design.cols() + importers.size())
#else
            dof (matrix_type::Constant (measurements.cols(), hypotheses.size(), std::numeric_limits<default_type>::signaling_NaN())),
            extra_column_data (matrix_type::Constant (measurements.rows(), importers.size(), std::numeric_limits<default_type>::signaling_NaN())),
            element_mask (measurements.rows()),
            permuted_mask (measurements.rows()),
            intermediate_shuffling_matrix (shuffle_matrix_type::Constant (measurements.rows(), measurements.rows(), std::numeric_limits<int8_t>::min())),
            shuffling_matrix_masked (shuffle_matrix_type::Constant (measurements.rows(), measurements.rows(), std::numeric_limits<int8_t>::min())),
            Mfull_masked (matrix_type::Constant (measurements.rows(), design.cols() + importers.size(), std::numeric_limits<default_type>::signaling_NaN())),
            pinvMfull_masked (matrix_type::Constant (design.cols() + importers.size(), measurements.rows(), std::numeric_limits<default_type>::signaling_NaN())),
            Rm (matrix_type::Constant (measurements.rows(), measurements.rows(), std::numeric_limits<default_type>::signaling_NaN())),
            y_masked (vector_type::Constant (measurements.rows(), std::numeric_limits<default_type>::signaling_NaN())),
            Sy (vector_type::Constant (measurements.rows(), std::numeric_limits<default_type>::signaling_NaN())),
            lambda (vector_type::Constant (design.cols() + importers.size(), std::numeric_limits<default_type>::signaling_NaN()))
#endif
        { }



        TestVariableBase::TestVariableBase (const TestVariableBase& that) :
            TestBase (that),
#ifdef NDEBUG
            dof (that.dof.rows(), that.dof.cols()),
            extra_column_data (that.extra_column_data.rows(), that.extra_column_data.cols()),
            element_mask (that.element_mask.size()),
            permuted_mask (that.permuted_mask.size()),
            intermediate_shuffling_matrix (that.intermediate_shuffling_matrix.rows(), that.intermediate_shuffling_matrix.cols()),
            shuffling_matrix_masked (that.shuffling_matrix_masked.rows(), that.shuffling_matrix_masked.cols()),
            Mfull_masked (that.Mfull_masked.rows(), that.Mfull_masked.cols()),
            pinvMfull_masked (that.pinvMfull_masked.rows(), that.pinvMfull_masked.cols()),
            Rm (that.Rm.rows(), that.Rm.cols()),
            y_masked (that.y_masked.size()),
            Sy (that.Sy.size()),
            lambda (that.lambda.size())
#else
            dof (matrix_type::Constant (that.dof.rows(), that.dof.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            extra_column_data (matrix_type::Constant (that.extra_column_data.rows(), that.extra_column_data.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            element_mask (that.element_mask.size()),
            permuted_mask (that.permuted_mask.size()),
            intermediate_shuffling_matrix (shuffle_matrix_type::Constant (that.intermediate_shuffling_matrix.rows(), that.intermediate_shuffling_matrix.cols(), std::numeric_limits<int8_t>::min())),
            shuffling_matrix_masked (shuffle_matrix_type::Constant (that.shuffling_matrix_masked.rows(), that.shuffling_matrix_masked.cols(), std::numeric_limits<int8_t>::min())),
            Mfull_masked (matrix_type::Constant (that.Mfull_masked.rows(), that.Mfull_masked.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            pinvMfull_masked (matrix_type::Constant (that.pinvMfull_masked.rows(), that.pinvMfull_masked.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            Rm (matrix_type::Constant (that.Rm.rows(), that.Rm.cols(), std::numeric_limits<default_type>::signaling_NaN())),
            y_masked (vector_type::Constant (that.y_masked.size(), std::numeric_limits<default_type>::signaling_NaN())),
            Sy (vector_type::Constant (that.Sy.size(), std::numeric_limits<default_type>::signaling_NaN())),
            lambda (vector_type::Constant (that.lambda.size(), std::numeric_limits<default_type>::signaling_NaN()))
#endif
        { }



        template <class SharedType>
        void TestVariableBase::set_mask (const SharedType& s, const index_type ie)
        {
          element_mask.clear (true);
          if (s.nans_in_data) {
            for (index_type row = 0; row != num_inputs(); ++row) {
              if (!std::isfinite (y (row, ie)))
                element_mask[row] = false;
            }
          }
          if (s.nans_in_columns) {
            for (index_type row = 0; row != extra_column_data.rows(); ++row) {
              if (!extra_column_data.row (row).allFinite())
                element_mask[row] = false;
            }
          }
        }
        template void TestVariableBase::set_mask (const TestVariableHomoscedastic::Shared&, const index_type);
        template void TestVariableBase::set_mask (const TestVariableHeteroscedastic::Shared&, const index_type);



        void TestVariableBase::apply_mask (const index_type ie, const shuffle_matrix_type& shuffling_matrix)
        {
          const index_type finite_count = element_mask.count();
          // Do we need to reduce the size of our matrices / vectors
          //   based on the presence of non-finite values?
          if (finite_count == num_inputs()) {

            Mfull_masked.block (0, 0, num_inputs(), M.cols()) = M;
            if (num_importers())
              Mfull_masked.block (0, M.cols(), num_inputs(), extra_column_data.cols()) = extra_column_data;
            shuffling_matrix_masked = shuffling_matrix;
            y_masked = y.col (ie);

          } else {

            permuted_mask.clear (true);
            index_type out_index = 0;
            for (index_type in_index = 0; in_index != num_inputs(); ++in_index) {
              if (element_mask[in_index]) {
                Mfull_masked.block (out_index, 0, 1, M.cols()) = M.row (in_index);
                if (num_importers())
                  Mfull_masked.block (out_index, M.cols(), 1, extra_column_data.cols()) = extra_column_data.row (in_index);
                y_masked[out_index++] = y(in_index, ie);
              } else {
                // Any row in the permutation matrix that contains a non-zero entry
                //   in the column corresponding to in_row needs to be removed
                //   from the permutation matrix
                for (index_type perm_row = 0; perm_row != shuffling_matrix.rows(); ++perm_row) {
                  if (shuffling_matrix (perm_row, in_index))
                    permuted_mask[perm_row] = false;
                }
              }
            }
            assert (out_index == finite_count);
            assert (permuted_mask.count() == finite_count);
            assert (Mfull_masked.topRows (finite_count).allFinite());
            assert (y_masked.head (finite_count).allFinite());
#ifndef NDEBUG
            Mfull_masked.bottomRows(num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
            y_masked.tail (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
#endif
            // Only after we've reduced the design matrix do we now reduce the shuffling matrix
            // Step 1: Remove rows that contain non-zero entries in columns to be removed
            out_index = 0;
            for (index_type in_index = 0; in_index != num_inputs(); ++in_index) {
              if (permuted_mask[in_index])
                intermediate_shuffling_matrix.row (out_index++) = shuffling_matrix.row (in_index);
            }
            assert (out_index == finite_count);
#ifndef NDEBUG
            intermediate_shuffling_matrix.bottomRows (num_inputs() - finite_count).fill (std::numeric_limits<int8_t>::min());
#endif
            // Step 2: Remove columns
            out_index = 0;
            for (index_type in_index = 0; in_index != num_inputs(); ++in_index) {
              if (element_mask[in_index])
                shuffling_matrix_masked.col (out_index++) = intermediate_shuffling_matrix.col (in_index);
            }
            assert (out_index == finite_count);
#ifndef NDEBUG
            shuffling_matrix_masked.rightCols (num_inputs() - finite_count).fill (std::numeric_limits<int8_t>::min());
#endif

          }
        }













        TestVariableHomoscedastic::Shared::Shared (const vector<CohortDataImport>& importers,
                                                   const bool nans_in_data,
                                                   const bool nans_in_columns) :
            SharedVariableBase (importers, nans_in_data, nans_in_columns)
#ifdef MRTRIX_USE_ZSTATISTIC_LOOKUP
            , SharedHomoscedasticBase ()
#endif
            { }



        TestVariableHomoscedastic::TestVariableHomoscedastic (const matrix_type& measurements,
                                                              const matrix_type& design,
                                                              const vector<Hypothesis>& hypotheses,
                                                              const vector<CohortDataImport>& importers,
                                                              const bool nans_in_data,
                                                              const bool nans_in_columns) :
            TestVariableBase (measurements, design, hypotheses, importers)
        {
          shared.reset (new Shared(importers, nans_in_data, nans_in_columns));
          for (index_type ih = 0; ih != hypotheses.size(); ++ih) {
#ifdef NDEBUG
            XtX.emplace_back (hypotheses[ih].matrix().rows(), hypotheses[ih].matrix().rows());
            beta.emplace_back (hypotheses[ih].matrix().rows(), 1);
#else
            XtX.emplace_back (matrix_type::Constant (hypotheses[ih].matrix().rows(), hypotheses[ih].matrix().rows(), std::numeric_limits<default_type>::signaling_NaN()));
            beta.emplace_back (matrix_type::Constant (hypotheses[ih].matrix().rows(), 1, std::numeric_limits<default_type>::signaling_NaN()));
#endif
          }
        }



        TestVariableHomoscedastic::TestVariableHomoscedastic (const TestVariableHomoscedastic& that) :
            TestVariableBase (that)
        {
          for (index_type ih = 0; ih != num_hypotheses(); ++ih) {
#ifdef NDEBUG
            XtX.emplace_back (c[ih].matrix().rows(), c[ih].matrix().rows());
            beta.emplace_back (c[ih].matrix().rows(), 1);
#else
            XtX.emplace_back (matrix_type::Constant (c[ih].matrix().rows(), c[ih].matrix().rows(), std::numeric_limits<default_type>::signaling_NaN()));
            beta.emplace_back (matrix_type::Constant (c[ih].matrix().rows(), 1, std::numeric_limits<default_type>::signaling_NaN()));
#endif
          }
        }




        std::unique_ptr<TestBase> TestVariableHomoscedastic::__clone() const
        {
          std::unique_ptr<TestBase> result (new TestVariableHomoscedastic (*this));
          return result;
        }






        void TestVariableHomoscedastic::operator() (const shuffle_matrix_type& shuffling_matrix,
                                                    matrix_type& stats,
                                                    matrix_type& zstats)
        {
          stats .resize (num_elements(), num_hypotheses());
          zstats.resize (num_elements(), num_hypotheses());

          // Let's loop over elements first, then hypotheses in the inner loop
          for (index_type ie = 0; ie != num_elements(); ++ie) {

            // For each element (row in y), need to load the additional data for that element
            //   for all subjects in order to construct the design matrix
            // Would it be preferable to pre-calculate and store these per-element design matrices,
            //   rather than re-generating them each time? (More RAM, less CPU)
            // No, most of the time that subject data will be memory-mapped, so pre-loading (in
            //   addition to the duplication of the fixed design matrix contents) would hurt bad
            for (index_type col = 0; col != num_importers(); ++col)
              extra_column_data.col (col) = S().importers[col] (ie);

            // What can we do here that's common across all hypotheses?
            // - Import the element-wise data
            // - Identify rows to be excluded based on NaNs in the design matrix
            // - Identify rows to be excluded based on NaNs in the input data
            //
            // Note that this is going to have to operate slightly differently to
            //   how it used to be done, i.e. via the permutation labelling vector,
            //   if we are to support taking the shuffling matrix as input to this functor
            // I think the approach will have to be:
            //   - Both NaNs in design matrix and NaNs in input data need to be removed
            //     in order to perform the initial regression against nuisance variables
            //   - Can then remove the corresponding _columns_ of the permutation matrix?
            //     No, don't think it's removal of columns; think it's removal of any rows
            //     that contain non-zero values in those columns
            //
            set_mask (S(), ie);
            const index_type finite_count = element_mask.count();
#ifdef GLM_TEST_DEBUG
            VAR (element_mask.size());
            VAR (finite_count);
#endif
            // Additional rejection here:
            // If the number of finite elements is _not_ equal to the number of inputs
            //   (i.e. at least one input has been removed), there needs to be a
            //   more stringent criterion met in order to proceed with the test.
            if (finite_count < std::min (num_inputs(), 2 * num_factors())) {
              stats.row (ie).setZero();
              zstats.row (ie).setZero();
              dof.row (ie).fill (std::numeric_limits<default_type>::quiet_NaN());
            } else {
              apply_mask (ie, shuffling_matrix);
              assert (Mfull_masked.topRows (finite_count).allFinite());

              // Test condition number of NaN-masked & data-filled design matrix;
              //   need to skip statistical testing if it is too poor
              // TODO Condition number testing may be quite slow;
              //   would a rank calculation with tolerance be faster?
              // TODO JacobiSVD refuses to run on an Eigen::Block due to member "Options" not being defined
              const default_type condition_number = Math::condition_number (Mfull_masked.topRows (finite_count).eval());
#ifdef GLM_TEST_DEBUG
              VAR (condition_number);
#endif
              if (!std::isfinite (condition_number) || condition_number > 1e5) {
                stats.row (ie).fill (0.0);
                zstats.row (ie).fill (0.0);
                dof.row (ie).fill (std::numeric_limits<default_type>::quiet_NaN());
              } else {

                pinvMfull_masked.leftCols (finite_count).noalias() = Math::pinv (Mfull_masked.topRows (finite_count));
#ifndef NDEBUG
                pinvMfull_masked.rightCols (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
#endif
                Rm.topLeftCorner (finite_count, finite_count).noalias() = matrix_type::Identity (finite_count, finite_count) - (Mfull_masked.topRows (finite_count) * pinvMfull_masked.leftCols (finite_count));
#ifndef NDEBUG
                Rm.bottomRows (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
                Rm.rightCols (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
#endif

                // We now have our permutation (shuffling) matrix and design matrix prepared,
                //   and can commence regressing the partitioned model of each hypothesis
                for (index_type ih = 0; ih != num_hypotheses(); ++ih) {

                  const auto partition = c[ih].partition (Mfull_masked.topRows (finite_count));
                  dof (ie, ih) = finite_count - partition.rank_x - partition.rank_z;
                  if (dof (ie, ih) < 1) {
                    stats (ie, ih) = zstats (ie, ih) = dof (ie, ih) = value_type(0);
                  } else {
                    XtX[ih].noalias() = partition.X.transpose() * partition.X;
#ifdef GLM_TEST_DEBUG
                    VAR (XtX[ih].rows());
                    VAR (XtX[ih].cols());
#endif
                    // Now that we have the individual hypothesis model partition for these data,
                    //   the rest of this function should proceed similarly to the fixed
                    //   design matrix case
                    Sy.head (finite_count) = shuffling_matrix_masked.topLeftCorner (finite_count, finite_count).cast<default_type>() * partition.Rz * y_masked.head (finite_count).matrix();
                    lambda = pinvMfull_masked.leftCols(finite_count) * Sy.head(finite_count).matrix();
                    beta[ih].noalias() = c[ih].matrix() * lambda.matrix();
#ifdef GLM_TEST_DEBUG
                    VAR (Sy.size());
                    VAR (lambda.size());
                    VAR (beta[ih].rows());
                    VAR (beta[ih].cols());
#endif
                    const default_type sse = (Rm.topLeftCorner (finite_count, finite_count) * Sy.head (finite_count).matrix()).squaredNorm();

                    const default_type F = ((beta[ih].transpose() * XtX[ih] * beta[ih]) (0, 0) / c[ih].rank()) /
                                            (sse / dof (ie, ih));

                    if (!std::isfinite (F)) {
                      stats  (ie, ih) = zstats (ie, ih) = value_type(0);
                    } else if (c[ih].is_F()) {
                      stats  (ie, ih) = F;
                      zstats (ie, ih) =
#ifdef MRTRIX_USE_ZSTATISTIC_LOOKUP
                      S().
#else
                      Math::
#endif
                      F2z (F, c[ih].rank(), dof (ie, ih));
                    } else {
                      assert (beta[ih].rows() == 1);
                      stats  (ie, ih) = std::sqrt (F) * (beta[ih].sum() > 0 ? 1.0 : -1.0);
                      zstats (ie, ih) =
#ifdef MRTRIX_USE_ZSTATISTIC_LOOKUP
                      S().
#else
                      Math::
#endif
                      t2z (stats (ie, ih), dof (ie, ih));
                    }

                  } // End checking for sufficient degrees of freedom

                } // End looping over hypotheses

              } // End checking for adequate condition number after NaN removal

            } // End checking for adequate number of remaining inputs after NaN removal

          } // End looping over elements

        } // End functor














        TestVariableHeteroscedastic::Shared::Shared (const vector<Hypothesis>& hypotheses,
                                                     const index_array_type& variance_groups,
                                                     const vector<CohortDataImport>& importers,
                                                     const bool nans_in_data,
                                                     const bool nans_in_columns) :
            SharedVariableBase (importers, nans_in_data, nans_in_columns),
            SharedHeteroscedasticBase (hypotheses, variance_groups) { }



        TestVariableHeteroscedastic::TestVariableHeteroscedastic (const matrix_type& measurements,
                                                                  const matrix_type& design,
                                                                  const vector<Hypothesis>& hypotheses,
                                                                  const index_array_type& variance_groups,
                                                                  const vector<CohortDataImport>& importers,
                                                                  const bool nans_in_data,
                                                                  const bool nans_in_columns) :
            TestVariableBase (measurements, design, hypotheses, importers),
#ifdef NDEBUG
            W (measurements.rows()),
            sq_residuals (measurements.rows()),
#else
            W (vector_type::Constant (measurements.rows(), std::numeric_limits<default_type>::signaling_NaN())),
            sq_residuals (vector_type::Constant (measurements.rows(), std::numeric_limits<default_type>::signaling_NaN())),
#endif
            VG_masked (measurements.rows())
        {
          shared.reset (new Shared (hypotheses, variance_groups, importers, nans_in_data, nans_in_columns));
          // Require shared to have been initialised
#ifdef NDEBUG
          sse_per_vg.resize (num_variance_groups());
          Rnn_sums.resize (num_variance_groups());
          Wterms.resize (num_variance_groups());
          VG_counts.resize (num_variance_groups());
#else
          sse_per_vg = vector_type::Constant (num_variance_groups(), std::numeric_limits<default_type>::signaling_NaN());
          Rnn_sums = vector_type::Constant (num_variance_groups(), std::numeric_limits<default_type>::signaling_NaN());
          Wterms = vector_type::Constant (num_variance_groups(), std::numeric_limits<default_type>::signaling_NaN());
          VG_counts = index_array_type::Zero (num_variance_groups());
#endif
        }



        TestVariableHeteroscedastic::TestVariableHeteroscedastic (const TestVariableHeteroscedastic& that) :
            TestVariableBase (that),
#ifdef NDEBUG
            W (num_inputs()),
            sq_residuals (num_inputs()),
            sse_per_vg (num_variance_groups()),
            Rnn_sums (num_variance_groups()),
            Wterms (num_variance_groups()),
            VG_masked (num_inputs()),
            VG_counts (num_variance_groups())
#else
            W (vector_type::Constant (num_inputs(), std::numeric_limits<default_type>::signaling_NaN())),
            sq_residuals (vector_type::Constant (num_inputs(), std::numeric_limits<default_type>::signaling_NaN())),
            sse_per_vg (vector_type::Constant (num_variance_groups(), std::numeric_limits<default_type>::signaling_NaN())),
            Rnn_sums (vector_type::Constant (num_variance_groups(), std::numeric_limits<default_type>::signaling_NaN())),
            Wterms (vector_type::Constant (num_variance_groups(), std::numeric_limits<default_type>::signaling_NaN())),
            VG_masked (num_inputs()),
            VG_counts (index_array_type::Zero (num_variance_groups()))
#endif
        { }



        std::unique_ptr<TestBase> TestVariableHeteroscedastic::__clone() const
        {
          std::unique_ptr<TestBase> result (new TestVariableHeteroscedastic (*this));
          return result;
        }







        void TestVariableHeteroscedastic::operator() (const shuffle_matrix_type& shuffling_matrix, matrix_type& stats, matrix_type& zstats)
        {
          stats.resize (num_elements(), num_hypotheses());
          zstats.resize (num_elements(), num_hypotheses());

          for (index_type ie = 0; ie != num_elements(); ++ie) {
            // Common ground to the TestVariableHomoscedastic case
            for (index_type col = 0; col != num_importers(); ++col)
              extra_column_data.col (col) = S().importers[col] (ie);
            set_mask (S(), ie);
            const index_type finite_count = element_mask.count();
            if (finite_count < std::min (num_inputs(), 2 * num_factors())) {
              stats.row (ie).setZero();
              zstats.row (ie).setZero();
            } else {
              apply_mask (ie, shuffling_matrix);
              const default_type condition_number = Math::condition_number (Mfull_masked.topRows (finite_count).eval());
              if (!std::isfinite (condition_number) || condition_number > 1e5) {
                stats.row (ie).fill (0.0);
                zstats.row (ie).fill (0.0);
              } else {
                VG_counts.setZero();
                index_type out_index = 0;
                for (index_type in_index = 0; in_index != num_inputs(); ++in_index) {
                  if (element_mask[in_index]) {
                    VG_masked[out_index++] = S().VG[in_index];
                    VG_counts[S().VG[in_index]]++;
                  }
                }
                assert (out_index == finite_count);
                if (VG_counts.minCoeff() <= 1) {
                  stats.row (ie).fill (0.0);
                  zstats.row (ie).fill (0.0);
                } else {
                  pinvMfull_masked.leftCols (finite_count) = Math::pinv (Mfull_masked.topRows (finite_count));
#ifndef NDEBUG
                  pinvMfull_masked.rightCols (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
#endif
                  Rm.topLeftCorner (finite_count, finite_count).noalias() = matrix_type::Identity (finite_count, finite_count) - (Mfull_masked.topRows (finite_count) * pinvMfull_masked.leftCols (finite_count));
#ifndef NDEBUG
                  Rm.bottomRows (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
                  Rm.rightCols (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
#endif
                  for (index_type ih = 0; ih != num_hypotheses(); ++ih) {
                    const auto partition = c[ih].partition (Mfull_masked.topRows (finite_count));

                    // At this point the implementation diverges from the TestVariableHomoscedastic case,
                    //   more closely mimicing the TestFixedHeteroscedastic case
                    Sy.head (finite_count) = shuffling_matrix_masked.topLeftCorner (finite_count, finite_count).cast<default_type>() * partition.Rz * y_masked.head (finite_count).matrix();
#ifndef NDEBUG
                    Sy.tail (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
#endif
                    lambda = pinvMfull_masked.leftCols (finite_count) * Sy.head (finite_count).matrix();
                    sq_residuals.head (finite_count) = (Rm.topLeftCorner (finite_count, finite_count) * Sy.head (finite_count).matrix()).array().square();
#ifndef NDEBUG
                    sq_residuals.tail (num_inputs() - finite_count).fill (std::numeric_limits<default_type>::signaling_NaN());
#endif
                    sse_per_vg.setZero();
                    Rnn_sums.setZero();
                    for (index_type input = 0; input != finite_count; ++input) {
                      sse_per_vg[VG_masked[input]] += sq_residuals[input];
                      Rnn_sums[VG_masked[input]] += Rm.diagonal()[input];
                    }
                    Wterms = sse_per_vg.inverse() * Rnn_sums;
                    for (index_type vg = 0; vg != num_variance_groups(); ++vg) {
                      if (!std::isfinite (Wterms[vg]))
                        Wterms[vg] = 0.0;
                    }
                    default_type W_trace (0.0);
                    for (index_type input = 0; input != finite_count; ++input) {
                      W[input] = Wterms[VG_masked[input]];
                      W_trace += W[input];
                    }

                    const default_type numerator = lambda.matrix().transpose() * c[ih].matrix().transpose() * (c[ih].matrix() * (Mfull_masked.topRows (finite_count).transpose() * W.matrix().head (finite_count).asDiagonal() * Mfull_masked.topRows (finite_count)).inverse() * c[ih].matrix().transpose()).inverse() * c[ih].matrix() * lambda.matrix();

                    default_type gamma (0.0);
                    for (index_type vg_index = 0; vg_index != num_variance_groups(); ++vg_index)
                      gamma += Math::pow2 (1.0 - ((Wterms[vg_index] * VG_counts[vg_index]) / W_trace)) / Rnn_sums[vg_index];
                    gamma = 1.0 + (S().gamma_weights[ih] * gamma);

                    const default_type denominator = gamma * c[ih].rank();
                    const default_type G = numerator / denominator;

                    if (!std::isfinite (G)) {
                      stats  (ie, ih) = zstats (ie, ih) = value_type(0);
                    } else {
                      stats  (ie, ih) = c[ih].is_F() ?
                                        G :
                                        std::sqrt (G) * ((c[ih].matrix() * lambda.matrix()).sum() > 0.0 ? 1.0 : -1.0);
                      if (c[ih].is_F() && c[ih].rank() > 1) {
                        const default_type dof = 2.0 * default_type(c[ih].rank() - 1) / (3.0 * (gamma - 1.0));
                        zstats (ie, ih) = Math::G2z (G, c[ih].rank(), dof);
                      } else {
                        const default_type dof = Math::welch_satterthwaite (Wterms.inverse(), VG_counts);
                        zstats (ie, ih) = c[ih].is_F() ?
                                          Math::G2z (G, c[ih].rank(), dof) :
                                          Math::v2z (stats (ie, ih), dof);
                      } // End switching for F-test with rank > 1

                    } // End checking for G being finite

                  } // End looping over hypotheses for this element

                } // End check for preservation of at least two elements in each VG

              } // End checking for adequate condition number after NaN removal

            } // End checking for adequate number of remaining inputs after NaN removal

          } // End looping over elements
        }






      }
    }
  }
}
