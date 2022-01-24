#ifndef __dwi_sdeconv_response_h__
#define __dwi_sdeconv_response_h__

#include "math/ZSH.h"
#include "math/math.h"
#include "mrtrix.h"
#include "debug.h"

namespace MR
{
  namespace DWI
  {
    namespace SDeconv
    {
      class Response
      {

      public:
        Response() {}
        Response(const std::string &filename) { load(filename); }

        Eigen::VectorXd coeffs(const double bval)
        {
          if (bval < original_bvals[0])
            throw Exception("bvalue out of bounds");

          if (bval >= original_bvals.back())
            return original_coefs.row(original_coefs.rows() - 1);

          int i = 0;
          while (bval > original_bvals[i + 1])
            i++;

          float ratio = (bval - original_bvals[i]) / (original_bvals[i + 1] - original_bvals[i]);
          return (1.0 - ratio) * original_coefs.row(i) + ratio * original_coefs.row(i + 1);
        }

        int lmax() const { return Math::ZSH::LforN(original_coefs.cols()); }

        const Eigen::VectorXd &bvalues() const { return original_bvals; }

        // LOAD MATRIX
        void load(const std::string &filename)
        {
          vector<std::string> comments;
          auto vec = load_matrix_2D_vector(filename, &comments);

          //load bvalues
          for (const auto &comment : comments)
          {
            auto n = comment.rfind("Shells:");
            if (n != comment.npos)
              original_bvals = parse_floats(comment.substr(n + 7));
          }

          if (original_bvals.size() != vec.size())
            throw Exception("Number of b-values does not match the number of rows in response");

          if (!std::is_sorted(original_bvals.begin(), original_bvals.end()))
            throw Exception("B-values are not sorted in ascending order");

          original_coefs.resize(vec.size(), vec[0].size());

          for (ssize_t r = 0; r < vec.size(); r++)
            for (ssize_t c = 0; c < vec[0].size(); c++)
              original_coefs(r, c) = vec[r][c];
        }

      private:
        Eigen::MatrixXd original_coefs;
        vector<double> original_bvals;
      };
    }
  }
}

#endif
