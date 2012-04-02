// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_AFW_MATH_LeastSquares_h_INCLUDED
#define LSST_AFW_MATH_LeastSquares_h_INCLUDED

#include "lsst/base.h"
#include "ndarray/eigen.h"

namespace lsst { namespace afw { namespace math {

/**
 *  @brief Solver for linear least-squares problems.
 *
 *  Linear least-squares problems are defined as finding the vector @f$x@f$ that minimizes 
 *  @f$\left|A x - b\right|_2@f$, with the number of rows of @f$A@f$ generally
 *  greater than the number of columns.  We call @f$A@f$ the design matrix, @f$b@f$
 *  the data vector, and @f$x@f$ the solution vector.  When the rank of @f$A@f$ is
 *  equal to the number of columns, we can obtain using the solution using the normal equations:
 *  @f[
 *      A^T A x = A^T b
 *  @f]
 *  If @f$A@f$ is not full-rank, the problem is underconstrained, and we usually wish to
 *  solve the minimum-norm least-squares problem, which also minimizes @f$|x|_2@f$.
 *  This can be done by computing a pseudo-inverse of @f$A@f$ using an SVD or complete
 *  orthogonal factorization of @f$A@f$, or by performing an Eigen decomposition of
 *  @f$A^T A@f$.
 *
 *  This class can be constructed from the design matrix and data vector, or from the two terms
 *  in the normal equations (below, we call the matrix @f$A^TA@f$ the Fisher matrix, and the
 *  vector @f$A^T b@f$ simply the "right-hand side" (RHS) vector).  If initialized with
 *  the design matrix and data vector, we can still use the normal equations to solve it.
 *  The solution via the normal equations is more susceptible to round-off error, but it is also
 *  faster, and if the normal equation terms can be computed directly it can be significantly
 *  less expensive in terms of memory.  The Fisher matrix is a symmetric matrix, and it should
 *  be exactly symmetric when provided as input, because which triangle will be used is an
 *  implementation detail that is subject to change.
 *
 *  The solver always operates in double precision, and returns all results in double precision.
 *  However, it can be initialized from single precision inputs.  It isn't usually a good idea
 *  to construct the normal equations in single precision, however, even when the data are
 *  single precision.
 */
class LeastSquares {
public:

    class Impl; ///< Private implementation; forward-declared publicly so we can inherit from it in .cc

    enum Factorization {
        NORMAL_EIGENSYSTEM,  /**<
                              *   @brief Use the normal equations with a symmetric Eigensystem decomposition.
                              *
                              *   This method is fully robust and computes the minimum-norm solution when
                              *   the problem does not have full rank.  It is affected slightly more by
                              *   round-off error than the SVD method, but because it can handle singular
                              *   matrices this usually isn't a problem.
                              */
        NORMAL_CHOLESKY,     /**<
                              *   @brief Use the normal equations with a Cholesky decomposition.
                              *
                              *   While this method uses a robust LDL^T decomposition that does not
                              *   require square roots, it is not appropriate for problems that do
                              *   not have full rank, and cannot be used to determine whether a
                              *   problem has full rank.  It is the fastest decomposition.
                              */
        DIRECT_SVD           /**<
                              *   @brief Use a thin singular value decomposition of the design matrix.
                              *
                              *   This method is the most robust and computes the minimum-norm solution
                              *   the problem does not have full rank.  However, it is also the slowest
                              *   and is not available when the solver is initialized with the
                              *   normal equations.
                              */
    };

    /// @brief Initialize from the design matrix and data vector given as ndarrays.
    template <typename T1, typename T2, int C1, int C2>
    static LeastSquares fromDesignMatrix(
        ndarray::Array<T1 const,2,C1> const & design,
        ndarray::Array<T2 const,1,C2> const & data,
        Factorization factorization = NORMAL_EIGENSYSTEM
    ) {
        LeastSquares r(factorization, design.template getSize<1>());
        r._getDesignMatrix() = design.asEigen();
        r._getDataVector() = data.asEigen();
        r._factor(false);
        return r;
    }

    /// @brief Initialize from the design matrix and data vector given as an Eigen objects.
    template <typename D1, typename D2>
    static LeastSquares fromDesignMatrix(
        Eigen::MatrixBase<D1> const & design,
        Eigen::MatrixBase<D2> const & data,
        Factorization factorization = NORMAL_EIGENSYSTEM
    ) {
        LeastSquares r(factorization, design.cols());
        r._getDesignMatrix() = design;
        r._getDataVector() = data;
        r._factor(false);
        return r;
    }

    /// @brief Reset the design matrix and data vector given as ndarrays; dimension must not change.
    template <typename T1, typename T2, int C1, int C2>
    void setDesignMatrix(
        ndarray::Array<T1 const,2,C1> const & design,
        ndarray::Array<T2 const,1,C2> const & data
    ) {
        _getDesignMatrix() = design.asEigen();
        _getDataVector() = data.asEigen();
        _factor(false);
    }

    /// @brief Reset the design matrix and data vector given as Eigen objects; dimension must not change.
    template <typename D1, typename D2>
    void setDesignMatrix(
        Eigen::MatrixBase<D1> const & design,
        Eigen::MatrixBase<D2> const & data
    ) {
        _getDesignMatrix() = design;
        _getDataVector() = data;
        _factor(false);
    }

    /// @brief Reset the design matrix given as an ndarray; dimension and data are not changed.
    template <typename T1, int C1>
    void setDesignMatrix(ndarray::Array<T1 const,2,C1> const & design) {
        _getDesignMatrix() = design.asEigen();
        _factor(false);
    }

    /// @brief Reset the design matrix given as an Eigen object; dimension and data are not changed.
    template <typename D1, typename D2>
    void setDesignMatrix(Eigen::MatrixBase<D1> const & design) {
        _getDesignMatrix() = design;
        _factor(false);
    }
    

    /// @brief Initialize from the terms in the normal equations, given as ndarrays.
    template <typename T1, typename T2, int C1, int C2>
    static LeastSquares fromNormalEquations(
        ndarray::Array<T1 const,2,C1> const & fisher,
        ndarray::Array<T2 const,1,C2> const & rhs,
        Factorization factorization = NORMAL_EIGENSYSTEM
    ) {
        LeastSquares r(factorization, fisher.template getSize<0>());
        if ((C1 > 0) == Eigen::MatrixXd::IsRowMajor)
            r._getFisherMatrix() = fisher.asEigen();
        else
            r._getFisherMatrix() = fisher.asEigen().transpose();
        r._getRhsVector() = rhs.asEigen();
        r._factor(true);
        return r;
    }

    /// @brief Initialize from the terms in the normal equations, given as Eigen objects.
    template <typename D1, typename D2>
    static LeastSquares fromNormalEquations(
        Eigen::MatrixBase<D1> const & fisher,
        Eigen::MatrixBase<D2> const & rhs,
        Factorization factorization = NORMAL_EIGENSYSTEM
    ) {
        LeastSquares r(factorization, fisher.rows());
        if (Eigen::MatrixBase<D1>::isRowMajor == Eigen::MatrixXd::IsRowMajor)
            r._getFisherMatrix() = fisher;
        else
            r._getFisherMatrix() = fisher.transpose();
        r._getRhsVector() = rhs;
        r._factor(true);
        return r;
    }

    /// @brief Reset the terms in the normal equations given as ndarrays; dimension must not change.
    template <typename T1, typename T2, int C1, int C2>
    void setNormalEquations(
        ndarray::Array<T1 const,2,C1> const & fisher,
        ndarray::Array<T2 const,1,C2> const & rhs
    ) {
        if ((C1 > 0) == Eigen::MatrixXd::IsRowMajor)
            _getFisherMatrix() = fisher.asEigen();
        else
            _getFisherMatrix() = fisher.asEigen().transpose();
        _getRhsVector() = rhs.asEigen();
        _factor(true);
    }

    /// @brief Reset the terms in the normal equations given as Eigen objects; dimension must not change.
    template <typename D1, typename D2>
    void setNormalEquations(
        Eigen::MatrixBase<D1> const & fisher,
        Eigen::MatrixBase<D2> const & rhs
    ) {
        if (Eigen::MatrixBase<D1>::isRowMajor == Eigen::MatrixXd::IsRowMajor)
            _getFisherMatrix() = fisher;
        else
            _getFisherMatrix() = fisher.transpose();
        _getRhsVector() = rhs;
        _factor(true);        
    }

    /**
     *  @brief Set the threshold used to determine when to truncate Eigenvalues.
     *
     *  The rank of the matrix is determined by comparing the product of this threshold
     *  and the first (largest) element of the array returned by getCondition() to all other
     *  elements of that array.  Elements smaller than this product are ignored and reduce
     *  the rank.
     *
     *  In the DIRECT_SVD case, the condition array contains the singular values of the design
     *  matrix, while in the NORMAL_EIGENSYSTEM case the condition array holds the Eigenvalues
     *  of the Fisher matrix, and the latter are the square of the former.  The default
     *  threshold is the same (std::numeric_limits<double>::epsilon()) in both cases,
     *  reflecting the fact that using the normal equations squares the condition number
     *  of the problem.
     *
     *  The NORMAL_CHOLESKY method does not use the threshold and assumes the problem is 
     *  full-rank.
     */
    void setThreshold(double threshold);

    /// @brief Get the threshold used to determine when to truncate Eigenvalues.
    double getThreshold() const;

    /**
     *  @brief Return the vector solution to the least squares problem.
     *
     *  The returned array is owned by the LeastSquares object and may be modified in-place
     *  by future calls to LeastSquares member functions, so it's best to promptly copy the
     *  result elsewhere.
     *
     *  If you want an Eigen object instead, just use solve().asEigen().
     */
    ndarray::Array<double const,1,1> solve();

    /**
     *  @brief Return the covariance matrix of the least squares problem.
     *
     *  The returned array is owned by the LeastSquares object and may be modified in-place
     *  by future calls to LeastSquares member functions, so it's best to promptly copy the
     *  result elsewhere.
     *
     *  If you want an Eigen object instead, just use computeCovariance().asEigen().
     */
    ndarray::Array<double const,2,2> computeCovariance();

    /**
     *  @brief Return the Fisher matrix (inverse of the covariance) of the parameters.
     *
     *  Note that the Fisher matrix is exactly the same as the matrix on the lhs of the
     *  normal equations.
     *
     *  The returned array is owned by the LeastSquares object and may be modified in-place
     *  by future calls to LeastSquares member functions, so it's best to promptly copy the
     *  result elsewhere.
     *
     *  If you want an Eigen object instead, just use computeFisherMatrix().asEigen().
     */
    ndarray::Array<double const,2,2> computeFisherMatrix();

    /**
     *  @brief Return a factorization-dependent vector that can be used to characterize
     *         the stability of the solution.
     *
     *  The returned array's size is always equal to getDimension().  When getRank() is
     *  less than getDimension(), some elements of the array were considered effectively
     *  zero (see setThreshold).
     *
     *  For the NORMAL_EIGENSYSTEM method, this is the vector of Eigenvalues of the Fisher
     *  matrix, including those rejected as being below the threshold.
     *
     *  For the DIRECT_SVD method, this is the vector of singular values of the design
     *  matrix, including those rejected as being below the threshold; with exact arithmetic,
     *  these would be the square roots of the Eigenvalues of the Fisher matrix.
     *
     *  For the NORMAL_CHOLESKY method, this is @f$D@f$ in the pivoted Cholesky factorization
     *  @f$P L D L^T P^T@f$ of the Fisher matrix.  This does not provide a reliable way to 
     *  test the stability of the problem, but it does provide a way to compute the determinant
     *  of the Fisher matrix.
     */
    ndarray::Array<double const,1,1> getCondition();

    /// @brief Return the number of parameters.
    int getDimension() const;

    /**
     *  @brief Return the rank of the problem (number of nonzero Eigenvalues).
     *
     *  The returned value is always the same as getDimension() when the factorization is NORMAL_CHOLESKY
     *  (which may be incorrect, because a Cholesky decomposition is not rank-revealing).
     */
    int getRank() const;

    // Need to define dtor in source file so it can see Impl declaration.
    ~LeastSquares();
    
private:

    // We want a column-major design matrix so the self-adjoint product is cache-friendly, hence '-2'...
    // so we always copy a (possibly row-major) design matrix into a col-major one.  This is an
    // unnecessarily and cache-unfriendly operation when solver is DIRECT_SVD, but right now it doesn't
    // seem to be worth special-casing the design for that case.  In other cases it's a cache-unfriendly
    // op that avoids an even worse one, and it can always be avoided by using a column-major design matrix.
    Eigen::MatrixXd & _getDesignMatrix();
    Eigen::VectorXd & _getDataVector();

    // Storage order matters less here, so we just go with what Eigen is most accustomed to.
    // Because the Fisher matrix is symmetric, we can just transpose it before copying to avoid doing
    // any expensive copies between different storage orders.
    Eigen::MatrixXd & _getFisherMatrix();
    Eigen::VectorXd & _getRhsVector();

    void _factor(bool haveNormalEquations);

    LeastSquares(Factorization factorization, int dimension);

    PTR(Impl) _impl;
};

}}} // namespace lsst::afw::math

#endif // !LSST_AFW_MATH_LeastSquares_h_INCLUDED
