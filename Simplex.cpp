// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <cmath>
#include <limits>

#include "Simplex.hpp"
#include "Factorial.hpp"
#include "lapack_names.h"


Simplex::Simplex(size_t npdims, size_t ndims) :
 mnpdims(npdims), 
 mndims(ndims) {
  // default coords are the canonical positions
  mvertices.push_back( std::vector<double>(ndims, 0.0) );
  for (size_t i = 0; i < npdims; ++i) {
    std::vector<double> pt(ndims, 0.0);
    pt[i] = 1.0; // canonical coordinate
    mvertices.push_back(pt);
  }
}

void 
Simplex::setPoint(size_t index, const std::vector<double>& pt) {
    mvertices[index] = pt;
}

bool
Simplex::getParametricCoordinates(const std::vector<double>& position, 
                                  double tol, std::vector<double>& xiCoords) {

  // set the right hand side
  std::vector<double> rhs(mndims);
    for (size_t i = 0; i < mndims; ++i) {
        rhs[i] = position[i] - mvertices[0][i];
    }

    int ier = this->solve(rhs);
    if (ier >= 0) {
    // success
    bool res = true;
    double sum = 0;
    std::copy(rhs.begin(), rhs.begin() + mnpdims, xiCoords.begin());
    for (size_t i = 0; i < mnpdims; ++i) {
      res &= xiCoords[i] > 0.0 - tol && xiCoords[i] <= 1.0 + tol;
      sum += xiCoords[i];
    }
    res &= sum < 1.0 + tol;
        return res;
    }
    else {
        // error
        return false;
    }

}

double
Simplex::getVolume() const {

  // See http://www.math.niu.edu/~rusin/known-math/97/volumes.polyh

  double res = 0.0;
  if (mndims == 0) return res;
  if (mnpdims == 0) return res;

  // all vertices have been set
  std::vector< std::vector<double> > wmat(mndims, std::vector<double>(mnpdims));
  for (size_t jcol = 0; jcol < mnpdims; ++jcol) {
    for (size_t irow = 0; irow < mndims; ++irow) {
      wmat[irow][jcol] = mvertices[jcol + 1][irow] - mvertices[0][irow];
    }
  }
  // transpose(wmat) * wmat
  std::vector<double> wmatTDotWmat(mnpdims * mnpdims, 0.0);
  for (size_t jcol = 0; jcol < mnpdims; ++jcol) {
    for (size_t irow = 0; irow < mnpdims; ++irow) {
      size_t k = irow + mnpdims*jcol;
      for (size_t el = 0; el < mndims; ++el) {
        wmatTDotWmat[k] += wmat[el][irow] * wmat[el][jcol];
      }
    }
  }

  double det = this->getDeterminant(wmatTDotWmat);
  res = det / factorial(mnpdims);
  return res;
}

double
Simplex::getDeterminant(const std::vector<double>& amat) const {

    double res = 1;
    size_t nSquare = amat.size();
    int n = (int) std::sqrt( (double) nSquare );
    std::vector<double> mat(amat); // make a copy

    // LU decomposition
      std::vector<int> ipiv((size_t) n);
      int err = 0;
      _GETRF_(&n, &n, &mat.front(), &n, &ipiv.front(), &err);
      // no error checking

      // compute the determinant
      for (size_t i = 0; i < (size_t)n; ++i) {
        double aii = mat[i + (size_t)n*i];
        // ipiv holds the fortran indices
        if (ipiv[i] != (int)i + 1) {
              res *= -aii;
        }
        else {
              res *= aii;
        }
      }
      return res;

}

int
Simplex::solve(std::vector<double>& rhs) const {

  // save the input matrix and rhs vector so we get a chance to
  // recover in case of a singular matrix.

  const double residualTol = 1.e-8;

  std::vector<double> mat(mnpdims * mndims);
  // build the matrix system
  for (size_t jcol = 0; jcol < mnpdims; ++jcol) {
    for (size_t irow = 0; irow < mndims; ++irow) {
      mat[irow + jcol*mndims] = mvertices[jcol + 1][irow] - mvertices[0][irow];
    }
  }

  std::vector<double> matCopy(mnpdims * mndims);
  std::vector<double> bCopy(rhs.size());
  std::copy(mat.begin(), mat.end(), matCopy.begin());
  std::copy(rhs.begin(), rhs.end(), bCopy.begin());

  int nrow = (int) mndims;
  int ncol = (int) mnpdims;
  char t = 'n';
  int one = 1;
  int mn = ncol;
  int nb = 1; // optimal block size
  int lwork = mn + mn*nb;
  std::vector<double> work((size_t) lwork);
  int errCode = 0;
  _GELS_(&t, &nrow, &ncol, &one,
     &matCopy.front(), &nrow,
     &rhs.front(), &nrow,
     &work.front(), &lwork, &errCode);

  if (!errCode) {
    // merrCode == 0 indicates everything was fine
    // merrCode < 0 indicates bad entry
    // in either case return
    return errCode;
  }
  else if (errCode > 0) {

    if (nrow <= 1) {
      return errCode;
    }

    for (size_t i = 0; i < mndims; ++i) {
      rhs[i] = bCopy[i];
    }
    errCode = 0;
    // relative accuracy in the matrix data
    double rcond = std::numeric_limits<double>::epsilon();
    int rank;
    std::vector<int> jpvt(ncol);
    lwork = mn + 3*ncol + 1 > 2*mn + nb*1? mn + 3*ncol + 1: 2*mn + nb*1;
    work.resize(lwork);
    _GELSY_(&nrow, &ncol, &one,
        &matCopy.front(), &nrow,
        &rhs.front(), &nrow,
        &jpvt.front(), &rcond, &rank, &work.front(), &lwork, &errCode);

    // check if this is a good solution
    double residualError = 0.0;
    for (size_t i = 0; i < mndims; ++i) {
      double rowSum = 0.0;
      for (size_t j = 0; j < mnpdims; ++j) {
        rowSum += mat[i + nrow*j]*rhs[j];
      }
      residualError += std::abs(rowSum - bCopy[i]);
    }
    if (residualError < residualTol) {
      // good enough
      errCode = 1;
      return errCode;
    }

    // some error
    errCode = 2;
    return errCode;
  }

  // we should never reach that point
  errCode = 0;
  return errCode;
}
