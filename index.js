'use strict';

const _Math = require('minimatrix');
const _Decomp = require('minimatrix-decomp');
const TOLERANCE = 1e-14;

const helpers = {
  vectorCoordMap: {
    0: 'x',
    1: 'y',
    2: 'z'
  },
  applyValuesToFunctionArray: (F, x, y) => {
    if (y) {
      for (let i = 0; i < F.length; ++i) {
        y.setComponent(i, (F[i])(x));
      }
      return y;
    } else {
      return F.map(f => f(x));
    }
  },
  getPartialDerivatives: (f, x, delta = 0.01) => {
    const del = x.clone().setScalar(0);
    const n = x.dimension;
    const d = delta;
    const d2 = 2 * delta;
    for (let i = 0; i < n; ++i) {
      const coord = helpers.vectorCoordMap[i];
      x[coord] += d;
      const dFdxip1 = f(x);
      x[coord] -= d2;
      const dFdxim1 = f(x);
      x[coord] += d;
      del.setComponent(i, (dFdxip1 - dFdxim1) / d2);
    }
    return del;
  },
  getSquareMatrix: (n) => {
    if (n === 1) {
      return undefined;
    } else if (n === 2) {
      return new _Math.Matrix2();
    } else if (n === 3) {
      return new _Math.Matrix3();
    } else {
      throw new Error('getSquareMatrix(): No matrix for bounds.');
    }
  },
  getJacobian: (F, x, JC) => {
    if (F.length !== x.dimension) {
      throw new Error('getJacobian(): Jacobian will not be square.');
    }
    const J = JC || helpers.getSquareMatrix(F.length);
    const n = F.length;
    for (let i = 0; i < n; ++i) {
      const row = helpers.getPartialDerivatives(F[i], x);
      J.setRow(i, row);
    }
    return J;
  },
  luDecompAndSolve: () => {
    throw new Error('Not implemented yet.');
  }
};

const RootFinders = {
  bisectionMethod: function bisectionMethod (f, options) {
    const imax = options.maxIterations;
    const es = options.rootTolerance || TOLERANCE;
    let xl = options.lowerBound;
    let xu = options.upperBound;
    let xr = options.initialValue;

    let iter = 0;
    let fl = f(xl);
    let fr = 0;
    do {
      xr = (xl + xu) / 2;
      fr = f(xr);
      iter++;
      const test = fl * fr;
      if (test < 0) {
        xu = xr;
      } else if (test > 0) {
        xl = xr;
        fl = fr;
      }
    } while (Math.abs(fr) > es && iter < imax);
    if (iter >= imax) {
      console.warn('bisectionMethod(): Iteration max reached. Solution may not be accurate.');
    }
    return {
      root: xr,
      iterations: iter
    };
  },
  regulaFalsiMethod: function regulaFalsiMethod (f, options) {
    const imax = options.maxIterations;
    const es = options.rootTolerance || TOLERANCE;
    let xl = options.lowerBound;
    let xu = options.upperBound;
    let xr = options.initialValue;
    let iter = 0;
    let fl = f(xl);
    let fu = f(xu);
    let ea = Number.POSITIVE_INFINITY;
    let iu = 0;
    let il = 0;
    let fr = 0;
    do {
      const xrold = xr;
      xr = xu - fu * (xl - xu) / (fl - fu);
      fr = f(xr);
      iter++;
      if (xr !== 0) {
        ea = Math.abs((xr - xrold) / xr) * 100;
      }
      const test = fl * fr;
      if (test < 0) {
        xu = xr;
        fu = f(xu);
        iu = 0;
        il++;
        if (il >= 2) {
          fl /= 2;
        }
      } else if (test > 0) {
        xl = xr;
        fl = f(xl);
        il = 0;
        iu++;
        if (iu >= 2) {
          fu /= 2;
        }
      } else {
        ea = 0;
      }
    } while (ea > es && Math.abs(fr) > es && iter < imax);
    return {
      root: xr,
      iterations: iter
    };
  },
  /**
   * Newton-Raphson method for nonlinear root finding
   *
   * @param {Function[]} functions - The array of functions for which we are simultaneously solving (equals the number of variables)
   * @param {Object} options - parameter list
   */
  newtonsMethod: function newtonsMethod (F, options) {
    const imax = options.maxIterations;
    const es = options.rootTolerance || TOLERANCE;
    const DF = options.derivativeMatrix;
    let x0 = options.initialValue;
    let iter = 0;

    let n = x0.dimension;
    if (n !== F.length) {
      throw new Error('Under- and over-constrained systems not yet supported.');
    }
    const haveAnalyticalDerivs = (Array.isArray(DF)) &&
      DF.reduce((acc, d) => acc && typeof d === 'function', true);

    if (n === 1) {
      // single function
    } else if (n === 2 || n === 3) {
      // multivariate function
      const xi = x0.clone();
      // const dx = new _Math.Vector2();
      const Fx = x0.clone().setScalar(0);
      const Jx = (n === 2 ? new _Math.Matrix2() : new _Math.Matrix3());
      helpers.applyValuesToFunctionArray(F, xi, Fx);

      // Begin Newton-Raphson
      while (Fx.length() > es && iter < imax) {
        ++iter;
        if (haveAnalyticalDerivs) {
          for (let i = 0; i < n; ++i) {
            for (let j = 0; j < n; ++j) {
              Jx.elements[j * n + i] = DF[i * n + j](xi);
            }
          }
        } else {
          helpers.getJacobian(F, xi, Jx);
        }
        const jInfo = _Decomp.luDecomposition(Jx);
        const dx = _Decomp.luSolve(jInfo.A, jInfo.P, Fx.negate());
        // helpers.luDecompAndSolve(Jx, Fx.negate(), dx);
        xi.add(dx);
        helpers.applyValuesToFunctionArray(F, xi, Fx);
      }
      return {
        root: xi,
        iterations: iter
      };
    } else {
      throw new Error('newtonsMethod(): Too many functions.');
    }
  }
};

module.exports = RootFinders;
