'use strict';

const TOLERANCE = 1e-14;

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
  }
};

try {
  const Minimatrix = require('minimatrix');
  Minimatrix.RootFinders = RootFinders;
} catch (e) {
  module.exports = RootFinders;
}
