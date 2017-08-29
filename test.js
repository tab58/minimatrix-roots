'use strict';
/* globals describe it */
const RootFinders = require('./index.js');
const _Math = require('minimatrix');
const expect = require('chai').expect;

describe('Root Finding', () => {
  describe('Bisection Method', () => {
    it('should find the root', () => {
      const f = (c) => 9.8 * 68.1 / c * (1 - Math.exp(c / 68.1 * -10)) - 40;
      const TOL = 1e-14;
      const rootInfo = RootFinders.bisectionMethod(f, {
        maxIterations: 1024,
        rootTolerance: TOL,
        lowerBound: 12,
        upperBound: 16,
        initialValue: 14
      });
      expect(Math.abs(f(rootInfo.root)) < TOL);
    });
  });
  describe('Regula Falsi Method', () => {
    const f = (c) => 9.8 * 68.1 / c * (1 - Math.exp(c / 68.1 * -10)) - 40;
    const TOL = 1e-14;
    it('should find the root', () => {
      // const imax = options.maxIterations;
      // const es = options.rootTolerance || TOLERANCE;
      // const DF = options.derivativeMatrix;
      // let x0 = options.initialValue;
      const rootOptions = {
        maxIterations: 1024,
        rootTolerance: TOL,
        lowerBound: 12,
        upperBound: 16,
        initialValue: 14
      };
      const bInfo = RootFinders.bisectionMethod(f, rootOptions);
      const rfInfo = RootFinders.regulaFalsiMethod(f, rootOptions);
      expect(Math.abs(f(rfInfo.root)) < TOL);
      expect(Math.abs(rfInfo.root - bInfo.root) < TOL);
    });
    it('should find the root quicker than bisection', () => {
      const rootOptions = {
        maxIterations: 1024,
        rootTolerance: TOL,
        lowerBound: 12,
        upperBound: 16,
        initialValue: 14
      };
      const bInfo = RootFinders.bisectionMethod(f, rootOptions);
      const rfInfo = RootFinders.regulaFalsiMethod(f, rootOptions);
      expect(bInfo.iterations).to.be.greaterThan(rfInfo.iterations);
    });
  });
  describe('Newton-Raphson Method', () => {
    const f1 = (v) => v.x * v.x - 2 * v.x + v.y * v.y - v.z + 1;
    const f2 = (v) => v.x * v.y * v.y - v.x - 3 * v.y + v.y * v.z + 2;
    const f3 = (v) => v.x * v.z * v.z - 3 * v.z + v.y * v.z * v.z + v.x * v.y;
    const F = [f1, f2, f3];
    const TOL = 1e-14;
    it('should find the root', () => {
      const rootOptions1 = {
        maxIterations: 1024,
        rootTolerance: TOL,
        initialValue: new _Math.Vector3(1, 2, 3),
        derivativeMatrix: null
      };
      const nInfo1 = RootFinders.newtonsMethod(F, rootOptions1);
      const r1 = nInfo1.root;
      expect(Math.abs(f1(r1)) < TOL);
      expect(Math.abs(f2(r1)) < TOL);
      expect(Math.abs(f3(r1)) < TOL);

      const f1x = (v) => 2 * v.x - 2;
      const f1y = (v) => 2 * v.y;
      const f1z = (v) => -1;
      const f2x = (v) => v.y * v.y - 1;
      const f2y = (v) => 2 * v.x * v.y - 3 + v.z;
      const f2z = (v) => v.y;
      const f3x = (v) => v.z * v.z + v.y;
      const f3y = (v) => v.z * v.z + v.x;
      const f3z = (v) => 2 * v.x * v.z - 3 + 2 * v.y * v.z;
      const rootOptions2 = {
        maxIterations: 1024,
        rootTolerance: TOL,
        initialValue: new _Math.Vector3(1, 2, 3),
        derivativeMatrix: [
          f1x, f1y, f1z,
          f2x, f2y, f2z,
          f3x, f3y, f3z
        ]
      };
      const nInfo2 = RootFinders.newtonsMethod(F, rootOptions2);
      const r2 = nInfo2.root;
      expect(Math.abs(f1(r2)) < TOL);
      expect(Math.abs(f2(r2)) < TOL);
      expect(Math.abs(f3(r2)) < TOL);
    });
  });
});
