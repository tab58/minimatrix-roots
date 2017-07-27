'use strict';
/* globals describe it */
const RootFinders = require('../index.js');
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
});
