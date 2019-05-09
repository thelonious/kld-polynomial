import assert from "assert";
import {Polynomial} from "../index.js";

/**
 * Assert that two polynomials are the same size and have the same valued coefficients
 * @param {Polynomial} poly1
 * @param {Polynomial} poly2
 */
function assertEqualPolynomials(poly1, poly2) {
    assert.strictEqual(poly1.getDegree(), poly2.getDegree());

    for (let i = 0; i <= poly1.getDegree(); i++) {
        assert.strictEqual(poly1.coefs[i], poly2.coefs[i]);
    }
}

/**
 * Test that two floats are equal withing a tolerance
 * @param {number} a
 * @param {number} b
 * @param {number} TOLERANCE
 */
function assertEqualWithinTolerance(a, b, TOLERANCE = 1e-12) {
    assert(Math.abs(a - b) < TOLERANCE, `${a} === ${b}`);
}

describe("Polynomial", () => {
    describe("Text Representation", () => {
        it("Default variable name", () => {
            const poly = new Polynomial(2, 1, 0);

            assert.strictEqual(poly.toString(), "2t^2 + t");
        });
        it("Variable name as 's'", () => {
            const poly = new Polynomial(2, 1, 0);

            poly._variable = "s";
            assert.strictEqual(poly.toString(), "2s^2 + s");
        });
    });
    describe("Operations", () => {
        it("clone", () => {
            const poly = new Polynomial(2, 1, 0);
            const copy = poly.clone();

            assertEqualPolynomials(poly, copy);
        });
        it("eval", () => {
            const poly = new Polynomial(2, 1, 0);

            for (let x = 0; x < 100; x++) {
                assert.strictEqual(2 * x * x + x, poly.eval(x));
            }
        });
        it("add", () => {
            const poly1 = new Polynomial(2, 1, 0);
            const poly2 = new Polynomial(3, 5, 7, 9);
            const expected = new Polynomial(3, 7, 8, 9);
            const result = poly1.add(poly2);

            assertEqualPolynomials(result, expected);
        });
        it("multiply", () => {
            const poly1 = new Polynomial(2, 1, 0);
            const poly2 = new Polynomial(1, 3, 5);
            const expected = new Polynomial(2, 7, 13, 5, 0);
            const result = poly1.multiply(poly2);

            assertEqualPolynomials(result, expected);
        });
        it("divideEqualsScalar", () => {
            const poly = new Polynomial(6, 4, 2);
            const expected = new Polynomial(3, 2, 1);

            poly.divideEqualsScalar(2);

            assertEqualPolynomials(poly, expected);
        });
        it("simplifyEquals", () => {
            const poly = new Polynomial(1e-13, 4, 2);
            const expected = new Polynomial(4, 2);

            poly.simplifyEquals();

            assertEqualPolynomials(poly, expected);
        });
        it("removeZeroEquals", () => {
            const poly = new Polynomial(1e-15, 2, 1e-16, 4, 1e-17);
            const expected = new Polynomial(1e-15, 2, 0, 4, 0);

            poly.removeZerosEquals();

            assertEqualPolynomials(poly, expected);
        });
        it("monicEquals", () => {
            const poly = new Polynomial(2, 4, 6, 8);
            const expected = new Polynomial(1, 2, 3, 4);

            poly.monicEquals();

            assertEqualPolynomials(poly, expected);
        });
        it("derivative", () => {
            const poly = new Polynomial(2, 4, 6, 8);
            const result = poly.getDerivative();
            const expected = new Polynomial(6, 8, 6);

            assertEqualPolynomials(result, expected);
        });
    });
    describe("Roots", () => {
        it("linear", () => {
            const poly = new Polynomial(10, -5);
            const roots = poly.getRoots();

            assert.strictEqual(roots.length, 1);
            assert.strictEqual(roots[0], 0.5);
        });
        it("quadratic", () => {
            const poly1 = new Polynomial(1, -5);
            const poly2 = new Polynomial(1, -6);
            const poly = poly1.multiply(poly2);
            const roots = poly.getRoots().sort();

            assert.strictEqual(roots.length, 2);
            assert.strictEqual(roots[0], 5);
            assert.strictEqual(roots[1], 6);
        });
        it("cubic", () => {
            const poly1 = new Polynomial(1, -5);
            const poly2 = new Polynomial(1, -6);
            const poly3 = new Polynomial(1, -7);
            const poly = poly1.multiply(poly2).multiply(poly3);
            const roots = poly.getRoots().sort();

            assert.strictEqual(roots.length, 3);
            assert.strictEqual(roots[0], 5);
            assert.strictEqual(roots[1], 6);
            assert.strictEqual(roots[2], 7);
        });
        it("quartic", () => {
            const poly1 = new Polynomial(1, -5);
            const poly2 = new Polynomial(1, -6);
            const poly3 = new Polynomial(1, -7);
            const poly4 = new Polynomial(1, -8);
            const poly = poly1.multiply(poly2).multiply(poly3).multiply(poly4);
            const roots = poly.getRoots().sort();

            assert.strictEqual(roots.length, 4);
            assertEqualWithinTolerance(roots[0], 5);
            assertEqualWithinTolerance(roots[1], 6);
            assertEqualWithinTolerance(roots[2], 7);
            assertEqualWithinTolerance(roots[3], 8);
        });
        it("bisection", () => {
            const poly = new Polynomial(1, -0.25);
            /* eslint-disable-next-line no-shadow */
            const root = poly.bisection(0, 1);

            assertEqualWithinTolerance(root, 0.25);
        });
        it("roots in interval", () => {
            const poly1 = new Polynomial(1, -0.25);
            const poly = new Polynomial(1, -0.75).multiply(poly1);
            const roots = poly.getRootsInInterval(0, 1).sort();

            assert.strictEqual(roots.length, 2);
            assertEqualWithinTolerance(roots[0], 0.25);
            assertEqualWithinTolerance(roots[1], 0.75);
        });
    });
});
