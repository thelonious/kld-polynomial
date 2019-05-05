"use strict";

const assert = require("assert");
const {SqrtPolynomial} = require("../index");

describe("SqrtPolynomial", () => {
    it("toString", () => {
        const poly = new SqrtPolynomial(2, 1, 0);

        assert.strictEqual(poly.toString(), "sqrt(2t^2 + t)");
    });
});
