import assert from "assert";
import {Polynomial} from "../index.js";

describe("Polynomial", () => {
    it("toString", () => {
        const poly = new Polynomial(2, 1, 0);

        assert.strictEqual(poly.toString(), "2t^2 + t");
    });
});
