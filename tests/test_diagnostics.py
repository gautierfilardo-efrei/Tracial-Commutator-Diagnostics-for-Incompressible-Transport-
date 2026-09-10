"""Independent identities and counterexamples, not tests of fitted data."""
import unittest
import numpy as np
from diagnostics import (cube_modes, mode_count, energy_coefficient,
                         cube_increment, cube_increment_direct,
                         cube_integrand, universal_bound, abc_flow, mean_se)


class MathematicalChecks(unittest.TestCase):
    def test_dirichlet_formula_against_explicit_sum(self):
        delta = np.random.default_rng(17).uniform(-1, 1, (31, 3))
        for k in (1, 2, 4):
            np.testing.assert_allclose(cube_increment(delta, k),
                                       cube_increment_direct(delta, k), rtol=2e-14, atol=1e-11)

    def test_identity_and_integer_translation(self):
        np.testing.assert_allclose(cube_increment(np.zeros((3, 3)), 4), 0)
        np.testing.assert_allclose(cube_increment(np.ones((3, 3)), 4), 0, atol=1e-20)

    def test_small_displacement_stability(self):
        d = np.array([[1e-12, -2e-12, 3e-12]])
        expected = 4 * np.pi**2 * energy_coefficient(4) * np.sum(d*d)
        self.assertAlmostEqual(float(cube_increment(d, 4)[0]) / expected, 1, places=12)

    def test_universal_bound(self):
        delta = np.random.default_rng(21).uniform(-0.5, 0.5, (501, 3))
        for k in (1, 2, 4, 8):
            s = cube_integrand(delta, k)
            self.assertTrue(np.all((s >= 0) & (s <= universal_bound(k) + 1e-12)))

    def test_isotropic_fourier_moment(self):
        for k in (1, 2, 4):
            modes = cube_modes(k)
            np.testing.assert_allclose(modes.T @ modes,
                                       energy_coefficient(k) * np.eye(3))
            self.assertEqual(len(modes), mode_count(k))

    def test_cat_map_closed_form(self):
        t = (np.arange(2**15) + 0.5) / 2**15
        value = np.mean(np.log(3 - 2 * np.cos(2 * np.pi * t)))
        self.assertAlmostEqual(value, np.log((3 + np.sqrt(5))/2), places=13)

    def test_shear_frequency_blindness(self):
        # n odd permutes this power-of-two midpoint grid; no fitted tolerance.
        t = (np.arange(2**12) + 0.5) / 2**12
        values = []
        for n in (1, 3, 17, 63):
            d = np.zeros((len(t), 3)); d[:, 0] = 0.2 * np.sin(2*np.pi*n*t)
            values.append(cube_integrand(d, 4).mean())
        np.testing.assert_allclose(values, values[0], atol=2e-12, rtol=0)

    def test_false_cocycle_counterexample(self):
        e11=np.array([[1,0],[0,0]]);e12=np.array([[0,1],[0,0]]);e21=e12.T
        phi=lambda a,b,c: np.trace(a@(b@c-c@b))/2
        a,b,c,d=e11,e11,e12,e21
        cob=phi(a@b,c,d)-phi(a,b@c,d)+phi(a,b,c@d)-phi(d@a,b,c)
        self.assertEqual(cob, -0.5)

    def test_abc_rk4_step_refinement(self):
        # Step refinement on a nontrivial particle ODE: checks integration error.
        x=np.random.default_rng(29).random((64,3))
        y1=abc_flow(x,0.2,16);y2=abc_flow(x,0.2,32);y3=abc_flow(x,0.2,64)
        coarse=np.max(np.abs(y1-y2));fine=np.max(np.abs(y2-y3))
        self.assertGreater(coarse/fine,10)
        self.assertLess(fine,1e-7)

    def test_invalid_fourier_parameters_are_rejected(self):
        for invalid in (0, -1, 2.5, 2.0, True, np.nan):
            for operation in (mode_count, cube_modes, energy_coefficient):
                with self.subTest(value=invalid, operation=operation.__name__):
                    with self.assertRaises(ValueError):
                        operation(invalid)
        with self.assertRaises(ValueError):
            mode_count(2, 2.5)
        with self.assertRaises(ValueError):
            cube_increment_direct(np.zeros((1,3)), 2, block_size=0)

    def test_nonfinite_or_malformed_inputs_are_rejected(self):
        for points in (np.array([1.,2.,3.]), np.empty((2,0)),
                       np.array([[np.nan,0,0]]), np.array([[np.inf,0,0]])):
            for operation in (cube_increment, cube_increment_direct):
                with self.subTest(operation=operation.__name__, shape=points.shape):
                    with self.assertRaises(ValueError):
                        operation(points, 2)
        for values in ([1.], [1., np.nan], [[1.,2.],[3.,4.]]):
            with self.assertRaises(ValueError):
                mean_se(values)

    def test_invalid_flow_parameters_are_rejected(self):
        points = np.zeros((2,3))
        for kwargs in ({'steps':1.5}, {'steps':True}, {'viscosity':-1},
                       {'lag':np.nan}, {'start_time':np.inf}):
            params = {'lag':0.2, **kwargs}
            with self.subTest(parameters=params):
                with self.assertRaises(ValueError):
                    abc_flow(points, **params)
        with self.assertRaises(ValueError):
            abc_flow(np.zeros((2,2)), 0.2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
