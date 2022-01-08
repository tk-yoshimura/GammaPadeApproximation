using MultiPrecision;
using System;
using System.Collections.Generic;
using System.IO;

namespace GammaPade {
    internal class Program {
        static void Main(string[] args) {
            using (StreamWriter sw = new("../../../../results/gamma_pade_table_e31.csv")) {
                for (MultiPrecision<Pow2.N16> x = 1; x <= 16; x += 1) {
                    bool is_e30 = false;

                    MultiPrecision<Pow2.N16> ddx = Math.ScaleB(1, -1);

                    MultiPrecision<Pow2.N32>[] ds = FiniteDifference<Pow2.N32>.Diff(
                        x.Convert<Pow2.N32>(), MultiPrecision<Pow2.N32>.Gamma, Math.ScaleB(1, -24)
                    );

                    List<MultiPrecision<Pow2.N16>> expecteds = new();

                    for (MultiPrecision<Pow2.N16> dx = -ddx, h = ddx / 4096; dx <= ddx; dx += h) {
                        MultiPrecision<Pow2.N16> expected = MultiPrecision<Pow2.N16>.Gamma(x + dx);

                        expecteds.Add(expected);
                    }

                    for (int n = 4; n <= 16; n += 1) {
                        MultiPrecision<Pow2.N16>[] cs = new MultiPrecision<Pow2.N16>[n * 2 + 1];
                        cs[0] = MultiPrecision<Pow2.N16>.Gamma(x);
                        for (int i = 0; i < n * 2; i++) {
                            cs[i + 1] = ds[i].Convert<Pow2.N16>() * MultiPrecision<Pow2.N16>.TaylorSequence[i + 1];
                        }

                        (MultiPrecision<Pow2.N16>[] ms, MultiPrecision<Pow2.N16>[] ns) = PadeSolver<Pow2.N16>.Solve(cs, n, n);

                        MultiPrecision<Pow2.N16> err = 0;

                        for ((MultiPrecision<Pow2.N16> dx, MultiPrecision<Pow2.N16> h, int i) = (-ddx, ddx / 4096, 0); i < expecteds.Count; dx += h, i++) {
                            MultiPrecision<Pow2.N16> expected = expecteds[i];
                            MultiPrecision<Pow2.N16> actual = PadeSolver<Pow2.N16>.Approx(dx, ms, ns);

                            err = MultiPrecision<Pow2.N16>.Max(err, MultiPrecision<Pow2.N16>.Abs(expected / actual - 1));
                        }

                        Console.WriteLine($"x={x}, n={n}, |dx| = {ddx}");
                        Console.WriteLine($"relative error = {err:e10}");

                        if (err < 2e-31) {
                            sw.WriteLine($"x={x}, n={n}, |dx| = {ddx}");

                            sw.WriteLine($"i,p_i,q_i");
                            for (int i = 0; i <= n; i++) {
                                sw.WriteLine($"{i},{ms[i]:e64},{ns[i]:e64}");
                            }

                            sw.WriteLine($"relative error = {err:e10}\n");

                            is_e30 = true;

                            sw.Flush();

                            break;
                        }
                    }

                    if (!is_e30) {
                        sw.WriteLine($"{x} not convergence");
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
