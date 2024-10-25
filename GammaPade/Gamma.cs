using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaPade {
    internal class Gamma {
        static void Main_(string[] args) {
            for (int exponent = -950; exponent <= -50; exponent += 25) {
                MultiPrecision<Pow2.N32> x = MultiPrecision<Pow2.N32>.Ldexp(1, exponent);
                MultiPrecision<Pow2.N32> y = 1 / MultiPrecision<Pow2.N32>.Gamma(x);

                Console.WriteLine($"{y:e40}");
            }

            for (int exponent = -950; exponent <= -50; exponent += 25) {
                MultiPrecision<Pow2.N32> x = MultiPrecision<Pow2.N32>.Ldexp(-1, exponent);
                MultiPrecision<Pow2.N32> y = 1 / MultiPrecision<Pow2.N32>.Gamma(x);

                Console.WriteLine($"{y:e40}");
            }


            List<(MultiPrecision<Pow2.N32> x, MultiPrecision<Pow2.N32> y)> expecteds = new();

            using StreamReader sr = new("../../../../results_disused/invgamma_e32.csv");

            sr.ReadLine();
            while (!sr.EndOfStream) {
                string? line = sr.ReadLine();
                if (string.IsNullOrWhiteSpace(line)) {
                    break;
                }

                string[] line_split = line.Split(",");
                MultiPrecision<Pow2.N32> x = line_split[1], y = line_split[0];

                expecteds.Add((x, y));
            }

            //using StreamWriter sw = new("../../../../results_disused/invgamma_e32_2.csv");

            //sw.WriteLine("pade approximant invgamma");

            //for (MultiPrecision<Pow2.N32> x = 2; x <= 180; x += 1 / 1024d) {
            //    MultiPrecision<Pow2.N32> y = f32(x);

            //    if (y > 1024) {
            //        break;
            //    }

            //    sw.WriteLine($"{y},{x}");

            //    expecteds.Add((x, y));
            //}

            //sw.Flush();

            using StreamWriter sw_result = new("../../../../results_disused/gamma_e32_pade_15.csv");

            for (double xmin = 2; xmin <= 3; xmin += 1) {
                double xmax = xmin + 1;

                sw_result.WriteLine($"\nrange x in [{xmin}, {xmax}]");

                sw_result.WriteLine("pade results");

                List<(MultiPrecision<Pow2.N32> x, MultiPrecision<Pow2.N32> y)> expecteds_range =
                    expecteds.Where((item) => item.x >= xmin && item.x <= xmax).ToList();

                int y0 = (int)MultiPrecision<Pow2.N32>.Floor(expecteds_range.Select(item => item.y).Min());

                Vector<Pow2.N32> xs = expecteds_range.Select(item => item.x - xmin).ToArray(), ys = expecteds_range.Select(item => item.y - y0).ToArray();

                for (int m = 4; m <= 32; m++) {
                    PadeFitter<Pow2.N32> pade = new(xs, ys, m, m, intercept: 0);

                    Vector<Pow2.N32> param = pade.Fit();
                    Vector<Pow2.N32> errs = pade.Error(param);

                    MultiPrecision<Pow2.N32> max_rateerr = 0;
                    for (int i = 0; i < errs.Dim; i++) {
                        if (ys[i] == 0) {
                            continue;
                        }

                        max_rateerr = MultiPrecision<Pow2.N32>.Max(errs[i] / ys[i], max_rateerr);
                    }

                    Console.WriteLine($"m={m},n={m}");
                    Console.WriteLine($"{max_rateerr:e20}");

                    if (max_rateerr < 2e-32) {
                        sw_result.WriteLine($"y0={y0}");
                        sw_result.WriteLine($"m={m},n={m}");
                        sw_result.WriteLine("numer");
                        foreach (var v in param[..m]) {
                            sw_result.WriteLine(v.val);
                        }
                        sw_result.WriteLine("denom");
                        foreach (var v in param[m..]) {
                            sw_result.WriteLine(v.val);
                        }
                        sw_result.WriteLine("hexcode");
                        for (int i = 0; i < m; i++) {
                            sw_result.WriteLine($"({ToFP128(param[i])}, {ToFP128(param[i + m])}),");
                        }
                        //sw_result.WriteLine($"({ToFP128(param[m - 1])}, ddouble.Zero),");

                        sw_result.WriteLine("relative err");
                        sw_result.WriteLine($"{max_rateerr:e20}");
                        sw_result.Flush();

                        break;
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        public static string ToFP128(MultiPrecision<Pow2.N32> x) {
            Sign sign = x.Sign;
            long exponent = x.Exponent;
            uint[] mantissa = x.Mantissa.Reverse().ToArray();

            string code = $"({(sign == Sign.Plus ? "+1" : "-1")}, {exponent}, 0x{mantissa[0]:X8}{mantissa[1]:X8}uL, 0x{mantissa[2]:X8}{mantissa[3]:X8}uL)";

            return code;
        }
    }
}
