using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaPade {
    internal class GammaExponent {
        static void Main_(string[] args) {
            static MultiPrecision<Pow2.N32> f32(MultiPrecision<Pow2.N32> x) {
                return MultiPrecision<Pow2.N32>.Log2(MultiPrecision<Pow2.N32>.Gamma(x));
            }

            using StreamWriter sw = new("../../../../results_disused/gamma_e32_pade_21.csv");

            List<(MultiPrecision<Pow2.N32> xmin, MultiPrecision<Pow2.N32> xmax, MultiPrecision<Pow2.N32> limit_range)> ranges = [];
            for (int x0 = 28; x0 < 36; x0 += 8) {
                ranges.Add((x0, x0 + 8, 1));
            }

            bool approximate(MultiPrecision<Pow2.N32> xmin, MultiPrecision<Pow2.N32> xmax) {
                Console.WriteLine($"[{xmin}, {xmax}]");

                MultiPrecision<Pow2.N32> xrange = xmax - xmin;

                long exponent_xmin = MultiPrecision<Pow2.N32>.Gamma(xmin).Exponent;
                long exponent_xmax = MultiPrecision<Pow2.N32>.Gamma(xmax).Exponent;

                long bias = exponent_xmax - exponent_xmin;

                List<(MultiPrecision<Pow2.N64> x, MultiPrecision<Pow2.N64> y)> expecteds_range = [];

                for (MultiPrecision<Pow2.N32> u = 0, h = 1 / 8192d; u <= 1; u += h) {
                    MultiPrecision<Pow2.N32> y = f32(xmin + u * xrange) - (exponent_xmin + bias * u);

                    expecteds_range.Add(((u * xrange).Convert<Pow2.N64>(), y.Convert<Pow2.N64>()));
                }

                if (bias > 1 && expecteds_range.Any(item => item.y < -1 || item.y >= 1)) {
                    Console.WriteLine("over one");
                    return false;
                }

                Vector<Pow2.N64> xs = expecteds_range.Select(item => item.x).ToArray();
                Vector<Pow2.N64> ys = expecteds_range.Select(item => item.y).ToArray();

                SumTable<Pow2.N64> sum_table = new(xs, ys);

                for (int coefs = 4; coefs <= ((bias > 1) ? 36 : 72); coefs++) {
                    foreach ((int m, int n) in CurveFittingUtils.EnumeratePadeDegree(coefs, 1)) {
                        PadeFitter<Pow2.N64> pade = new(sum_table, m, n);

                        Vector<Pow2.N64> param = pade.Fit();
                        Vector<Pow2.N64> errs = pade.Error(param);

                        MultiPrecision<Pow2.N64> max_abserr = CurveFittingUtils.MaxAbsoluteError(ys, pade.Regress(xs, param));

                        Console.WriteLine($"m={m},n={n}");
                        Console.WriteLine($"{max_abserr:e20}");

                        if (max_abserr > "1e-20") {
                            break;
                        }

                        if (max_abserr < "1e-40") {
                            return false;
                        }

                        if (max_abserr < "1e-31" &&
                            !CurveFittingUtils.HasLossDigitsPolynomialCoef(param[m..], 0, xrange.Convert<Pow2.N64>())) {

                            sw.WriteLine($"x=[{xmin},{xmax}]");
                            sw.WriteLine($"samples={expecteds_range.Count}");
                            sw.WriteLine($"m={m},n={n}");
                            sw.WriteLine("numer");
                            foreach (var (_, val) in param[..m]) {
                                sw.WriteLine($"{val:e38}");
                            }
                            sw.WriteLine("denom");
                            foreach (var (_, val) in param[m..]) {
                                sw.WriteLine($"{val:e38}");
                            }

                            sw.WriteLine($"exponent: {exponent_xmin}");
                            sw.WriteLine($"bias: {bias}");

                            sw.WriteLine("coef");
                            foreach ((var numer, var denom) in CurveFittingUtils.EnumeratePadeCoef(param, m, n)) {
                                sw.WriteLine($"({ToFP128(numer)}, {ToFP128(denom)}),");
                            }

                            sw.WriteLine("absolute err");
                            sw.WriteLine($"{max_abserr:e20}");
                            sw.Flush();

                            return true;
                        }
                    }
                }

                return false;
            }

            Segmenter<Pow2.N32> segmenter = new(ranges, approximate);
            segmenter.Execute();

            foreach ((var xmin, var xmax, bool is_successs) in segmenter.ApproximatedRanges) {
                sw.WriteLine($"[{xmin},{xmax}],{(is_successs ? "OK" : "NG")}");
            }

            Console.WriteLine("END");
            Console.Read();
        }

        public static string ToFP128(MultiPrecision<Pow2.N64> x) {
            Sign sign = x.Sign;
            long exponent = x.Exponent;
            uint[] mantissa = x.Mantissa.Reverse().ToArray();

            string code = $"({(sign == Sign.Plus ? "+1" : "-1")}, {exponent}, 0x{mantissa[0]:X8}{mantissa[1]:X8}uL, 0x{mantissa[2]:X8}{mantissa[3]:X8}uL)";

            return code;
        }
    }
}