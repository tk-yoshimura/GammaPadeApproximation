using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaPade {
    internal class GammaNearZero {
        static void Main(string[] args) {
            using StreamWriter sw = new("../../../../results_disused/gamma_e32_pade_nearzero_5.csv");

            List<(MultiPrecision<Pow2.N32> xmin, MultiPrecision<Pow2.N32> xmax, MultiPrecision<Pow2.N32> limit_range)> ranges = [];
            for (double x0 = 1.5; x0 < 2; x0 += 0.5) {
                ranges.Add((x0, x0 + 0.5, 0.5));
            }

            bool approximate(MultiPrecision<Pow2.N32> xmin, MultiPrecision<Pow2.N32> xmax) {
                Console.WriteLine($"[{xmin}, {xmax}]");

                MultiPrecision<Pow2.N32> xrange = xmax - xmin;

                List<(MultiPrecision<Pow2.N64> x, MultiPrecision<Pow2.N64> y)> expecteds_range = [];

                for (MultiPrecision<Pow2.N32> u = 0, h = 1 / 8192d; u <= 1; u += h) {
                    MultiPrecision<Pow2.N32> y = MultiPrecision<Pow2.N32>.Gamma(xmin + u * xrange);

                    expecteds_range.Add(((u * xrange).Convert<Pow2.N64>(), y.Convert<Pow2.N64>()));
                }

                Vector<Pow2.N64> xs = expecteds_range.Select(item => item.x).ToArray();
                Vector<Pow2.N64> ys = expecteds_range.Select(item => item.y).ToArray();

                SumTable<Pow2.N64> sum_table = new(xs, ys);

                for (int coefs = 4; coefs <= 72; coefs++) {
                    foreach ((int m, int n) in CurveFittingUtils.EnumeratePadeDegree(coefs, 1)) {
                        PadeFitter<Pow2.N64> pade = new(sum_table, m, n, intercept: ys[0]);

                        Vector<Pow2.N64> param = pade.Fit();
                        Vector<Pow2.N64> errs = pade.Error(param);

                        MultiPrecision<Pow2.N64> max_rateerr = CurveFittingUtils.MaxRelativeError(ys, pade.Regress(xs, param));

                        Console.WriteLine($"m={m},n={n}");
                        Console.WriteLine($"{max_rateerr:e20}");

                        if (max_rateerr > "1e-20") {
                            break;
                        }

                        if (max_rateerr < "1e-40") {
                            return false;
                        }

                        if (max_rateerr < "1e-34" &&
                            !CurveFittingUtils.HasLossDigitsPolynomialCoef(param[..m], 0, xrange.Convert<Pow2.N64>()) &&
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
                            sw.WriteLine("coef");
                            foreach ((var numer, var denom) in CurveFittingUtils.EnumeratePadeCoef(param, m, n)) {
                                sw.WriteLine($"({ToFP128(numer)}, {ToFP128(denom)}),");
                            }

                            sw.WriteLine("relative err");
                            sw.WriteLine($"{max_rateerr:e20}");
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