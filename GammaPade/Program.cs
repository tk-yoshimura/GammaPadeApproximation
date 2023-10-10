using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaPade {
    internal class Program {
        static void Main(string[] args) {
            static MultiPrecision<Pow2.N32> f32(MultiPrecision<Pow2.N32> x) {
                return MultiPrecision<Pow2.N32>.Log2(MultiPrecision<Pow2.N32>.Gamma(x));
            }

            List<(MultiPrecision<Pow2.N32> x, MultiPrecision<Pow2.N32> y)> expecteds = new();

            //using StreamReader sr = new("../../../../results_disused/invgamma_e32.csv");

            //sr.ReadLine();
            //while (!sr.EndOfStream) {
            //    string? line = sr.ReadLine();
            //    if (string.IsNullOrWhiteSpace(line)) {
            //        break;
            //    }

            //    string[] line_split = line.Split(",");
            //    MultiPrecision<Pow2.N32> x = line_split[1], y = line_split[0];

            //    expecteds.Add((x, y));
            //}

            using StreamWriter sw = new("../../../../results_disused/invgamma_e32_2.csv");

            sw.WriteLine("pade approximant invgamma");

            for (MultiPrecision<Pow2.N32> x = 2; x <= 180; x += 1 / 1024d) {
                MultiPrecision<Pow2.N32> y = f32(x);

                if (y > 1024) {
                    break;
                }

                sw.WriteLine($"{y},{x}");

                expecteds.Add((x, y));
            }

            sw.Flush();

            using StreamWriter sw_result = new("../../../../results_disused/invgamma_e32_pade_2.csv");

            foreach ((double ymin, double ymax) in new (double, double)[] { (0, 0.5), (0.5, 1), (1, 2), (2, 4), (4, 8), (8, 16), (16, 32), (32, 64), (64, 128), (128, 256), (256, 512), (512, 1024) }) {
                sw_result.WriteLine($"\nrange y in [{ymin}, {ymax}]");

                sw_result.WriteLine("pade results");

                List<(MultiPrecision<Pow2.N32> x, MultiPrecision<Pow2.N32> y)> expecteds_range =
                    expecteds.Where((item) => item.y >= ymin && item.y <= ymax).ToList();
                
                Vector<Pow2.N32> xs = expecteds_range.Select(item => item.x).ToArray(), ys = expecteds_range.Select(item => item.y - ymin).ToArray();

                for (int m = 2; m <= 32; m++) {
                    PadeFitter<Pow2.N32> pade = new(ys, xs, m, m);

                    Vector<Pow2.N32> param = pade.ExecuteFitting();
                    Vector<Pow2.N32> errs = pade.Error(param);

                    MultiPrecision<Pow2.N32> max_rateerr = 0;
                    for (int i = 0; i < errs.Dim; i++) {
                        if (xs[i] == 0) {
                            continue;
                        }

                        max_rateerr = MultiPrecision<Pow2.N32>.Max(errs[i] / xs[i], max_rateerr);
                    }

                    Console.WriteLine($"m={m},n={m}");
                    Console.WriteLine($"{max_rateerr:e20}");

                    if (max_rateerr < 2e-32) {
                        sw_result.WriteLine($"m={m},n={m}");
                        sw_result.WriteLine("numer");
                        foreach(var v in param[..m]) {
                            sw_result.WriteLine(v.val);
                        }
                        sw_result.WriteLine("denom");
                        foreach(var v in param[m..]) {
                            sw_result.WriteLine(v.val);
                        }
                        sw_result.WriteLine("numer(hexcode)");
                        foreach(var v in param[..m]) {
                            sw_result.WriteLine(ToFP128(v.val));
                        }
                        sw_result.WriteLine("denom(hexcode)");
                        foreach(var v in param[m..]) {
                            sw_result.WriteLine(ToFP128(v.val));
                        }
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
