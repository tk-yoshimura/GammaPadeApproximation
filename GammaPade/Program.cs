using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;
using System;
using System.Collections.Generic;
using System.Diagnostics.Metrics;
using System.IO;
using System.Linq;

namespace GammaPade {
    internal class Program {
        static void Main(string[] args) {
            static MultiPrecision<Pow2.N32> f32(MultiPrecision<Pow2.N32> x) {
                return MultiPrecision<Pow2.N32>.Log2(MultiPrecision<Pow2.N32>.Gamma(x));
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

            using StreamWriter sw_result = new("../../../../results_disused/gamma_e32_pade_7.csv");

            foreach ((double xmin, double xmax) in new (double, double)[] { 
                (36, 44), (44, 52), (52, 60), (60, 68), 
                (68, 76), (76, 84), (84, 92), (92, 100),
                (100, 108), (108, 116), (116, 124), (124, 132),
                (132, 140), (140, 148), (148, 156), (156, 164),
                (164, 172)
            }) {
                
                sw_result.WriteLine($"\nrange x in [{xmin}, {xmax}]");

                sw_result.WriteLine("pade results");

                List<(MultiPrecision<Pow2.N32> x, MultiPrecision<Pow2.N32> y)> expecteds_range =
                    expecteds.Where((item) => item.x >= xmin && item.x <= xmax).ToList();
                
                Vector<Pow2.N32> xs = expecteds_range.Select(item => item.x - xmin).ToArray(), ys = expecteds_range.Select(item => item.y / item.x).ToArray();

                for (int m = 2; m <= 32; m++) {
                    PadeFitter<Pow2.N32> pade = new(xs, ys, m, m);

                    Vector<Pow2.N32> param = pade.ExecuteFitting();
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
                        sw_result.WriteLine($"m={m},n={m}");
                        sw_result.WriteLine("numer");
                        foreach(var v in param[..m]) {
                            sw_result.WriteLine(v.val);
                        }
                        sw_result.WriteLine("denom");
                        foreach(var v in param[m..]) {
                            sw_result.WriteLine(v.val);
                        }
                        sw_result.WriteLine("hexcode");
                        for(int i = 0; i < m; i++) {
                            sw_result.WriteLine($"({ToFP128(param[i])}, {ToFP128(param[i + m])}),");
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
