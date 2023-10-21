using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionCurveFitting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaPade {
    internal class InvGammaPolar {
        static void Main() {
            const string filepath_expected = "../../../../results_disused/invgamma_polar_expected.csv";

            List<(MultiPrecision<Pow2.N32> u, MultiPrecision<Pow2.N32> v)> expecteds = new();

            if (!Path.Exists(filepath_expected)) {
                using StreamWriter sw = new(filepath_expected);

                sw.WriteLine($"x,y:=invgamma(x),u:=x-sqrt(pi)/2,v:=(y-3/2)^2");

                MultiPrecision<Pow2.N32> c = MultiPrecision<Pow2.N32>.Sqrt(MultiPrecision<Pow2.N32>.PI) / 2;

                for (MultiPrecision<Pow2.N32> y = 1.5; y <= 2; y += 1 / 8192d) {
                    MultiPrecision<Pow2.N32> x = MultiPrecision<Pow2.N32>.Gamma(y);

                    MultiPrecision<Pow2.N32> u = x - c;
                    MultiPrecision<Pow2.N32> v = MultiPrecision<Pow2.N32>.Square(y - 1.5);

                    sw.WriteLine($"{x},{y},{u},{v}");

                    expecteds.Add((u, v));
                }
            }
            else {
                using StreamReader sr = new(filepath_expected);

                sr.ReadLine();

                while (!sr.EndOfStream) {
                    string? line = sr.ReadLine();

                    if (string.IsNullOrWhiteSpace(line)) {
                        continue;
                    }

                    string[] line_split = line.Split(',');

                    MultiPrecision<Pow2.N32> u = line_split[2], v = line_split[3];

                    expecteds.Add((u, v));
                }
            }

            using StreamWriter sw_result = new("../../../../results_disused/invgamma_polar_e32_pade.csv");

            foreach ((double umin, double umax) in new[] { (0, 1 / 256d), (1 / 256d, 1 / 128d), (1 / 128d, 1 / 32d), (1 / 32d, 1 / 8d) }) {
                List<(MultiPrecision<Pow2.N32> u, MultiPrecision<Pow2.N32> v)> expecteds_range =
                    expecteds.Where((item) => item.u >= umin && item.u <= umax).ToList();

                Vector<Pow2.N32> us = expecteds_range.Select(item => item.u - umin).ToArray(), vs = expecteds_range.Select(item => item.v).ToArray();

                for (int m = 4; m <= 50; m++) {
                    PadeFitter<Pow2.N32> pade = new(us, vs, m, m, intercept: umin == 0 ? 0 : null);

                    Vector<Pow2.N32> param = pade.ExecuteFitting();
                    Vector<Pow2.N32> errs = pade.Error(param);

                    MultiPrecision<Pow2.N32> max_rateerr = 0;
                    for (int i = 0; i < errs.Dim; i++) {
                        if (vs[i] == 0) {
                            continue;
                        }

                        max_rateerr = MultiPrecision<Pow2.N32>.Max(errs[i] / vs[i], max_rateerr);
                    }

                    Console.WriteLine($"m={m},n={m}");
                    Console.WriteLine($"{max_rateerr:e20}");

                    if (max_rateerr < 2e-32) {
                        sw_result.WriteLine($"u=[{umin},{umax}]");
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
