using MultiPrecision;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GammaPade {
    internal class Program {
        static void Main(string[] args) {
            static MultiPrecision<Pow2.N256> f256(MultiPrecision<Pow2.N256> x) {
                if (x == 0) {
                    return -1;
                }
                return LambertWPrototype<Pow2.N256>.LambertW(x * x - 1 / MultiPrecision<Pow2.N256>.E);
            }
            
            static MultiPrecision<Pow2.N128> f128(MultiPrecision<Pow2.N128> x) {
                if (x == 0) {
                    return -1;
                }                
                return LambertWPrototype<Pow2.N128>.LambertW(x * x - 1 / MultiPrecision<Pow2.N128>.E);
            }
            
            using StreamWriter sw = new("../../../../results_disused/lambert_w_pade_table_e32.csv");
            
            MultiPrecision<Pow2.N128>[] xs = new MultiPrecision<Pow2.N128>[]{
                0, 1d / 4
            };
            
            sw.WriteLine("pade approximant lambert_w(x^2-1/e)");
            
            for (int j = 0; j < xs.Length - 1; j++) {
                MultiPrecision<Pow2.N128> x0 = xs[j], x1 = xs[j + 1];
                MultiPrecision<Pow2.N128> dx = (x1 - x0) / 1024;
            
                sw.WriteLine($"\nrange x in [{x0}, {x1}]");
            
                sw.WriteLine("expected");
                List<MultiPrecision<Pow2.N128>> expecteds = new();
                for (MultiPrecision<Pow2.N128> x = x0; x <= x1; x += dx) {
                    MultiPrecision<Pow2.N128> y = f128(x);
            
                    expecteds.Add(y);
            
                    sw.WriteLine($"{x},{y:e40}");
                }
            
                sw.WriteLine($"diffs x = {x0}");
                MultiPrecision<Pow2.N256>[] diffs = ForwardFiniteDifference<Pow2.N256>.Diff(x0.Convert<Pow2.N256>(), f256, Math.ScaleB(1, -24));
            
                sw.Flush();
            
                MultiPrecision<Pow2.N128>[] cs = new MultiPrecision<Pow2.N128>[diffs.Length + 1];
                cs[0] = f128(x0);
                for (int i = 0; i < diffs.Length; i++) {
                    cs[i + 1] = diffs[i].Convert<Pow2.N128>() * MultiPrecision<Pow2.N128>.TaylorSequence[i + 1];
                }
            
                sw.WriteLine("pade results");
            
                for (int m = 2; m <= 128; m++) {
                    (MultiPrecision<Pow2.N128>[] ms, MultiPrecision<Pow2.N128>[] ns) =
                        PadeSolver<Pow2.N128>.Solve(cs.Take(m + m + 1).ToArray(), m, m);
            
                    MultiPrecision<Pow2.N128> err = 0;
                    for ((int i, MultiPrecision<Pow2.N128> x) = (0, x0); i < expecteds.Count; i++, x += dx) {
                        MultiPrecision<Pow2.N128> expected = expecteds[i];
                        MultiPrecision<Pow2.N128> actual = PadeSolver<Pow2.N128>.Approx(x - x0, ms, ns);
            
                        err = MultiPrecision<Pow2.N128>.Max(err, MultiPrecision<Pow2.N128>.Abs(expected / actual - 1));
                    }
            
                    Console.WriteLine($"m={m},n={m}");
                    Console.WriteLine($"{err:e20}");
            
                    if (err < 2e-32) {
                        sw.WriteLine($"m={m},n={m}");
                        for (int i = 0; i <= m; i++) {
                            sw.WriteLine($"({ToFP128(ms[i])}, {ToFP128(ns[i])}), ");
                        }
                        sw.WriteLine("relative err");
                        sw.WriteLine($"{err:e20}");
                        sw.Flush();
            
                        break;
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        public static string ToFP128(MultiPrecision<Pow2.N128> x) {
            Sign sign = x.Sign;
            long exponent = x.Exponent;
            uint[] mantissa = x.Mantissa.Reverse().ToArray();

            string code = $"({(sign == Sign.Plus ? "+1" : "-1")}, {exponent}, 0x{mantissa[0]:X8}{mantissa[1]:X8}uL, 0x{mantissa[2]:X8}{mantissa[3]:X8}uL)";

            return code;
        }
    }
}
