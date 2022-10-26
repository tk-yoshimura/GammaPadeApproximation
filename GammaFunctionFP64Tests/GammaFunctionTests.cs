using GammaFunctionFP64;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;

namespace GammaFunctionFP64Tests {
    [TestClass()]
    public class GammaTests {
        [TestMethod()]
        public void GammaPositiveValueTest() {
            for ((uint i, ulong f) = (1, 1); i <= 21; f *= i, i++) {
                Assert.AreEqual(f, GammaFunction.Gamma(i), $"gamma({i})");
            }

            for ((double x, int i) = (1d / 1024, 1); x <= 171.625; x += 1d / 1024, i++) {
                MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.Gamma(x);
                double actual = GammaFunction.Gamma(x);

                if ((i % 256) == 0) {
                    Console.WriteLine(x);
                    Console.WriteLine(actual);
                    Console.WriteLine(expected);
                }

                if (expected > double.MaxValue) {
                    Console.WriteLine($"expected > maxvalue: {x}");

                    break;
                }

                if (MultiPrecision<Pow2.N8>.Abs(actual - expected) > MultiPrecision<Pow2.N8>.Abs(expected) * 8e-16 || !double.IsFinite(actual)) {
                    if (!double.IsNormal(actual) && !double.IsNormal((double)expected)) {
                        Console.WriteLine("skip checking abnormal value");
                        continue;
                    }

                    Assert.Fail($"{x}\n{actual}\n{expected}");
                }
            }
        }

        [TestMethod()]
        public void GammaNegativeValueTest() {
            for ((double x, int i) = (-1d / 1024, 1); x >= -32; x -= 1d / 1024, i++) {
                MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.Gamma(x);
                double actual = GammaFunction.Gamma(x);

                if ((i % 256) == 0) {
                    Console.WriteLine(x);
                    Console.WriteLine(actual);
                    Console.WriteLine(expected);
                }

                if (!expected.IsFinite) {
                    continue;
                }

                if (MultiPrecision<Pow2.N8>.Abs(actual - expected) > MultiPrecision<Pow2.N8>.Abs(expected) * 8e-16 || !double.IsFinite(actual)) {
                    if (!double.IsNormal(actual) && !double.IsNormal((double)expected)) {
                        Console.WriteLine("skip checking abnormal value");
                        continue;
                    }

                    Assert.Fail($"{x}\n{actual}\n{expected}");
                }
            }
        }
    }
}