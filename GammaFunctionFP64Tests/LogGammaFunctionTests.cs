using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;

namespace GammaFunctionFP64Tests {
    [TestClass()]
    public class LogGammaTests {
        [TestMethod()]
        public void LogGammaValueTest() {
            for ((double x, int i) = (1d / 1024, 1); x <= 128; x += 1d / 1024, i++) {
                MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.LogGamma(x);
                double actual = GammaFunction.LogGamma(x);

                if ((i % 256) == 0) {
                    Console.WriteLine(x);
                    Console.WriteLine(actual);
                    Console.WriteLine(expected);
                }

                if (MultiPrecision<Pow2.N8>.Abs(actual - expected) > MultiPrecision<Pow2.N8>.Abs(expected) * 8e-16 || !double.IsFinite(actual)) {
                    if (!double.IsNormal(actual) && !double.IsNormal((double)expected)) {
                        Console.WriteLine("skip checking abnormal value");
                        continue;
                    }

                    Assert.Fail($"{x}\n{actual}\n{expected}");
                }
            }

            for (double x = 256; x < double.MaxValue; x *= 2) {
                MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.LogGamma(x);
                double actual = GammaFunction.LogGamma(x);

                Console.WriteLine(x);
                Console.WriteLine(actual);
                Console.WriteLine(expected);

                if (MultiPrecision<Pow2.N8>.Abs(actual - expected) > MultiPrecision<Pow2.N8>.Abs(expected) * 8e-16 || !double.IsFinite(actual)) {
                    if (!double.IsNormal(actual) && !double.IsNormal((double)expected)) {
                        Console.WriteLine("skip checking abnormal value");
                        break;
                    }

                    Assert.Fail($"{x}\n{actual}\n{expected}");
                }
            }
        }
    }
}