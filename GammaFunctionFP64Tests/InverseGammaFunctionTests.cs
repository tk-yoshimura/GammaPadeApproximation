using GammaFunctionFP64;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace GammaFunctionFP64Tests {
    [TestClass()]
    public class InverseGammaTests {
        [TestMethod()]
        public void InverseGammaPositiveValueTest() {
            for ((double x, int i) = (2, 0); x <= 171.625; x += 1d / 1024, i++) {
                double y = GammaFunction.Gamma(x);
                double z = GammaFunction.InverseGamma(y);

                if ((i % 256) == 0) {
                    Console.WriteLine(x);
                    Console.WriteLine(y);
                    Console.WriteLine(z);
                }

                if (y > double.MaxValue) {
                    Console.WriteLine($"expected > maxvalue: {x}");

                    break;
                }

                Assert.AreEqual(x, z, x * 5.0e-16);
            }
        }
    }
}