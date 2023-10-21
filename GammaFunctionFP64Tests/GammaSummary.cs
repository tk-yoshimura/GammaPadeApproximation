using DoubleDouble;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;

namespace GammaFunctionFP64Tests {
    [TestClass()]
    public class GammaSummary {
        [TestMethod()]
        public void PlotGamma() {
            using StreamWriter sw = new("../../../../results/gamma_approx.csv");

            sw.WriteLine("x,y_expected,y_actual,error(abs),error(rate)");

            MultiPrecision<Pow2.N8> maxrateerror = 0;

            for (double x = 1d / 32; x <= 128; x += 1d / 32) {
                double actual = GammaFunction.Gamma((double)x);
                MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.Gamma(x);
                MultiPrecision<Pow2.N8> error = MultiPrecision<Pow2.N8>.Abs(expected - actual);
                MultiPrecision<Pow2.N8> rateerror = (error != 0) ? error / MultiPrecision<Pow2.N8>.Abs(expected) : 0;

                if (double.IsNormal(actual)) {
                    maxrateerror = MultiPrecision<Pow2.N8>.Max(rateerror, maxrateerror);
                }

                sw.WriteLine($"{x},{expected:e20},{actual},{error:e8},{rateerror:e8}");
            }

            Console.WriteLine($"max error (rate): {maxrateerror:e8}");
        }

        [TestMethod()]
        public void PlotLogGamma() {
            using StreamWriter sw = new("../../../../results/loggamma_approx.csv");

            sw.WriteLine("x,y_expected,y_actual,error(abs),error(rate)");

            MultiPrecision<Pow2.N8> maxrateerror = 0;

            for (double x = 1d / 32; x <= 256; x += 1d / 32) {
                double actual = GammaFunction.LogGamma((double)x);
                MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.LogGamma(x);
                MultiPrecision<Pow2.N8> error = MultiPrecision<Pow2.N8>.Abs(expected - actual);
                MultiPrecision<Pow2.N8> rateerror = (error != 0) ? error / MultiPrecision<Pow2.N8>.Abs(expected) : 0;

                if (double.IsNormal(actual)) {
                    maxrateerror = MultiPrecision<Pow2.N8>.Max(rateerror, maxrateerror);
                }

                sw.WriteLine($"{x},{expected:e20},{actual},{error:e8},{rateerror:e8}");
            }

            Console.WriteLine($"max error (rate): {maxrateerror:e8}");
        }

        [TestMethod()]
        public void PlotDigamma() {
            using StreamWriter sw = new("../../../../results/digamma_approx.csv");

            sw.WriteLine("x,y_expected,y_actual,error(abs),error(rate)");

            MultiPrecision<Pow2.N8> maxrateerror = 0;

            for (double x = 1d / 32; x <= 256; x += 1d / 32) {
                double actual = GammaFunction.Digamma((double)x);
                MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.Digamma(x);
                MultiPrecision<Pow2.N8> error = MultiPrecision<Pow2.N8>.Abs(expected - actual);
                MultiPrecision<Pow2.N8> rateerror = (error != 0) ? error / MultiPrecision<Pow2.N8>.Abs(expected) : 0;

                if (double.IsNormal(actual)) {
                    maxrateerror = MultiPrecision<Pow2.N8>.Max(rateerror, maxrateerror);
                }

                sw.WriteLine($"{x},{expected:e20},{actual},{error:e8},{rateerror:e8}");
            }

            Console.WriteLine($"max error (rate): {maxrateerror:e8}");
        }

        [TestMethod()]
        public void PlotInverseGamma() {
            using StreamWriter sw = new("../../../../results/invgamma_approx.csv");

            sw.WriteLine("x,y_expected,y_actual,error(abs),error(rate)");

            ddouble maxrateerror = 0;

            for (double x = 1; x <= 64; x += 1d / 1024) {
                double actual = GammaFunction.InverseGamma((double)x);
                ddouble expected = ddouble.InverseGamma(x);
                ddouble error = ddouble.Abs(expected - actual);
                ddouble rateerror = (error != 0) ? error / ddouble.Abs(expected) : 0;

                if (double.IsNormal(actual)) {
                    maxrateerror = ddouble.Max(rateerror, maxrateerror);
                }

                sw.WriteLine($"{x},{expected:e20},{actual},{error:e8},{rateerror:e8}");
            }

            Console.WriteLine($"max error (rate): {maxrateerror:e8}");
        }

        [TestMethod()]
        public void PlotPolygamma() {
            for (int n = 1; n <= 16; n++) {
                using StreamWriter sw = new($"../../../../results/polygamma{n}_approx.csv");

                sw.WriteLine("x,y_expected,y_actual,error(abs),error(rate)");

                ddouble maxrateerror = 0;

                for (double x = 1d / 32; x <= 256; x += 1d / 32) {
                    double actual = GammaFunction.Polygamma(n, (double)x);
                    ddouble expected = ddouble.Polygamma(n, x);
                    ddouble error = ddouble.Abs(expected - actual);
                    ddouble rateerror = (error != 0) ? error / ddouble.Abs(expected) : 0;

                    if (double.IsNormal(actual)) {
                        maxrateerror = ddouble.Max(rateerror, maxrateerror);
                    }

                    sw.WriteLine($"{x},{expected:e20},{actual},{error:e8},{rateerror:e8}");
                }

                Console.WriteLine($"n: {n} max error (rate): {maxrateerror:e8}");
            }
        }
    }
}