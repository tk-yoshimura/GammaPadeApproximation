// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

namespace GammaFunctionFP64 {
    public static partial class GammaFunction {
        public static double InverseGamma(double x) {
            if (!(x >= 1.0)) {
                return double.NaN;
            }
            if (double.IsPositiveInfinity(x)) {
                return double.PositiveInfinity;
            }

            static double lambert_w(double x) {
                double y;

                if (x < 8) {
                    y = x * (30.0 + x * (21.0 + x)) / (30.0 + x * (36.0 + x * 9.0));
                }
                else {
                    double lnx = Math.Log(x), lnlnx = Math.Log(lnx);
                    y = lnx - lnlnx + lnlnx / (lnx + lnx);
                }

                double exp_y = Math.Exp(y);
                double d = y * exp_y - x;

                y -= d / (exp_y * (y + 1.0) - (y + 2.0) * d / (2.0 * y + 2.0));

                return y;
            }

            double c = 0.0365338144849004166003358471688;
            double s = 0.3989422804014326779399460599344;

            double lnx = Math.Log(x);
            double l = Math.Log((x + c) * s);
            double y = l / lambert_w(l / Math.E) + 0.5;

            for (int i = 0; i < 4; i++) {
                double lng = LogGamma(y), psi = Digamma(y);
                double delta = (lnx - lng) / psi;

                y += delta;

                if (Math.Abs(delta) < y * 5.0e-16) {
                    break;
                }
            }

            for (int i = 0; i < 4; i++) {
                double g = Gamma(y);
                if (x == g || !double.IsFinite(g)) {
                    break;
                }

                double psi = Digamma(y);
                double delta = (x / g - 1.0) / psi;

                y += delta;

                if (Math.Abs(delta) < y * 2.5e-16) {
                    break;
                }
            }

            return y;
        }
    }
}