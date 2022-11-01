// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

using System.Runtime.CompilerServices;

namespace GammaFunctionFP64 {
    public static partial class GammaFunction {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Fma(double z, double x, double y) {
            return Math.FusedMultiplyAdd(x, y, z);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Theta(double x) {
            double r = Math.Round(x);
            double t = (((int)r & 1) == 0) ? (x - r) : (r - x);

            return t;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double SinPI(double x) {
            double t = Theta(x);
            double y = Math.Sin(t * Math.PI);

            return y;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double CosPI(double x) {
            double y = SinPI(x + 0.5);

            return y;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double TanPI(double x) {
            double t = x - Math.Round(x);
            double y = Math.Tan(t * Math.PI);

            return y;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Square(double x) => x * x;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Pow(double x, int n) {
            if (double.IsNaN(x)) {
                return double.NaN;
            }

            if (n == 0) {
                return 1d;
            }

            int n_abs = Math.Abs(n);
            double y = 1d, z = x;

            while (n_abs > 0) {
                if ((n_abs & 1) == 1) {
                    y *= z;
                }

                z *= z;
                n_abs >>= 1;
            }

            return (n > 0) ? y : 1.0 / y;
        }
    }
}