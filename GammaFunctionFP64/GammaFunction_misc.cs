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
    }
}