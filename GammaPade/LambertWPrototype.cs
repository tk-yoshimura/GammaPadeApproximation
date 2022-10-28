﻿using MultiPrecision;
using System;

namespace GammaPade {

    internal struct Plus1<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 1);
    }
    public static class LambertWPrototype<N> where N : struct, IConstant {

        public static MultiPrecision<N> LambertW(MultiPrecision<N> x) {
            if (x.IsZero) {
                return 0;
            }
            if (x.IsNaN) {
                return MultiPrecision<N>.NaN;
            }
            if (x > 0 && !x.IsFinite) {
                return MultiPrecision<N>.PositiveInfinity;
            }
            if (x <= -1 / MultiPrecision<N>.E) {
                return (x < -1 / MultiPrecision<N>.E) ? MultiPrecision<N>.NaN : -1;
            }

            MultiPrecision<Plus1<N>> y;
            if (x >= -0.365) {
                double xd = (double)x, yd;

                if (x < 8) {
                    yd = xd * (60.0 + xd * (114.0 + xd * 17.0)) / (60.0 + xd * (174.0 + xd * 101.0));
                }
                else {
                    double logx = Math.Log(xd), loglogx = Math.Log(logx);

                    yd = logx - loglogx + loglogx / (logx + logx);
                }

                double exp_y, d, dy;

                for (int i = 0; i < 4; i++) {
                    exp_y = Math.Exp(yd);
                    d = yd * exp_y - xd;
                    dy = d / (exp_y * (yd + 1d) - (yd + 2d) * d / (yd + yd + 2d));

                    if (!double.IsFinite(dy)) {
                        break;
                    }

                    yd -= dy;

                    if (Math.Abs(dy) <= Math.Abs(yd) * 1e-15) {
                        break;
                    }
                }

                y = yd;
            }
            else { 
                MultiPrecision<Plus1<N>> v = MultiPrecision<Plus1<N>>.Log(x.Convert<Plus1<N>>() + 1 / MultiPrecision<Plus1<N>>.E);

                y = MultiPrecision<Plus1<N>>.Exp(
                    ((v + MultiPrecision<Plus1<N>>.Log(2 * MultiPrecision<Plus1<N>>.E)) * 2
                    - (MultiPrecision<Plus1<N>>.PI * MultiPrecision<Plus1<N>>.Exp(v / 2))) / 4
                    ) - 1.0;
            }

            {
                MultiPrecision<Plus1<N>> exp_y, d, dy;

                for (int i = 0; i < 65536; i++) {
                    exp_y = MultiPrecision<Plus1<N>>.Exp(y);
                    d = y * exp_y - x.Convert<Plus1<N>>();
                    dy = d / (exp_y * (y + 1) - (y + 2) * d / (y + y + 2));
                    y -= dy;

                    if (y.Exponent - dy.Exponent > MultiPrecision<N>.Bits) {
                        return y.Convert<N>();
                    }
                }
            }

            return MultiPrecision<N>.NaN;
        }
    }
}
