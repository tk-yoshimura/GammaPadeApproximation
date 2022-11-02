// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

using System.Collections.ObjectModel;

namespace GammaFunctionFP64 {
    public static partial class GammaFunction {
        public static double Polygamma(int n, double x) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (n > 16) {
                throw new ArgumentOutOfRangeException(
                    nameof(n),
                    "In the calculation of the Polygamma function, n greater than 16 is not supported."
                );
            }
            if (n == 0) {
                return Digamma(x);
            }
            if (double.IsNaN(x) || double.IsNegativeInfinity(x)) {
                return double.NaN;
            }
            if (double.IsPositiveInfinity(x)) {
                return 0;
            }

            if (x >= 0) {
                double y = PolygammaPlusX.Polygamma(n, x);

                return ((n & 1) == 1) ? y : -y;
            }
            else {
                return PolygammaMinusX.Polygamma(n, x);
            }
        }

        private static class PolygammaPlusX {

            public static double Polygamma(int n, double x) {
                if (x <= Math.ScaleB(1, -96)) {
                    return ((n & 1) == 1) ? double.PositiveInfinity : double.NaN;
                }

                if (n == 1 && x <= 8) {
                    return PolygammaN1Recurrence(x);
                }

                if ((n == 2 && x <= 16) || (n >= 3 && x <= 32)) {
                    return PolygammaNearZero(n, x);
                }

                return PolygammaLimit(n, x);
            }

            public static double PolygammaN1Recurrence(double x) {
                int k = (int)Math.Floor(x);
                double f = x - k;

                double y = PolygammaLimit(1, 8.0 + f);
                for (int i = k; i < 8; i++) {
                    y += 1 / Square(i + f);
                }

                return y;
            }

            public static double PolygammaNearZero(int n, double x) {
#if DEBUG
                if (n <= 1) {
                    throw new ArgumentOutOfRangeException(nameof(n));
                }
#endif

                if (x < 1) {
                    double v = PolygammaNearZero(n, x + 1d);
                    double y = v + factorial[n] / Pow(x, n + 1);

                    return y;
                }

                double scale = (n == 2) ? 2.0 : (n == 3) ? 1.375 : 1, r = scale * n / x;
                double ir = 0, it = 0;

                double polygamma_ir(double t) {
                    double y = Pow(t, n) * Math.Exp(-x * t) / (1d - Math.Exp(-t));

                    return y;
                }

                double polygamma_it(double u) {
                    double v = (u + scale * n) / x;
                    double y = Pow(v, n) / (1d - Math.Exp(-v));

                    return y;
                }

                foreach ((double t, double w) in gles) {
                    ir += w * polygamma_ir(t * r);
                }
                foreach ((double t, double w) in glas) {
                    it += w * polygamma_it(t);
                }

                ir *= r;
                it *= Math.Exp(-scale * n) / x;

                double i = ir + it;

                return i;
            }

            public static double PolygammaLimit(int n, double x) {
                double inv_x2 = 1d / (x * x), c = Math.Pow(x, -n);
                double v = c * factorial[n - 1] * (2 * x + n) / (2 * x);
                double u = c * factorial[n + 1] / 2 * inv_x2;
                double dv = bernoulli_sequence[1] * u;

                v += dv;

                for (int k = 2; k <= 10; k++) {
                    u *= inv_x2 * ((n + 2 * k - 2) * (n + 2 * k - 1)) / ((2 * k) * (2 * k - 1));
                    dv = bernoulli_sequence[k] * u;
                    double next_v = v + dv;

                    if (v == next_v) {
                        break;
                    }
                    if (double.IsNaN(next_v)) {
                        return 0;
                    }

                    v = next_v;
                }

                return v;
            }

            static readonly ReadOnlyCollection<(double x, double w)> gles = new(new (double, double)[] {
                (3.43570040745253760694e-3, 8.80700356957605915593e-3),
                (1.80140363610431043662e-2, 2.03007149001934706655e-2),
                (4.38827858743370470661e-2, 3.13360241670545317848e-2),
                (8.04415140888905883027e-2, 4.16383707883523743624e-2),
                (1.26834046769924603693e-1, 5.09650599086202175184e-2),
                (1.81973159636742487274e-1, 5.90972659807592086562e-2),
                (2.44566499024586450998e-1, 6.58443192245883134492e-2),
                (3.13146955642290219664e-1, 7.10480546591910256646e-2),
                (3.86107074429177460960e-1, 7.45864932363018733939e-2),
                (4.61736739433251333123e-1, 7.63766935653629253490e-2),
                (5.38263260566748666877e-1, 7.63766935653629253490e-2),
                (6.13892925570822539040e-1, 7.45864932363018733939e-2),
                (6.86853044357709780336e-1, 7.10480546591910256646e-2),
                (7.55433500975413549002e-1, 6.58443192245883134492e-2),
                (8.18026840363257512726e-1, 5.90972659807592086562e-2),
                (8.73165953230075396307e-1, 5.09650599086202175184e-2),
                (9.19558485911109411697e-1, 4.16383707883523743624e-2),
                (9.56117214125662952934e-1, 3.13360241670545317848e-2),
                (9.81985963638956895634e-1, 2.03007149001934706655e-2),
                (9.96564299592547462393e-1, 8.80700356957605915593e-3),
            });

            static readonly ReadOnlyCollection<(double x, double w)> glas = new(new (double, double)[] {
                (7.05398896919887533667e-2, 1.68746801851113862149e-1),
                (3.72126818001611443794e-1, 2.91254362006068281717e-1),
                (9.16582102483273564668e-1, 2.66686102867001288550e-1),
                (1.70730653102834388069e0, 1.66002453269506840031e-1),
                (2.74919925530943212965e0, 7.48260646687923705401e-2),
                (4.04892531385088692237e0, 2.49644173092832210728e-2),
                (5.61517497086161651410e0, 6.20255084457223684745e-3),
                (7.45901745367106330977e0, 1.14496238647690824204e-3),
                (9.59439286958109677247e0, 1.55741773027811974780e-4),
                (1.20388025469643163096e1, 1.54014408652249156894e-5),
                (1.48142934426307399785e1, 1.08648636651798235148e-6),
                (1.79488955205193760174e1, 5.33012090955671475093e-8),
                (2.14787882402850109757e1, 1.75798117905058200358e-9),
                (2.54517027931869055035e1, 3.72550240251232087263e-11),
                (2.99325546317006120067e1, 4.76752925157819052449e-13),
                (3.50134342404790000063e1, 3.37284424336243841237e-15),
                (4.08330570567285710620e1, 1.15501433950039883096e-17),
                (4.76199940473465021399e1, 1.53952214058234355346e-20),
                (5.58107957500638988908e1, 5.28644272556915782880e-24),
                (6.65244165256157538186e1, 1.65645661249902329591e-28),
            });

            static readonly ReadOnlyCollection<double> factorial = new(new double[] {
                1,
                1,
                2,
                6,
                24,
                120,
                720,
                5040,
                40320,
                362880,
                3628800,
                39916800,
                479001600,
                6227020800,
                87178291200,
                1307674368000,
                20922789888000,
                355687428096000
            });

            static readonly ReadOnlyCollection<double> bernoulli_sequence = new(new double[] {
                1,
                1d/6,
                -1d/30,
                1d/42,
                -1d/30,
                5d/66,
                -691d/2730,
                7d/6,
                -3617d/510,
                43867d/798,
                -174611d/330,
            });
        }

        private static class PolygammaMinusX {
            public static double Polygamma(int n, double x) {
                if (Math.Abs(x - Math.Round(x)) <= Math.ScaleB(1, -96)) {
                    return ((n & 1) == 1) ? double.PositiveInfinity : double.NaN;
                }

                double p = PolygammaPlusX.Polygamma(n, 1d - x);
                double g = Math.Pow(Math.PI, n + 1) * reflecs[n](x);

                double y = -(g + p);

                return y;
            }

            static readonly ReadOnlyCollection<Func<double, double>> reflecs = new(new Func<double, double>[] {
               (x) => 1d / TanPI(x),
               (x) => -1d / Square(SinPI(x)),
               (x) => 2d / (TanPI(x) * Square(SinPI(x))),
               (x) => -2d * (2d + CosPI(2 * x)) / Pow(SinPI(x), 4),
               (x) => 4d * (5d + CosPI(2 * x)) / (TanPI(x) * Pow(SinPI(x), 4)),
               (x) => -2d * (33d + 26d * CosPI(2 * x) + CosPI(4 * x)) / Pow(SinPI(x), 6),
               (x) => 4d * (123d + 56d * CosPI(2 * x) + CosPI(4 * x)) / (TanPI(x) * Pow(SinPI(x), 6)),
               (x) => -2d * (1208d + 1191d * CosPI(2 * x) + 120d * CosPI(4 * x) + CosPI(6 * x)) / Pow(SinPI(x), 8),
               (x) => 4d * (5786d + 4047d * CosPI(2 * x) + 246d * CosPI(4 * x) + CosPI(6 * x)) / (TanPI(x) * Pow(SinPI(x), 8)),
               (x) => -2d * (78095d + 88234d * CosPI(2 * x) + 14608d * CosPI(4 * x) + 502d * CosPI(6 * x) + CosPI(8 * x)) / Pow(SinPI(x), 10),
               (x) => 4d * (450995d + 408364d * CosPI(2 * x) + 46828d * CosPI(4 * x) + 1012d * CosPI(6 * x) + CosPI(8 * x)) / (TanPI(x) * Pow(SinPI(x), 10)),
               (x) => -2d * (7862124d + 9738114d * CosPI(2 * x) + 2203488d * CosPI(4 * x) + 152637d * CosPI(6 * x) + 2036d * CosPI(8 * x) + CosPI(10 * x)) / Pow(SinPI(x), 12),
               (x) => 4d * (52953654d + 56604978d * CosPI(2 * x) + 9713496d * CosPI(4 * x) + 474189d * CosPI(6 * x) + 4082d * CosPI(8 * x) + CosPI(10 * x)) / (TanPI(x) * Pow(SinPI(x), 12)),
               (x) => -2d * (1137586002d + 1505621508d * CosPI(2 * x) + 423281535d * CosPI(4 * x) + 45533450d * CosPI(6 * x) + 1479726d * CosPI(8 * x) + 8178d * CosPI(10 * x) + CosPI(12 * x)) / Pow(SinPI(x), 14),
               (x) => 4d * (8752882782d + 10465410528d * CosPI(2 * x) + 2377852335d * CosPI(4 * x) + 193889840d * CosPI(6 * x) + 4520946d * CosPI(8 * x) + 16368d * CosPI(10 * x) + CosPI(12 * x)) / (TanPI(x) * Pow(SinPI(x), 14)),
               (x) => -2d * (223769408736d + 311387598411d * CosPI(2 * x) + 102776998928d * CosPI(4 * x) + 15041229521d * CosPI(6 * x) + 848090912d * CosPI(8 * x) + 13824739d * CosPI(10 * x) + 32752d * CosPI(12 * x) + CosPI(14 * x)) / Pow(SinPI(x), 16),
               (x) => 4d * (1937789122548d + 2507220680379d * CosPI(2 * x) + 700262497778d * CosPI(4 * x) + 81853020521d * CosPI(6 * x) + 3530218028d * CosPI(8 * x) + 41867227d * CosPI(10 * x) + 65518d * CosPI(12 * x) + CosPI(14 * x)) / (TanPI(x) * Pow(SinPI(x), 16)),
            });
        }
    }
}