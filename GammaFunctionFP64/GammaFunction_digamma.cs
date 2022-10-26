// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

namespace GammaFunctionFP64 {
    public static partial class GammaFunction {
        public static double DigammaZero = 1.46163214496836234126266;

        public static double Digamma(double x) {
            if (double.IsNaN(x) || double.IsNegativeInfinity(x)) {
                return double.NaN;
            }
            
            if (x < 4.5) {
                if (Math.Abs(x - DigammaZero) < 0.25) { 
                    return DigammaNearRoot(x);
                }
                if (x < 0.25) {
                    double y = Digamma(1d - x) - Math.PI / TanPI(x);

                    return y;
                }
                if (x < 0.75) {
                    return DigammaNear0p5(x);
                }
                if (x < 1.5) {
                    return DigammaNear1(x);
                }
                if (x < 2.5) {
                    return DigammaNear2(x);
                }
                if (x < 3.5) {
                    return DigammaNear3(x);
                }
                else {
                    return DigammaNear4(x);
                }
            }
            if (x < 8.5) {
                if (x < 5.5) {
                    return DigammaNear5(x);
                }
                if (x < 6.5) {
                    return DigammaNear6(x);
                }
                if (x < 7.5) {
                    return DigammaNear7(x);
                }
                else {
                    return DigammaNear8(x);
                }
            }
            if (x < 12.5) {
                if (x < 9.5) {
                    return DigammaNear9(x);
                }
                if (x < 10.5) {
                    return DigammaNear10(x);
                }
                if (x < 11.5) {
                    return DigammaNear11(x);
                }
                else {
                    return DigammaNear12(x);
                }
            }
            if (x < 16.5) {
                if (x < 13.5) {
                    return DigammaNear13(x);
                }
                if (x < 14.5) {
                    return DigammaNear14(x);
                }
                if (x < 15.5) {
                    return DigammaNear15(x);
                }
                else {
                    return DigammaNear16(x);
                }
            }
            if (x < 20.5) {
                if (x < 17.5) {
                    return DigammaNear17(x);
                }
                if (x < 18.5) {
                    return DigammaNear18(x);
                }
                if (x < 19.5) {
                    return DigammaNear19(x);
                }
                else {
                    return DigammaNear20(x);
                }
            }
            if (x < 24.5) {
                if (x < 21.5) {
                    return DigammaNear21(x);
                }
                if (x < 22.5) {
                    return DigammaNear22(x);
                }
                if (x < 23.5) {
                    return DigammaNear23(x);
                }
                else {
                    return DigammaNear24(x);
                }
            }

            return DigammaLimit(x);
        }

        private static double DigammaNearRoot(double x) {
            double w = x - DigammaZero;

#if DEBUG
            if (w < -0.25 || w > 0.25) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                -9.24126552172942751679e-17, w, Fma(
                9.67672245447621103443e-1, w, Fma(
                1.14417380943934137200e0, w, Fma(
                4.71507845373433295698e-1, w, Fma(
                8.15984636391829491530e-2, w, Fma(
                5.50124017486102931828e-3, w,
                9.61671107604462865050e-5))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.63995297569876651027e0, w, Fma(
                9.70492711080937200640e-1, w, Fma(
                2.59707918978674903314e-1, w, Fma(
                3.16765151172736884890e-2, w, Fma(
                1.51426838914381237592e-3, w,
                1.70010400090620010203e-5))))));

            return y;
        }

        private static double DigammaNear0p5(double x) {
            double w = x - 0.5;

#if DEBUG
            if (w < -0.25 || w > 0.25) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                -1.96351002602142347944e0, w, Fma(
                -1.83958166615166356305e0, w, Fma(
                1.51730075678507209003e0, w, Fma(
                1.95525330344012965161e0, w, Fma(
                6.23462683530581086926e-1, w, Fma(
                6.81954718140085611071e-2, w,
                1.84784443198280621560e-3))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                3.45013968704961880129e0, w, Fma(
                3.61294709317636867840e0, w, Fma(
                1.56756259632741560316e0, w, Fma(
                2.94934802021506528650e-1, w, Fma(
                2.11363118706076315056e-2, w,
                3.50681082886986691533e-4))))));

            return y;
        }

        private static double DigammaNear1(double x) {
            double w = x - 1.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                -5.77215664901532860607e-1, w, Fma(
                2.52598982937360917596e-1, w, Fma(
                1.50936509640295219406e0, w, Fma(
                1.20532800369376795987e0, w, Fma(
                4.14909624982080297149e-1, w, Fma(
                7.06351062770466725398e-2, w, Fma(
                5.83881602630548368500e-3, w, Fma(
                2.01906837891786847234e-4, w,
                1.86350366774983395496e-6))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.41215748042524758184e0, w, Fma(
                2.17668731288860525169e0, w, Fma(
                9.66626214321677728020e-1, w, Fma(
                2.29430855276794470932e-1, w, Fma(
                2.91279029559269458753e-2, w, Fma(
                1.84097788336962884652e-3, w, Fma(
                4.78470201787889622613e-5, w,
                2.98674272351684867351e-7))))))));

            return y;
        }

        private static double DigammaNear2(double x) {
            double w = x - 2.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                4.22784335098467139394e-1, w, Fma(
                1.24525313770042489660e0, w, Fma(
                1.04116036862522883922e0, w, Fma(
                3.82396012007726909696e-1, w, Fma(
                6.88531823973700970273e-2, w, Fma(
                5.97118109946465816766e-3, w, Fma(
                2.15616636061591103405e-4, w,
                2.07574374014576233016e-6)))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.41991796056583600804e0, w, Fma(
                7.74541012287014034409e-1, w, Fma(
                2.06840963333788942014e-1, w, Fma(
                2.83626262046562062390e-2, w, Fma(
                1.89553320611896293198e-3, w, Fma(
                5.15231403725752668348e-5, w,
                3.34817161748632685229e-7)))))));

            return y;
        }

        private static double DigammaNear3(double x) {
            double w = x - 3.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                9.22784335098467139394e-1, w, Fma(
                1.25959780951528628163e0, w, Fma(
                5.95967665149065610155e-1, w, Fma(
                1.26557108168070125725e-1, w, Fma(
                1.24081435962089094404e-2, w, Fma(
                4.96037060109409052872e-4, w,
                5.24629987240131189040e-6))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                9.37016060827251206087e-1, w, Fma(
                3.28316154902781381034e-1, w, Fma(
                5.33973044114325275595e-2, w, Fma(
                4.03325218827777764994e-3, w, Fma(
                1.20935289093861540223e-4, w,
                8.58227475963243961442e-7))))));

            return y;
        }

        private static double DigammaNear4(double x) {
            double w = x - 4.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.25611766843180047273e0, w, Fma(
                1.20191946359942727820e0, w, Fma(
                4.19815823173294349734e-1, w, Fma(
                6.76151558963431902711e-2, w, Fma(
                5.10776014410146728221e-3, w, Fma(
                1.58796871955403379019e-4, w,
                1.31130878316236508930e-6))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                7.30900082799176861489e-1, w, Fma(
                2.00928204252060567133e-1, w, Fma(
                2.57619876319143148605e-2, w, Fma(
                1.53983978508570895794e-3, w, Fma(
                3.66424371382534354707e-5, w,
                2.06787487914784715751e-7))))));

            return y;
        }

        private static double DigammaNear5(double x) {
            double w = x - 5.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.50611766843180047273e0, w, Fma(
                9.90922754268388791598e-1, w, Fma(
                2.27953048185686112723e-1, w, Fma(
                2.21603033954495091703e-2, w, Fma(
                8.45478096630384806391e-4, w,
                8.41501154687051438227e-6)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                5.10982517941374400924e-1, w, Fma(
                9.24600753488674383585e-2, w, Fma(
                7.03186203236275637340e-3, w, Fma(
                2.03044431464879243795e-4, w,
                1.36386463158749740654e-6)))));

            return y;
        }

        private static double DigammaNear6(double x) {
            double w = x - 6.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.70611766843180047273e0, w, Fma(
                9.11734027593969346007e-1, w, Fma(
                1.72108504441107160117e-1, w, Fma(
                1.38186403423134952982e-2, w, Fma(
                4.37150342456275851005e-4, w,
                3.61238066985781471390e-6)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                4.28112952214029029972e-1, w, Fma(
                6.49877008663698653471e-2, w, Fma(
                4.15119560848517915381e-3, w, Fma(
                1.00774509864664398772e-4, w,
                5.69572846851354444315e-7)))));

            return y;
        }

        private static double DigammaNear7(double x) {
            double w = x - 7.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.87278433509846713939e0, w, Fma(
                8.42276462690944111007e-1, w, Fma(
                1.34581989881132360555e-1, w, Fma(
                9.18104098725102739134e-3, w, Fma(
                2.47344467161404618635e-4, w,
                1.74149025669570436177e-6)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                3.67757926966745208474e-1, w, Fma(
                4.79925896151135482275e-2, w, Fma(
                2.63728684819489886031e-3, w, Fma(
                5.51117069303840907463e-5, w,
                2.68278285032553478834e-7)))));

            return y;
        }

        private static double DigammaNear8(double x) {
            double w = x - 8.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.01564147795560999654e0, w, Fma(
                7.82113573071722430720e-1, w, Fma(
                1.08199870786626819451e-1, w, Fma(
                6.40613285597556769326e-3, w, Fma(
                1.50000288436937885735e-4, w,
                9.17952400717247806941e-7)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                3.21970234029875061635e-1, w, Fma(
                3.68039159824573819816e-2, w, Fma(
                1.77230424037595407238e-3, w, Fma(
                3.24685535906226539999e-5, w,
                1.38612396300550350445e-7)))));

            return y;
        }

        private static double DigammaNear9(double x) {
            double w = x - 9.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.14064147795560999654e0, w, Fma(
                6.11944274355184465426e-1, w, Fma(
                5.69115163179519793673e-2, w, Fma(
                1.89781157648830297366e-3, w,
                1.59764778917420568520e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.30973876173488763659e-1, w, Fma(
                1.71284966504604905512e-2, w, Fma(
                4.38601993092376073435e-4, w,
                2.52429287897219229788e-6))));

            return y;
        }

        private static double DigammaNear10(double x) {
            double w = x - 10.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.25175258906672110765e0, w, Fma(
                5.72164346817399228148e-1, w, Fma(
                4.73886466382085169272e-2, w, Fma(
                1.40888732880089833384e-3, w,
                1.05741545509924276769e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.07393127203761364415e-1, w, Fma(
                1.38127026215235017860e-2, w, Fma(
                3.17723940528208437900e-4, w,
                1.64295810260619186065e-6))));

            return y;
        }

        private static double DigammaNear11(double x) {
            double w = x - 11.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.35175258906672110765e0, w, Fma(
                5.37611202952677022145e-1, w, Fma(
                4.01145530501526731856e-2, w, Fma(
                1.07528517449443164539e-3, w,
                7.27516132204947010093e-6))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.88134104466564174020e-1, w, Fma(
                1.13683035031015421396e-2, w, Fma(
                2.37289615965861400940e-4, w,
                1.11360783536413717592e-6))));

            return y;
        }

        private static double DigammaNear12(double x) {
            double w = x - 12.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.44266167997581201674e0, w, Fma(
                5.07327255831244412602e-1, w, Fma(
                3.44286905128829142369e-2, w, Fma(
                8.39814807260406637238e-4, w,
                5.16934815223207606372e-6))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.72117729772401064816e-1, w, Fma(
                9.51623399846104186844e-3, w, Fma(
                1.81765370547135472298e-4, w,
                7.80687940754933719941e-7))));

            return y;
        }

        private static double DigammaNear13(double x) {
            double w = x - 13.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.52599501330914535007e0, w, Fma(
                4.80565842822420569743e-1, w, Fma(
                2.98968199068732640898e-2, w, Fma(
                6.68810518340851205988e-4, w,
                3.77432733273561102451e-6))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.58594301368111185824e-1, w, Fma(
                8.08035102187623323385e-3, w, Fma(
                1.42239250887779257393e-4, w,
                5.63077092529842557697e-7))));

            return y;
        }

        private static double DigammaNear14(double x) {
            double w = x - 14.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.60291809023222227315e0, w, Fma(
                4.56740482682184314643e-1, w, Fma(
                2.62239663588130540079e-2, w, Fma(
                5.41577291719578059431e-4, w,
                2.82056830329849669889e-6))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.47027374950562100315e-1, w, Fma(
                6.94518699205500904545e-3, w, Fma(
                1.13356899997680659101e-4, w,
                4.16104158480291521609e-7))));

            return y;
        }

        private static double DigammaNear15(double x) {
            double w = x - 15.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.67434666166079370172e0, w, Fma(
                4.35386107019977009524e-1, w, Fma(
                2.32040360810724857473e-2, w, Fma(
                4.44915735469093712141e-4, w,
                2.15057448454754357251e-6))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.37023327762873390378e-1, w, Fma(
                6.03257321904729416617e-3, w, Fma(
                9.17729049148988939280e-5, w,
                3.14008486976913988637e-7))));

            return y;
        }

        private static double DigammaNear16(double x) {
            double w = x - 16.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.74101332832746036839e0, w, Fma(
                4.16130160005833134206e-1, w, Fma(
                2.06894757728644074611e-2, w, Fma(
                3.70135719065377457587e-4, w,
                1.66870646309874250632e-6))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.28286999909321442079e-1, w, Fma(
                5.28810520526920944931e-3, w, Fma(
                7.53255195002857058099e-5, w,
                2.41334628546542932513e-7))));

            return y;
        }

        private static double DigammaNear17(double x) {
            double w = x - 17.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.80351332832746036839e0, w, Fma(
                3.98671170250165009237e-1, w, Fma(
                1.85723628103509406266e-2, w, Fma(
                3.11355764349987080808e-4, w,
                1.31491657491845978479e-6))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.20592840929571022656e-1, w, Fma(
                4.67299525097470331501e-3, w, Fma(
                6.25764587337145270927e-5, w,
                1.88485962021280533412e-7))));

            return y;
        }

        private static double DigammaNear18(double x) {
            double w = x - 18.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.86233685773922507427e0, w, Fma(
                3.82762813960556502181e-1, w, Fma(
                1.67722594550886627795e-2, w, Fma(
                2.64497708741952020901e-4, w,
                1.05038555518959310885e-6))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.13765606339909383701e-1, w, Fma(
                4.15899944775983371548e-3, w, Fma(
                5.25439799305198623281e-5, w,
                1.49321943304884325582e-7))));

            return y;
        }

        private static double DigammaNear19(double x) {
            double w = x - 19.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.91789241329478062982e0, w, Fma(
                3.68201978781885967270e-1, w, Fma(
                1.52282014622783492813e-2, w, Fma(
                2.26672583097258135792e-4, w,
                8.49361886287754105096e-7))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.07667120046228925572e-1, w, Fma(
                3.72516400661054536562e-3, w, Fma(
                4.45426532079604458533e-5, w,
                1.19807864167214172389e-7))));

            return y;
        }

        private static double DigammaNear20(double x) {
            double w = x - 20.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.97052399224214905088e0, w, Fma(
                3.54819745024736225627e-1, w, Fma(
                1.38932874147960162547e-2, w, Fma(
                1.95794445516301226649e-4, w,
                6.94366415224113279569e-7))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.02186995588079607001e-1, w, Fma(
                3.35568313365231898554e-3, w, Fma(
                3.80842143164442368235e-5, w,
                9.72294541686021407096e-8))));

            return y;
        }

        private static double DigammaNear21(double x) {
            double w = x - 21.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.02052399224214905088e0, w, Fma(
                3.42474506187406621194e-1, w, Fma(
                1.27309464788206261092e-2, w, Fma(
                1.70331218563754440403e-4, w,
                5.73288918252885152945e-7))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                9.72360040862267368018e-2, w, Fma(
                3.03845550154005963494e-3, w, Fma(
                3.28145380782015835875e-5, w,
                7.97218376252323008711e-8))));

            return y;
        }

        private static double DigammaNear22(double x) {
            double w = x - 22.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.06814303986119666992e0, w, Fma(
                3.31046670541034551198e-1, w, Fma(
                1.17123145332018294833e-2, w, Fma(
                1.49140172212490159197e-4, w,
                4.77582537275256747261e-7))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                9.27412501976600642601e-2, w, Fma(
                2.76409099956524999002e-3, w, Fma(
                2.84725957414889543024e-5, w,
                6.59791928483540076830e-8))));

            return y;
        }

        private static double DigammaNear23(double x) {
            double w = x - 23.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.11359758531574212447e0, w, Fma(
                3.20434541634982588226e-1, w, Fma(
                1.08143575557367396864e-2, w, Fma(
                1.31357016413867465871e-4, w,
                4.01113630645120428681e-7))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                8.86426073170325017455e-2, w, Fma(
                2.52521435835779052301e-3, w, Fma(
                2.48631423269347265946e-5, w,
                5.50713389424288255884e-8))));

            return y;
        }

        private static double DigammaNear24(double x) {
            double w = x - 24.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.15707584618530734186e0, w, Fma(
                3.10551086042388441977e-1, w, Fma(
                1.00185083524552990159e-2, w, Fma(
                1.16319714871342597970e-4, w,
                3.39413028873614935376e-7))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                8.48900453240238698495e-2, w, Fma(
                2.31596855012054448778e-3, w, Fma(
                2.18381809749460847191e-5, w,
                4.63252681252297658389e-8))));

            return y;
        }

        private static double DigammaLimit(double x) {
            double s = 1.0 / x, t = s * s;

            double v = t * Fma(
                8.33333333333333333333e-2, t, Fma(
                -8.33333333333333333334e-3, t, Fma(
                3.96825396825396825397e-3, t, Fma(
                -4.16666666666666666667e-3, t, Fma(
                7.57575757575757575758e-3, t,
                -2.10927960927960927961e-2)))));

            double a = Math.Log(x);
            double b = s * 0.5;

            double y = a - b - v;

            return y;
        }
    }
}