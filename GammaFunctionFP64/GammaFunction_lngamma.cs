// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

namespace GammaFunctionFP64 {
    public static partial class GammaFunction {
        public static double LogGamma(double x) {
            if (double.IsNaN(x) || x < 0) {
                return double.NaN;
            }
            if (x == 0) {
                return double.PositiveInfinity;
            }

            if (x < 4.5) {
                if (x < 1.875) {
                    if (x < 0.375) {
                        return Math.Log(Gamma(x));
                    }
                    if (x < 0.625) {
                        return LogGammaNear0p50(x);
                    }
                    if (x < 0.875) {
                        return LogGammaNear0p75(x);
                    }
                    if (x < 1.125) {
                        return LogGammaNear1p00(x);
                    }
                    if (x < 1.375) {
                        return LogGammaNear1p25(x);
                    }
                    if (x < 1.625) {
                        return LogGammaNear1p50(x);
                    }
                    else {
                        return LogGammaNear1p75(x);
                    }
                }
                if (x < 2.5) {
                    return LogGammaNear2(x);
                }
                if (x < 3.5) {
                    return LogGammaNear3(x);
                }
                else {
                    return LogGammaNear4(x);
                }
            }
            if (x < 8.5) {
                if (x < 5.5) {
                    return LogGammaNear5(x);
                }
                if (x < 6.5) {
                    return LogGammaNear6(x);
                }
                if (x < 7.5) {
                    return LogGammaNear7(x);
                }
                else {
                    return LogGammaNear8(x);
                }
            }
            if (x < 12.5) {
                if (x < 9.5) {
                    return LogGammaNear9(x);
                }
                if (x < 10.5) {
                    return LogGammaNear10(x);
                }
                if (x < 11.5) {
                    return LogGammaNear11(x);
                }
                else {
                    return LogGammaNear12(x);
                }
            }
            if (x < 16.5) {
                if (x < 13.5) {
                    return LogGammaNear13(x);
                }
                if (x < 14.5) {
                    return LogGammaNear14(x);
                }
                if (x < 15.5) {
                    return LogGammaNear15(x);
                }
                else {
                    return LogGammaNear16(x);
                }
            }
            if (x < 20.5) {
                if (x < 17.5) {
                    return LogGammaNear17(x);
                }
                if (x < 18.5) {
                    return LogGammaNear18(x);
                }
                if (x < 19.5) {
                    return LogGammaNear19(x);
                }
                else {
                    return LogGammaNear20(x);
                }
            }
            if (x < 24.5) {
                if (x < 21.5) {
                    return LogGammaNear21(x);
                }
                if (x < 22.5) {
                    return LogGammaNear22(x);
                }
                if (x < 23.5) {
                    return LogGammaNear23(x);
                }
                else {
                    return LogGammaNear24(x);
                }
            }

            return LogGammaLimit(x);
        }

        private static double LogGammaNear0p50(double x) {
            double w = x - 0.50;

#if DEBUG
            if (w < -0.125 || w > 0.125) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                5.72364942924700087072e-1, w, Fma(
                1.76274975930175720731e0, w, Fma(
                -9.42016134520741232048e-1, w, Fma(
                -7.50677018009491334622e0, w, Fma(
                -6.08342991220481330994e0, w, Fma(
                1.63115195313417291601e0, w, Fma(
                2.86558568287474980506e0, w,
                5.85802144954420781179e-1)))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                6.51028654250301401939e0, w, Fma(
                1.63769563100532498015e1, w, Fma(
                1.99013918197879721582e1, w, Fma(
                1.18560167772445248741e1, w, Fma(
                3.04961852869671560670e0, w, Fma(
                2.13627592464765158113e-1, w,
                -3.81357795097451111612e-3)))))));

            return y;
        }

        private static double LogGammaNear0p75(double x) {
            double w = x - 0.75;

#if DEBUG
            if (w < -0.125 || w > 0.125) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.03280951431295371481e-1, w, Fma(
                -3.48698679927728318292e-1, w, Fma(
                -1.69187431336532354929e0, w, Fma(
                -9.16163360173121547804e-1, w, Fma(
                8.23625465661105641747e-1, w, Fma(
                7.21569152525773283703e-1, w,
                1.18666395649547975799e-1))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                3.62632206642287859321e0, w, Fma(
                4.79566395852590364602e0, w, Fma(
                2.78528605586054934412e0, w, Fma(
                6.61299595315117680597e-1, w, Fma(
                4.17717990592655621144e-2, w,
                -6.20683971151943902869e-4))))));

            return y;
        }

        private static double LogGammaNear1p00(double x) {
            double w = x - 1.0;

#if DEBUG
            if (w < -0.125 || w > 0.125) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = w * Fma(
                -5.77215664901532860607e-1, w, Fma(
                -9.08286876739912372625e-1, w, Fma(
                1.02780689563411686947e-1, w, Fma(
                8.23805420232929772298e-1, w, Fma(
                4.67903059609853915943e-1, w, Fma(
                8.69877581020716087047e-2, w,
                4.02560838654105430412e-3))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.99845277147714738108e0, w, Fma(
                3.40022343677035696557e0, w, Fma(
                1.80506210009889151132e0, w, Fma(
                4.47351944512267540725e-1, w, Fma(
                4.40712554538182199479e-2, w,
                1.02438390986121967437e-3))))));

            return y;
        }

        private static double LogGammaNear1p25(double x) {
            double w = x - 1.25;

#if DEBUG
            if (w < -0.125 || w > 0.125) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                -9.82718364218131614639e-2, w, Fma(
                -4.38099071604208491053e-1, w, Fma(
                -5.30406091078421511993e-2, w, Fma(
                6.25972810090628916993e-1, w, Fma(
                5.04146932102961014759e-1, w, Fma(
                1.31923101431277783301e-1, w,
                9.98993052006533186353e-3))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.14349854340553061181e0, w, Fma(
                1.67045691679466194265e0, w, Fma(
                5.70077080203670826693e-1, w, Fma(
                7.95579140468213426621e-2, w, Fma(
                3.01096497841799463410e-3, w,
                -2.11406891560420621612e-5))))));

            return y;
        }

        private static double LogGammaNear1p50(double x) {
            double w = x - 1.50;

#if DEBUG
            if (w < -0.125 || w > 0.125) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                -1.20782237635245222346e-1, w, Fma(
                -1.37566760698731677312e-1, w, Fma(
                4.36651936882752200682e-1, w, Fma(
                5.45924617983819471147e-1, w, Fma(
                1.85908053825366729369e-1, w,
                1.73697014788622832547e-2)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.44107890435801175677e0, w, Fma(
                6.89953231056627737169e-1, w, Fma(
                1.21548078495995400765e-1, w, Fma(
                5.49775770062710077104e-3, w,
                -4.81168528646484739527e-5)))));

            return y;
        }

        private static double LogGammaNear1p75(double x) {
            double w = x - 1.75;

#if DEBUG
            if (w < -0.125 || w > 0.125) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                -8.44011210204855559578e-2, w, Fma(
                1.43819901728416524274e-1, w, Fma(
                6.43709329080319467315e-1, w, Fma(
                4.93117199598983456474e-1, w, Fma(
                1.27989681501814207218e-1, w,
                9.62210177411197904907e-3)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.22809449169859297249e0, w, Fma(
                5.00718028429664201101e-1, w, Fma(
                7.51485546489315584334e-2, w, Fma(
                2.91105609847390276606e-3, w,
                -2.04476157549167520614e-5)))));

            return y;
        }

        private static double LogGammaNear2(double x) {
            double w = x - 2.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = w * Fma(
                4.22784335098467139394e-1, w, Fma(
                1.05000085049473750916e0, w, Fma(
                9.81253367349466482875e-1, w, Fma(
                4.48612936190413778215e-1, w, Fma(
                1.06617723221536780904e-1, w, Fma(
                1.26787174098271901045e-2, w, Fma(
                6.46123281924455599896e-4, w,
                9.04485505477592573373e-6)))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.72081545287428995176e0, w, Fma(
                1.16773345949285709047e0, w, Fma(
                3.95893248149539042506e-1, w, Fma(
                6.99517764734187788468e-2, w, Fma(
                6.08267140325837670731e-3, w, Fma(
                2.16059117399485166600e-4, w,
                1.83240727523092522095e-6)))))));

            return y;
        }

        private static double LogGammaNear3(double x) {
            double w = x - 3.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                6.93147180559945309417e-1, w, Fma(
                1.64822499587290081976e0, w, Fma(
                1.45513143714130124409e0, w, Fma(
                6.26006619924649985446e-1, w, Fma(
                1.41688236402187754650e-1, w, Fma(
                1.66223939905079946789e-2, w, Fma(
                9.03454429692092654296e-4, w,
                1.62987450796249325198e-5)))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.04658964375848823112e0, w, Fma(
                4.21105190018624902482e-1, w, Fma(
                8.14199796965950217686e-2, w, Fma(
                7.68562955400308296033e-3, w, Fma(
                3.14018129862454594998e-4, w, Fma(
                3.67081281704172084269e-6, w,
                -5.48424001476313559479e-9)))))));

            return y;
        }

        private static double LogGammaNear4(double x) {
            double w = x - 4.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.79175946922805500081e0, w, Fma(
                2.42489606547589785701e0, w, Fma(
                1.23783039591498011500e0, w, Fma(
                3.01723753281378246993e-1, w, Fma(
                3.63500625453952220588e-2, w, Fma(
                1.98030224608904033096e-3, w,
                3.54652303015663785912e-5))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                6.52307643473843622322e-1, w, Fma(
                1.54342012181444967112e-1, w, Fma(
                1.59742638196502398073e-2, w, Fma(
                6.77529922848484073284e-4, w, Fma(
                7.99122190377047366703e-6, w,
                -1.18779889415294899335e-8))))));

            return y;
        }

        private static double LogGammaNear5(double x) {
            double w = x - 5.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.17805383034794561965e0, w, Fma(
                2.83792038450225676077e0, w, Fma(
                9.26586780390846139775e-1, w, Fma(
                1.35248527089471094724e-1, w, Fma(
                8.52112732316600627607e-3, w,
                1.73149238073185636393e-4)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                4.19062352988730017189e-1, w, Fma(
                5.81387535815053033767e-2, w, Fma(
                2.97109059313005100048e-3, w, Fma(
                4.00954142801902376854e-5, w,
                -6.88018736340780488262e-8)))));

            return y;
        }

        private static double LogGammaNear6(double x) {
            double w = x - 6.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                4.78749174278204599425e0, w, Fma(
                3.37107426552932052381e0, w, Fma(
                8.75725698119438005654e-1, w, Fma(
                1.02526083667651869460e-1, w, Fma(
                5.21419526307899512647e-3, w,
                8.59111298512121122400e-5)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                3.47772212789238516699e-1, w, Fma(
                4.00468374092953735074e-2, w, Fma(
                1.69961220553540175398e-3, w, Fma(
                1.91091368374610057764e-5, w,
                -2.59988510697487303277e-8)))));

            return y;
        }

        private static double LogGammaNear7(double x) {
            double w = x - 7.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                6.57925121201010099506e0, w, Fma(
                3.82654221657386515370e0, w, Fma(
                8.25048723562245785762e-1, w, Fma(
                8.05363594317671370051e-2, w, Fma(
                3.42781032081579712940e-3, w,
                4.73928710952802790584e-5)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.96957483232880460536e-1, w, Fma(
                2.92037506400474236389e-2, w, Fma(
                1.05902376659721286209e-3, w, Fma(
                1.01996999733867004179e-5, w,
                -1.14214833227813368168e-8)))));

            return y;
        }

        private static double LogGammaNear8(double x) {
            double w = x - 8.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                8.52516136106541430017e0, w, Fma(
                4.22322775473684592725e0, w, Fma(
                7.77860351464628842457e-1, w, Fma(
                6.50458416149321199921e-2, w, Fma(
                2.37727678406673479093e-3, w,
                2.82710118685253538332e-5)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.58949500576414594665e-1, w, Fma(
                2.22098423760789731157e-2, w, Fma(
                7.02720842080342449481e-4, w, Fma(
                5.91756866530098086969e-6, w,
                -5.60351451948032738879e-9)))));

            return y;
        }

        private static double LogGammaNear9(double x) {
            double w = x - 9.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.06046029027452502284e1, w, Fma(
                4.57415859034055476749e0, w, Fma(
                7.34974614315377781697e-1, w, Fma(
                5.37157548803008546411e-2, w, Fma(
                1.71854723852444726337e-3, w,
                1.79100573113048319769e-5)))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.29477438684193613699e-1, w, Fma(
                1.74442819936306207028e-2, w, Fma(
                4.89363892122983333449e-4, w, Fma(
                3.66013928344683552706e-6, w,
                -2.99176983877679129477e-9)))));

            return y;
        }

        private static double LogGammaNear10(double x) {
            double w = x - 10.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.28018274800814696112e1, w, Fma(
                4.22307301587780553648e0, w, Fma(
                4.84590943965715231487e-1, w, Fma(
                2.21414367367920507759e-2, w,
                3.18482878516380123110e-4))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.53987423270489123427e-1, w, Fma(
                6.66047071986419161219e-3, w, Fma(
                6.93804996908471561939e-5, w,
                -8.17677227917308108757e-8))));

            return y;
        }

        private static double LogGammaNear11(double x) {
            double w = x - 11.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.51044125730755152952e1, w, Fma(
                4.46004855038640382151e0, w, Fma(
                4.58531520457467715476e-1, w, Fma(
                1.87900477089924616293e-2, w,
                2.42596457160591350959e-4))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.39581460127608113813e-1, w, Fma(
                5.47437988572981767285e-3, w, Fma(
                5.17887351893501754155e-5, w,
                -5.41031514943323769012e-8))));

            return y;
        }

        private static double LogGammaNear12(double x) {
            double w = x - 12.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.75023078458738858393e1, w, Fma(
                4.67645911051682693223e0, w, Fma(
                4.35335515085076958813e-1, w, Fma(
                1.61652443068101723709e-2, w,
                1.89225591185571939156e-4))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.27628736176505183066e-1, w, Fma(
                4.57829654424639161963e-3, w, Fma(
                3.96687250178703033204e-5, w,
                -3.71451726892869419271e-8))));

            return y;
        }

        private static double LogGammaNear13(double x) {
            double w = x - 13.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.99872144956618861495e1, w, Fma(
                4.87557489985071113322e0, w, Fma(
                4.14571154190595978255e-1, w, Fma(
                1.40689220203003117681e-2, w,
                1.50566975428393339138e-4))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.17554143777840030808e-1, w, Fma(
                3.88504656399907874843e-3, w, Fma(
                3.10506789980880924800e-5, w,
                -2.63036780698309113797e-8))));

            return y;
        }

        private static double LogGammaNear14(double x) {
            double w = x - 14.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.25521638531234228856e1, w, Fma(
                5.05994785399323505489e0, w, Fma(
                3.95880450438922587000e-1, w, Fma(
                1.23665213458643378515e-2, w,
                1.21857841986620818691e-4))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.08948736793641192168e-1, w, Fma(
                3.33784725443870569872e-3, w, Fma(
                2.47570813517491200337e-5, w,
                -1.91224076950813541146e-8))));

            return y;
        }

        private static double LogGammaNear15(double x) {
            double w = x - 15.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.51912211827386815001e1, w, Fma(
                5.23160631775986209144e0, w, Fma(
                3.78967936881866425678e-1, w, Fma(
                1.09639848164520790522e-2, w,
                1.00076935122792864797e-4))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.01513921756652770979e-1, w, Fma(
                2.89844643363522911991e-3, w, Fma(
                2.00548617524927994828e-5, w,
                -1.42193634257410935042e-8))));

            return y;
        }

        private static double LogGammaNear16(double x) {
            double w = x - 16.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.78992713838408915661e1, w, Fma(
                5.39218839168732518627e0, w, Fma(
                3.63589445354946675726e-1, w, Fma(
                9.79391686149728330412e-3, w,
                8.32434746803225803754e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                9.50266774671187722553e-2, w, Fma(
                2.54032312149660812618e-3, w, Fma(
                1.64715426139192602800e-5, w,
                -1.07830443141362245533e-8))));

            return y;
        }

        private static double LogGammaNear17(double x) {
            double w = x - 17.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.06718601060806728038e1, w, Fma(
                5.54303559958914465739e0, w, Fma(
                3.49542383463865503441e-1, w, Fma(
                8.80697167817834039644e-3, w,
                7.00212609033063876554e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                8.93171220065188056302e-2, w, Fma(
                2.24462665554058099492e-3, w, Fma(
                1.36933447552971871476e-5, w,
                -8.31907461858460799889e-9))));

            return y;
        }

        private static double LogGammaNear18(double x) {
            double w = x - 18.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.35050734501368888840e1, w, Fma(
                5.68526016955897586439e0, w, Fma(
                3.36657749824160710481e-1, w, Fma(
                7.96632640286458886004e-3, w,
                5.94861796708135602139e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                8.42536076221679617836e-2, w, Fma(
                1.99766404135628237764e-3, w, Fma(
                1.15063652538934225923e-5, w,
                -6.51646398723553245217e-9))));

            return y;
        }

        private static double LogGammaNear19(double x) {
            double w = x - 19.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.63954452080330535762e1, w, Fma(
                5.81979433060015902094e0, w, Fma(
                3.24793724089157793309e-1, w, Fma(
                7.24403678543767838053e-3, w,
                5.09849329043915355198e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                7.97325572119909797741e-2, w, Fma(
                1.78929662535013184569e-3, w, Fma(
                9.76139372652218715543e-6, w,
                -5.17399323286538274656e-9))));

            return y;
        }

        private static double LogGammaNear20(double x) {
            double w = x - 20.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.93398841871994940362e1, w, Fma(
                5.94742707810639485729e0, w, Fma(
                3.13830566576490591130e-1, w, Fma(
                6.61857575360235237572e-3, w,
                4.40465888436927865092e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                7.56713739089470344558e-2, w, Fma(
                1.61188891875457742839e-3, w, Fma(
                8.35220154589895508774e-6, w,
                -4.15817951617463639600e-9))));

            return y;
        }

        private static double LogGammaNear21(double x) {
            double w = x - 21.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                4.23356164607534850297e1, w, Fma(
                6.06883203299853415644e0, w, Fma(
                3.03666567643001645590e-1, w, Fma(
                6.07313493274965546252e-3, w,
                3.83257132392134191709e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                7.20033932559424841941e-2, w, Fma(
                1.45960268182851620360e-3, w, Fma(
                7.20177006369346803046e-6, w,
                -3.37848232901374596572e-9))));

            return y;
        }

        private static double LogGammaNear22(double x) {
            double w = x - 22.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                4.53801388984769080262e1, w, Fma(
                6.18458885657127674530e0, w, Fma(
                2.94214824915709802967e-1, w, Fma(
                5.59442989816328981732e-3, w,
                3.35649290078279919669e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                6.86742238423311047759e-2, w, Fma(
                1.32791260360124884730e-3, w, Fma(
                6.25333723950446977703e-6, w,
                -2.77225036445119036978e-9))));

            return y;
        }

        private static double LogGammaNear23(double x) {
            double w = x - 23.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                4.84711813518352238796e1, w, Fma(
                6.29519992641136212814e0, w, Fma(
                2.85400669307642946462e-1, w, Fma(
                5.17184507228078013623e-3, w,
                2.95697369780526193087e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                6.56390509239188916525e-2, w, Fma(
                1.21326756311863869086e-3, w, Fma(
                5.46444694644260946067e-6, w,
                -2.29534279223649900615e-9))));

            return y;
        }

        private static double LogGammaNear24(double x) {
            double w = x - 24.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                5.16066755677643735704e1, w, Fma(
                6.40110347584193822235e0, w, Fma(
                2.77159599014788155838e-1, w, Fma(
                4.79681202174865229543e-3, w,
                2.61912540909426876192e-5))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                6.28606201419992716702e-2, w, Fma(
                1.11284956204304871514e-3, w, Fma(
                4.80291440460830058975e-6, w,
                -1.91614807585740017884e-9))));

            return y;
        }

        private static double LogGammaLimit(double x) {
            double s = 1.0 / x, t = s * s;

            double v = s * Fma(
                8.3333333333333333333e-2, t, Fma(
                -2.7777777777777777778e-3, t, Fma(
                7.9365079365079365079e-4, t, Fma(
                -5.952380952380952381e-4, t, Fma(
                8.4175084175084175084e-4, t,
                -1.9175269175269175269e-3)))));

            double a = 9.18938533204672741781e-1;
            double b = (x - 0.5) * Math.Log(x) - x;

            double y = a + b + v;

            return y;
        }
    }
}