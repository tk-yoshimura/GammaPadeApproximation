// Copyright (c) T.Yoshimura 2022
// https://github.com/tk-yoshimura

using System.Runtime.CompilerServices;

namespace GammaFunctionFP64 {
    public static partial class GammaFunction {
        public static double Gamma(double x) {
            if (double.IsNaN(x) || double.IsNegativeInfinity(x)) {
                return double.NaN;
            }
            if (x == 0 || x >= 171.625) {
                return double.PositiveInfinity;
            }

            if (x < 0.5) {
                if (Math.Truncate(x) == x) {
                    return double.NaN;
                }

                double t = Theta(x);
                double y = Math.PI / (Math.Sin(t * Math.PI) * Gamma(1d - x));

                return y;
            }
            if (x < 4.5) {
                if (x < 1.5) {
                    return GammaNear1(x);
                }
                if (x < 2.5) {
                    return GammaNear2(x);
                }
                if (x < 3.5) {
                    return GammaNear3(x);
                }
                else {
                    return GammaNear4(x);
                }
            }
            if (x < 8.5) {
                if (x < 5.5) {
                    return GammaNear5(x);
                }
                if (x < 6.5) {
                    return GammaNear6(x);
                }
                if (x < 7.5) {
                    return GammaNear7(x);
                }
                else {
                    return GammaNear8(x);
                }
            }
            if (x < 12.5) {
                if (x < 9.5) {
                    return GammaNear9(x);
                }
                if (x < 10.5) {
                    return GammaNear10(x);
                }
                if (x < 11.5) {
                    return GammaNear11(x);
                }
                else {
                    return GammaNear12(x);
                }
            }
            if (x < 16.5) {
                if (x < 13.5) {
                    return GammaNear13(x);
                }
                if (x < 14.5) {
                    return GammaNear14(x);
                }
                if (x < 15.5) {
                    return GammaNear15(x);
                }
                else {
                    return GammaNear16(x);
                }
            }
            if (x < 20.5) {
                if (x < 17.5) {
                    return GammaNear17(x);
                }
                if (x < 18.5) {
                    return GammaNear18(x);
                }
                if (x < 19.5) {
                    return GammaNear19(x);
                }
                else {
                    return GammaNear20(x);
                }
            }
            if (x < 24.5) {
                if (x < 21.5) {
                    return GammaNear21(x);
                }
                if (x < 22.5) {
                    return GammaNear22(x);
                }
                if (x < 23.5) {
                    return GammaNear23(x);
                }
                else {
                    return GammaNear24(x);
                }
            }

            if (x <= 170.625) {
                return GammaLimit(x);
            }
            else {
                return GammaExtremeLarge(x);
            }
        }

        private static double GammaNear1(double x) {
            double w = x - 1.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.00000000000000000000e0, w, Fma(
                7.44911371245221012548e-1, w, Fma(
                3.38201312691908083967e-1, w, Fma(
                1.04119303551766298426e-1, w, Fma(
                2.47080687424251963882e-2, w, Fma(
                4.50527222853504011495e-3, w, Fma(
                6.55736910244071441396e-4, w, Fma(
                7.15522156225798597022e-5, w, Fma(
                5.75266978105177154926e-6, w,
                2.44426356208484719181e-7)))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                1.32212703614675387315e0, w, Fma(
                1.12297753617677035891e-1, w, Fma(
                -2.31239269532121445724e-1, w, Fma(
                -1.76109205971347748613e-3, w, Fma(
                1.81310130440256331024e-2, w, Fma(
                -2.05448182537120808469e-3, w, Fma(
                -4.22677168398699062544e-4, w, Fma(
                1.07127608695013124298e-4, w,
                -6.79514988445133153244e-6)))))))));

            return y;
        }

        private static double GammaNear2(double x) {
            double w = x - 2.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.00000000000000000000e0, w, Fma(
                7.08947211249356502985e-1, w, Fma(
                3.10190819869412358300e-1, w, Fma(
                9.06517080384026018481e-2, w, Fma(
                2.03635393210163369744e-2, w, Fma(
                3.41890297997255512848e-3, w, Fma(
                4.51252372517573980218e-4, w, Fma(
                4.07513512581266609766e-5, w,
                2.54996497810593499350e-6))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                2.86162876150889363592e-1, w, Fma(
                -2.22634691880346095549e-1, w, Fma(
                -1.46521645019620545759e-2, w, Fma(
                2.06548934233331053419e-2, w, Fma(
                -2.09778625130634035372e-3, w, Fma(
                -5.20315883991121309839e-4, w, Fma(
                1.28964433154387247686e-4, w,
                -8.40791722737135184029e-6))))))));

            return y;
        }

        private static double GammaNear3(double x) {
            double w = x - 3.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.00000000000000000000e0, w, Fma(
                1.39337734236528765608e0, w, Fma(
                6.04095151272401423469e-1, w, Fma(
                1.74685164820137351773e-1, w, Fma(
                3.89656989452367941910e-2, w, Fma(
                6.48796329187210944188e-3, w, Fma(
                8.53001289741259962951e-4, w, Fma(
                7.65257696964918857420e-5, w,
                4.82476144303238622349e-6))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -2.26095663915823311352e-1, w, Fma(
                -1.12547385444263049153e-1, w, Fma(
                4.46126275977445932252e-2, w, Fma(
                -1.57742244324347490814e-3, w, Fma(
                -1.59544186166538522063e-3, w, Fma(
                3.15788005495380004019e-4, w, Fma(
                -2.18891230350873854047e-5, w,
                3.32603623763660750489e-7))))))));

            return y;
        }

        private static double GammaNear4(double x) {
            double w = x - 4.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                6.00000000000000000000e0, w, Fma(
                4.52493387757083287269e0, w, Fma(
                1.93093531773874771960e0, w, Fma(
                5.50195487362800293975e-1, w, Fma(
                1.15951899406011003567e-1, w, Fma(
                1.81980021551047021503e-2, w, Fma(
                2.14106760669559715374e-3, w, Fma(
                1.72269523193409352773e-4, w,
                8.30893379258166840405e-6))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -5.01962022169994993946e-1, w, Fma(
                2.15186415441148254811e-2, w, Fma(
                3.66679935209228677511e-2, w, Fma(
                -9.04166878383391116360e-3, w, Fma(
                2.47213236307537973364e-4, w, Fma(
                2.08240558241630435061e-4, w, Fma(
                -3.45857239575982434490e-5, w,
                1.86408074957083462636e-6))))))));

            return y;
        }

        private static double GammaNear5(double x) {
            double w = x - 5.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.40000000000000000000e1, w, Fma(
                1.68904920732362593792e1, w, Fma(
                6.82927283507664422108e0, w, Fma(
                1.84155651418729151358e0, w, Fma(
                3.68690726763262773945e-1, w, Fma(
                5.48742576644790623885e-2, w, Fma(
                6.14346915041939487417e-3, w, Fma(
                4.68929204503566430166e-4, w,
                2.17697460586845221793e-5))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -8.02347165380289665259e-1, w, Fma(
                2.48125583340507169187e-1, w, Fma(
                -2.61156338781585557437e-2, w, Fma(
                -4.79853635551869913993e-3, w, Fma(
                2.01013052872769920112e-3, w, Fma(
                -3.00243344934956063355e-4, w, Fma(
                2.25855888344283596750e-5, w,
                -7.15707222831102189665e-7))))))));

            return y;
        }

        private static double GammaNear6(double x) {
            double w = x - 6.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.20000000000000000000e2, w, Fma(
                9.34217024749973640151e1, w, Fma(
                3.88690883606056799009e1, w, Fma(
                1.06623739336809535875e1, w, Fma(
                2.10808642800666379641e0, w, Fma(
                3.05862893269188548891e-1, w, Fma(
                3.22126327458339575025e-2, w, Fma(
                2.26916742357825729829e-3, w,
                8.69755084189571533853e-5))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -9.27603481140155772601e-1, w, Fma(
                3.60429531007638978397e-1, w, Fma(
                -6.88520823768772222173e-2, w, Fma(
                3.70825546561125303162e-3, w, Fma(
                1.20544243062654516411e-3, w, Fma(
                -3.03432903458283237686e-4, w, Fma(
                3.04718094245288708385e-5, w,
                -1.25586753838294830173e-6))))))));

            return y;
        }

        private static double GammaNear7(double x) {
            double w = x - 7.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                7.20000000000000000000e2, w, Fma(
                5.38314647777797956117e2, w, Fma(
                2.16609495021229984064e2, w, Fma(
                5.75277010650204209679e1, w, Fma(
                1.10319923652426023563e1, w, Fma(
                1.55320720570140127463e0, w, Fma(
                1.58959996154975655173e-1, w, Fma(
                1.08851776430491155353e-2, w,
                4.07607315804536161501e-4))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.12512510207374775590e0, w, Fma(
                5.77530015177825482093e-1, w, Fma(
                -1.76822717400539237849e-1, w, Fma(
                3.53320973601991119156e-2, w, Fma(
                -4.71370988036609138665e-3, w, Fma(
                4.08520705445062789404e-4, w, Fma(
                -2.08200365856950514597e-5, w,
                4.65396969700973760907e-7))))))));

            return y;
        }

        private static double GammaNear8(double x) {
            double w = x - 8.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                5.04000000000000000000e3, w, Fma(
                3.96761969935028866504e3, w, Fma(
                1.60077994774981071359e3, w, Fma(
                4.18492924955670126697e2, w, Fma(
                7.71610279488244398559e1, w, Fma(
                1.02429803101487709388e1, w, Fma(
                9.59830420463809004899e-1, w, Fma(
                5.83218690006924504499e-2, w,
                1.79268976531152659771e-3))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.22841534713214002332e0, w, Fma(
                6.95686203660205110194e-1, w, Fma(
                -2.38126353480685766258e-1, w, Fma(
                5.41144443578088126154e-2, w, Fma(
                -8.40249964512875461132e-3, w, Fma(
                8.75720023183541034815e-4, w, Fma(
                -5.64089289156022730943e-5, w,
                1.73433256021146924798e-6))))))));

            return y;
        }

        private static double GammaNear9(double x) {
            double w = x - 9.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                4.03200000000000000000e4, w, Fma(
                3.53526615657470659354e4, w, Fma(
                1.53954425189805330906e4, w, Fma(
                4.29025124033213891181e3, w, Fma(
                8.34225064618917970411e2, w, Fma(
                1.15918665718504920725e2, w, Fma(
                1.12758958507126338519e1, w, Fma(
                7.05695818214480450306e-1, w,
                2.19361214870971739667e-2))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.26383935578926411520e0, w, Fma(
                7.37329382237256295566e-1, w, Fma(
                -2.60356717879036380312e-1, w, Fma(
                6.11315934187382868564e-2, w, Fma(
                -9.82445338133739523786e-3, w, Fma(
                1.06184673074230406735e-3, w, Fma(
                -7.10890072584460888255e-5, w,
                2.27742019784540119920e-6))))))));

            return y;
        }

        private static double GammaNear10(double x) {
            double w = x - 10.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.62880000000000000000e5, w, Fma(
                3.46991264844610396185e5, w, Fma(
                1.61879721320849231440e5, w, Fma(
                4.78696450948412799425e4, w, Fma(
                9.80899998726203844692e3, w, Fma(
                1.42902532503554222795e3, w, Fma(
                1.45131083745670389052e2, w, Fma(
                9.45087434570684235823e0, w,
                3.04625698988808970695e-1))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.29553768374096494532e0, w, Fma(
                7.75549420204888687827e-1, w, Fma(
                -2.81305394686737534637e-1, w, Fma(
                6.79288088722608934399e-2, w, Fma(
                -1.12420865935930794315e-2, w, Fma(
                1.25309678193093716918e-3, w, Fma(
                -8.66598140482971532878e-5, w,
                2.87308254669950868044e-6))))))));

            return y;
        }

        private static double GammaNear11(double x) {
            double w = x - 11.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.62880000000000000000e6, w, Fma(
                3.72205135538709301571e6, w, Fma(
                1.84209294255439190227e6, w, Fma(
                5.73822262125472497405e5, w, Fma(
                1.23217470018720114792e5, w, Fma(
                1.87351435011684832049e4, w, Fma(
                1.97932142633655716669e3, w, Fma(
                1.33714798046177452580e2, w,
                4.46072229567500135928e0))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.32605501538200632157e0, w, Fma(
                8.13231460319916988031e-1, w, Fma(
                -3.02475989822813967417e-1, w, Fma(
                7.49769818488199399928e-2, w, Fma(
                -1.27520115799734191935e-2, w, Fma(
                1.46258683110372092539e-3, w, Fma(
                -1.04223794136230785698e-4, w,
                3.56603994281946089009e-6))))))));

            return y;
        }

        private static double GammaNear12(double x) {
            double w = x - 12.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.99168000000000000000e7, w, Fma(
                4.33978610005730310459e7, w, Fma(
                2.26003461127440430190e7, w, Fma(
                7.36966606479502762528e6, w, Fma(
                1.64995434047099322745e6, w, Fma(
                2.60728454749954706115e5, w, Fma(
                2.85509856912645840126e4, w, Fma(
                1.99469778654288073699e3, w,
                6.86832415643108482973e1))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.35545376249312224587e0, w, Fma(
                8.50352307138933137746e-1, w, Fma(
                -3.23819703692710450503e-1, w, Fma(
                8.22555955746614568012e-2, w, Fma(
                -1.43507520225035583541e-2, w, Fma(
                1.69024860921016484924e-3, w, Fma(
                -1.23837140824545917593e-4, w,
                4.36217480676931854632e-6))))))));

            return y;
        }

        private static double GammaNear13(double x) {
            double w = x - 13.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                4.79001600000000000000e8, w, Fma(
                5.47142773272545831293e8, w, Fma(
                2.97848708905929606408e8, w, Fma(
                1.01134539634934094332e8, w, Fma(
                2.35044176245543991930e7, w, Fma(
                3.84573156629393096971e6, w, Fma(
                4.35093685527969460848e5, w, Fma(
                3.13472250844357578824e4, w,
                1.11128488676519449865e3))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.38373834178123013790e0, w, Fma(
                8.86823544716992188798e-1, w, Fma(
                -3.45249820677872607179e-1, w, Fma(
                8.97297083462729021668e-2, w, Fma(
                -1.60311310977587409995e-2, w, Fma(
                1.93539796675520058374e-3, w, Fma(
                -1.45495829881901616053e-4, w,
                5.26473504470165023218e-6))))))));

            return y;
        }

        private static double GammaNear14(double x) {
            double w = x - 14.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                6.22702080000000000000e9, w, Fma(
                7.42254472095116624489e9, w, Fma(
                4.20113542632415537598e9, w, Fma(
                1.47882409512659789144e9, w, Fma(
                3.55433810361630327684e8, w, Fma(
                6.00188382397389676556e7, w, Fma(
                6.99553727094738378415e6, w, Fma(
                5.18435113015743875678e5, w,
                1.88791862972795056438e4))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.41092837968698589866e0, w, Fma(
                9.22581756705847091928e-1, w, Fma(
                -3.66692968765934428079e-1, w, Fma(
                9.73672908453891038212e-2, w, Fma(
                -1.77860288070750548355e-2, w, Fma(
                2.19725536441069309769e-3, w, Fma(
                -1.69178226274081660021e-4, w,
                6.27592190605702592758e-6))))))));

            return y;
        }

        private static double GammaNear15(double x) {
            double w = x - 15.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                8.71782912000000000000e10, w, Fma(
                1.07864846386103848040e11, w, Fma(
                6.31985167844824396530e10, w, Fma(
                2.29765730883816774869e10, w, Fma(
                5.69271536192035138372e9, w, Fma(
                9.89275560019858223490e8, w, Fma(
                1.18491270322853103445e8, w, Fma(
                9.01227318469855575539e6, w,
                3.36431439704738902021e5))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.43705644982759998066e0, w, Fma(
                9.57587080111246153509e-1, w, Fma(
                -3.88089789908359059826e-1, w, Fma(
                1.05140086838055122124e-1, w, Fma(
                -1.96087090804403828732e-2, w, Fma(
                2.47500798352576537484e-3, w, Fma(
                -1.94851071804162615274e-4, w,
                7.39711730735237276457e-6))))))));

            return y;
        }

        private static double GammaNear16(double x) {
            double w = x - 16.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.30767436800000000000e12, w, Fma(
                1.67231875709778968332e12, w, Fma(
                1.01061248362274323401e12, w, Fma(
                3.78290861418336239265e11, w, Fma(
                9.63495805472670573322e10, w, Fma(
                1.71888127808712045536e10, w, Fma(
                2.11099809104946039689e9, w, Fma(
                1.64451566387072745128e8, w,
                6.28170505289827573633e6))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.46216379359543931273e0, w, Fma(
                9.91818434875349214687e-1, w, Fma(
                -4.09393234407971699003e-1, w, Fma(
                1.13023535378532450824e-1, w, Fma(
                -2.14929260490782136285e-2, w, Fma(
                2.76784096028999772451e-3, w, Fma(
                -2.22473127485349958462e-4, w,
                8.62903947955886892704e-6))))))));

            return y;
        }

        private static double GammaNear17(double x) {
            double w = x - 17.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.09227898880000000000e13, w, Fma(
                2.75598444459628886135e13, w, Fma(
                1.71264219777905583788e13, w, Fma(
                6.58279446066095812922e12, w, Fma(
                1.71944413814137324906e12, w, Fma(
                3.14230937145680959084e11, w, Fma(
                3.94926501929973244624e10, w, Fma(
                3.14553196484873534932e9, w,
                1.22743665706652073552e8))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.48629681019144033259e0, w, Fma(
                1.02526911574332657769e0, w, Fma(
                -4.30566646938896568595e-1, w, Fma(
                1.20996485720091402993e-1, w, Fma(
                -2.34329621230194012219e-2, w, Fma(
                3.07495747007261550370e-3, w, Fma(
                -2.51997935659983490885e-4, w,
                9.97187552177053227677e-6))))))));

            return y;
        }

        private static double GammaNear18(double x) {
            double w = x - 18.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                3.55687428096000000000e14, w, Fma(
                4.81185497574559452747e14, w, Fma(
                3.06703152246804482204e14, w, Fma(
                1.20773486823625845475e14, w, Fma(
                3.22854071157410410633e13, w, Fma(
                6.03275503748335467743e12, w, Fma(
                7.74570397120133215195e11, w, Fma(
                6.29764407555497150605e10, w,
                2.50675480187999861318e9))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.50950439989737091463e0, w, Fma(
                1.05794307657907074785e0, w, Fma(
                -4.51581942634519924038e-1, w, Fma(
                1.29040829325865791816e-1, w, Fma(
                -2.54236258265314943922e-2, w, Fma(
                3.39559103291099647485e-3, w, Fma(
                -2.83375881654932344952e-4, w,
                1.14253918380288900208e-5))))))));

            return y;
        }

        private static double GammaNear19(double x) {
            double w = x - 19.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                6.40237370572800000000e15, w, Fma(
                8.87405081494282699435e15, w, Fma(
                5.78894736980493590404e15, w, Fma(
                2.33081205634101717618e15, w, Fma(
                6.36528970610807573194e14, w, Fma(
                1.21411426525825897481e14, w, Fma(
                1.59008905998445498055e13, w, Fma(
                1.31784778437351976002e12, w,
                5.34388998416577323146e10))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.53183605001134777656e0, w, Fma(
                1.08985193893144990971e0, w, Fma(
                -4.72417987401404560903e-1, w, Fma(
                1.37141113398803384456e-1, w, Fma(
                -2.74602263963405050736e-2, w, Fma(
                3.72901236037712584348e-3, w, Fma(
                -3.16555698115322052699e-4, w,
                1.29890241561688366532e-5))))))));

            return y;
        }

        private static double GammaNear20(double x) {
            double w = x - 20.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.21645100408832000000e17, w, Fma(
                1.71391363767209224593e17, w, Fma(
                1.14535301640829165100e17, w, Fma(
                4.77958762157648550217e16, w, Fma(
                1.37751014671046262026e16, w, Fma(
                2.85528598814762755265e15, w, Fma(
                4.27716481469730599436e14, w, Fma(
                4.49041630863876491247e13, w, Fma(
                3.01834179829296074561e12, w,
                9.97727666403493192747e10)))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.56157810629041683850e0, w, Fma(
                1.14261636243637449729e0, w, Fma(
                -5.15902100282936456456e-1, w, Fma(
                1.59005342761034537970e-1, w, Fma(
                -3.48425365243671694301e-2, w, Fma(
                5.45540868184180760676e-3, w, Fma(
                -5.91955009525368716219e-4, w, Fma(
                4.06687281488037794482e-5, w,
                -1.35888236041166244558e-6)))))))));

            return y;
        }

        private static double GammaNear21(double x) {
            double w = x - 21.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.43290200817664000000e18, w, Fma(
                3.50023414704946580088e18, w, Fma(
                2.38664046336652409344e18, w, Fma(
                1.01547179430516378433e18, w, Fma(
                2.98206752676021023080e17, w, Fma(
                6.29440715713435768276e16, w, Fma(
                9.59622269719188352746e15, w, Fma(
                1.02481068014614375587e15, w, Fma(
                7.00368979574632420465e13, w,
                2.35274172237795189833e12)))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.58181658220850476775e0, w, Fma(
                1.17273196827297128893e0, w, Fma(
                -5.36648909855536883885e-1, w, Fma(
                1.67680911221815468549e-1, w, Fma(
                -3.72615552048946433362e-2, w, Fma(
                5.91826850919210717096e-3, w, Fma(
                -6.51652579729245318498e-4, w, Fma(
                4.54464063317504256019e-5, w,
                -1.54202747270471068025e-6)))))))));

            return y;
        }

        private static double GammaNear22(double x) {
            double w = x - 22.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                5.10909421717094400000e19, w, Fma(
                7.49389795711444594561e19, w, Fma(
                5.20598693731364973251e19, w, Fma(
                2.25541100816968369277e19, w, Fma(
                6.74020349155570789729e18, w, Fma(
                1.44704258414048752923e18, w, Fma(
                2.24277830389678882807e17, w, Fma(
                2.43383732444871670153e16, w, Fma(
                1.68947330085700712294e15, w,
                5.76236782771151484304e13)))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.60136680936450443215e0, w, Fma(
                1.20218470428080315043e0, w, Fma(
                -5.57195960861611526759e-1, w, Fma(
                1.76383910234217846351e-1, w, Fma(
                -3.97202046211630689959e-2, w, Fma(
                6.39504199551052837291e-3, w, Fma(
                -7.13990019860395389231e-4, w, Fma(
                5.05054620699229824520e-5, w,
                -1.73874600996591505612e-6)))))))));

            return y;
        }

        private static double GammaNear23(double x) {
            double w = x - 23.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                1.12400072777760768000e21, w, Fma(
                1.67850618239934726550e21, w, Fma(
                1.18649726864952956536e21, w, Fma(
                5.22770842428239346934e20, w, Fma(
                1.58807006719055360600e20, w, Fma(
                3.46411086171322626617e19, w, Fma(
                5.45287904041321551431e18, w, Fma(
                6.00739948165864201208e17, w, Fma(
                4.23192785802745685063e16, w,
                1.46428478023566904376e15)))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.62026564974118365392e0, w, Fma(
                1.23099355488261706902e0, w, Fma(
                -5.77536448613077450421e-1, w, Fma(
                1.85105405953306928266e-1, w, Fma(
                -4.22149810851205969192e-2, w, Fma(
                6.88500970236322343506e-3, w, Fma(
                -7.78888843271381080050e-4, w, Fma(
                5.58425846179929139126e-5, w,
                -1.94910243112866506150e-6)))))))));

            return y;
        }

        private static double GammaNear24(double x) {
            double w = x - 24.0;

#if DEBUG
            if (w < -0.5 || w > 0.5) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double y = Fma(
                2.58520167388849766400e22, w, Fma(
                3.92570013023396516407e22, w, Fma(
                2.82044198209349808293e22, w, Fma(
                1.26247219958341765604e22, w, Fma(
                3.89454192974468055238e21, w, Fma(
                8.62347221104336773597e20, w, Fma(
                1.37739159608066574271e20, w, Fma(
                1.53923862151612663415e19, w, Fma(
                1.09951365087528528667e18, w,
                3.85650189872096289763e16)))))))))
                / Fma(
                1.00000000000000000000e0, w, Fma(
                -1.63854823192410578411e0, w, Fma(
                1.25917866366258507458e0, w, Fma(
                -5.97665941015605634035e-1, w, Fma(
                1.93837729459897352378e-1, w, Fma(
                -4.47427142721710133417e-2, w, Fma(
                7.38749696405145903748e-3, w, Fma(
                -8.46272387088990782284e-4, w, Fma(
                6.14542205244802881877e-5, w,
                -2.17313646012790655241e-6)))))))));

            return y;
        }

        private static double GammaLimit(double x) {
            double s = 1.0 / x, t = s * s;

            double v = Fma(
                8.3333333333333333333e-2, t, Fma(
                -2.7777777777777777778e-3, t, Fma(
                7.9365079365079365079e-4, t, Fma(
                -5.952380952380952381e-4, t, Fma(
                8.4175084175084175084e-4, t, 
                -1.9175269175269175269e-3)))));

            double a = Math.Sqrt(2.0 * Math.PI / x);
            double b = Math.Exp(v * s);
            double c = (x < 128.0)
                ? (Math.Pow(x, x) / Math.Exp(x))
                : Math.Pow(Math.Pow(x, 128.0) / Math.Exp(128.0), Math.ScaleB(x, -7));

            double y = a * b * c;

            return y;
        }

        private static double GammaExtremeLarge(double x) {
#if DEBUG
            if (x <= 170.625) {
                throw new ArgumentOutOfRangeException(nameof(x));
            }
#endif

            double u = x - 1.0; 

            double y = u * GammaLimit(u);

            return y;
        }
    }
}