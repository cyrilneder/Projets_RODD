using PlotlyJS


# N varie de 8 à 50, G = 5/8 * N, A = 2, T = 50, theta1 = 0.001
data = [0.041114091873168945, 0.03421499729156494, 0.020683431625366212, 0.01063852310180664, 0.011694622039794923, 0.009366822242736817, 0.011053323745727539, 0.011786365509033203, 0.013962984085083008, 0.018955230712890625, 0.009188151359558106, 0.014823102951049804, 0.01615908145904541, 0.014908742904663087, 0.030535340309143066, 0.014059925079345703, 0.018848943710327148, 0.013434624671936036, 0.019196844100952147, 0.027320003509521483, 0.03193578720092773, 0.020642733573913573, 0.032025480270385744, 0.03303048610687256, 0.02104794979095459, 0.021877002716064454, 0.036464381217956546, 0.024904847145080566, 0.030553698539733887, 0.045178675651550294, 0.03212959766387939, 0.04408564567565918, 0.05748105049133301, 0.02812614440917969, 0.029541993141174318, 0.10859699249267578, 0.046123409271240236, 0.03764169216156006, 0.04674539566040039, 0.03791742324829102, 0.0454211950302124, 0.045184659957885745, 0.05377511978149414, 0.04576032161712647, 0.048656272888183597, 0.06500833034515381, 0.047956228256225586, 0.05276763439178467, 0.05524475574493408, 0.05553083419799805, 0.057883405685424806, 0.0602299690246582, 0.07495276927947998, 0.06877779960632324, 0.06316671371459961, 0.0656196117401123, 0.06916484832763672, 0.08131146430969238, 0.11749186515808105, 0.07607941627502442, 0.1317366600036621, 0.1005007266998291, 0.10896089076995849, 0.09355173110961915, 0.13814697265625, 0.09609789848327636, 0.10790016651153564, 0.12509431838989257, 0.13317627906799318, 0.12347025871276855, 0.09248402118682861, 0.10060315132141114, 0.12785160541534424, 0.1036022424697876, 0.12677581310272218, 0.10859494209289551, 0.11395003795623779, 0.20540060997009277, 0.12685091495513917, 0.1300602674484253, 0.13051972389221192, 0.15479559898376466, 0.20969147682189943, 0.1667940378189087, 0.18419127464294432, 0.2099278211593628, 0.1585092544555664, 0.1744229793548584, 0.20967655181884765, 0.21609852313995362, 0.1527625560760498, 0.15508558750152587, 0.16013460159301757]

[0.2112689733505249, 0.9184395551681519,2.182102179527283,4.203507351875305,5.63663375377655,10.03566746711731]

plot(scatter(x=8:50, y=data, mode="markers"))