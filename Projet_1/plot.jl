using Plots


data = [0.6439661979675293 0.7172893762588501 0.38086321353912356 0.4747746706008911 0.47377514839172363 0.849184226989746; 0.5288702011108398 0.4986548900604248 0.6071965217590332 0.7613292932510376 0.6650416851043701 0.8369988679885865; 0.3993539571762085 0.4698573112487793 0.5295150756835938 0.5699986457824707 0.6023004055023193 0.8266937255859375; 0.4993598937988281 0.4335439920425415 0.713708472251892 0.7672348499298096 0.868030595779419 0.7675543069839478; 0.5826488971710205 0.6243920087814331 0.8643115043640137 0.8304750204086304 1.1744812965393066 1.6382641553878785; 0.7048612833023071 0.6152467489242553 0.6293207883834839 0.9148398160934448 0.9966335296630859 1.011458396911621]

# n m varient, p = 6, q = 3
data = [0.3120940804481506 0.32439026832580564 0.2997978925704956 0.3525524616241455 0.37066586017608644 0.4009990692138672 0.5218647718429565 0.5381778240203857 0.5950982332229614 0.6291195392608643 0.44167802333831785; 0.31489388942718505 0.36685919761657715 0.4699180841445923 0.4657161712646484 0.46457111835479736 0.4342073917388916 0.46590800285339357 0.6303716897964478 0.799160099029541 0.8771427631378174 0.7636188745498658; 0.29558260440826417 0.4655693769454956 0.4116262197494507 0.4976963520050049 0.35755798816680906 0.5295346975326538 0.5439676761627197 0.6758466482162475 0.5862280368804932 0.925901222229004 0.7698856830596924; 0.3564456939697266 0.5053547859191895 0.4999713897705078 0.4739309072494507 0.5959487676620483 0.46453559398651123 0.7577349185943604 0.9161274194717407 0.7921461820602417 0.8371710538864136 0.9983999967575073; 0.39372959136962893 0.47513551712036134 0.5042593717575073 0.5048408269882202 0.502909779548645 0.7469016075134277 0.867728877067566 0.8566803932189941 0.7807393074035645 0.7622183561325073 0.8488499879837036; 0.3503655195236206 0.4506392478942871 0.6438525199890137 0.761597752571106 0.6150080442428589 0.8269150018692016 0.9719278573989868 0.6889304876327514 1.0078853845596314 0.953293251991272 0.9623032093048096; 0.4768620729446411 0.4945739030838013 0.5028961420059204 0.6718421220779419 0.7667878627777099 0.7972378253936767 0.9872243165969848 0.9404699563980102 1.0085751056671142 1.2282402992248536 1.3140836954116821; 0.4961317777633667 0.4665502071380615 0.6932224988937378 0.7115610361099243 0.7986215114593506 0.8474962472915649 1.0657920837402344 1.4019457817077636 1.2786837100982666 1.4059460639953614 1.6100718259811402; 0.6647876739501953 0.4711137056350708 0.5719128370285034 0.7116880893707276 0.960002326965332 1.718245768547058 1.33927161693573 0.9755719184875489 1.2265551328659057 1.2055525541305543 1.200948214530945; 0.6470071077346802 0.7042775869369506 0.6679670095443726 0.968664002418518 1.0827741861343383 1.3789743661880494 1.290052890777588 1.568287181854248 1.3253049850463867 1.663319706916809 1.4050142765045166; 0.5242993593215942 0.8502145051956177 0.6056270122528076 0.9125085592269897 0.9432740926742553 1.1758940935134887 2.141046929359436 1.0887176036834716 1.269672703742981 1.3262391328811645 1.5961167097091675]



# n=m varient de 10 à 25, p = int(0.6*n), q = p/2
n_data = [0.406596302986145, 0.3954850912094116, 0.532216477394104, 0.8308518886566162, 1.7647210597991942, 1.4780856370925903, 2.1815401554107665, 5.537090802192688, 4.834748053550721, 3.905217504501343, 8.010083198547363, 12.713053369522095, 31.756101107597352, 18.42375361919403, 39.85250599384308, 52.690928769111636]




#title!("Log-log plot")
#xlabel!("x")
#ylabel!("y")

# p varie de 3 à 20, n = 10, q = p/2
p_data = [0.16821639537811278, 0.2602193593978882, 0.2941837549209595, 0.27802343368530275, 0.47490177154541013, 0.4379671812057495, 0.4452111005783081, 0.5077829837799073, 0.586903738975525, 0.8800825834274292, 0.8301659345626831, 0.7863914728164673, 0.9015446901321411, 0.9221567153930664, 0.604929780960083, 0.7257765769958496, 0.8068867206573487, 0.7362793684005737]

#plot(heatmap(z=data, x = 10:20, y = 10:20))
#plot(10:25, n_data, ls=:dot)


x = 3:20
y = p_data

p = plot(x, y)
xlims!(p, 3, 20)
ylims!(p, 0.1, 1)
plot!(p, xscale=:log10, yscale=:log10, minorgrid=true)
display(p)