#######################
## Network inference ##
#######################

using CSV, DataFrames, PhyloNetworks, RCall, PhyloPlots
cd("/Users/kevinsanchez/proj/Hdar_postdoc1/analisis/phylonetworks_K8_Hwil")
cftable = readTableCF("CF.csv")
spptree = readTopology("((Pa1,Pa2),(WRN,((SEL,WCH),(SNQ,(CRN_NEC,(wil,MZ_NNQ))))));") # Snapper tree

net0 = snaq!(spptree, cftable, hmax = 0, filename = "net0", runs = 30, seed = 233342); rootonedge!(net0, 15)
net0 = readSnaqNetwork("net0.out"); rootonedge!(net0, 15)
net1 = snaq!(net0, cftable, hmax = 1, filename = "net1", runs = 30, seed = 2323)
net1 = readSnaqNetwork("net1.out"); rootonedge!(net1, 18)
# networks 1-8 ~loglik
net2 = snaq!(net1, cftable, hmax = 2, filename = "net2", runs = 30, seed = 56764) # all h1 networks, let's try with hmax = 3
rootonedge!(net1, 7)
net2 = readTopology("(Pa1,(WRN,(((((wil,(CRN_NEC,#H10:0.0292637052888941::0.4769371432898904):0.2321688665904056):0.05615355550762497,MZ_NNQ):0.07569174128530191,SNQ):0.1280272959132389,((WCH,SEL):0.06159291975958376)#H10:0.22867408015552557::0.5230628567101097):0.0)#H11:0.0::0.9761675334937375):1.5119302107622745,(Pa2,#H11:0.11090396943249527::0.02383246650626249):0.004025397193045896);") # best network with 2 reticulations, inferred in iteration 23, -loglik = -17920.131861206788
rootonedge!(net2, 19)
net3 = snaq!(net2, cftable, hmax = 3, filename = "net3", runs = 30, seed = 123)
net3 = readTopology("(((WCH,(#H18:0.039582827000212956::0.07318655925097312,SEL):0.7003049452383624):0.1321574193394584,(((MZ_NNQ)#H16:::0.7311037714738124,(CRN_NEC,(wil,#H16:::0.26889622852618755):0.3866893578583273):0.07530926209044952):0.009033116712668684)#H18:0.0::0.9268134407490269):0.026231492110830104,(SNQ,(WRN,#H10:0.5055646005636519::0.06775393028856581):6.183277854364634):0.0,((Pa1,Pa2):0.873141429545499)#H10:0.7331963212449063::0.9322460697114342);")  # best network with 3 reticulations, inferred in iteration 29 (compatible with root), -loglik = 20848.426709607684
rootonedge!(net3, 24)
plot(net0, style=:majortree, showedgenumber = true, arrowlen = .2, showgamma = true)
plot(net1, style=:majortree, showedgenumber = true, arrowlen = .2, showgamma = true)
plot(net2, style=:majortree, showedgenumber = true, arrowlen = .2, showgamma = true)
plot(net3, style=:majortree, showedgenumber = true, arrowlen = .2, showgamma = true)

# plot changes in pseudo-lik values across iterations
y_h1 = [17639.53161563339, 17639.544114599506, 17639.79891201193, 17640.04557236954, 17640.06278388722, 17640.362758049818, 17642.040105323525, 17644.03178682189, 17783.68315207952, 17783.688346498955, 17783.73619519972, 17783.802009055322, 17784.20179536901, 18052.56625862218, 18053.816466504042, 18157.0256771539, 18173.384978749527, 18561.85776251252, 18561.857762541622, 18561.85776283924, 18561.857763000862, 18561.857764756976, 18561.857769963834, 18720.626280968125, 18796.26721626872, 20013.953526011275, 20013.953550308484, 20434.86138089757, 20447.282739746086, 20484.603613181836]
y_h2 = [17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.532450387273, 17639.532949949356, 17639.538113028448, 17639.540721140158, 17639.58007592538, 17639.58819449439, 17639.591746644925, 17639.624410020802, 17642.644853839607, 17920.131861206788, 17933.434259893493, 18561.939029284364, 18582.685399817812, 18583.33293401178, 19470.902172641185, 19481.651491884724, 20276.639762201637]
y_h3 = [17785.329875984156, 17918.740733783095, 17919.743651391313, 17919.743651391313, 17919.885472113652, 17920.64053219712, 17926.193174748216, 17936.769489124188, 17939.980993914716, 17941.058924349254, 17941.05892447903, 17941.17976658675, 17941.179766655823, 17941.179766655823, 17943.884960778676, 17946.30714486394, 17952.94681527914, 17973.0241128015, 17979.454053855992, 18179.84581341461, 18588.869626465428, 18720.801428965635, 18796.082797590712, 18869.7833249198, 19158.459299111317, 19412.442858020953, 20061.790990985974, 20331.237656136433, 20848.426709607684, 21349.295541566556]
x = 1:30
@rput x y_h1 y_h2 y_h3

R"""
ymin <- min(c(y_h1, y_h2, y_h3))
plot(x, y_h1, type = "b", col = "blue", ylab = "network score", xlab = "it", ylim = c(ymin, ymax))
lines(x, y_h2, type = "b", col = "red")
lines(x, y_h3, type = "b", col = "green")
legend("bottomright", legend = c("h1", "h2", "h3"), col = c("blue", "red", "green"), lty = 1, pch = 1)
"""


# pseudo-loglik scores
net2_score = 17920.131861206788
net3_score = 20848.426709607684
scores = [net0.loglik, net1.loglik, net2_score, net3_score]
hmax = collect(0:3)
R"plot"(hmax, scores, type = "b", ylab = "network score", xlab = "hmax", col = "blue")



# Boootstrap

# cftable = CSV.read("CF.csv", DataFrame)
# bootnet0 = bootsnaq(net0, cftable, nrep = 100, hmax = 0, runs = 20, filename = "boot_h0")
# bootnet1 = bootsnaq(net1, cftable, nrep = 100, hmax = 1, runs = 20, filename = "boot_h1")

# net = readSnaqNetwork("net1.out")
# bootnet = readMultiTopology("boot_h1.out")
# rootatnode!(net, "Pa")
# BSe_tree, tree1 = treeEdgesBootstrap(bootnet, net)
# # support for hybrid branches
# BSn, BSe, BSc, BSgam, BSedgenum = hybridBootstrapSupport(bootnet, net)
# # plot backbone tree
# plot(net, :R, edgeLabel = BSe_tree[BSe_tree[!, :proportion] .< 100.0, :])
# # plot major and minor hybrid edges 
# plot(net, :R, edgeLabel = BSe[!, [:edge,:BS_hybrid_edge]]);


#######################


using CSV, DataFrames, SNaQ, PhyloNetworks, RCall, PhyloPlots
cd("/Users/kevinsanchez/proj/Hdar_postdoc1/analisis/phylonetworks_MO")
cftable = readTableCF("CF.csv")
spptree = readTopology("((Pa1,Pa2),(WRN,((SEL,WCH),(SNQ,(CRN_NEC,(Hwil,MZ_NNQ))))));") # Snapper tree

net0 = snaq!(spptree, cftable, hmax = 0, filename = "net0", runs = 30, seed = 233342); rootonedge!(net0, 15)
net0 = readTopology("((WCH,SEL):0.06135098505308582,((((Pa1,Pa2):1.3668286516642203,WRN):0.08359977235362863,SNQ):0.038386028639381016,(wil,MZ_NNQ):0.08517916811229667):0.06951089982772049,CRN_SEC);")
rootonedge!(net0, 6)
net1 = snaq!(net0, cftable, hmax = 1, filename = "net1", runs = 30, seed = 2323)
net1 = readTopology("(SEL,(CRN_SEC)#H10:::0.5033217682848337,(WCH,(((Pa1,Pa2):1.379710585213527,WRN):0.08495834916165623,(SNQ,((wil,#H10:::0.4966782317151663):0.0630262188908492,MZ_NNQ):0.10907616887296148):0.017517228146696966):0.15655189339986558):0.1030906016112008);")
rootonedge!(net1, 6)
net2 = snaq!(net1, cftable, hmax = 2, filename = "net2", runs = 30, seed = 56764)
net2 = readTopology("(Pa2,(((SNQ,((wil,#H10:::0.4210290762512214):0.17199834084206828,MZ_NNQ):0.0894376947205943):0.04779139762982771,(WCH,(SEL,(CRN_SEC)#H10:::0.5789709237487786):0.08525878815880393):0.14304243262781985):0.07722804092269135,(WRN,#H11:::0.008722896182078962):1.399943865606343):1.399943865606343,(Pa1)#H11:::0.991277103817921);")
rootonedge!(net2, 19)
net3 = snaq!(net2, cftable, hmax = 3, filename = "net3", runs = 30, seed = 123) # todas colapsan en net2
plot(net0, style=:majortree, showedgenumber = true, arrowlen = .2, showgamma = true)
plot(net1, style=:majortree, showedgenumber = true, arrowlen = .2, showgamma = true)
plot(net2, style=:majortree, showedgenumber = true, arrowlen = .2, showgamma = true)

# plot changes in pseudo-lik values across iterations
y_h0 = [17639.53161563339, 17639.544114599506, 17639.79891201193, 17640.04557236954, 17640.06278388722, 17640.362758049818, 17642.040105323525, 17644.03178682189, 17783.68315207952, 17783.688346498955, 17783.73619519972, 17783.802009055322, 17784.20179536901, 18052.56625862218, 18053.816466504042, 18157.0256771539, 18173.384978749527, 18561.85776251252, 18561.857762541622, 18561.85776283924, 18561.857763000862, 18561.857764756976, 18561.857769963834, 18720.626280968125, 18796.26721626872, 20013.953526011275, 20013.953550308484, 20434.86138089757, 20447.282739746086, 20484.603613181836]
y_h1 = [17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.531615277665, 17639.532450387273, 17639.532949949356, 17639.538113028448, 17639.540721140158, 17639.58007592538, 17639.58819449439, 17639.591746644925, 17639.624410020802, 17642.644853839607, 17920.131861206788, 17933.434259893493, 18561.939029284364, 18582.685399817812, 18583.33293401178, 19470.902172641185, 19481.651491884724, 20276.639762201637]
y_h2 = [17785.329875984156, 17918.740733783095, 17919.743651391313, 17919.743651391313, 17919.885472113652, 17920.64053219712, 17926.193174748216, 17936.769489124188, 17939.980993914716, 17941.058924349254, 17941.05892447903, 17941.17976658675, 17941.179766655823, 17941.179766655823, 17943.884960778676, 17946.30714486394, 17952.94681527914, 17973.0241128015, 17979.454053855992, 18179.84581341461, 18588.869626465428, 18720.801428965635, 18796.082797590712, 18869.7833249198, 19158.459299111317, 19412.442858020953, 20061.790990985974, 20331.237656136433, 20848.426709607684, 21349.295541566556]
x = 1:30
@rput x y_h0 y_h1 y_h2

R"""
ymin <- min(c(y_h0, y_h1, y_h2))
plot(x, y_h0, type = "b", col = "blue", ylab = "network score", xlab = "it", ylim = c(ymin, ymax))
lines(x, y_h1, type = "b", col = "red")
lines(x, y_h2, type = "b", col = "green")
legend("bottomright", legend = c("h0", "h1", "h2"), col = c("blue", "red", "green"), lty = 1, pch = 1)
"""


# pseudo-loglik scores
net0_score = 8028.860825550925
net1_score = 6203.894948605068
net2_score = 5821.1808766639915
scores = [net0_score, net1_score, net2_score]
hmax = collect(0:2)
R"plot"(hmax, scores, type = "b", ylab = "network score", xlab = "hmax", col = "blue")

