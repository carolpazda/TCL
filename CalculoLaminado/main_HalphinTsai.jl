#
# CÁLCULO POR MEIO DE Halphin e Tsai
# O presente código é de autoria exclusiva do grupo: 
# - DIANA MARIA VICENTIN RIBEIRO
# - LEONARDO DA SILVA
# - VERÔNICA CAROLINE HERBST PAZDA
#
# Esse 'main' deve ser utilizado apenas se o usuário desejar calcular pelo método de Halphin e Tsai!!!

# Códigos auxiliares
include("calculo_laminado.jl")
include("plot.jl")
include("CriteriosFalha.jl")
include("HalphinTsai.jl")

# Bibliotecas que vamos utilizar
using LinearAlgebra, Plots

function main()

    # MÉTODO = 1 -> Temos a fração volumétrica
    # MÉTODO = 2 -> Temos a fração mássica
    # Mude conforme no 'main' conforme a sua necessidade.

    # Digite a fração volumétrica da fibra e da matriz
    fracaoVolumetrica_fibra  = 0.4
    fracaoVolumetrica_matriz = 0.6

    # Digite a fração mássica da fibra e da matriz
    fracaoMassica_fibra = 0.7
    fracaoMassica_matriz = 0.3

    # Digite a densidade da fibra e da matriz (kg/m^3)
    df = 1200 # -> Densidade da fibra
    dm = 2600 # -> Densidade da matriz

    # Digite o módulo de elasticidade da fibra e da matriz, em Pa
    Ef = 85.0E9 # -> Módulo de elasticidade da fibra
    Em = 3.4E9  # -> Módulo de elasticidade da matriz

    # Digite os coeficientes de poisson da fibra e da matriz
    vf = 0.2 # -> Coeficiente de poisson da fibra
    vm = 0.3 # -> Coeficiente de poisson da matriz

    # Digite aqui o seu ξ1 -> Para arranjo quadrado e fibra circular ξ1 = 2
    ξ1 = 2

    # Digite aqui o seu ξ2 -> Para fibra circular e arrando quadrado ξ1 = 1
    ξ2 = 1

    # Primeiro vamos informar ao programa o número de camadas do laminado
    n_camadas = 4 

    # Digite aqui o ângulo de orientação das fibras [°] (em ordem, por favor!)
    angulo = [45.0; 0.0; 0.0; 45.0]

    # Digite aqui a espessura do laminado [m]
    espessura = 3E-3

    # Espessura total do laminado
    v = espessura * n_camadas

    # E quais são as propriedades do material? Vamos puxar do outro código!
    # Halphin-Tsai -> Fração volumétrica
    E11, E22, v12, G12 = FracaoVolumetrica_HalphinTsai(Ef, Em, vf, vm, fracaoVolumetrica_fibra, fracaoVolumetrica_matriz, ξ1, ξ2)
    # @show E11, E22, v12, G12 -> Verificado tudo certo!

    # Halphin-Tsai -> Fração mássica
    #E11, E22, v12, G12 = FracaoMassica_HalphinTsai(fracaoMassica_fibra, fracaoMassica_matriz, df, dm, vf, vm, Ef, Em, ξ1, ξ2)

    # Carregamento do Laminado? [N e Nm]
    Nx = 1000000.0
    Ny = 200000.0
    Nxy = 0.0 
    Mx = 0.0
    My = 0.0
    Mxy = 0.0

    # De forma que o vetor de carregamento vai ser dado por
    N =[Nx;Ny;Nxy;Mx;My;Mxy] 

    # Propriedades de resistência do material -> Todas em Pa!!!
    Xt = 1447.0E6      # -> Resistência à tração 
    Xc = -1447.0E6     # -> Resistência à compressão 
    Yt = 51.7E6        # -> Resistência à tração 
    Yc = -206.0E6      # -> Resistência à compressão 
    S12 = 93.0E6       # -> Resistência ao cisalhamento no plano 1-2 

    # Calculamos a matriz de rigidez do nosso problema em questão
    Q = matriz_Q(E11,E22,G12,v12)

    # Calculamos o h e o zm
    #h,zm =  planoMedio(espessura,n_camadas) -> Não estou precisando utilizar
	
    # Calculamos a matriz ABBD
    ABBD= ABBD_laminado(v,n_camadas,E11, E22, G12, v12, angulo) # -> Verificado

    # Deformação no plano médio do laminado
    Defo_PlanoMedio = ABBD\N # Verificado

    # Calcula as Tensões e Deformações globais # -> Verificado também
    tensao_global, tensao_local, deformacao_global, deformacao_local = Defor_Tensoes(E11, E22, G12, v12, 
	                                                n_camadas, v, angulo, Defo_PlanoMedio[1:3], Defo_PlanoMedio[4:6])

    # CÁLCULO DOS CRITÉRIOS DE FALHA
    # 1 -> Máxima Tensão
    SF_MaxTensao = MaximaTensao(n_camadas, tensao_local, Xt, Xc, Yt, Yc, S12)

    # 2 -> Tsai-Hill
    SF_Hill, Msh = Tsai_Hill(n_camadas, tensao_local, Xt, Xc, Yt, Yc, S12)

    # 3 -> Tsai-Wu
    SF_Wu, Msw = Wu(Xt, Xc, Yt, Yc, S12, n_camadas, tensao_local)

    # Plotando os gráficos
    # Primeiro vamos chamar a função que está nos retornando os pontos de interesse
    sig1, sig2, sig3, sigt, sigc, plot1 = Plot_CriteriosFalha(Xt,Xc,Yt,Yc,S12,0)

    # Plotando a Máxima Tensão
    Plots.plot([Xt,Xc,Xc,Xt, Xt], [Yt,Yt,Yc,Yc,Yt], label="Máxima Tensão", linecolor=:tomato)

    # Plotando Tsai-Hill
    Plots.plot!(sigt, plot1[1:250], label="Tsai-Hill", linecolor=:green)
    Plots.plot!(sigc, plot1[251:500], label = false, linecolor=:green)
    Plots.plot!(sigc, plot1[501:750], label = false, linecolor=:green)
    Plots.plot!(sigt, plot1[751:1000], label = false, linecolor=:green)

    # Vamos plotar os pontos σ1 e σ2 dos critérios de falha utilizados
    tens1 = tensao_local[1,:]
    tens2 = tensao_local[2,:]
    #@show tensao_local
    #@show tens1
    #@show tens2
    Plots.plot!(tens1,tens2, seriestype =:scatter, label = false)

    # Plotando Tsai-Wu
    Plots.plot!(sig1,sig3, label="Tsai-Wu", linecolor=:darkmagenta)
    display(Plots.plot!(sig1, sig2, label=false, linecolor=:darkmagenta, title="Critérios de falha", 
	                    xlabel = "σ1 [Pa]", ylabel = "σ2 [Pa]"))
	
    # Aqui tive que fazer algumas adaptações para plotar a Distribuição de Tensão e Deformação em relação a espessura					
    # Distribuição de Tensão nas Lâminas
    h, z = Auxiliar(espessura, n_camadas, v) # -> Verificado
    deformacao_global_p, tensao_global_p = Deformacao_e_TensaoGlobal(Defo_PlanoMedio, n_camadas, z, E11, E22, G12, v12, angulo) # -> Certo, foi verificado

    # Pontos que vamos utilizar para o Plot
    t1 = tensao_global_p[1,:]
    passo1 = (maximum(t1) - minimum(t1))/1
    t2 = tensao_global_p[2,:]
    passo2 = (maximum(t2) - minimum(t2))/1
    t3 = tensao_global_p[3,:]
    passo3 = (maximum(t3) - minimum(t3))/1

    # Plots -> Distribuição de Tensão nas Lâminas
    p1 = Plots.plot(t1,z, label=false, xticks=minimum(t1):passo1:maximum(t1))
    p2 = Plots.plot(t2,z, label=false, xticks=minimum(t2):passo2:maximum(t2), title="Distribuição de Tensão nas Lâminas")
    p3 = Plots.plot(t3,z, label=false, xticks=minimum(t3):passo3:maximum(t3))
    display(Plots.plot(p1,p2,p3, layout=(1,3)))

    # Distribuição de Deformação nas Lâminas
    d1 = deformacao_global_p[1,:]
    d2 = deformacao_global_p[2,:]
    d3 = deformacao_global_p[3,:]

    # Plots -> Distribuição de Deformação na Lâmina
    r1 = Plots.plot(d1,z, label=false)
    r2 = Plots.plot(d2,z, label=false, title="Distribuição de Deformação nas Lâminas")
    r3 = Plots.plot(d3,z, label=false)
    #display(Plots.plot(r1,r2,r3, layout=(1,3)))

    # Vamos retornar tudo que é de interesse para o usuário
    return ABBD, Defo_PlanoMedio, tensao_global, tensao_local, deformacao_global, deformacao_local, SF_MaxTensao, SF_Hill, Msh, SF_Wu, Msw

end
