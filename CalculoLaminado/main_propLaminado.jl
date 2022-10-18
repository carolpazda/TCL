#
# PROGRAMA DE RELAÇÕES CONSTITUTIVAS
# O presente código é de autoria exclusiva do grupo: 
# - DIANA MARIA VICENTIN RIBEIRO
# - LEONARDO DA SILVA
# - VERÔNICA CAROLINE HERBST PAZDA
#
# Fica estritamente proibida a utilização desse código para outros fins,
# não sendo o didático.
# Cálculo realizado por meio das Propriedades do Laminado -> Já temos E11, E22, G12, v12.
#
# Esse 'main' deve ser utilizado apenas se o usuário tiver as propriedades em!!!

# Códigos auxiliares
include("calculo_laminado.jl")
include("plot.jl")
include("CriteriosFalha.jl")

# Bibliotecas que vamos utilizar
using LinearAlgebra, Plots

function main()

    # Primeiro vamos informar ao programa o número de camadas do laminado
    n_camadas = 4 

    # Digite aqui o ângulo de orientação das fibras [°] (em ordem, por favor!)
    angulo = [45.0; 0.0; 0.0; 45.0]

    # Digite aqui a espessura do laminado [m]
    espessura = 3E-3

    # Espessura total do laminado
    v = espessura * n_camadas

    # E quais são as propriedades do material (Em Pa meu amigo!)?
    E11 = 19.76E9
    E22 = 1.97E9
    G12 = 0.7E9
    v12 = 0.35

    # Carregamento do Laminado? [N e N/m]
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

    # Vamos retornar tudo que é de interesse para o usuário
    return ABBD, Defo_PlanoMedio, tensao_global, tensao_local, deformacao_global, deformacao_local, SF_MaxTensao, SF_Hill, Msh, SF_Wu, Msw

end
