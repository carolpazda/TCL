#
# CÁLCULO REALIZADO POR MEIO DA REGRA DAS MISTURAS
# O presente código é de autoria exclusiva do grupo: 
# - DIANA MARIA VICENTIN RIBEIRO
# - LEONARDO DA SILVA
# - VERÔNICA CAROLINE HERBST PAZDA
#
# Para esse caso não temos todas as propriedades -> Realizamos os ensaios
# Então vamos encontrar o E11, E22, G12, v12 por meio da regra das misturas!
#
# Esse 'main' deve ser utilizado apenas se o usuário quiser utilizar a regra das misturas!!!

# Códigos auxiliares
include("calculo_laminado.jl")
include("plot.jl")
include("CriteriosFalha.jl")
include("RegraMisturas.jl")

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

    # Primeiro vamos informar ao programa o número de camadas do laminado
    n_camadas = 4 

    # Digite aqui o ângulo de orientação das fibras [°] (em ordem, por favor!)
    angulo = [45.0; 0.0; 0.0; 45.0]

    # Digite aqui a espessura do laminado [m]
    espessura = 3E-3

    # Espessura total do laminado
    v = espessura * n_camadas

    # E quais são as propriedades do material? Vamos puxar do outro código!
	# RegraMisturas -> Fração volumétrica
    E11, E22, v12, G12 = FracaoVolumetrica(Ef, Em, vf, vm, fracaoVolumetrica_fibra, fracaoVolumetrica_matriz)
    # @show E11, E22, v12, G12 -> Verificado tudo certo!

    # RegraMisturas -> Fração mássica
    #E11, E22, v12, G12 = FracaoMassica(fracaoMassica_fibra, fracaoMassica_matriz, df, dm, vf, vm, Ef, Em)

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
    MaxTensao = MaximaTensao(n_camadas, tensao_local, Xt, Xc, Yt, Yc, S12)

    # 2 -> Tsai-Hill
    TsaiHill = Tsai_Hill(n_camadas, tensao_local, Xt, Xc, Yt, Yc, S12)

    # 3 -> Tsai-Wu
    TsaiWu = Wu(Xt, Xc, Yt, Yc, S12, n_camadas, tensao_local)

    # Vamos retornar tudo que é de interesse para o usuário
    return ABBD, Defo_PlanoMedio, tensao_global, tensao_local, deformacao_global, deformacao_local, MaxTensao, TsaiHill, TsaiWu

end
