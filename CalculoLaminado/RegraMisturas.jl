#
# Cálculos das regras das misturas -> Separados para caso você tenha a fração volumétrica ou a fração mássica.
#

# Regra das misturas -> Fração volumétrica
function FracaoVolumetrica(Ef, Em, vf, vm, fvf, fvm)

    # Vamos calcular o módulo de elasticidade longitudinal
    E1 = (Ef * fvf) + (Em * fvm)

    # Módulo de elasticidade transversal
    E2 = 1/ ((fvf/Ef) + (fvm/Em))

    # Coeficiente de poisson v12
    v12 = (vf * fvf) + (vm * fvm)

    # Módulos de cisalhamento
    # Módulo de cisalhamento da fibra
    Gf = Ef / (2*(1+vf))
    #@show Gf
    # Módulo de cisalhamento da matriz
    Gm = Em / (2*(1+vm))
    #@show Gm
    # Módulo de cisalhamento G12
    G12 = 1 / ((fvf/Gf) + (fvm/Gm))

    # Retornando os valores de interesse
    return E1, E2, v12, G12

end

# Regra das misturas -> Fração mássica
function FracaoMassica(fmf, fmm, df, dm, vf, vm, Ef, Em)

    # Primeiramente, vamos encontrar a densidade do compósito
    pc = 1 / ((fmf/df) + (fmm/dm))

    # Calculando a fração volumétrica da fibra
    fvf = (fmf * pc) / df

    # Calculando a fração volumétrica da matriz
    fvm = (fmm * pc) / dm

    # Vamos calcular o módulo de elasticidade longitudinal
    E1 = (Ef * fvf) + (Em * fvm)

    # Módulo de elasticidade transversal
    E2 = 1/ ((fvf/Ef) + (fvm/Em))
    
    # Coeficiente de poisson v12
    v12 = (vf * fvf) + (vm * fvm)
    
    # Módulos de cisalhamento
    # Módulo de cisalhamento da fibra
    Gf = Ef / (2*(1+vf))

    # Módulo de cisalhamento da matriz
    Gm = Em / (2*(1+vm))
        
    # Módulo de cisalhamento G12
    G12 = 1 / ((fvf/Gf) + (fvm/Gm))
    
    # Retornando os valores de interesse
    return E1, E2, v12, G12

#end