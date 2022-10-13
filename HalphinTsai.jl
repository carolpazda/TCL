#
# Método de Halphin e Tsai
#
#
# Cálculos realizados por meio de Halphin-Tsai -> Separados para caso você tenha a fração volumétrica ou a fração mássica.
#

# Halphin e Tsai -> Fração volumétrica
function FracaoVolumetrica_HalphinTsai(Ef, Em, vf, vm, fvf, fvm, ξ1, ξ2)

    # Vamos calcular o módulo de elasticidade longitudinal
    E1 = (Ef * fvf) + (Em * fvm)

    # Vamos calcular n
    n1 = ((Ef/Em) - 1)/ (((Ef/Em) + ξ1))

    # Calculando o módulo de elasticidade transversal
    E2 = (((1 + (ξ1*n1*vf))) / (1 - (n1*vf))) * Em

    # Coeficiente de poisson v12
    v12 = (vf * fvf) + (vm * fvm)

    # Módulos de cisalhamento
    # Módulo de cisalhamento da fibra
    Gf = Ef / (2*(1+vf))

    # Módulo de cisalhamento da matriz
    Gm = Em / (2*(1+vm))

    # Calculando o n
    n2 = ((Gf/Gm) - 1)/ (((Gf/Gm) + ξ2))

    # Calculando o módulo de cisalhamento
    G12 = (((1 + (ξ2*n2*vf))) / (1 - (n2*vf))) * Gm

    # Retornando os valores de interesse
    return E1, E2, v12, G12

end

# Halphin e Tsai -> Fração mássica
function FracaoMassica_HalphinTsai(fmf, fmm, df, dm, vf, vm, Ef, Em, ξ1, ξ2)

    # Primeiramente, vamos encontrar a densidade do compósito
    pc = 1 / ((fmf/df) + (fmm/dm))

    # Calculando a fração volumétrica da fibra
    fvf = (fmf * pc) / df

    # Calculando a fração volumétrica da matriz
    fvm = (fmm * pc) / dm

    # Vamos calcular o módulo de elasticidade longitudinal
    E1 = (Ef * fvf) + (Em * fvm)

    # Vamos calcular n
    n1 = ((Ef/Em) - 1)/ (((Ef/Em) + ξ1))

    # Calculando o módulo de elasticidade transversal
    E2 = (((1 + (ξ1*n1*vf))) / (1 - (n1*vf))) * Em

    # Coeficiente de poisson v12
    v12 = (vf * fvf) + (vm * fvm)

    # Módulos de cisalhamento
    # Módulo de cisalhamento da fibra
    Gf = Ef / (2*(1+vf))

    # Módulo de cisalhamento da matriz
    Gm = Em / (2*(1+vm))

    # Calculando o n
    n2 = ((Gf/Gm) - 1)/ (((Gf/Gm) + ξ2))

    # Calculando o módulo de cisalhamento
    G12 = (((1 + (ξ2*n2*vf))) / (1 - (n2*vf))) * Gm

    # Retornando os valores de interesse
    return E1, E2, v12, G12

end