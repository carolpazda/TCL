#
# PROGRAMA PARA CÁLCULO DE LAMINADO - DISCIPLINA DE RELAÇÕES CONSTITUTIVAS
#
# Cálculo realizado por meio da regra das misturas

# Bibliotecas utilizadas
using LinearAlgebra

# Calcula a matriz de rotação para cada uma das lâminas
# ang -> Ângulos que serão considerados
function matriz_T(ang::Float64)
	
    # Definição do cosseno e do seno
	c = cosd(ang)
	s = sind(ang)

    # Retornando a matriz T, de transformação
	return [c^2 s^2 2*s*c; s^2 c^2 -2*s*c; -s*c (s*c) ((c^2)-(s^2))]

end #function

# Função que irá realizar a inversão da matriz T, de transformação
function inverse_T(ang::Float64)

    # Da mesma forma que a função anterior, irá fazer o cálculo do seno e do cosseno
	c = cosd(ang)
	s = sind(ang)

    # Retornando a inversa da matriz T
	return [ c^2/(s^4+2*c^2*s^2+c^4)      s^2/(s^4+2*c^2*s^2+c^4)      -(2*c*s)/(s^4+2*c^2*s^2+c^4);
             s^2/(s^4+2*c^2*s^2+c^4)      c^2/(s^4+2*c^2*s^2+c^4)       (2*c*s)/(s^4+2*c^2*s^2+c^4);
            (c*s)/(s^4+2*c^2*s^2+c^4)   -(c*s)/(s^4+2*c^2*s^2+c^4)    -(s^2-c^2)/(s^4+2*c^2*s^2+c^4)]
    
end #function

# Cálculo da matriz de rigidez Q
# Vamos receber os módulos de Elasticidade E11 e E22;
# O módulo de elasticidade transversal G12;
# e o coeficiente de poisson v12
function matriz_Q(E11::Float64, E22::Float64, G12::Float64, v12::Float64) 

	# Cálculo dos termos da matriz de rigidez
    v21 = E22*v12/E11 
	Q11 = E11/(1-v12*v21)
	Q22 = E22/(1-v12*v21)
	Q12 = E22*v12/(1-v12*v21)
	Q66 = G12
 
    # Matriz de rigidez do laminado
 	Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66]

end #function

# Código para a matriz ABBD arcaico, feio... mas que funciona.
function ABBD_laminado(h, n, E11, E22, G12, v12, angulo)

    # Primeiramente vamos armazenar o local dessas matrizes
    A = zeros(3,3)
    B = zeros(3,3)
    D = zeros(3,3)

    # Agora vamos calcular a matriz de rigidez
    Q = matriz_Q(E11, E22, G12, v12) 

    # Calculamos a matriz de Reuter
    R = [1 0 0;0 1 0; 0 0 2] 

    # Matriz de Reuter invertida
    inversa_R= inv(R)

    for i=1:n
        # Calculando a matriz de transformação (rotação) para dada orientação do laminado
        T = matriz_T(angulo[i]) 
	    # E sua inversa também
	    inversa_T = inverse_T(angulo[i])
        Qg = inversa_T*Q*R*T*inversa_R

        # Matriz ABBD
        A = A .+ Qg*(((((n/2)-(i))/n)*-h) - (((n/2 - (i-1))/n)*-h))
        #A = A .+ Qg*(h/n)
        B = B .+ 0.5*Qg*(((((n/2)-(i))/n)*h)^2 - (((n/2 - (i-1))/n)h)^2)
        D = D .+ (1/3)*Qg*(((((n/2)-(i))/n)*-h)^3 - (((n/2 - (i-1))/n)*-h)^3)
    end

    # De forma que como queremos a matriz ABBD
    ABBD = [A B ; B D]

    # Basta retornarmos a mesma
    return ABBD
    
end

# Deformação e tensão em cada lâmina --> Global
function Defor_Tensoes(E11, E22, G12, v12, n, h, angulo, Defor_i_global, K_global)

    tensao_global = zeros(3,n)
    deformacao_global = zeros(3,n)
    tensao_local = zeros(3,n)
    deformacao_local = zeros(3,n)
    z = zeros(n+1)

    # Agora vamos calcular a matriz de rigidez
    Q = matriz_Q(E11, E22, G12, v12) 

    # Calculamos a matriz de Reuter
    R = [1 0 0;0 1 0; 0 0 2] 
    
    # Matriz de Reuter invertida
    inversa_R= inv(R)

    # Delimitação do for para encontrarmos a deformação e tensão global
    for j=1:n

        # Vamos definir o z
        z[j] = 0.5*(((((n/2)-(j))/n)*-h) - (((n/2 - (j-1))/n)*-h))

        # Calculando a matriz de transformação (rotação) para dada orientação do laminado
        T = matriz_T(angulo[j]) 

        # E sua inversa também
        inversa_T = inverse_T(angulo[j])

        # Definição da Qg
        Qg = inversa_T*Q*R*T*inversa_R

        # De forma que temos o seguinte:
        # Deformação Global -> Defo_PlanoMedio[1:3]
        # Rigidez Global -> Defo_PlanoMedio[4:6]

        # E a tensão global vai ser dada por
        tensao_global[:,j] = Qg*(Defor_i_global .+ z[j]*K_global)

        # Tensão no sistema local
        tensao_local[:,j] = T * tensao_global[:,j]

        # Deformação global
        deformacao_global[:,j] = Defor_i_global .+ z[j]*K_global

        # Deformação local
        deformacao_local[:,j] = T * deformacao_global[:,j]

    end #for

    return tensao_global, tensao_local, deformacao_global, deformacao_local

end #function

