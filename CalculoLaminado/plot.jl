#
# Plot do critério de falha
#
using LinearAlgebra
using Plots

function Plot_CriteriosFalha(Xt,Xc,Yt,Yc,S12,tau12)

    # Definindo o plot
    plot1 = collect(LinRange(0.0,0.0,1000))
    sigc = collect(LinRange(Xc, 0.0, 250))
    sigt = collect(LinRange(0.0, Xt, 250))
        
    for k=1:250

        plot1[k] = (Yt*sqrt(-4*Xt^4*tau12^2+(S12^2*Yt^2-4*S12^2*Xt^2)*sigt[k]^2+4*S12^2*Xt^4)-S12*Yt^2*sigt[k])/(2*S12*Xt^2)
        plot1[k+250] = (Yt*sqrt(-4*Xc^4*tau12^2+(S12^2*Yt^2-4*S12^2*Xc^2)*sigc[k]^2+4*S12^2*Xc^4)+S12*Yt^2*sigc[k])/(2*S12*Xc^2)
        plot1[k+500] = (Yc*sqrt(-4*Xc^4*tau12^2+(S12^2*Yc^2-4*S12^2*Xc^2)*sigc[k]^2+4*S12^2*Xc^4)-S12*Yc^2*sigc[k])/(2*S12*Xc^2)
        plot1[k+750] = (Yc*sqrt(-4*Xt^4*tau12^2+(S12^2*Yc^2-4*S12^2*Xt^2)*sigt[k]^2+4*S12^2*Xt^4)+S12*Yc^2*sigt[k])/(2*S12*Xt^2)
    
    end #for

    x_inter1 = -(2*(Yc*Yt*(S12^2*Xc*Xt*Yc^2 + S12^2*Xc*Xt*Yt^2 + S12^2*Xc^2*Yc*Yt + S12^2*Xt^2*Yc*Yt + S12^2*Xc*Xt*Yc*Yt + 3*Xc*Xt*Yc*Yt*tau12^2
               - S12^2*Xc*Xt^2*Yc*Yt^2*(1/(Xc*Xt*Yc*Yt))^(1/2) - S12^2*Xc*Xt^2*Yc^2*Yt*(1/(Xc*Xt*Yc*Yt))^(1/2) - S12^2*Xc^2*Xt*Yc*Yt^2*(1/(Xc*Xt*Yc*Yt))^(1/2) 
               - S12^2*Xc^2*Xt*Yc^2*Yt*(1/(Xc*Xt*Yc*Yt))^(1/2)))^(1/2) - 2*S12*Xc*Yc*Yt - 2*S12*Xt*Yc*Yt + S12*Xc*Xt*Yc*Yt^2*(1/(Xc*Xt*Yc*Yt))^(1/2) + 
                 S12*Xc*Xt*Yc^2*Yt*(1/(Xc*Xt*Yc*Yt))^(1/2))/(3*S12*Yc*Yt)

    x_inter2 = (2*(Yc*Yt*(S12^2*Xc*Xt*Yc^2 + S12^2*Xc*Xt*Yt^2 + S12^2*Xc^2*Yc*Yt + S12^2*Xt^2*Yc*Yt + S12^2*Xc*Xt*Yc*Yt + 3*Xc*Xt*Yc*Yt*tau12^2 
               - S12^2*Xc*Xt^2*Yc*Yt^2*(1/(Xc*Xt*Yc*Yt))^(1/2) - S12^2*Xc*Xt^2*Yc^2*Yt*(1/(Xc*Xt*Yc*Yt))^(1/2) - S12^2*Xc^2*Xt*Yc*Yt^2*(1/(Xc*Xt*Yc*Yt))^(1/2) 
               - S12^2*Xc^2*Xt*Yc^2*Yt*(1/(Xc*Xt*Yc*Yt))^(1/2)))^(1/2) + 2*S12*Xc*Yc*Yt + 2*S12*Xt*Yc*Yt - S12*Xc*Xt*Yc*Yt^2*(1/(Xc*Xt*Yc*Yt))^(1/2) 
               - S12*Xc*Xt*Yc^2*Yt*(1/(Xc*Xt*Yc*Yt))^(1/2))/(3*S12*Yc*Yt)

    sig1 = collect(LinRange(x_inter1,x_inter2, 1000))
    sig2 = collect(LinRange(0.0,0.0,1000))
    sig3 = collect(LinRange(0.0,0.0,1000))

    for k=1:1000

        sig2[k] = (sqrt(4*Xc^2*Xt^2*Yc*Yt*tau12^2-3*S12^2*Xc*Xt*Yc*Yt*sig1[k]^2+(((4*S12^2*Xc*Xt^2+4*S12^2*Xc^2*Xt)*Yc-2*S12^2*Xc^2*Xt^2*Yc^2*sqrt(1/(Xc*Xt*Yc*Yt)))*Yt
                   -2*S12^2*Xc^2*Xt^2*Yc*sqrt(1/(Xc*Xt*Yc*Yt))*Yt^2)*sig1[k]+S12^2*Xc^2*Xt^2*Yt^2-2*S12^2*Xc^2*Xt^2*Yc*Yt+S12^2*Xc^2*Xt^2*Yc^2)
                    -S12*Xc*Xt*Yc*sqrt(1/(Xc*Xt*Yc*Yt))*Yt*sig1[k]+S12*Xc*Xt*Yt+S12*Xc*Xt*Yc)/(2*S12*Xc*Xt)

        sig3[k] = -(sqrt(4*Xc^2*Xt^2*Yc*Yt*tau12^2-3*S12^2*Xc*Xt*Yc*Yt*sig1[k]^2+(((4*S12^2*Xc*Xt^2+4*S12^2*Xc^2*Xt)*Yc-2*S12^2*Xc^2*Xt^2*Yc^2*sqrt(1/(Xc*Xt*Yc*Yt)))*Yt
                  -2*S12^2*Xc^2*Xt^2*Yc*sqrt(1/(Xc*Xt*Yc*Yt))*Yt^2)*sig1[k]+S12^2*Xc^2*Xt^2*Yt^2-2*S12^2*Xc^2*Xt^2*Yc*Yt+S12^2*Xc^2*Xt^2*Yc^2)
                  +S12*Xc*Xt*Yc*sqrt(1/(Xc*Xt*Yc*Yt))*Yt*sig1[k]-S12*Xc*Xt*Yt-S12*Xc*Xt*Yc)/(2*S12*Xc*Xt)
    
    end # for

    #@show sig2
    #@show sig3

    return sig1, sig2, sig3, sigt, sigc, plot1
    
end # function

# Função auxiliar que vai calcular Z
function Auxiliar(espessura, n_camadas, v)

    # Armazenando o h e o z
    h = zeros(n_camadas+1,1)
    z = zeros(3*n_camadas,1)
    h[1] = -v/2

    for j=1:n_camadas

        h[j+1] = h[j] + espessura
        z[3*j-2] = -v/2
        z[3*j-1] = z[3*j-2] + (espessura/2)
        z[3*j] = z[3*j-1] + (espessura/2)
        v = v - (2*espessura)

    end #for

    return h, z

end

# Função auxiliar para teste
function Deformacao_e_TensaoGlobal(Defo_PlanoMedio, n_camadas, z, E11, E22, G12, v12, angulo)

    # Definindo a Deformação e Tensão global
    tensao_global_p = zeros(3, n_camadas*3)
    deformacao_global_p = zeros(3, n_camadas*3)

    # Agora vamos calcular a matriz de rigidez
    Q = matriz_Q(E11, E22, G12, v12) 

    # Calculamos a matriz de Reuter e a sua inversa
    R = [1 0 0;0 1 0; 0 0 2] 
    Rinv = inv(R)


    # Loop que vai calcular
    for j=1:n_camadas

        # Primeiro as deformações correspondentes
        deformacao_global_p[:,(3*j)-2] = Defo_PlanoMedio[1:3] + (z[(3*j)-2]*Defo_PlanoMedio[4:6])
        deformacao_global_p[:,(3*j)-1] = Defo_PlanoMedio[1:3] + (z[(3*j)-1]*Defo_PlanoMedio[4:6])
        deformacao_global_p[:,3*j] = Defo_PlanoMedio[1:3] + (z[3*j]*Defo_PlanoMedio[4:6])

        # E agora as tensões correspondentes

        # Calculando a matriz de transformação (rotação) para dada orientação do laminado
        T = matriz_T(angulo[j]) 

        # E sua inversa também
        inversa_T = inverse_T(angulo[j])

        tensao_global_p[:,(3*j)-2] = inversa_T * Q * R * T * Rinv * deformacao_global_p[:,(3*j)-2]
        tensao_global_p[:,(3*j)-1] = inversa_T * Q * R * T * Rinv * deformacao_global_p[:,(3*j)-1]
        tensao_global_p[:,3*j] = inversa_T * Q * R * T * Rinv * deformacao_global_p[:,3*j]

    end # for

    # Retornando os valores que queremos
    return deformacao_global_p, tensao_global_p

end #function
