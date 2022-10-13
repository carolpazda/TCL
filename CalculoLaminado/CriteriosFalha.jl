#
# CRITÉRIOS DE FALHA
# 1 -> Critério da máxima tensão
# 2 -> Critério Tsai-Hill
# 3 -> Critério Tsai-Wu
#
# 1 - Critério da Máxima Tensão
function MaximaTensao(n, stress, XT, XC, YT, YC, S12)

    # primeiro vamos delimitar as entradas
    Cs1 = zeros(n)
    Cs2 = zeros(n)
    Cs12 = zeros(n)
    Ms = zeros(n)

    # For que vai delimitar o nosso Critério
    for i=1:n

        # Pegar as posições corretas
        tensao_local = stress[:,i]

        # Primeiro 
        if tensao_local[1]>0
            Cs1[i] = XT/tensao_local[1]
        else
            Cs1[i] = (-XC)/tensao_local[1]
        end

        # Segundo
        if tensao_local[2]>0
            Cs2[i] = YT/tensao_local[2]
        else
            Cs2[i] = (-YC)/tensao_local[2]
        end

        # Terceiro
        if tensao_local[3]>0
            Cs12[i] = S12/abs(tensao_local[3])
        else
            Cs12[i] = S12/abs(tensao_local[3])
        end

    end #for

    # De forma que temos o seguinte
    Coef = minimum([minimum(Cs1), minimum(Cs2), minimum(Cs12)])
    return Coef

end

# 2 - Critério de Tsai-Hill
function Tsai_Hill(n, stress, XT, XC, YT, YC, S12)

    # Delimitando as entradas
    FS = zeros(n)
    MS = zeros(n)

    # Loop para o critério de falha
    for i=1:n

        # Pegando as posições corretas
        tensao_local = stress[:,i]

        if tensao_local[1]>0 && tensao_local[2]>0
            FS[i] = sqrt(((tensao_local[1]/XT)^2) + ((tensao_local[2]/YT)^2)- (tensao_local[1]*tensao_local[2]/(XT^2)) + ((tensao_local[3]/S12)^2))
            MS[i] = (1/FS[i]) - 1
        
        elseif tensao_local[1]<0 && tensao_local[2]>0
            FS[i] = sqrt(((tensao_local[1]/XC)^2) + ((tensao_local[2]/YT)^2)- (tensao_local[1]*tensao_local[2]/(XC^2)) + ((tensao_local[3]/S12)^2))
            MS[i] = (1/FS[i]) - 1

        elseif tensao_local[1]<0 && tensao_local[2]<0
            FS[i] = sqrt(((tensao_local[1]/XC)^2) + ((tensao_local[2]/YC)^2)- (tensao_local[1]*tensao_local[2]/(XC^2)) + ((tensao_local[3]/S12)^2))
            MS[i] = (1/FS[i]) - 1
        
        else tensao_local[1]>0 && tensao_local[2]<0
            FS[i] = sqrt(((tensao_local[1]/XT)^2) + ((tensao_local[2]/YC)^2)- (tensao_local[1]*tensao_local[2]/(XT^2)) + ((tensao_local[3]/S12)^2))
            MS[i] = (1/FS[i]) - 1

        end #Fim dos if's

    end
            # Retornando o que nos interessa
            return FS, MS

end #function

# 3 - Critério de Tsai-Wu 
function Wu(Xt::Float64, Xc::Float64, Yt::Float64, Yc::Float64, S12::Float64, n::Int64, stress::Array{Float64,2})

    # Delimitando as entradas
	MSw = zeros(n)
	SF = zeros(n)

    # Vamos delimitar o que vamos precisar
	F1 = 1/Xt + 1/Xc 
	F2 = 1/Yt + 1/Yc 
	F11 = -1/(Xt*Xc)
	F22 = -1/(Yt*Yc)
	F66 = (1/S12)^2
	F12 = (-1/2)*sqrt(F11*F22)

    # For para calcular o critério de falha
	for j=1:n 

		σm = stress[:,j]
		A1 = F11*σm[1]^2 + F22*σm[2]^2 + F66*σm[3]^2 + 2*F12*σm[1]*σm[2]
		B1 = F1*σm[1] + F2*σm[2]
		Sf1=((-B1 + sqrt(B1^2 +4*A1))/(2*A1))
        Sf2=abs((-B1 - sqrt(B1^2 +4*A1))/(2*A1))

        if abs(Sf1) > abs(Sf2) 
            Sf=Sf2
        else 
            Sf=Sf1
        end

        SF[j] = Sf
        MSw[j] = Sf-1

	end

    # Retornando o que queremos
	return SF
	
end