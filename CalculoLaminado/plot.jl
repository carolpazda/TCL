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
