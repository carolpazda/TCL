# TCL- Teoria Clássica de Laminados
_Código feito para a disciplina de Relações Constitutivas, utilizando a Teoria Clássica dos Laminados_.

_De autoria: Verônica Caroline Herbst Pazda_

O código está separado em três códigos principais:

**Propriedades Laminado** -> Você deve utilizar esse main caso você já possua as propriedades do seu laminado: E11, E22, v12 e G12. Você deve entrar com as suas propriedades da seguinte forma:

```
function main() 

  # Informe ao programa o número de camadas do seu laminado
  n_camadas = 4 

  # Digite os ângulos de orientação das fibras [°] (em ordem, por favor!)
  angulo = [45.0; 0.0; 0.0; 45.0]

  # Digite aqui a espessura do laminado [m]
  espessura = 3E-3
  
  # Informe as propriedades do material (Em Pa meu amigo!)?
  E11 = 19.76E9
  E22 = 1.97E9
  G12 = 0.7E9
  v12 = 0.35
  
  # Informe o carregamento do laminado [N e N/m]
  Nx = 1000000.0
  Ny = 200000.0
  Nxy = 0.0 
  Mx = 0.0
  My = 0.0
  Mxy = 0.0
  
  # E por fim, informe as propriedades de resistência do material -> Todas em Pa!!!
  Xt = 1447.0E6      # -> Resistência à tração 
  Xc = -1447.0E6     # -> Resistência à compressão 
  Yt = 51.7E6        # -> Resistência à tração 
  Yc = -206.0E6      # -> Resistência à compressão 
  S12 = 93.0E6       # -> Resistência ao cisalhamento no plano 1-2 

```

**(2) main_RegraMisturas** -> Você deve utilizar esse main caso você possua as frações volumétricas ou mássicas do laminado (sua escolha); e os módulos de elasticidade da fibra e da matriz, bem como os coeficientes de poisson da fibra e da matriz.

**(3) main_HalphinTsai** -> Você deve utilizar esse main caso você possua as frações volumétricas ou mássicas do laminado (sua escolha); e os módulos de elasticidade da fibra e da matriz, bem como os coeficientes de poisson da fibra e da matriz.

Os critério de falha implementados são:

**a)** Máxima tensão;

**b)** Tsai-Hill;

**c)** Tsai-Wu.
