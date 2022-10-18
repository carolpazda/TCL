# TCL- Teoria Clássica de Laminados
_Código feito para a disciplina de Relações Constitutivas, utilizando a Teoria Clássica dos Laminados_.

_De autoria: Verônica Caroline Herbst Pazda_

O código está separado em três códigos principais:

**I - Propriedades Laminado**: Você deve utilizar esse main caso você já possua as propriedades do seu laminado: E11, E22, v12 e G12. Primeiro entre no _main_propLaminado.jl_ (é lá que você irá alterar sua entrada). Você deve entrar com as suas propriedades da seguinte forma:

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
  
  # Informe o carregamento do laminado [N e Nm]
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

**II - Regra das Misturas**: Você deve utilizar esse main caso você possua as frações volumétricas ou mássicas do laminado (sua escolha); e os módulos de elasticidade da fibra e da matriz, bem como os coeficientes de poisson da fibra e da matriz. Entre no _main_RegraMisturas.jl_; e nesse caso forneça os seguintes valores:

```
function main() 

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

```

O código irá calcular as propriedades E11, E22, v12, G12 com a função auxiliar presente em _RegraMisturas.jl_. É importante ressaltar que você pode entrar no seu main tanto com a fração volumétrica, tanto com a mássica (o que você possuir).

**III - Halphin-Tsai**: Você deve utilizar esse main caso você possua as frações volumétricas ou mássicas do laminado (sua escolha); e os módulos de elasticidade da fibra e da matriz, bem como os coeficientes de poisson da fibra e da matriz. Entre no _main_HalphinTsai.jl_; e nesse caso forneça os seguintes valores:

```
function main() 

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

```

**E claro, modifique suas entradas conforme a sua necessidade.**

Todos os _mains_ irão retornar a matriz ABBD, calculada pela função:

```
ABBD = ABBD_laminado(v,n_camadas,E11, E22, G12, v12, angulo)
```

Irá calcular as tensões e deformações locais e globais utilizando:

```
tensao_global, tensao_local, deformacao_global, deformacao_local = Defor_Tensoes(E11, E22, G12, v12, n_camadas, v, angulo, Defo_PlanoMedio[1:3], Defo_PlanoMedio[4:6])
```

**CRITÉRIOS DE FALHA:**

Independente do _main_ que você esteja utilizando, o seu código irá retornar os valores de FS e Ms para os seguintes critérios de falha:

**a)** Máxima tensão:
```
SF_MaxTensao = MaximaTensao(n_camadas, tensao_local, Xt, Xc, Yt, Yc, S12)
```

**b)** Tsai-Hill:
```
SF_Hill, Msh = Tsai_Hill(n_camadas, tensao_local, Xt, Xc, Yt, Yc, S12)
```

**c)** Tsai-Wu:
```
SF_Wu, Msw = Wu(Xt, Xc, Yt, Yc, S12, n_camadas, tensao_local)
```
