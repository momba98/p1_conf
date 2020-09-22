```python
from handcalcs import render 
import forallpeople as si

si.environment('default', top_level=True)

from scipy.optimize import fsolve

import sympy as sp
import numpy as np

from numpy import log as ln
from numpy import pi, arctan, sin, e

from math import ceil as Arredondar
from math import fabs as Absoluto
from math import sqrt

import pandas as pd
```

# PONTIFÍCIA UNIVERSIDADE CATÓLICA DO RIO GRANDE DO SUL

## Curso de Engenharia Mecânica - Conformação Mecânica

---

> **Professor: José Fazzi**
>
> **1ª Tarefa Conformação Mecânica 2020/02**
>
> **Nome: André Mombach dos Santos**

---

Essa tarefa foi desenvolvida pelo autor/aluno (21-22/09/2020) em linguagem Python no ambiente [Jupyter Notebook](https://jupyter.org/), com auxilio das bibliotecas [Numpy](https://numpy.org/) e [Handcalcs](https://github.com/connorferster/handcalcs).

## Questão I  

Uma chapa de latão ($\sigma_e = 26 \cdot |\phi| \cdot 0.50 \cdot \frac{kgf}{mm^2}$), de espessura 15 mm e largura 500 mm, deve ser transformada por laminação a frio em um chapa de 5 mm de espessura. 

Considerando o coeficiente de atrito $\mu=0.10$ e que os cilindros de laminação em cada passe são 15% maior que o mínimo necessário:

1. Programe um processo de laminação para obtenção da chapa desejada e calcule o comprimento final da chapa, sabendo que o comprimento inicial é de 5 m; 
2. Calcule força e potência ideais para o primeiro e o último passe, sabendo que cada par de cilindros gira a 150 RPM;
3. Calcule a força e a potência reais de laminação para o primeiro e o último passe. 

---

## Resposta I 

### 1. Determinação das etapas do projeto e do comprimento final

O primeiro passo é descobrir quantos passes serão necessários para que o processo seja feito de forma homogênea. Para tanto, calcula-se o $\phi_h$ para o processo inteiro, e então o comparamos com o coeficiente $n$, que indica qual a deformação verdadeira máxima homogênea:


```python
%%render 
#short

h_i = (15*mm)

h_f = (5*mm)

phi_h_total = ln(h_f/h_i)
```


\[
\begin{aligned}
h_{i} &= 15.000\ \text{mm}\;
\\[10pt]
h_{f} &= 5.000\ \text{mm}\;
\\[10pt]
\phi_{h_{total}} &= \operatorname{ln} \left( \frac{ h_{f} }{ h_{i} } \right) = \operatorname{ln} \left( \frac{ 5.000\ \text{mm} }{ 15.000\ \text{mm} } \right) &= -1.099
\end{aligned}
\]



```python
%%render
#short

n = 0.5 #retirado da equação de Hollomon

N_passes = Arredondar((Absoluto(phi_h_total))/n)
```


\[
\begin{aligned}
n &= 0.5\;\;\textrm{(retirado da equação de Hollomon)}
\\[10pt]
N_{passes} &= \operatorname{Arredondar} \left( \frac{ \operatorname{Absoluto} \left( \phi_{h_{total}} \right) }{ n } \right) = \operatorname{Arredondar} \left( \frac{ \operatorname{Absoluto} \left( -1.099 \right) }{ 0.5 } \right) &= 3
\end{aligned}
\]


Como demonstrado, serão necessários 3 passes.

Logo, cada passe apresentará a deformação de:


```python
%%render
#sgort

phi_h = phi_h_total/N_passes
```


\[
\begin{aligned}
\phi_{h} &= \frac{ \phi_{h_{total}} }{ N_{passes} } = \frac{ -1.099 }{ 3 } &= -0.366
\end{aligned}
\]


Com $\phi_h$, torna-se possível a determinação das espessuras em cada etapa:


```python
%%render
#short

h_i_1 = h_i
h_f_1 = h_i_1*e**phi_h
Delta_h_1 = h_i_1 - h_f_1

h_i_2 = h_f_1
h_f_2 = h_i_2*e**phi_h
Delta_h_2 = h_i_2 - h_f_2

h_i_3 = h_f_2
h_f_3 = h_i_3*e**phi_h
Delta_h_3 = h_i_3 - h_f_3

```


\[
\begin{aligned}
h_{i_{1}} &= 15.000\ \text{mm}\;
\\[10pt]
h_{f_{1}} &= h_{i_{1}} \cdot \left( e \right) ^{ \phi_{h} } = 15.000\ \text{mm} \cdot \left( 2.718 \right) ^{ -0.366 } &= 10.400\ \text{mm}
\\[10pt]
\Delta_{h_{1}} &= h_{i_{1}} - h_{f_{1}} = 15.000\ \text{mm} - 10.400\ \text{mm} &= 4.600\ \text{mm}
\\[10pt]
h_{i_{2}} &= 10.400\ \text{mm}\;
\\[10pt]
h_{f_{2}} &= h_{i_{2}} \cdot \left( e \right) ^{ \phi_{h} } = 10.400\ \text{mm} \cdot \left( 2.718 \right) ^{ -0.366 } &= 7.211\ \text{mm}
\\[10pt]
\Delta_{h_{2}} &= h_{i_{2}} - h_{f_{2}} = 10.400\ \text{mm} - 7.211\ \text{mm} &= 3.189\ \text{mm}
\\[10pt]
h_{i_{3}} &= 7.211\ \text{mm}\;
\\[10pt]
h_{f_{3}} &= h_{i_{3}} \cdot \left( e \right) ^{ \phi_{h} } = 7.211\ \text{mm} \cdot \left( 2.718 \right) ^{ -0.366 } &= 5.000\ \text{mm}
\\[10pt]
\Delta_{h_{3}} &= h_{i_{3}} - h_{f_{3}} = 7.211\ \text{mm} - 5.000\ \text{mm} &= 2.211\ \text{mm}
\end{aligned}
\]


Com todos estes dados, já se torna possível o cálculo do comprimento final. 

De acordo com a conservação de massa, o volume da peça não deve mudar. Como não se considera variação na largura no processo ed conformação a frio, conclui-se que a deformação verdadeira no comprimento é exatamente a deformação verdadeira na espessura com o sinal invertido.


```python
%%render
#long

phi_l = -phi_h

l_i = (5*m) #comprimento inicial da peça

l_i_1 = l_i #comprimento da peça antes de passar pelo primeiro rolo
l_f_1 = l_i_1*e**phi_l #comprimento da peça depois de passar pelo primeiro rolo
Delta_l_1 = l_f_1 - l_i_1 #variação de comprimento na etapa do primeiro rolo

l_i_2 = l_f_1
l_f_2 = l_i_2*e**phi_l
Delta_l_2 = l_f_2 - l_i_2

l_i_3 = l_f_2
l_f_3 = l_i_3*e**phi_l
Delta_l_3 = l_f_3 - l_i_3


```


\[
\begin{aligned}
\phi_{l} &= 0.366\;
\\[10pt]
l_{i} &= 5.000\ \text{m}\;\;\textrm{(comprimento inicial da peça)}
\\[10pt]
l_{i_{1}} &= 5.000\ \text{m}\;\;\textrm{(comprimento da peça antes de passar pelo primeiro rolo)}
\\[10pt]
l_{f_{1}} &= l_{i_{1}} \cdot \left( e \right) ^{ \phi_{l} } \\&= 5.000\ \text{m} \cdot \left( 2.718 \right) ^{ 0.366 } \\&= 7.211\ \text{m}\;\;\textrm{(comprimento da peça depois de passar pelo primeiro rolo)}\\
\\[10pt]
\Delta_{l_{1}} &= l_{f_{1}} - l_{i_{1}} \\&= 7.211\ \text{m} - 5.000\ \text{m} \\&= 2.211\ \text{m}\;\;\textrm{(variação de comprimento na etapa do primeiro rolo)}\\
\\[10pt]
l_{i_{2}} &= 7.211\ \text{m}\;
\\[10pt]
l_{f_{2}} &= l_{i_{2}} \cdot \left( e \right) ^{ \phi_{l} } \\&= 7.211\ \text{m} \cdot \left( 2.718 \right) ^{ 0.366 } \\&= 10.400\ \text{m}\\
\\[10pt]
\Delta_{l_{2}} &= l_{f_{2}} - l_{i_{2}} \\&= 10.400\ \text{m} - 7.211\ \text{m} \\&= 3.189\ \text{m}\\
\\[10pt]
l_{i_{3}} &= 10.400\ \text{m}\;
\\[10pt]
l_{f_{3}} &= l_{i_{3}} \cdot \left( e \right) ^{ \phi_{l} } \\&= 10.400\ \text{m} \cdot \left( 2.718 \right) ^{ 0.366 } \\&= 15.000\ \text{m}\\
\\[10pt]
\Delta_{l_{3}} &= l_{f_{3}} - l_{i_{3}} \\&= 15.000\ \text{m} - 10.400\ \text{m} \\&= 4.600\ \text{m}\\
\end{aligned}
\]


De modo óbvio, se a conformação tenta diminuir a espessura por um fator de 3, a largura será aumentada também por um fator de 3.

Continuando com o projeto do processo:

Com todos $\Delta_h$, calcula-se o raio mínimo requirido pelo rolo de compressão de cada etapa:


```python
%%render
#short

mu = 0.10 #identificado pelo enunciado

CS = 1.25 #25% de segurança

R_1 = CS*Delta_h_1/((sin(arctan(mu)))**2)
R_2 = CS*Delta_h_2/((sin(arctan(mu)))**2)
R_3 = CS*Delta_h_3/((sin(arctan(mu)))**2)
```


\[
\begin{aligned}
\mu &= 0.1\;\;\textrm{(identificado pelo enunciado)}
\\[10pt]
CS &= 1.25\;\;\textrm{(25% de segurança)}
\\[10pt]
R_{1} &= CS \cdot \frac{ \Delta_{h_{1}} }{ \left( \sin{ \left( \operatorname{arctan} \left( \mu \right) \right) } \right) ^{ 2 } } = 1.25 \cdot \frac{ 4.600\ \text{mm} }{ \left( \sin{ \left( \operatorname{arctan} \left( 0.1 \right) \right) } \right) ^{ 2 } } &= 580.697\ \text{mm}
\\[10pt]
R_{2} &= CS \cdot \frac{ \Delta_{h_{2}} }{ \left( \sin{ \left( \operatorname{arctan} \left( \mu \right) \right) } \right) ^{ 2 } } = 1.25 \cdot \frac{ 3.189\ \text{mm} }{ \left( \sin{ \left( \operatorname{arctan} \left( 0.1 \right) \right) } \right) ^{ 2 } } &= 402.633\ \text{mm}
\\[10pt]
R_{3} &= CS \cdot \frac{ \Delta_{h_{3}} }{ \left( \sin{ \left( \operatorname{arctan} \left( \mu \right) \right) } \right) ^{ 2 } } = 1.25 \cdot \frac{ 2.211\ \text{mm} }{ \left( \sin{ \left( \operatorname{arctan} \left( 0.1 \right) \right) } \right) ^{ 2 } } &= 279.170\ \text{mm}
\end{aligned}
\]


E tambem o parâmetro $l_d$:


```python
%%render
#short

l_d_1 = (R_1*Delta_h_1)**0.5
l_d_2 = (R_2*Delta_h_2)**0.5
l_d_3 = (R_3*Delta_h_3)**0.5
```


\[
\begin{aligned}
l_{d_{1}} &= \left( R_{1} \cdot \Delta_{h_{1}} \right) ^{ 0.5 } = \left( 580.697\ \text{mm} \cdot 4.600\ \text{mm} \right) ^{ 0.5 } &= 51.681\ \text{mm}
\\[10pt]
l_{d_{2}} &= \left( R_{2} \cdot \Delta_{h_{2}} \right) ^{ 0.5 } = \left( 402.633\ \text{mm} \cdot 3.189\ \text{mm} \right) ^{ 0.5 } &= 35.834\ \text{mm}
\\[10pt]
l_{d_{3}} &= \left( R_{3} \cdot \Delta_{h_{3}} \right) ^{ 0.5 } = \left( 279.170\ \text{mm} \cdot 2.211\ \text{mm} \right) ^{ 0.5 } &= 24.846\ \text{mm}
\end{aligned}
\]


Como o exercício menciona um processo a frio, não se deve preocupar com a variação na espessura. Pode-se estimar que a variação na entrada, na saída e consequentemente a média são todos o mesmo valor.


```python
%%render
#short

b_m = (500*mm)
```


\[
\begin{aligned}
b_{m} &= 500.000\ \text{mm}\;
\end{aligned}
\]


A área de contato $A_c$, no entanto, varia entre etapas:


```python
%%render
#short

A_c_1 = l_d_1 * b_m
A_c_2 = l_d_2 * b_m
A_c_3 = l_d_3 * b_m

```


\[
\begin{aligned}
A_{c_{1}} &= l_{d_{1}} \cdot b_{m} = 51.681\ \text{mm} \cdot 500.000\ \text{mm} &= 25840.681\ \text{mm}^{2.0}
\\[10pt]
A_{c_{2}} &= l_{d_{2}} \cdot b_{m} = 35.834\ \text{mm} \cdot 500.000\ \text{mm} &= 17916.928\ \text{mm}^{2.0}
\\[10pt]
A_{c_{3}} &= l_{d_{3}} \cdot b_{m} = 24.846\ \text{mm} \cdot 500.000\ \text{mm} &= 12422.904\ \text{mm}^{2.0}
\end{aligned}
\]


Durante um processo de laminação a frio, se faz necessário a realização de recozimento entre as etapas. Caso não seja feito, o material não voltará às suas características iniciais e o encruamento passará do seu ponto de homogeniedade, acarretando em deformações heterogêneas ao longo da lamina. Com certeza, não é algo desejável.

Dito isto, o cálculo da tensão de escoamento $\sigma_e$ para cada etapa será feito com os valores de $\sigma_e$ antes de passar ser comprimido rolo e $\sigma_e$ depois de ser comprimido pelo rolo. Como o material será recozido em cada etapa, a tensão de escoamento inicial sempre será a mesma e correspondente ao início da região elástica, com $\phi=0.002$.

Como o material será conformado sempre a mesma quantidade ($\phi_{h}=cte$), a $\sigma_e$ de saída também terá sempre o mesmo valor. 


```python
%%render
#long

phi_plastico = 0.002 #deformação verdadeira no início da deformação plástica

sigma_o = 26*9.81*MPa #transformar de kgf/mm² para MPa

sigma_e_i = sigma_o*phi_plastico**n #como o material é recozido, a tensão de entrada é sempre a mesma

sigma_e_f = sigma_o*Absoluto(phi_h)**n #como $\phi$ é constate durante as etapas, a tensão de saída também é sempre a mesma

sigma_medio = (sigma_e_i + sigma_e_f)/2
```


\[
\begin{aligned}
\phi_{plastico} &= 0.002\;\;\textrm{(deformação verdadeira no início da deformação plástica)}
\\[10pt]
\sigma_{o} &= 26 \cdot 9.81 \cdot MPa \\&= 26 \cdot 9.81 \cdot MPa \\&= 255.060\ \text{MPa}\;\;\textrm{(transformar de kgf/mm² para MPa)}\\
\\[10pt]
\sigma_{e_{i}} &= \sigma_{o} \cdot \left( \phi_{plastico} \right) ^{ n } \\&= 255.060\ \text{MPa} \cdot \left( 0.002 \right) ^{ 0.5 } \\&= 11.407\ \text{MPa}\;\;\textrm{(como o material é recozido, a tensão de entrada é sempre a mesma)}\\
\\[10pt]
\sigma_{e_{f}} &= \sigma_{o} \cdot \operatorname{Absoluto} \left( \phi_{h} \right) ^{ n } \\&= 255.060\ \text{MPa} \cdot \operatorname{Absoluto} \left( -0.366 \right) ^{ 0.5 } \\&= 154.349\ \text{MPa}\;\;\textrm{(como $\phi$ é constate durante as etapas, a tensão de saída também é sempre a mesma)}\\
\\[10pt]
\sigma_{medio} &= \frac{ \sigma_{e_{i}} + \sigma_{e_{f}} }{ 2 } \\&= \frac{ 11.407\ \text{MPa} + 154.349\ \text{MPa} }{ 2 } \\&= 82.878\ \text{MPa}\\
\end{aligned}
\]


### 2. Obtenção das grandezas ideais

As forças ideais $F_{id}$, os momentos ideais $M_{id}$ e as potências $P_{id}$ de cada etapa:


```python
%%render
#long

omega_angular = 150*(2*pi/60)*(1/s) #transformação da rotação de RPM para rad/s

F_i_1 = sigma_medio*A_c_1
M_i_1 = F_i_1*l_d_1
P_i_1 = M_i_1*omega_angular

F_i_2 = sigma_medio*A_c_2
M_i_2 = F_i_2*l_d_2
P_i_2 = M_i_2*omega_angular

F_i_3 = sigma_medio*A_c_3
M_i_3 = F_i_3*l_d_3
P_i_3 = M_i_3*omega_angular
```


\[
\begin{aligned}
\omega_{angular} &= 150 \cdot \left( 2 \cdot \frac{ \pi }{ 60 } \right) \cdot \left( \frac{ 1 }{ s } \right) \\&= 150 \cdot \left( 2 \cdot \frac{ 3.142 }{ 60 } \right) \cdot \left( \frac{ 1 }{ s } \right) \\&= 15.708\ \text{Hz}\;\;\textrm{(transformação da rotação de RPM para rad/s)}\\
\\[10pt]
F_{i_{1}} &= \sigma_{medio} \cdot A_{c_{1}} \\&= 82.878\ \text{MPa} \cdot 25840.681\ \text{mm}^{2.0} \\&= 2.142\ \text{MN}\\
\\[10pt]
M_{i_{1}} &= F_{i_{1}} \cdot l_{d_{1}} \\&= 2.142\ \text{MN} \cdot 51.681\ \text{mm} \\&= 110.682\ \text{kN} \cdot \text{m}\\
\\[10pt]
P_{i_{1}} &= M_{i_{1}} \cdot \omega_{angular} \\&= 110.682\ \text{kN} \cdot \text{m} \cdot 15.708\ \text{Hz} \\&= 1.739\ \text{MW}\\
\\[10pt]
F_{i_{2}} &= \sigma_{medio} \cdot A_{c_{2}} \\&= 82.878\ \text{MPa} \cdot 17916.928\ \text{mm}^{2.0} \\&= 1.485\ \text{MN}\\
\\[10pt]
M_{i_{2}} &= F_{i_{2}} \cdot l_{d_{2}} \\&= 1.485\ \text{MN} \cdot 35.834\ \text{mm} \\&= 53.210\ \text{kN} \cdot \text{m}\\
\\[10pt]
P_{i_{2}} &= M_{i_{2}} \cdot \omega_{angular} \\&= 53.210\ \text{kN} \cdot \text{m} \cdot 15.708\ \text{Hz} \\&= 835.825\ \text{kW}\\
\\[10pt]
F_{i_{3}} &= \sigma_{medio} \cdot A_{c_{3}} \\&= 82.878\ \text{MPa} \cdot 12422.904\ \text{mm}^{2.0} \\&= 1.030\ \text{MN}\\
\\[10pt]
M_{i_{3}} &= F_{i_{3}} \cdot l_{d_{3}} \\&= 1.030\ \text{MN} \cdot 24.846\ \text{mm} \\&= 25.581\ \text{kN} \cdot \text{m}\\
\\[10pt]
P_{i_{3}} &= M_{i_{3}} \cdot \omega_{angular} \\&= 25.581\ \text{kN} \cdot \text{m} \cdot 15.708\ \text{Hz} \\&= 401.823\ \text{kW}\\
\end{aligned}
\]


Potências ideais relativamente altas. Provavelmente se deve ao raio ser grande, bem como a área de contato.

### 3. Finalmente, obtenção dos parâmetros não-ideais ou reais

#### 3.1 Primeiro Passo

Alguns fatores serão utilizados para correção de forças, momentos, etc.

O primeiro passo para obter a força real é o cálculo do raio corrigido. Marca, também, o início do processo iterativo. Considera-se os rolos feitos de aço. A força referência é a força ideal $F_{i_1}$, por enquanto:


```python
%%render
#long

c = (((4.6*10**(-4))/9.81)*(mm**2)/N) # constante do aço transformada para de kgf/mm² para Pa

Rlinha_1_1 = R_1*(1+((c*F_i_1)/(b_m*Delta_h_1)))
```


\[
\begin{aligned}
c &= 46.891\ \text{pPa}^{-1}\;\;\textrm{(constante do aço transformada para de kgf/mm² para Pa)}
\\[10pt]
Rlinha_{1_{1}} &= R_{1} \cdot \left( 1 + \left( \frac{ c \cdot F_{i_{1}} }{ b_{m} \cdot \Delta_{h_{1}} } \right) \right) \\&= 580.697\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 2.142\ \text{MN} }{ 500.000\ \text{mm} \cdot 4.600\ \text{mm} } \right) \right) \\&= 606.054\ \text{mm}\\
\end{aligned}
\]


As constantes que serão utilizados para obter o primeiro fator de correção:


```python
%%render
#long

Curva_empírica = mu * sqrt(Rlinha_1_1/h_f_1)

epsilon = Delta_h_1/h_i_1
```


\[
\begin{aligned}
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{1_{1}} }{ h_{f_{1}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 606.054\ \text{mm} }{ 10.400\ \text{mm} } } \\&= 0.763\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{1}} }{ h_{i_{1}} } \\&= \frac{ 4.600\ \text{mm} }{ 15.000\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]


Do cruzamento deste dois valores do Gráfico 1 da apostila (pg. 35), conclui-se que:


```python
%%render
#short

f_1_Passo1_1 = 1.1 #aproximadamente
```


\[
\begin{aligned}
f_{1_{Passo1_{1}}} &= 1.1\;\;\textrm{(aproximadamente)}
\end{aligned}
\]


E, consequentemente, a primeira força corrigida do primeiro passo $F_{c_{1_{1}}}$


```python
%%render
#long

F_c_1_1 = b_m*((Rlinha_1_1*Delta_h_1)**0.5)*1.15*sigma_medio*f_1_Passo1_1
```


\[
\begin{aligned}
F_{c_{1_{1}}} &= b_{m} \cdot \left( \left( Rlinha_{1_{1}} \cdot \Delta_{h_{1}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo1_{1}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 606.054\ \text{mm} \cdot 4.600\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.1 \\&= 2.768\ \text{MN}\\
\end{aligned}
\]


A força ideal $F_{i_1}$ foi encontrada como valendo $2.142 MN.$ Agora, a primeira força corrigida $F_{c_{1_{1}}}$ foi encontrada como $2.768 MN$. A diferença é gritante. Mais um processo iterativo é necessário. Dessa vez, a força de referência mudou:


```python
%%render
#long

Rlinha_1_2 = R_1*(1+((c*F_c_1_1)/(b_m*Delta_h_1)))
```


\[
\begin{aligned}
Rlinha_{1_{2}} &= R_{1} \cdot \left( 1 + \left( \frac{ c \cdot F_{c_{1_{1}}} }{ b_{m} \cdot \Delta_{h_{1}} } \right) \right) \\&= 580.697\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 2.768\ \text{MN} }{ 500.000\ \text{mm} \cdot 4.600\ \text{mm} } \right) \right) \\&= 613.466\ \text{mm}\\
\end{aligned}
\]


A discrepância entre os raios foi bem menor nesta etapa. Bom indicador.

Continuando:


```python
%%render
#long

Curva_empírica = mu * sqrt(Rlinha_1_2/h_f_1)

epsilon = Delta_h_1/h_i_1
```


\[
\begin{aligned}
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{1_{2}} }{ h_{f_{1}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 613.466\ \text{mm} }{ 10.400\ \text{mm} } } \\&= 0.768\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{1}} }{ h_{i_{1}} } \\&= \frac{ 4.600\ \text{mm} }{ 15.000\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]


Valores bem próximos aos da iteração breve. Portanto:


```python
%%render
#long

f_1_Passo1_2 = 1.12 #aproximadamente

F_c_1_2 = b_m*((Rlinha_1_2*Delta_h_1)**0.5)*1.15*sigma_medio*f_1_Passo1_2
```


\[
\begin{aligned}
f_{1_{Passo1_{2}}} &= 1.12\;\;\textrm{(aproximadamente)}
\\[10pt]
F_{c_{1_{2}}} &= b_{m} \cdot \left( \left( Rlinha_{1_{2}} \cdot \Delta_{h_{1}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo1_{2}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 613.466\ \text{mm} \cdot 4.600\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.12 \\&= 2.835\ \text{MN}\\
\end{aligned}
\]


Bem parecido. Mais uma última iteração para confirmar:


```python
%%render
#long

Rlinha_1_3 = R_1*(1+((c*F_c_1_2)/(b_m*Delta_h_1)))

Curva_empírica = mu * sqrt(Rlinha_1_3/h_f_1)

epsilon = Delta_h_1/h_i_1
```


\[
\begin{aligned}
Rlinha_{1_{3}} &= R_{1} \cdot \left( 1 + \left( \frac{ c \cdot F_{c_{1_{2}}} }{ b_{m} \cdot \Delta_{h_{1}} } \right) \right) \\&= 580.697\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 2.835\ \text{MN} }{ 500.000\ \text{mm} \cdot 4.600\ \text{mm} } \right) \right) \\&= 614.265\ \text{mm}\\
\\[10pt]
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{1_{3}} }{ h_{f_{1}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 614.265\ \text{mm} }{ 10.400\ \text{mm} } } \\&= 0.769\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{1}} }{ h_{i_{1}} } \\&= \frac{ 4.600\ \text{mm} }{ 15.000\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]



```python
%%render
#long

f_1_Passo1_3 = 1.12 #aproximadamente

F_c_1_3 = b_m*((Rlinha_1_3*Delta_h_1)**0.5)*1.15*sigma_medio*f_1_Passo1_3
```


\[
\begin{aligned}
f_{1_{Passo1_{3}}} &= 1.12\;\;\textrm{(aproximadamente)}
\\[10pt]
F_{c_{1_{3}}} &= b_{m} \cdot \left( \left( Rlinha_{1_{3}} \cdot \Delta_{h_{1}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo1_{3}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 614.265\ \text{mm} \cdot 4.600\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.12 \\&= 2.837\ \text{MN}\\
\end{aligned}
\]


Houve convergência da grandeza. Calcula-se agora o momento real $M_{c_1}$, com as mesmas grandezas utilizadas no Gráfico 1, porém agora olhando para o Gráfico 2. 


```python
%%render
#long

f_2_Passo1 = 0.115

M_c_1 = 2*R_1*(((h_i_1)**2)/h_f_1)*sigma_medio*b_m*f_2_Passo1
```


\[
\begin{aligned}
f_{2_{Passo1}} &= 0.115\;
\\[10pt]
M_{c_{1}} &= 2 \cdot R_{1} \cdot \left( \frac{ \left( h_{i_{1}} \right) ^{ 2 } }{ h_{f_{1}} } \right) \cdot \sigma_{medio} \cdot b_{m} \cdot f_{2_{Passo1}} \\&= 2 \cdot 580.697\ \text{mm} \cdot \left( \frac{ \left( 15.000\ \text{mm} \right) ^{ 2 } }{ 10.400\ \text{mm} } \right) \cdot 82.878\ \text{MPa} \cdot 500.000\ \text{mm} \cdot 0.115 \\&= 119.734\ \text{kN} \cdot \text{m}\\
\end{aligned}
\]


E finalmente a potência corrigida $P_{c_1}$ para o primeiro passo:


```python
%%render
#short

P_c_1 = M_c_1*omega_angular
```


\[
\begin{aligned}
P_{c_{1}} &= M_{c_{1}} \cdot \omega_{angular} = 119.734\ \text{kN} \cdot \text{m} \cdot 15.708\ \text{Hz} &= 1.881\ \text{MW}
\end{aligned}
\]


#### 3.2 Segundo Passo

Repete-se todo o processo do primeiro passo. Não será comentado para economizar espaço.


```python
%%render
#long

Rlinha_2_1 = R_2*(1+((c*F_i_2)/(b_m*Delta_h_2)))
```


\[
\begin{aligned}
Rlinha_{2_{1}} &= R_{2} \cdot \left( 1 + \left( \frac{ c \cdot F_{i_{2}} }{ b_{m} \cdot \Delta_{h_{2}} } \right) \right) \\&= 402.633\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 1.485\ \text{MN} }{ 500.000\ \text{mm} \cdot 3.189\ \text{mm} } \right) \right) \\&= 420.214\ \text{mm}\\
\end{aligned}
\]



```python
%%render
#long

Curva_empírica = mu * sqrt(Rlinha_2_1/h_f_2)

epsilon = Delta_h_2/h_i_2
```


\[
\begin{aligned}
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{2_{1}} }{ h_{f_{2}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 420.214\ \text{mm} }{ 7.211\ \text{mm} } } \\&= 0.763\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{2}} }{ h_{i_{2}} } \\&= \frac{ 3.189\ \text{mm} }{ 10.400\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]



```python
%%render
#short

f_1_Passo2_1 = 1.1 #aproximadamente
```


\[
\begin{aligned}
f_{1_{Passo2_{1}}} &= 1.1\;\;\textrm{(aproximadamente)}
\end{aligned}
\]



```python
%%render
#long

F_c_2_1 = b_m*((Rlinha_2_1*Delta_h_2)**0.5)*1.15*sigma_medio*f_1_Passo2_1
```


\[
\begin{aligned}
F_{c_{2_{1}}} &= b_{m} \cdot \left( \left( Rlinha_{2_{1}} \cdot \Delta_{h_{2}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo2_{1}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 420.214\ \text{mm} \cdot 3.189\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.1 \\&= 1.919\ \text{MN}\\
\end{aligned}
\]



```python
%%render
#long

Rlinha_2_2 = R_2*(1+((c*F_c_2_1)/(b_m*Delta_h_2)))
```


\[
\begin{aligned}
Rlinha_{2_{2}} &= R_{2} \cdot \left( 1 + \left( \frac{ c \cdot F_{c_{2_{1}}} }{ b_{m} \cdot \Delta_{h_{2}} } \right) \right) \\&= 402.633\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 1.919\ \text{MN} }{ 500.000\ \text{mm} \cdot 3.189\ \text{mm} } \right) \right) \\&= 425.354\ \text{mm}\\
\end{aligned}
\]



```python
%%render
#long

Curva_empírica = mu * sqrt(Rlinha_2_2/h_f_2)

epsilon = Delta_h_2/h_i_2
```


\[
\begin{aligned}
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{2_{2}} }{ h_{f_{2}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 425.354\ \text{mm} }{ 7.211\ \text{mm} } } \\&= 0.768\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{2}} }{ h_{i_{2}} } \\&= \frac{ 3.189\ \text{mm} }{ 10.400\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]



```python
%%render
#long

f_1_Passo2_2 = 1.11 #aproximadamente

F_c_2_2 = b_m*((Rlinha_2_2*Delta_h_2)**0.5)*1.15*sigma_medio*f_1_Passo2_2
```


\[
\begin{aligned}
f_{1_{Passo2_{2}}} &= 1.11\;\;\textrm{(aproximadamente)}
\\[10pt]
F_{c_{2_{2}}} &= b_{m} \cdot \left( \left( Rlinha_{2_{2}} \cdot \Delta_{h_{2}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo2_{2}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 425.354\ \text{mm} \cdot 3.189\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.11 \\&= 1.948\ \text{MN}\\
\end{aligned}
\]



```python
%%render
#long

Rlinha_2_3 = R_2*(1+((c*F_c_2_2)/(b_m*Delta_h_2)))

Curva_empírica = mu * sqrt(Rlinha_2_3/h_f_2)

epsilon = Delta_h_2/h_i_2
```


\[
\begin{aligned}
Rlinha_{2_{3}} &= R_{2} \cdot \left( 1 + \left( \frac{ c \cdot F_{c_{2_{2}}} }{ b_{m} \cdot \Delta_{h_{2}} } \right) \right) \\&= 402.633\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 1.948\ \text{MN} }{ 500.000\ \text{mm} \cdot 3.189\ \text{mm} } \right) \right) \\&= 425.700\ \text{mm}\\
\\[10pt]
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{2_{3}} }{ h_{f_{2}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 425.700\ \text{mm} }{ 7.211\ \text{mm} } } \\&= 0.768\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{2}} }{ h_{i_{2}} } \\&= \frac{ 3.189\ \text{mm} }{ 10.400\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]



```python
%%render
#long

f_1_Passo2_3 = 1.12 #aproximadamente

F_c_2_3 = b_m*((Rlinha_2_3*Delta_h_2)**0.5)*1.15*sigma_medio*f_1_Passo2_3
```


\[
\begin{aligned}
f_{1_{Passo2_{3}}} &= 1.12\;\;\textrm{(aproximadamente)}
\\[10pt]
F_{c_{2_{3}}} &= b_{m} \cdot \left( \left( Rlinha_{2_{3}} \cdot \Delta_{h_{2}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo2_{3}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 425.700\ \text{mm} \cdot 3.189\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.12 \\&= 1.967\ \text{MN}\\
\end{aligned}
\]



```python
%%render
#long

f_2_Passo2 = 0.112

M_c_2 = 2*R_2*(((h_i_2)**2)/h_f_2)*sigma_medio*b_m*f_2_Passo2
```


\[
\begin{aligned}
f_{2_{Passo2}} &= 0.112\;
\\[10pt]
M_{c_{2}} &= 2 \cdot R_{2} \cdot \left( \frac{ \left( h_{i_{2}} \right) ^{ 2 } }{ h_{f_{2}} } \right) \cdot \sigma_{medio} \cdot b_{m} \cdot f_{2_{Passo2}} \\&= 2 \cdot 402.633\ \text{mm} \cdot \left( \frac{ \left( 10.400\ \text{mm} \right) ^{ 2 } }{ 7.211\ \text{mm} } \right) \cdot 82.878\ \text{MPa} \cdot 500.000\ \text{mm} \cdot 0.112 \\&= 56.060\ \text{kN} \cdot \text{m}\\
\end{aligned}
\]


E finalmente a potência corrigida $P_{c_1}$ para o primeiro passo:


```python
%%render
#short

P_c_2 = M_c_2*omega_angular
```


\[
\begin{aligned}
P_{c_{2}} &= M_{c_{2}} \cdot \omega_{angular} = 56.060\ \text{kN} \cdot \text{m} \cdot 15.708\ \text{Hz} &= 880.596\ \text{kW}
\end{aligned}
\]


#### 3.3 Terceiro (e último!) Passo


```python
%%render
#long

Rlinha_3_1 = R_3*(1+((c*F_i_3)/(b_m*Delta_h_3)))
```


\[
\begin{aligned}
Rlinha_{3_{1}} &= R_{3} \cdot \left( 1 + \left( \frac{ c \cdot F_{i_{3}} }{ b_{m} \cdot \Delta_{h_{3}} } \right) \right) \\&= 279.170\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 1.030\ \text{MN} }{ 500.000\ \text{mm} \cdot 2.211\ \text{mm} } \right) \right) \\&= 291.360\ \text{mm}\\
\end{aligned}
\]



```python
%%render
#long

Curva_empírica = mu * sqrt(Rlinha_3_1/h_f_3)

epsilon = Delta_h_3/h_i_3
```


\[
\begin{aligned}
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{3_{1}} }{ h_{f_{3}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 291.360\ \text{mm} }{ 5.000\ \text{mm} } } \\&= 0.763\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{3}} }{ h_{i_{3}} } \\&= \frac{ 2.211\ \text{mm} }{ 7.211\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]



```python
%%render
#short

f_1_Passo3_1 = 1.1 #aproximadamente
```


\[
\begin{aligned}
f_{1_{Passo3_{1}}} &= 1.1\;\;\textrm{(aproximadamente)}
\end{aligned}
\]



```python
%%render
#long

F_c_3_1 = b_m*((Rlinha_3_1*Delta_h_3)**0.5)*1.15*sigma_medio*f_1_Passo3_1
```


\[
\begin{aligned}
F_{c_{3_{1}}} &= b_{m} \cdot \left( \left( Rlinha_{3_{1}} \cdot \Delta_{h_{3}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo3_{1}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 291.360\ \text{mm} \cdot 2.211\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.1 \\&= 1.331\ \text{MN}\\
\end{aligned}
\]



```python
%%render
#long

Rlinha_3_2 = R_3*(1+((c*F_c_3_1)/(b_m*Delta_h_3)))
```


\[
\begin{aligned}
Rlinha_{3_{2}} &= R_{3} \cdot \left( 1 + \left( \frac{ c \cdot F_{c_{3_{1}}} }{ b_{m} \cdot \Delta_{h_{3}} } \right) \right) \\&= 279.170\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 1.331\ \text{MN} }{ 500.000\ \text{mm} \cdot 2.211\ \text{mm} } \right) \right) \\&= 294.924\ \text{mm}\\
\end{aligned}
\]



```python
%%render
#long

Curva_empírica = mu * sqrt(Rlinha_3_2/h_f_3)

epsilon = Delta_h_3/h_i_3
```


\[
\begin{aligned}
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{3_{2}} }{ h_{f_{3}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 294.924\ \text{mm} }{ 5.000\ \text{mm} } } \\&= 0.768\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{3}} }{ h_{i_{3}} } \\&= \frac{ 2.211\ \text{mm} }{ 7.211\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]



```python
%%render
#long

f_1_Passo3_2 = 1.11 #aproximadamente

F_c_3_2 = b_m*((Rlinha_3_2*Delta_h_3)**0.5)*1.15*sigma_medio*f_1_Passo3_2
```


\[
\begin{aligned}
f_{1_{Passo3_{2}}} &= 1.11\;\;\textrm{(aproximadamente)}
\\[10pt]
F_{c_{3_{2}}} &= b_{m} \cdot \left( \left( Rlinha_{3_{2}} \cdot \Delta_{h_{3}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo3_{2}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 294.924\ \text{mm} \cdot 2.211\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.11 \\&= 1.351\ \text{MN}\\
\end{aligned}
\]



```python
%%render
#long

Rlinha_3_3 = R_3*(1+((c*F_c_3_2)/(b_m*Delta_h_3)))

Curva_empírica = mu * sqrt(Rlinha_3_3/h_f_3)

epsilon = Delta_h_3/h_i_3
```


\[
\begin{aligned}
Rlinha_{3_{3}} &= R_{3} \cdot \left( 1 + \left( \frac{ c \cdot F_{c_{3_{2}}} }{ b_{m} \cdot \Delta_{h_{3}} } \right) \right) \\&= 279.170\ \text{mm} \cdot \left( 1 + \left( \frac{ 46.891\ \text{pPa}^{-1} \cdot 1.351\ \text{MN} }{ 500.000\ \text{mm} \cdot 2.211\ \text{mm} } \right) \right) \\&= 295.164\ \text{mm}\\
\\[10pt]
Curva_{empírica} &= \mu \cdot \sqrt{ \frac{ Rlinha_{3_{3}} }{ h_{f_{3}} } } \\&= 0.1 \cdot \sqrt{ \frac{ 295.164\ \text{mm} }{ 5.000\ \text{mm} } } \\&= 0.768\\
\\[10pt]
\epsilon &= \frac{ \Delta_{h_{3}} }{ h_{i_{3}} } \\&= \frac{ 2.211\ \text{mm} }{ 7.211\ \text{mm} } \\&= 0.307\\
\end{aligned}
\]



```python
%%render
#long

f_1_Passo3_3 = 1.105 #aproximadamente

F_c_3_3 = b_m*((Rlinha_3_3*Delta_h_3)**0.5)*1.15*sigma_medio*f_1_Passo3_3
```


\[
\begin{aligned}
f_{1_{Passo3_{3}}} &= 1.105\;\;\textrm{(aproximadamente)}
\\[10pt]
F_{c_{3_{3}}} &= b_{m} \cdot \left( \left( Rlinha_{3_{3}} \cdot \Delta_{h_{3}} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot \sigma_{medio} \cdot f_{1_{Passo3_{3}}} \\&= 500.000\ \text{mm} \cdot \left( \left( 295.164\ \text{mm} \cdot 2.211\ \text{mm} \right) ^{ 0.5 } \right) \cdot 1.15 \cdot 82.878\ \text{MPa} \cdot 1.105 \\&= 1.345\ \text{MN}\\
\end{aligned}
\]



```python
%%render
#long

f_2_Passo3 = 0.115

M_c_3 = 2*R_3*(((h_i_3)**2)/h_f_3)*sigma_medio*b_m*f_2_Passo3
```


\[
\begin{aligned}
f_{2_{Passo3}} &= 0.115\;
\\[10pt]
M_{c_{3}} &= 2 \cdot R_{3} \cdot \left( \frac{ \left( h_{i_{3}} \right) ^{ 2 } }{ h_{f_{3}} } \right) \cdot \sigma_{medio} \cdot b_{m} \cdot f_{2_{Passo3}} \\&= 2 \cdot 279.170\ \text{mm} \cdot \left( \frac{ \left( 7.211\ \text{mm} \right) ^{ 2 } }{ 5.000\ \text{mm} } \right) \cdot 82.878\ \text{MPa} \cdot 500.000\ \text{mm} \cdot 0.115 \\&= 27.673\ \text{kN} \cdot \text{m}\\
\end{aligned}
\]


E finalmente a potência corrigida $P_{c_1}$ para o primeiro passo:


```python
%%render
#short

P_c_3 = M_c_3*omega_angular
```


\[
\begin{aligned}
P_{c_{3}} &= M_{c_{3}} \cdot \omega_{angular} = 27.673\ \text{kN} \cdot \text{m} \cdot 15.708\ \text{Hz} &= 434.686\ \text{kW}
\end{aligned}
\]


### Síntese

Como conclusão, plota-se uma tabela para comparação de resultados:


```python
tabela = pd.DataFrame(index=['Força','Momento','Potência'])

tabela['Etapa 1'] = [(F_i_1,F_c_1_3),(M_i_1,M_c_1),(P_i_1,P_c_1)]
tabela['Etapa 2'] = [(F_i_2,F_c_2_3),(M_i_2,M_c_2),(P_i_2,P_c_2)]
tabela['Etapa 3'] = [(F_i_3,F_c_3_3),(M_i_3,M_c_3),(P_i_3,P_c_3)]
```

Na tabela, serão expostos as grandezas de cada etapa em forma de pares. O valor da esquerda é referente ao processo ideal, enquanto o da direta é referente ao processo real: 


```python
tabela
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Etapa 1</th>
      <th>Etapa 2</th>
      <th>Etapa 3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Força</th>
      <td>(2.142 MN, 2.837 MN)</td>
      <td>(1.485 MN, 1.967 MN)</td>
      <td>(1.030 MN, 1.345 MN)</td>
    </tr>
    <tr>
      <th>Momento</th>
      <td>(110.682 kN·m, 119.734 kN·m)</td>
      <td>(53.210 kN·m, 56.060 kN·m)</td>
      <td>(25.581 kN·m, 27.673 kN·m)</td>
    </tr>
    <tr>
      <th>Potência</th>
      <td>(1.739 MW, 1.881 MW)</td>
      <td>(835.825 kW, 880.596 kW)</td>
      <td>(401.823 kW, 434.686 kW)</td>
    </tr>
  </tbody>
</table>
</div>



Nota-se que, como esperado, todas as grandezas reais requerem maior intensidade de energia do que os ideais.

Por último, plota-se outra tabela, agora com os valores sendo exibidos ao decorrer da correção. O primeiro valor é sempre ideal, enquanto os outros 3 são os corrigidos, lembrando que quanto mais iterado, mais corrigido - logo, quanto mais para a direita, mais corrigido está o valor.


```python
tabela2 = pd.DataFrame(index=['Raios', 'Forças'])

tabela2['Etapa 1'] = [(R_1,Rlinha_1_1,Rlinha_1_2,Rlinha_1_3), (F_i_1,F_c_1_1,F_c_1_2,F_c_1_3)]
tabela2['Etapa 2'] = [(R_2,Rlinha_2_1,Rlinha_2_2,Rlinha_2_3), (F_i_2,F_c_2_1,F_c_2_2,F_c_2_3)]
tabela2['Etapa 3'] = [(R_3,Rlinha_3_1,Rlinha_3_2,Rlinha_3_3), (F_i_3,F_c_3_1,F_c_3_2,F_c_3_3)]

tabela2
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Etapa 1</th>
      <th>Etapa 2</th>
      <th>Etapa 3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Raios</th>
      <td>(580.697 mm, 606.054 mm, 613.466 mm, 614.265 mm)</td>
      <td>(402.633 mm, 420.214 mm, 425.354 mm, 425.700 mm)</td>
      <td>(279.170 mm, 291.360 mm, 294.924 mm, 295.164 mm)</td>
    </tr>
    <tr>
      <th>Forças</th>
      <td>(2.142 MN, 2.768 MN, 2.835 MN, 2.837 MN)</td>
      <td>(1.485 MN, 1.919 MN, 1.948 MN, 1.967 MN)</td>
      <td>(1.030 MN, 1.331 MN, 1.351 MN, 1.345 MN)</td>
    </tr>
  </tbody>
</table>
</div>



---

## Questão II

Uma laminadora deverá realizar a deformação de uma barra de seção retangular de 5.0" x 5.0", entre cilindros de ferro fundido, com espaçamento entre eles de 60 mm, sendo a rotação dos cilindros da ordem de 80 RPM e diâmetro de 500 mm em uma gaiola com potência de motor de 700 cv e rendimento de $\eta_{mec}=0.85$.

1. Determine o comprimento final da barra sabendo que seu comprimento inicial é de 5 m; 
2. Será possível realizar o processo? 


O material será laminado a 1000 °C, e sua composição química é: 
- 0,40%C; 
- 0,87%Mn; 
- 0,95%Cr; 
- 0,25%Si; 
- 0,20%Mo; 
- 0,04%S; 
- 0,03%P. 

--- 

## Resposta II

### 1. Determinando o comprimento final da barra e outras grandezas importantes

O primeiro passo consta na determinação da deformação máxima na espessura que a chapa sofrerá


```python
%%render
#long

h_i = 5.0*25.4*mm #transformação de $inch$ para $mm$

h_f = (60*mm) #definido pelo enunciado

Delta_h = h_i - h_f

phi_h = ln(h_f/h_i)
```


\[
\begin{aligned}
h_{i} &= 5.0 \cdot 25.4 \cdot mm \\&= 5.0 \cdot 25.4 \cdot mm \\&= 127.000\ \text{mm}\;\;\textrm{(transformação de $inch$ para $mm$)}\\
\\[10pt]
h_{f} &= 60.000\ \text{mm}\;\;\textrm{(definido pelo enunciado)}
\\[10pt]
\Delta_{h} &= h_{i} - h_{f} \\&= 127.000\ \text{mm} - 60.000\ \text{mm} \\&= 67.000\ \text{mm}\\
\\[10pt]
\phi_{h} &= \operatorname{ln} \left( \frac{ h_{f} }{ h_{i} } \right) \\&= \operatorname{ln} \left( \frac{ 60.000\ \text{mm} }{ 127.000\ \text{mm} } \right) \\&= -0.75\\
\end{aligned}
\]


Como o processo é feito a quente, não se deve preocupar com deformação homogênea: a tensão de escoamento é constante.

Com tais informações, calcula-se o raio


```python
%%render
#parameters

Temp_processo = 1000*Celsius

C = 0.0005/Celsius
```


\[
\begin{aligned}
Temp_{processo} &= 1.000\ \text{k°C} &C &= 500.000\ \text{μ°C}^{-1}
\end{aligned}
\]



```python
%%render
#long

mu = 1.05-C*Temp_processo #específico para rolos de ferro fundido, de acordo com a apostila

R = ((500/2)*mm) #configuração da máquina

CritSeg = 1.25 #agindo com segurança

R_Minimo = CritSeg * Delta_h/((sin(arctan(mu)))**2)
```


\[
\begin{aligned}
\mu &= 1.05 - C \cdot Temp_{processo} \\&= 1.05 - 500.000\ \text{μ°C}^{-1} \cdot 1.000\ \text{k°C} \\&= 0.55\;\;\textrm{(específico para rolos de ferro fundido, de acordo com a apostila)}\\
\\[10pt]
R &= 250.000\ \text{mm}\;\;\textrm{(configuração da máquina)}
\\[10pt]
CritSeg &= 1.25\;\;\textrm{(agindo com segurança)}
\\[10pt]
R_{Minimo} &= CritSeg \cdot \frac{ \Delta_{h} }{ \left( \sin{ \left( \operatorname{arctan} \left( \mu \right) \right) } \right) ^{ 2 } } \\&= 1.25 \cdot \frac{ 67.000\ \text{mm} }{ \left( \sin{ \left( \operatorname{arctan} \left( 0.55 \right) \right) } \right) ^{ 2 } } \\&= 360.610\ \text{mm}\\
\end{aligned}
\]


Percebe-se que o raio atual da máquina não é adequado para a tarefa, uma vez que o raio mínimo que o processo requer é maior que este valor. Para que o projeto se torne acessível, o rolo deverá ser alterado.


```python
%%render
#short

R = R_Minimo
```


\[
\begin{aligned}
R &= 360.610\ \text{mm}\;
\end{aligned}
\]


Continuando, definindo alguns parâmetros de contato:


```python
%%render
#long

l_d = (R*Delta_h)**0.5
```


\[
\begin{aligned}
l_{d} &= \left( R \cdot \Delta_{h} \right) ^{ 0.5 } \\&= \left( 360.610\ \text{mm} \cdot 67.000\ \text{mm} \right) ^{ 0.5 } \\&= 155.438\ \text{mm}\\
\end{aligned}
\]


Como o processo é feito a quente, a deformação verdadeira $\phi_b$ na espessura **não** deverá ser desprezada.


```python
%%render
#long

b_i = 5*25.4*mm #transformação de $inch$ para $mm$

Temp_ref = 1000 #K
Temp_proc = (1000 + 273) #K

C_bmu = Temp_ref/Temp_proc
```


\[
\begin{aligned}
b_{i} &= 5 \cdot 25.4 \cdot mm \\&= 5 \cdot 25.4 \cdot mm \\&= 127.000\ \text{mm}\;\;\textrm{(transformação de $inch$ para $mm$)}\\
\\[10pt]
Temp_{ref} &= 1000\;\;\textrm{(K)}
\\[10pt]
Temp_{proc} &= 1273\;\;\textrm{(K)}
\\[10pt]
C_{b\mu} &= \frac{ Temp_{ref} }{ Temp_{proc} } \\&= \frac{ 1000 }{ 1273 } \\&= 0.786\\
\end{aligned}
\]



```python
%%render
#long

phi_b = -phi_h*e**(-C_bmu*(b_i/l_d))

b_f = b_i*e**phi_b #variação bastante considerável, como se percebe com o cálculo
```


\[
\begin{aligned}
\phi_{b} &= - \phi_{h} \cdot \left( e \right) ^{ \left( - C_{b\mu} \cdot \left( \frac{ b_{i} }{ l_{d} } \right) \right) } \\&= - -0.75 \cdot \left( 2.718 \right) ^{ \left( - 0.786 \cdot \left( \frac{ 127.000\ \text{mm} }{ 155.438\ \text{mm} } \right) \right) } \\&= 0.395\\
\\[10pt]
b_{f} &= b_{i} \cdot \left( e \right) ^{ \phi_{b} } \\&= 127.000\ \text{mm} \cdot \left( 2.718 \right) ^{ 0.395 } \\&= 188.453\ \text{mm}\;\;\textrm{(variação bastante considerável, como se percebe com o cálculo)}\\
\end{aligned}
\]


O importante, neste caso, é conhecer a variação média na espessura, ou seja, entre seu valor de entrada e de saída:


```python
%%render
#short

b_m = (b_i+b_f)/2
```


\[
\begin{aligned}
b_{m} &= \frac{ b_{i} + b_{f} }{ 2 } = \frac{ 127.000\ \text{mm} + 188.453\ \text{mm} }{ 2 } &= 157.727\ \text{mm}
\end{aligned}
\]


Com tal valor, torna-se possível a validação da área de contato $A_c$:


```python
%%render
#short

A_c = b_m*l_d
```


\[
\begin{aligned}
A_{c} &= b_{m} \cdot l_{d} = 157.727\ \text{mm} \cdot 155.438\ \text{mm} &= 24516.656\ \text{mm}^{2.0}
\end{aligned}
\]


De acordo com a conservação de massa, é possível obter o valor da deformação no comprimento $\phi_l$ com a deformação na largura e na espessura:


```python
%%render 
#long 

phi_l = - phi_b - phi_h

Conservação_massa = phi_l + phi_b + phi_h #a variação de massa ou volume sempre deve ser 0 entre os processos de conformação 
```


\[
\begin{aligned}
\phi_{l} &= - \phi_{b} - \phi_{h} \\&= - 0.395 - -0.75 \\&= 0.355\\
\\[10pt]
Conservação_{massa} &= \phi_{l} + \phi_{b} + \phi_{h} \\&= 0.355 + 0.395 + -0.75 \\&= 0.0\;\;\textrm{(a variação de massa ou volume sempre deve ser 0 entre os processos de conformação)}\\
\end{aligned}
\]


E com todas as deformação verdadeiras definidas, fica fácil determinar o comprimento final da peça: 


```python
%%render
#short

l_i = (5*m)

l_f = l_i*e**phi_l
```


\[
\begin{aligned}
l_{i} &= 5.000\ \text{m}\;
\\[10pt]
l_{f} &= l_{i} \cdot \left( e \right) ^{ \phi_{l} } = 5.000\ \text{m} \cdot \left( 2.718 \right) ^{ 0.355 } &= 7.132\ \text{m}
\end{aligned}
\]


A tensão de escoamento $\sigma_e$ é constante para o processo a quente e segue uma equação empírica guiada por sua composição:


```python
%%render
#long

C_perc = 0.40 #%

Mn_perc = 0.87 #%

Cr_perc = 0.95 #%

C_2 = (0.01*(1/Celsius))

sigma_e_kgfmm2 = (14 - C_2*Temp_processo)*(1.4+C_perc+Mn_perc+0.3*Cr_perc) #kgf/mm²

sigma_e = sigma_e_kgfmm2*9.81*MPa #transformado para MPa para manter a consistência dimensional
```


\[
\begin{aligned}
C_{perc} &= 0.4\;\;\textrm{(%)}
\\[10pt]
Mn_{perc} &= 0.87\;\;\textrm{(%)}
\\[10pt]
Cr_{perc} &= 0.95\;\;\textrm{(%)}
\\[10pt]
C_{2} &= 10.000\ \text{m°C}^{-1}\;
\\[10pt]
\sigma_{e_{kgfmm2}} &= \left( 14 - C_{2} \cdot Temp_{processo} \right) \cdot \left( 1.4 + C_{perc} + Mn_{perc} + 0.3 \cdot Cr_{perc} \right) \\&= \left( 14 - 10.000\ \text{m°C}^{-1} \cdot 1.000\ \text{k°C} \right) \cdot \left( 1.4 + 0.4 + 0.87 + 0.3 \cdot 0.95 \right) \\&= 11.82\;\;\textrm{(kgf/mm²)}\\
\\[10pt]
\sigma_{e} &= \sigma_{e_{kgfmm2}} \cdot 9.81 \cdot MPa \\&= 11.82 \cdot 9.81 \cdot MPa \\&= 115.954\ \text{MPa}\;\;\textrm{(transformado para MPa para manter a consistência dimensional)}\\
\end{aligned}
\]


A força $F$, na laminação a quente, é calculada pela área de contato e pela pressão específica $K_w$.

$K_w$ depende de dois novos fatores, a velocidade de deformação $\dot{\phi}$ e o coeficiente de plasticidade $\eta$:


```python
%%render
#long

eta = 0.01*(14-C_2*Temp_processo)

omega_angular = 80*(2*pi/60)*1/s
```


\[
\begin{aligned}
\eta &= 0.01 \cdot \left( 14 - C_{2} \cdot Temp_{processo} \right) \\&= 0.01 \cdot \left( 14 - 10.000\ \text{m°C}^{-1} \cdot 1.000\ \text{k°C} \right) \\&= 0.04\\
\\[10pt]
\omega_{angular} &= 80 \cdot \left( 2 \cdot \frac{ \pi }{ 60 } \right) \cdot \frac{ 1 }{ s } \\&= 80 \cdot \left( 2 \cdot \frac{ 3.142 }{ 60 } \right) \cdot \frac{ 1 }{ s } \\&= 8.378\ \text{Hz}\\
\end{aligned}
\]



```python
%%render
#long

V_linear = omega_angular*R

phi_dot = (2*V_linear*(Delta_h/R)**0.5)/(h_i+h_f)
```


\[
\begin{aligned}
V_{linear} &= \omega_{angular} \cdot R \\&= 8.378\ \text{Hz} \cdot 360.610\ \text{mm} \\&= 3.021\ \text{m} \cdot \text{s}^{-1}\\
\\[10pt]
\phi_{dot} &= \frac{ 2 \cdot V_{linear} \cdot \left( \frac{ \Delta_{h} }{ R } \right) ^{ 0.5 } }{ h_{i} + h_{f} } \\&= \frac{ 2 \cdot 3.021\ \text{m} \cdot \text{s}^{-1} \cdot \left( \frac{ 67.000\ \text{mm} }{ 360.610\ \text{mm} } \right) ^{ 0.5 } }{ 127.000\ \text{mm} + 60.000\ \text{mm} } \\&= 13.927\ \text{Hz}\\
\end{aligned}
\]


Finalmente, $K_w$:


```python
%%render
#long

K_w = (1+(1.6*mu*l_d-1.2*Delta_h)/(h_i+h_f))*(sigma_e_kgfmm2+eta*phi_dot*s)*(9.81*MPa) #cálculo é empírico e envolve kgf/mm², alguma complicação no dimensional acontece
```


\[
\begin{aligned}
K_{w} &= \left( 1 + \frac{ 1.6 \cdot \mu \cdot l_{d} - 1.2 \cdot \Delta_{h} }{ h_{i} + h_{f} } \right) \cdot \left( \sigma_{e_{kgfmm2}} + \eta \cdot \phi_{dot} \cdot s \right) \cdot \left( 9.81 \cdot MPa \right) \\&= \left( 1 + \frac{ 1.6 \cdot 0.55 \cdot 155.438\ \text{mm} - 1.2 \cdot 67.000\ \text{mm} }{ 127.000\ \text{mm} + 60.000\ \text{mm} } \right) \cdot \left( 11.82 + 0.04 \cdot 13.927\ \text{Hz} \cdot s \right) \cdot \left( 9.81 \cdot MPa \right) \\&= 158.030\ \text{MPa}\;\;\textrm{(cálculo é empírico e envolve kgf/mm², alguma complicação no dimensional acontece)}\\
\end{aligned}
\]


Perto do final, $F$:


```python
%%render
#long

F = K_w*A_c #valor elevado, uma vez que a área de contato é relativamente grande
```


\[
\begin{aligned}
F &= K_{w} \cdot A_{c} \\&= 158.030\ \text{MPa} \cdot 24516.656\ \text{mm}^{2.0} \\&= 3.874\ \text{MN}\;\;\textrm{(valor elevado, uma vez que a área de contato é relativamente grande)}\\
\end{aligned}
\]


O cálculo do momento $M$ deve ser feito com a força $F$ e $l_d$:


```python
%%render
#long

M = F * l_d
```


\[
\begin{aligned}
M &= F \cdot l_{d} \\&= 3.874\ \text{MN} \cdot 155.438\ \text{mm} \\&= 602.222\ \text{kN} \cdot \text{m}\\
\end{aligned}
\]


### 2. Determinando se o processo é viável (Síntese)

Enfim, a potência $\dot W$:


```python
%%render
#long

eta_mecanica = 0.85

W_dot = (M * omega_angular)/eta_mecanica #ou 8100 cv
```


\[
\begin{aligned}
\eta_{mecanica} &= 0.85\;
\\[10pt]
W_{dot} &= \frac{ M \cdot \omega_{angular} }{ \eta_{mecanica} } \\&= \frac{ 602.222\ \text{kN} \cdot \text{m} \cdot 8.378\ \text{Hz} }{ 0.85 } \\&= 5.935\ \text{MW}\;\;\textrm{(ou 8100 cv)}\\
\end{aligned}
\]


Como se percebe, é um valor muito alto para ser praticado. Se deve ao tamanho da área de contato e também à variação de espessura em um único passe (~67mm) muito alta.

Para resolver tal problema, sugere-se a divisão do processo em mais de uma gaiola.
