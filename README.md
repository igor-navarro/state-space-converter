# Conversor de Espaço de Estados para Função de Transferência

Um conversor simbólico, escrito em Python, que transforma a representação de
um sistema dinâmico linear no espaço de estados em sua função de
transferência no domínio da frequência:

```
G(s) = C (sI - A)⁻¹ B + D
```

Além da conversão em si, o programa verifica a **controlabilidade** e a
**observabilidade** do sistema (critério de Kalman) e detecta
**cancelamentos polo-zero** que podem ocorrer durante a simplificação
algébrica — situações em que a função de transferência final deixa de
revelar um modo interno do sistema.

## Por que isso é interessante

A conversão parece trivial à primeira vista, mas o cálculo manual de
`(sI - A)⁻¹` se torna rapidamente inviável para sistemas de 3ª ordem ou
mais — é exatamente o tipo de tarefa mecânica e algebricamente exigente que
compensa resolver com computação simbólica em vez de fazer na mão.

Ao validar o programa com um sistema clássico em **forma canônica de
Jordan** (autovalor duplo), usado como um dos benchmarks de teste, o
resultado calculado pelo programa divergiu de um valor de referência que eu
estava usando para conferência. Em vez de simplesmente confiar no valor de
referência, refiz a conta manualmente, por substituição direta nas
matrizes originais, e confirmei que o **programa estava correto** — a
referência é que continha um erro no numerador. Esse episódio acabou virando
o motivo pelo qual o projeto passou a comparar sistematicamente os polos
*antes* e *depois* da simplificação simbólica (veja `detectar_polos_zeros_e_cancelamentos`
em `conversor.py`), em vez de simplesmente confiar cegamente em um único
caso de referência.

## Funcionalidades

- Leitura interativa das matrizes `A`, `B`, `C`, `D`, com validação de
  dimensões em camadas (`A` precisa ser quadrada; `B`, `C`, `D` precisam
  ser compatíveis com `A`) — o programa não avança enquanto a entrada não
  for válida, sem valores padrão escondidos.
- Leitura segura da entrada com `ast.literal_eval()` em vez de `eval()`,
  evitando a execução de código arbitrário a partir de texto digitado
  pelo usuário.
- Conversão simbólica exata via [SymPy](https://www.sympy.org/), com
  coeficientes racionais/inteiros em vez de aproximações numéricas.
- Verificação cruzada de controlabilidade/observabilidade: posto
  calculado simbolicamente (SymPy) e numericamente (NumPy), como duas
  checagens independentes do mesmo resultado.
- Detecção de cancelamento polo-zero por comparação estrutural (polos da
  expressão bruta vs. polos após a simplificação), não apenas pela
  aparência do resultado final.
- Suporte nativo a sistemas MIMO (múltiplas entradas e saídas), analisando
  cada canal entrada-saída individualmente.
- Testado contra 5 casos de validação, incluindo um sistema de 4ª ordem
  fora do conjunto de exemplos padrão (ver seção de testes abaixo).

## Estrutura

```
.
├── conversor.py   # biblioteca principal (toda a lógica)
├── main.py        # interface de linha de comando interativa
├── exemplos.py     # 5 casos de teste, executáveis sem interação
└── requirements.txt
```

## Como usar

```bash
pip install -r requirements.txt

# Modo interativo — digite suas próprias matrizes
python main.py

# Rodar os 5 casos de teste de uma vez
python exemplos.py
```

Exemplo de entrada no modo interativo (sistema de 2ª ordem):

```
A = [[0, 1], [-2, -3]]
B = [[0], [1]]
C = [[1, 0]]
D = [[0]]
```

Saída esperada:

```
Funcao de Transferencia:

      1
--------------
 s^2 + 3s + 2

Polos: [-2.0, -1.0]
Zeros: []
Nenhum cancelamento polo-zero detectado.
```

## Casos de teste

| Caso | Sistema | G(s) | Controlável | Observável | Cancelamento |
|---|---|---|:---:|:---:|:---:|
| 1 | Circuito RC (1ª ordem) | `1/(s+1)` | ✅ | ✅ | Não |
| 2 | Sistema de 2ª ordem | `1/(s²+3s+2)` | ✅ | ✅ | Não |
| 3 | Cancelamento polo-zero (proposital) | `1/(s+2)` | ✅ | ❌ | Sim, em `s=-1` |
| 4 | MIMO (2 entradas / 2 saídas) | matriz 2×2 de transferências | ✅ | ✅ | Não |
| 5 | Forma de Jordan, polo duplo | `(s²+5s+7)/((s+2)²(s+3))` | ✅ | ✅ | Não |

O Caso 3 foi construído deliberadamente — a partir do Caso 2, alterando a
matriz `C` — para forçar um cancelamento e confirmar que a detecção
funciona na prática, não só na teoria.

## Tecnologias

- Python 3
- [SymPy](https://www.sympy.org/) — álgebra simbólica
- [NumPy](https://numpy.org/) — verificação numérica auxiliar

## Referências

- OGATA, Katsuhiko. *Modern Control Engineering*. 5. ed. Prentice Hall, 2010.
- DORF, Richard C.; BISHOP, Robert H. *Modern Control Systems*. 13. ed. Pearson, 2016.
- MEURER, Aaron et al. SymPy: symbolic mathematics in Python. *PeerJ Computer Science*, v. 3, e103, 2017.
