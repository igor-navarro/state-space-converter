"""
exemplos.py

Cinco casos de validacao do conversor, cobrindo: um sistema de 1a ordem,
um de 2a ordem, um sistema construido deliberadamente para apresentar
cancelamento polo-zero, um sistema MIMO (multiplas entradas/saidas) e um
sistema de 3a ordem em forma canonica de Jordan (autovalor duplo) — um
bom teste de estresse para a simplificacao simbolica, por envolver um
fator quadratico no denominador.

Uso:
    python exemplos.py
"""

import sympy as sp

from conversor import (
    analisar_kalman,
    calcular_funcao_transferencia,
    detectar_polos_zeros_e_cancelamentos,
    formatar_saida_ascii,
)


CASOS = [
    {
        "nome": "Circuito RC (1a ordem)",
        "A": [[-1]], "B": [[1]], "C": [[1]], "D": [[0]],
    },
    {
        "nome": "Sistema de 2a ordem",
        "A": [[0, 1], [-2, -3]], "B": [[0], [1]], "C": [[1, 0]], "D": [[0]],
    },
    {
        "nome": "Cancelamento polo-zero (construido deliberadamente)",
        "A": [[0, 1], [-2, -3]], "B": [[0], [1]], "C": [[1, 1]], "D": [[0]],
    },
    {
        "nome": "Sistema MIMO (2 entradas / 2 saidas)",
        "A": [[-1, 0], [0, -2]], "B": [[1, 0], [0, 1]],
        "C": [[1, 1], [0, 1]], "D": [[0, 0], [0, 0]],
    },
    {
        "nome": "Forma de Jordan, 3a ordem, polo duplo em s=-2",
        "A": [[-2, 1, 0], [0, -2, 0], [0, 0, -3]],
        "B": [[0], [1], [1]], "C": [[1, 0, 1]], "D": [[0]],
    },
]


def rodar_caso(caso):
    A = sp.Matrix(caso["A"])
    B = sp.Matrix(caso["B"])
    C = sp.Matrix(caso["C"])
    D = sp.Matrix(caso["D"])
    n, p, q = A.rows, B.cols, C.rows

    print("\n" + "=" * 70)
    print(f" {caso['nome']}")
    print("=" * 70)

    kalman = analisar_kalman(A, B, C)
    print(f"Ordem n={n} | Controlavel: {kalman['controlavel']} "
          f"(posto {kalman['posto_Co']}/{n}) | "
          f"Observavel: {kalman['observavel']} (posto {kalman['posto_Ob']}/{n})")

    G, s = calcular_funcao_transferencia(A, B, C, D)
    G_bruto = C * (s * sp.eye(n) - A).inv() * B + D

    for i in range(q):
        for j in range(p):
            info = detectar_polos_zeros_e_cancelamentos(G_bruto[i, j], G[i, j], s)
            print(f"\n  Canal y{i+1} <- u{j+1}:")
            for linha in formatar_saida_ascii(info["num"], info["den"]).split("\n"):
                print(f"  {linha}")
            print(f"  Polos: {info['polos']}  Zeros: {info['zeros']}")
            if info["cancelados"]:
                print(f"  >> cancelamento detectado em s = {info['cancelados']}")


if __name__ == "__main__":
    for caso in CASOS:
        rodar_caso(caso)
