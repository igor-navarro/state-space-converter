"""
main.py

Interface de linha de comando: le um sistema em espaco de estados do
usuario, analisa controlabilidade/observabilidade e imprime a funcao de
transferencia correspondente, ja simplificada, com deteccao de
cancelamento polo-zero.

Uso:
    python main.py
"""

import sympy as sp

from conversor import (
    ler_sistema_interativo,
    analisar_kalman,
    calcular_funcao_transferencia,
    detectar_polos_zeros_e_cancelamentos,
    formatar_saida_ascii,
)


def main():
    print("=" * 60)
    print(" CONVERSOR: ESPACO DE ESTADOS -> FUNCAO DE TRANSFERENCIA")
    print("=" * 60)

    A, B, C, D = ler_sistema_interativo()
    n, p, q = A.rows, B.cols, C.rows

    print(f"\n[+] Sistema carregado: ordem n={n}, entradas p={p}, saidas q={q}")
    print(f"    A = {A.tolist()}")
    print(f"    B = {B.tolist()}")
    print(f"    C = {C.tolist()}")
    print(f"    D = {D.tolist()}")

    kalman = analisar_kalman(A, B, C)
    print("\n" + "-" * 50)
    print(" CONTROLABILIDADE E OBSERVABILIDADE (KALMAN)")
    print("-" * 50)
    print(f"Posto(Co) = {kalman['posto_Co']} de {n} ->",
          "CONTROLAVEL" if kalman["controlavel"] else "NAO CONTROLAVEL")
    print(f"Posto(Ob) = {kalman['posto_Ob']} de {n} ->",
          "OBSERVAVEL" if kalman["observavel"] else "NAO OBSERVAVEL")
    print(f"[checagem cruzada NumPy] posto(Co)={kalman['posto_Co_numpy']}, "
          f"posto(Ob)={kalman['posto_Ob_numpy']}")

    G, s = calcular_funcao_transferencia(A, B, C, D)
    G_bruto = C * (s * sp.eye(n) - A).inv() * B + D

    print("\n" + "-" * 50)
    print(" FUNCAO DE TRANSFERENCIA")
    print("-" * 50)

    for i in range(q):
        for j in range(p):
            info = detectar_polos_zeros_e_cancelamentos(G_bruto[i, j], G[i, j], s)
            print(f"\nCanal saida y{i+1} <- entrada u{j+1}:")
            print(formatar_saida_ascii(info["num"], info["den"]))
            print(f"  Polos: {info['polos']}")
            print(f"  Zeros: {info['zeros']}")
            if info["cancelados"]:
                print(f"  ALERTA: cancelamento polo-zero detectado em s = {info['cancelados']}")
            else:
                print("  Nenhum cancelamento polo-zero detectado.")


if __name__ == "__main__":
    main()
