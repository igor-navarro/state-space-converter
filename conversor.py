"""
conversor.py

Conversao simbolica de sistemas lineares invariantes no tempo, da
representacao em espaco de estados para a funcao de transferencia no
dominio da frequencia:

    G(s) = C (sI - A)^-1 B + D

O modulo tambem verifica controlabilidade e observabilidade (criterio de
Kalman) e detecta cancelamentos polo-zero que possam ocorrer durante a
simplificacao algebrica.

Dependencias: sympy, numpy
"""

import ast
import sympy as sp
import numpy as np


class MatrizInvalidaError(Exception):
    """Levantada quando a entrada nao pode ser interpretada como uma
    matriz numerica valida."""
    pass


# ---------------------------------------------------------------------
# Leitura e validacao de matrizes
# ---------------------------------------------------------------------

def interpretar_matriz(entrada: str):
    """Converte uma string no formato de lista Python (ex.: "[[0, 1], [-2, -3]]")
    em uma lista de listas numerica. Usa ast.literal_eval em vez de eval()
    para nao correr o risco de executar codigo arbitrario a partir da
    entrada do usuario, e valida que todo elemento e realmente um numero.
    """
    try:
        lista = ast.literal_eval(entrada)
        for linha in lista:
            for valor in linha:
                if not isinstance(valor, (int, float)):
                    raise ValueError(f"valor nao numerico encontrado: {valor!r}")
        return lista
    except Exception as exc:
        raise MatrizInvalidaError(str(exc)) from exc


def ler_matriz_usuario(nome_matriz: str, exemplo: str, validador=None):
    """Pede repetidamente uma matriz ao usuario ate que uma entrada
    numericamente valida (e estruturalmente valida, se um `validador`
    for passado) seja fornecida. Nao ha valor padrao: a entrada e
    obrigatoria.

    validador: funcao opcional que recebe a sp.Matrix ja convertida e
    devolve (True, None) se ela for aceitavel, ou (False, "motivo")
    caso contrario — usada para checagens estruturais como "deve ser
    quadrada" ou "deve ter N linhas".
    """
    while True:
        entrada = input(
            f"Digite a matriz {nome_matriz} no formato de lista Python "
            f"(ex: {exemplo}): "
        ).strip()

        if not entrada:
            print("[-] Entrada obrigatoria. Por favor, digite a matriz.")
            continue

        try:
            lista = interpretar_matriz(entrada)
        except MatrizInvalidaError:
            print("[-] Erro: entrada invalida ou nao numerica. Digite novamente.")
            continue

        matriz = sp.Matrix(lista)

        if validador is not None:
            ok, motivo = validador(matriz)
            if not ok:
                print(f"[-] Erro de dimensao: {motivo}. Digite {nome_matriz} novamente.")
                continue

        return matriz


def ler_sistema_interativo():
    """Le A, B, C e D do usuario, validando dimensoes de forma cumulativa:
    A precisa ser quadrada; B precisa ter as mesmas linhas que A; C
    precisa ter as mesmas colunas que A; D precisa ter dimensao
    (linhas de C) x (colunas de B).
    """
    def validar_A(m):
        if not m.is_square:
            return False, f"a matriz A deve ser quadrada (recebido {m.rows}x{m.cols})"
        return True, None

    A = ler_matriz_usuario("A", "[[0, 1], [-2, -3]]", validador=validar_A)
    n = A.rows

    def validar_B(m):
        if m.rows != n:
            return False, f"B deve ter {n} linha(s), o mesmo numero de linhas que A"
        return True, None

    B = ler_matriz_usuario("B", "[[0], [1]]", validador=validar_B)
    p = B.cols

    def validar_C(m):
        if m.cols != n:
            return False, f"C deve ter {n} coluna(s), o mesmo numero de colunas que A"
        return True, None

    C = ler_matriz_usuario("C", "[[1, 0]]", validador=validar_C)
    q = C.rows

    def validar_D(m):
        if m.rows != q or m.cols != p:
            return False, f"D deve ter dimensao {q}x{p} (saidas x entradas)"
        return True, None

    D = ler_matriz_usuario("D", "[[0]]", validador=validar_D)

    return A, B, C, D


# ---------------------------------------------------------------------
# Analise de controlabilidade e observabilidade (criterio de Kalman)
# ---------------------------------------------------------------------

def matriz_controlabilidade(A, B):
    n = A.rows
    return sp.Matrix.hstack(*[A**i * B for i in range(n)])


def matriz_observabilidade(A, C):
    n = A.rows
    return sp.Matrix.vstack(*[C * A**i for i in range(n)])


def analisar_kalman(A, B, C):
    """Retorna um dicionario com as matrizes de Kalman, seus postos e se
    o sistema e completamente controlavel/observavel. Inclui tambem uma
    checagem numerica cruzada via NumPy (posto calculado a partir da
    versao numerica das matrizes), como segunda verificacao independente
    do resultado simbolico do SymPy.
    """
    n = A.rows
    Co = matriz_controlabilidade(A, B)
    Ob = matriz_observabilidade(A, C)

    co_rank = Co.rank()
    ob_rank = Ob.rank()

    Co_np = np.array(Co.evalf(), dtype=float)
    Ob_np = np.array(Ob.evalf(), dtype=float)

    return {
        "Co": Co,
        "Ob": Ob,
        "posto_Co": co_rank,
        "posto_Ob": ob_rank,
        "controlavel": co_rank == n,
        "observavel": ob_rank == n,
        "posto_Co_numpy": int(np.linalg.matrix_rank(Co_np)),
        "posto_Ob_numpy": int(np.linalg.matrix_rank(Ob_np)),
    }


# ---------------------------------------------------------------------
# Conversao para funcao de transferencia
# ---------------------------------------------------------------------

def calcular_funcao_transferencia(A, B, C, D):
    """Calcula G(s) = C (sI - A)^-1 B + D, ja simplificada."""
    s = sp.Symbol("s")
    n = A.rows
    G_bruto = C * (s * sp.eye(n) - A).inv() * B + D
    return sp.simplify(G_bruto), s


def extrair_numerador_denominador(elemento):
    """Isola numerador e denominador de uma expressao racional. Usa
    sympy.together() antes de sympy.fraction() — passo necessario porque,
    para alguns sistemas (percebido em sistemas de ordem mais alta),
    sympy.simplify() pode devolver uma soma de fracoes parciais em vez
    de uma unica fracao, o que faria sympy.fraction() sozinho extrair um
    denominador incorreto (igual a 1).
    """
    return sp.fraction(sp.together(elemento))


def _numero_para_exibicao(v):
    """Converte uma raiz simbolica para float quando possivel, ou complex
    quando a raiz nao e real (evita o TypeError do float() com numeros
    complexos)."""
    if not v.is_number:
        return v
    v_ev = v.evalf()
    try:
        return float(v_ev)
    except TypeError:
        return complex(v_ev)


def detectar_polos_zeros_e_cancelamentos(G_bruto_elemento, elemento_simplificado, s):
    """Compara os polos da expressao bruta (antes da simplificacao) com
    os polos remanescentes apos a simplificacao. A diferenca entre os
    dois conjuntos e exatamente o que foi cancelado durante a
    simplificacao algebrica — nao bastaria comparar numerador e
    denominador ja simplificados, pois isso so confirmaria que o
    resultado final esta em forma minima, sem indicar se um
    cancelamento de fato ocorreu.
    """
    num, den = extrair_numerador_denominador(elemento_simplificado)
    polos = sp.solve(den, s)
    zeros = sp.solve(num, s)

    _, den_bruto = sp.fraction(sp.together(G_bruto_elemento))
    polos_brutos = sp.solve(den_bruto, s)

    cancelados = [p for p in polos_brutos if p not in polos]

    return {
        "num": num,
        "den": den,
        "polos": [_numero_para_exibicao(p) for p in polos],
        "zeros": [_numero_para_exibicao(z) for z in zeros],
        "cancelados": [_numero_para_exibicao(c) for c in cancelados],
    }


def formatar_saida_ascii(num_expr, den_expr):
    """Formata numerador/denominador como texto simples, numerador sobre
    uma linha de tracos, denominador logo abaixo — util para logs e
    terminais sem suporte a LaTeX."""
    num_str = str(sp.expand(num_expr)).replace("**", "^").replace("*", "")
    den_str = str(sp.expand(den_expr)).replace("**", "^").replace("*", "")
    largura = max(len(num_str), len(den_str)) + 2
    linhas = [
        "Funcao de Transferencia:",
        "",
        num_str.center(largura),
        "-" * largura,
        den_str.center(largura),
    ]
    return "\n".join(linhas)
