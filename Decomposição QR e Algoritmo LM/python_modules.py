# %%
import numpy as np
import graphs as gr

# %%
def QRdecomposition_Gram_Schmidt(A):
    """
    entrada:
        A - matriz a decompor
    saída:
        Q - matriz ortogonal
        R - matriz triangular superior
    """
    q = [] #sera a lista de colunas de Q
    R = np.zeros((len(A[0]), len(A[0])), dtype=float) # matriz R, preenchida de zeros

    # percorrer as colunas da matriz A
    for i in range(len(np.transpose(A))):
        ai = np.transpose(A)[i] #ai = coluna i de A

        projecoes = 0 #a subtrair no calculo de vi
        # percorrer as colunas Q ja calculadas
        for n in range(len(q)):
            qn = q[n] #qn = coluna n de Q, ja calculada, pelo que n < i
            R[n][i] = np.inner(qn,ai) #calcular r_{n,i}
            projecoes += R[n][i]*qn #somar as projecoes a subtrair a ai
        
        # calcular vi
        vi = ai - projecoes

        # calcular r_{i,i}: norma 2 de vi
        R[i][i] = np.sqrt(np.inner(vi,vi))

        # calcular qi: normalizacao de vi
        qi = vi/R[i][i]
        q.append(qi)

    # converter q numa matriz (Q') e transpor
    Q = np.transpose(np.vstack(q))

    return Q,R

# %%
def QR_exemplo():
    x = list(range(19 + 1))
    y = [np.log(value) for value in [175.5, 188.1, 200.2, 213.7, 231.6, 254.7, 283.1,
                                    313.5, 342.9, 363.7, 389.1, 422.2, 481.1, 518.8,
                                    555.1, 591.6, 643.1, 693.0, 767.3, 846.1]]
    # Ax = b; x = [ln(b0), b1]
    A = np.column_stack((np.ones(len(x)), x))
    b = np.array(y).reshape(-1, 1)
    Q, R = QRdecomposition_Gram_Schmidt(A)  # QRx = b
    Rx = np.dot(np.transpose(Q), b)  # Rx = Qt*b

    x = np.zeros(len(R))
    for i in range(len(R) - 1, -1, -1):
        x[i] = (Rx[i, 0] - np.dot(R[i, i + 1:], x[i + 1:])) / R[i, i]

    b0, b1 = np.exp(x[0]), x[1]
    return b0, b1

# %%
def realidadeVSQR():
    x = list(range(19+1))
    y = [175.5, 188.1, 200.2, 213.7, 231.6, 254.7, 283.1, 313.5, 342.9, 363.7, 389.1, 422.2, 481.1, 518.8, 555.1, 591.6, 643.1, 693.0, 767.3, 846.1]
    b0,b1 = 171.00724094657875, 0.08356097087277026
    col_pontos = (84/255, 168/255, 217/255)
    col_linha = (50/255, 71/255, 143/255)

    gr.expetativaVSrealidade(x,y,b0,b1,
                        "realidadevsQR.png", 12, 5, "Modelo ajustado (QR)",
                        col_linha, col_pontos)

# %%
def LM_iteration(X, Y, ig, miu):
    """
    entrada:
        X,Y - conjunto de pontos
        ig, miu - aproximacao inicial (b0, b1), parametro mu
    saida:
        b0, b1 - coeficientes atualizados
    """
    matrix_r = np.zeros((len(X), 1))  # matriz dos residuos (R)
    grad_r = np.zeros((len(X), 2))  # gradiente de R

    b0, b1 = ig
    # atualizar matrizes
    for row in range(len(X)):
        matrix_r[row, 0] = Y[row] - b0 * np.exp(b1 * X[row])

        grad_r[row, 0] = -np.exp(b1 * X[row])  # r gradiente: col. 1
        grad_r[row, 1] = -b0 * X[row]* np.exp(b1 * X[row])  # r gradiente: col. 2

    # calcular A(s_LM) = b
    A = np.dot(grad_r.T, grad_r) + np.eye(2) * miu
    b = -np.dot(grad_r.T, matrix_r)

    # resolver s_LM atraves Cholesky
    L = np.linalg.cholesky(A)
    s_LM_y = np.linalg.solve(L, b)
    s_LM = np.linalg.solve(L.T, s_LM_y)

    # atualizar betas
    b0, b1 = float(b0 + s_LM[0]), float(b1 + s_LM[1])

    return b0, b1, s_LM

# %%
def LM_exemplo():
    import warnings
    warnings.filterwarnings("ignore")
    x = list(range(19 + 1))
    y = [175.5, 188.1, 200.2, 213.7, 231.6, 254.7, 283.1, 313.5, 342.9, 363.7, 389.1, 422.2,
    481.1, 518.8, 555.1, 591.6, 643.1, 693.0, 767.3, 846.1]

    ig = (170, 0.0836) #obtido atraves do exemplo da seccao 2.3.3
    miu = 0.1

    iteracoes = 1
    while True:
        b0, b1, s_LM = LM_iteration(x, y, ig, miu)

        # condicao de paragem
        if max(abs(s_LM)) < 0.5*10**(-6):
                break
    
        # atualizar parametros
        ig = (b0, b1)
        miu = 0.1
        iteracoes += 1

    return iteracoes, b0, b1

# %%
def realidadeVSLM():
    x = list(range(19+1))
    y = [175.5, 188.1, 200.2, 213.7, 231.6, 254.7, 283.1, 313.5, 342.9, 363.7, 389.1, 422.2, 481.1, 518.8, 555.1, 591.6, 643.1, 693.0, 767.3, 846.1]
    b0,b1 = 171.5072875301047, 0.08333471235500428
    col_pontos = (84/255, 168/255, 217/255)
    col_linha = (50/255, 71/255, 143/255)

    gr.expetativaVSrealidade(x,y,b0,b1,
                          "realidadevsLM.png", 12, 5, "Modelo ajustado (LM)",
                          col_linha, col_pontos)

# %%
def miu_grafico(width = 9, height = 4):
    x = list(range(19 + 1))
    y = [175.5, 188.1, 200.2, 213.7, 231.6, 254.7, 283.1, 313.5, 342.9, 363.7, 389.1, 422.2,
        481.1, 518.8, 555.1, 591.6, 643.1, 693.0, 767.3, 846.1]
    
    ig = (168,0.0825)
    miu = [100, 10**(-2)]

    gr.LM_miu(x, y, ig, miu, width, height)

# %%
def funcao_objetivo(y,x,beta):
    r_i = y - beta[0]*np.exp(beta[1]*x)
    valor_otimo = 0.5 * np.sum(r_i**2)
    return valor_otimo

# %%
def funcao_obj_conclusao():
    x = np.array(list(range(19 + 1)))
    y = np.array([175.5, 188.1, 200.2, 213.7, 231.6, 254.7, 283.1, 313.5, 342.9, 363.7, 389.1, 422.2,
        481.1, 518.8, 555.1, 591.6, 643.1, 693.0, 767.3, 846.1])
    beta_QR, beta_LM = (171.0, 0.08356), (171.5, 0.08333)

    return funcao_objetivo(y,x,beta_QR), funcao_objetivo(y,x,beta_LM)


