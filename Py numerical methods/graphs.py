# %%
from sympy import exp, zeros, eye, Matrix
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# LM_miu(x, y, ig, miu):
# expetativaVSrealidade(x,y,b0,b1,nome, largo, esticado, label_linha, cor_graph, cor_ponto):

# %%
def LM_method_alterado(X, Y, ig, miu):
    """LM_method mas retorna b0,b1 a cada iteracao"""
    b0_arr, b1_arr = np.array([]), np.array([])
    matrix_r, grad_r = zeros(len(X), 1), zeros(len(X), 2)
    b0, b1, s_LM = ig[0], ig[1], Matrix([[0], [0]])

    while True:
        b0, b1 = b0 + s_LM[0], b1 + s_LM[1]
        b0_arr, b1_arr = np.append(b0_arr, b0), np.append(b1_arr, b1)

        for row in range(len(X)):
            matrix_r[row, 0] = Y[row] - b0 * exp(b1 * X[row])
            grad_r[row, 0], grad_r[row, 1] = -exp(b1 * X[row]), -b0 * X[row] * exp(b1 * X[row])

        A = grad_r.transpose() * grad_r + eye(2) * miu
        b = -grad_r.transpose() * matrix_r
        s_LM = A.cholesky_solve(b)

        if max(abs(s_LM)) < 0.5 * 10**(-3):
            break

    return b0_arr, b1_arr

# %%
def plot_residuals(x, y, ig, miu, ax, PONTOS=100):
    b0, b1 = LM_method_alterado(x, y, ig, miu=miu)

    b0_values = np.linspace(float(min(b0)) * 0.995, float(max(b0)) * 1.005, PONTOS)
    b1_values = np.linspace(float(min(b1)) * 0.995, float(max(b1)) * 1.005, PONTOS)
    X, Y = np.meshgrid(b0_values, b1_values)
    Z = np.sum(np.square(y - X.reshape(-1, 1) * np.exp(Y.reshape(-1, 1) * x)), axis=1).reshape(PONTOS, PONTOS)

    x, y, z = X.flatten(), Y.flatten(), Z.flatten()

    tricontour = ax.tricontourf(x, y, z, levels=np.logspace(np.log10(np.min(z)), np.log10(np.max(z)), 10), cmap="GnBu", norm=LogNorm())
    ax.scatter(b0, b1, color='darkgreen', s=10, alpha=0.5) #marcar pontos
    ax.text(b0[0], b1[0], '1', color='black', fontsize=8, ha='right', va='bottom', fontweight='bold') #ponto 1
    ax.text(b0[-1], b1[-1], str(len(b0)), color='black', fontsize=8, ha='left', va='top', fontweight='bold') #ponto fim
    ax.set_xlim(float(min(x)), float(max(x)))
    ax.set_ylim(float(min(y)), float(max(y)))

    return tricontour, z  # Return tricontour and z for later use

# %%
def LM_miu(x, y, ig, miu, width, height):
    # len(miu) graficos/linha & mudar tamanho (9,4) da imagem
    fig, axs = plt.subplots(1, len(miu), figsize=(width, height), constrained_layout=True)
    fig.supxlabel(r'Valores para $\beta_0$')
    fig.supylabel(r'Valores para $\beta_1$')

    for i in range(len(miu)):
        tricontour, z = plot_residuals(x, y, ig, miu=miu[i], ax=axs[i])
        axs[i].set_title(r'$\mu = {:g}$'.format(miu[i]))

    # Create a single colorbar for the entire figure
    cbar = plt.colorbar(tricontour, ax=axs, label='Resíduos', format='%.3g')
    cbar.set_ticks(np.logspace(np.log10(min(min(z), min(z))), np.log10(max(max(z), max(z))), 10))

    plt.savefig('LM_miu.png')
    #plt.show()

# %%
def grafico_QRvsLM(b_LM, b_QR, pontos):
    """
    b_LM = betas do LM (b0, b1)
    b_QR = betas do QR (b0, b1)
    pontos = (X,Y)
    """
    X,Y = pontos
    x = np.linspace(0.9*min(X), 1.1*max(X), num=50)
    LM =  b_LM[0]*np.exp(x*b_LM[1])
    QR =  b_QR[0]*np.exp(x*b_QR[1])
    
    plt.plot(x, LM, label='LM', color="red")
    plt.plot(x, LM, label='QR', color="blue")

    plt.scatter(X, Y, color="black")

    plt.savefig('LMvsQR.png')
    #plt.show()

#%%
def expetativaVSrealidade(x,y,b0,b1,nome, largo, esticado, label_linha, cor_graph, cor_ponto):

    plt.figure(figsize=(largo, esticado))

    x_e = np.linspace(0.98*min(x), 1.02*max(x), num=95)
    expetativa = b0*np.exp(x_e*b1)
    
    plt.plot(x_e, expetativa, label=label_linha, color=cor_graph)
    plt.scatter(x, y, color=cor_ponto, label="Valores observados")

    plt.xlabel("Anos (2000 a 2019)")
    plt.ylabel("Rendimento Nacional Bruto (RNB)")
    plt.legend()

    plt.savefig(nome)
    #plt.show()

# %%
# X = [-1, 0, 1, 2, 3]
# Y = [0.430, 1.78, 8.00, 38.0, 159]
# b_LM = (2.09076203689172, 1.44393202455175)
# b_QR = (1.8497359616687743, 1.4886721340451912)
# grafico_QRvsLM(b_LM, b_QR, (X,Y))

# x = list(range(19 + 1))
# y = [175.5, 188.1, 200.2, 213.7, 231.6, 254.7, 283.1, 313.5, 342.9, 363.7, 389.1, 422.2,
#     481.1, 518.8, 555.1, 591.6, 643.1, 693.0, 767.3, 846.1]
# ig = (170,0.0836)
#
# miu = [150, 10**(-2)]
# width, height = 9, 4 #image sizes
#
# LM_miu(x, y, ig, miu, width, height)

#x = list(range(19 + 1))
#y = [175.5, 188.1, 200.2, 213.7, 231.6, 254.7, 283.1, 313.5, 342.9, 363.7, 389.1, 422.2,
#    481.1, 518.8, 555.1, 591.6, 643.1, 693.0, 767.3, 846.1]
#b0,b1 = 170,0.0836
#expetativaVSrealidade(x,y,b0,b1, "QRvsrealidade.png",12,5, "Modelo ajustado (LM)", (50/255, 71/255, 143/255), "black")

# %%
