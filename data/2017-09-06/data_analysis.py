import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, fsolve
from scipy.stats import chi2
chi2_scipy = chi2
import webbrowser

# analyse functies
# data is van de vorm [x-waarde, gemeten y-waarde, foutop y-waarde]
def Chi2(data, model,fout_model=None):
    # Chi kwadraat test is op Sum ( gemeten - verwacht)^2 / var
    # data is array als volgt
    # data = np.array([1/golflengtes^2, n van die golflengtes, fout (standaardeviatie) op n])//
    return lambda theta : afstand(data, model, theta, fout_model)

def afstand(data, model, theta, fout_model=None):
    x = data[0]
    n_model = (model(theta))(x)
    n_exp = data[1]

    if fout_model is None:
        var = data[2]**2
    elif len(data)>= 4:
        sdx = fout_model(theta)(np.array([data[0], data[3]]))
        #print(theta, sdx)
        var = data[2]**2 + sdx**2
    else:
        raise KeyboardInterrupt('data is niet van de juiste vorm. Moet zijn als np.array([X,Y,dY,dX])')

    return (np.sum((n_exp-n_model)**2/ var))**0.5

def optimalisatie(data, model, beginwaarden=np.array([1, 1]), fout_model=None, plot=False):
    # optimaliseren
    fun = Chi2(data, model, fout_model)
    mini = minimize(fun, beginwaarden)
    if plot:
        plot_fit(data, model, mini['x'])
    return mini

def kwaliteit_fit(data, mini):
    # theta heeft aantal parameters p
    # chi^2 reduced is delen door n - p

    n = np.size(data[0])
    p = np.size(mini["x"])

    chi = mini["fun"]
    chi_red = chi / (n-p)


    alfa = chi2_scipy.sf(chi_red, n-p)
    print("de p-waarde is ",alfa,"\nChi^2 reduced=", chi_red,"\naantal vrijheidgraden=",n-p,"\nchi^2=",chi)
    print("we behouden de fit tot op een betrouwbaarheidsniveau van ", alfa,"\nOftewel verwerp als p-waarde <\alpha-niveau ")
    chi2_plot(n-p, end=3*n, t=chi_red)
    return alfa


def chi2_plot(df,end,t=0):
    X=np.linspace(0,end,500)
    Y=chi2_scipy.pdf(X,df)
    h = max(Y)
    fig, ax = plt.subplots(1,1)
    ax.plot(X,Y)
    ax.plot([t,t],[0,h])
    plt.show()

def plot_fit(data, model, theta0, titel="Data met fit ",labelx=" $x[]$",
             labely=" $y$  []", figname=None , error=True, ylim=None):
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=120, figsize=(8, 4))

    # We genereren een linspace om het model te plotten
    x = np.linspace(np.min(data[0]), np.max(data[0]), 1000)

    if  error:
        if len(data)<=3:
            # We schalen ook de stroom voor een beter leesbare x-as
            ax.errorbar(data[0], data[1], yerr=data[2], label="data",

                    # De errorbars wat mooier maken :)
                    marker="o", markersize=4, fmt=" ", color="black", ecolor="black", capsize=2, capthick=0.6, linewidth=0.6)
        else:
            # We schalen ook de stroom voor een beter leesbare x-as
            ax.errorbar(data[0], data[1], yerr=data[2],xerr=data[3] , label="data",

                    # De errorbars wat mooier maken :)
                    marker="o", markersize=4, fmt=" ", color="black", ecolor="black", capsize=2, capthick=0.6, linewidth=0.6)


    elif not error:
        # gewoon punten plotten als de errorbars het model zouden overstemmen bijvoorbeeld
        ax.plot(data[0], data[1], label='data', marker='.', markersize=4, color="black")

    ax.plot(x, model(theta0)(x), 'r', label="fit")

    ax.set_title(titel)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.legend()
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    try:
        if figname is not None:
            fig.savefig(figname+'.png')
            print("De figuur is succesvol opgeslagen.")
    except:
        print("De grafiek kon niet worden opgeslagen. Fout waarschijnlijk in de naam.")

    plt.tight_layout() ; plt.show()

def plot_data(data,titel="titel ",labelx="x", labely="y", figname=None, error=True, grijs_gebied=False, vlines=None,
             x_domain=None, x_as=False):

    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=120, figsize=(8, 4))

    # We genereren een linspace om het model te plotten
    x = np.linspace(np.min(data[0]), np.max(data[0]), 1000)

    if  error:
        if len(data)<=3:
            # We schalen ook de stroom voor een beter leesbare x-as
            ax.errorbar(data[0], data[1], yerr=data[2], label="data",

                    # De errorbars wat mooier maken :)
                    marker="o", markersize=4, fmt=" ", color="black", ecolor="black", capsize=0.5, capthick=0.5, linewidth=0.5)
        else:
            # We schalen ook de stroom voor een beter leesbare x-as
            ax.errorbar(data[0], data[1], yerr=data[2],xerr=data[3] , label="data",

                    # De errorbars wat mooier maken :)
                    marker="o", markersize=4, fmt=" ", color="black", ecolor="black", capsize=1, capthick=0.5, linewidth=0.5)


    elif not error:

        if grijs_gebied:
            ax.plot(data[0], data[1], label='data', marker='.', linestyle=' ', markersize=1, color="black")

            y1 = data[1]-data[2]
            y2 = data[1]+data[2]
            ax.fill_between(data[0], y1, y2, color="grey", label="error", alpha=0.5)

        else:
            ax.plot(data[0], data[1], label='data', marker='.', linestyle=' ', markersize=4, color="black")


    if vlines is not None:
        ax.vlines(vlines, ymin=np.min(data[1]), ymax=np.max(data[1]))
    ax.set_title(titel)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    if x_domain is not None:
        ax.set_xlim(x_domain[0], x_domain[1])
    if x_as:
        ax.axhline(y=0, xmin=0, xmax=90, color='black', linewidth=0.5)

    ax.legend()
    try:
        if figname is not None:
            fig.savefig(figname+'.png')
            print("De figuur is succesvol opgeslagen.")
    except:
        print("De grafiek kon niet worden opgeslagen. Fout waarschijnlijk in de naam.")

    plt.tight_layout() ; plt.show()



# deze functie is de partiele chi2 functie die slechts in het i-de element van theta veranderd.
def part_chi2(i, data, model, theta0, fout_model):
    return lambda theta_i: afstand(data, model, subs(theta0, i, theta_i), fout_model)

# een substitutie functie die x op de i-de plaats in  vector v zet (overschrijft).
# hierbij wordt met een nieuwe V gewerkt zodat de oude v niet veranderd
def subs(v, i, x):
    V = np.copy(v)
    V[i]=x
    return V

def print_fout(waarde, fout, end='\n'):
    print_onzekerheid([waarde-fout, waarde, waarde+fout], end)

def print_onzekerheid(oz = ['ondergrens', 'verwachte waarde', 'bovengrens'], end="\n"):
    pas_op=[False, False, False] # om op te passen dat we geen 0'en vergeten bij X, T1 of T2

    if oz[1]<0:
        Negatief=True
        oz = [-i for i in oz[::-1]]
    else:
        Negatief = False

    T1 = oz[1]-oz[0]
    T2 = oz[2]-oz[1]

    if not (oz[0]<= oz[1]<=oz[2]):
        T1 = max(abs(T1), abs(T2))
        T2 = T1

    T = min(T1, T2)

    if T> 0:
         e = int(-np.floor(np.log10(T)) )
    else:
        e = 0

    T1 = T1 * 10**e
    T2 = T2 *10**e

    if T1 < 2:
        T1 = round(T1, 1)
        pas_op[1]=True
    else:
        T1 = int(round(T1))

    if T2 < 2:
        T2 = round(T2, 1)
        pas_op[2]=True
    else:
        T2 = int(round(T2))
    T = min(T1, T2)

    if T<2:
        X = round(oz[1], e+1)*10**e
    else:
        X = int(round(oz[1], e)*10**e)
    if str(X)[-1]=='0':
        pas_op[0]=True

    # schalen
    if X == 0:
        orde = 0
    else:
        orde= int(-np.floor(np.log10(X)) )
    if orde != 0:
        X=X*10**orde
        T=T*10**orde
        T1=T1*10**orde
        T2=T2*10**orde
        e=e+orde

    # foutjes in afronden opsporen:
    if pas_op[1] or pas_op[2]:
        k=1
    else:
        k=0
    if X != round(X, -orde+k):
        #print("Fout:", X,"=/=", round(X,-orde+k), e, orde, k)
        X=round(X, -orde+k)
    if T1 != round(T1, -orde+k):
        #print("Fout:", T1,"=/=", round(T1,-orde+k), e, orde, k)
        T1=round(T1, -orde+k)
    if T2 != round(T2, -orde+k):
        #print("Fout:", T2,"=/=", round(T2,-orde+k), e, orde, k)
        T2=round(T2, -orde+k)

    if Negatief:
        X = -X
        T1_nieuw = T2
        T2 = T1
        T1 = T1_nieuw

    # werkelijk weergeven waarachtig wünderbar
    if T1== T2:
        print('$', end="")
        if e!=0:
            print('\\left(', end='')
        print(str(X) + '\\pm'+str(T1), end='')
        if e!=0:
            print('\\right)e'+str(-e), end='')

        print('$', end)
    else:
        print('$', end="")
        if e!=0:
            print('\\left(', end='')
        print(str(X) + '_{-'+str(T1)+'}^{+'+str(T2)+'}', end='')
        if e!=0:
            print('\\right)e'+str(-e), end='')

        print('$', end)

def print_onzekerheid_gem_cl(gem, cl):
    #cl van de vorm [[np.array([]), np.array([])]]
    oz =[]

    oz.append(cl[0][0][0])
    oz.append(cl[0][1][0])
    oz.insert(1, gem)
    print_onzekerheid(oz)

def onzekerheidsplots(mini, namen_variabelen, data, model, fout_model=None, chi_in_punt=None):
    theta0 = mini['x']
    P = len(theta0)
    N = len(data[0])

    fig, ax = plt.subplots(nrows=P, ncols=1, dpi=100, figsize=(8, 2.5*P), sharex=False)

    # bepalen van de x-assen want $\theta = (x_0, \gamma, A, y_0)$
    X_as = namen_variabelen
    list_CL = []

    ####################################################
    # HET AANTAL VRIJHEIDSGRADEN???      N-P ??        #
    ####################################################

    # bepaal 1 sigma hoogte:
    t = chi2_scipy.ppf(0.68, df=N-P)
    h = mini["fun"] + t

    # genereer een 1 sigma hoogte lijn
    H = h*np.ones(500)

    # we overlopen alle variablen in theta
    for i in range(P):

        # bepaal de variable en zijn partieelfunctie
        a = theta0[i]
        f = part_chi2(i, data, model, theta0, fout_model)

        # bepaal het 1-sigma betrouwbaarheidsinterval
        CL = [fsolve(lambda x: f(x) - h, a*0.5), fsolve(lambda x:f(x) - h, a*1.5)]
        print_onzekerheid_gem_cl(a, [CL])
        CL.sort()
        list_CL.append(CL)

        # We genereren een linspace om het model te plotten
        D = CL[1]-CL[0]
        p=0.75
        X = np.linspace(theta0[i]-p*D, theta0[i]+p*D, 500)

        # genereer de chi2 waarden voor X, voor de partieelfuntie van chi2 voor theta[i]
        Y = []
        for j in X:
            Y.append(f(j)) # f = part_chi2(i, data, model, theta0, fout_model)
        Y = np.array(Y)

        MaxY = np.max(Y)
        if not MaxY is float:
            MaxY = np.quantile(Y, q=0.75)

        if P>1:
            # plot de gevonden waarden
            ax[i].plot(X, Y, label="$\chi_"+str(i+1)+"^2($"+X_as[i]+")")
            ax[i].plot(X, H, '--',  label="$\chi_{min}^2$+chi2.ppf(0.68, df=4)")
            ax[i].set_xlim(theta0[i]-p*D, theta0[i]+p*D)
            ax[i].set_ylim(np.min(Y), MaxY)
            ax[i].legend(fontsize='small')
            ax[i].set_xlabel(X_as[i])
            ax[i].set_title("de partiële functie $\chi^2_"+str(i+1)+"$ met als variabele "+X_as[i], size=10)
        elif P==1:
            # plot de gevonden waarden
            ax.plot(X, Y, label="$\chi_"+str(i+1)+"^2($"+X_as[i]+")")
            ax.plot(X, H, '--',  label="$\chi_{min}^2$+chi2.ppf(0.68, df=4)")
            ax.set_xlim(theta0[i]-p*D, theta0[i]+p*D)
            ax.set_ylim(np.min(Y), MaxY)
            ax.legend(fontsize='small')
            ax.set_xlabel(X_as[i])
            ax.set_title("de partiële functie $\chi^2_"+str(i+1)+"$ met als variabele "+X_as[i], size=10)

    fig.suptitle('De partiële functies van $\chi^2$')
    fig.savefig('onzekerheidsplot.png')

    plt.tight_layout() ; plt.show()
    if not chi_in_punt is None:
        functiewaarde = f(chi_in_punt)
        distance = functiewaarde - mini["fun"]
        kans = chi2_scipy.sf(distance, df=N-P)
        print('De functiewaarde in het punt', chi_in_punt, 'is', functiewaarde, '\nDe afstand tot het minimum is', distance,
              '\nDe kans om punten verder dan dit punt te berijken met de chikwadraat verdeling met vrijheidsgraden N-P=',N, '-',
              P,'=', N-p,'is', kans)
    return list_CL

def local_minima(data, threshold=0.1):
    minima=[]
    laatste_maximum=0
    laatste_minima= []
    huidig_minimum=0
    nieuw_maximum=0
    p = punt(data)
    for i in range(len(data[0])):
       # print("punt", p(i), "minimum", p(huidig_minimum), "laatste max:", p(laatste_maximum)
        #      ,"nieuw max:", p(nieuw_maximum), "beslissing: ", end="")

        if data[1][i]<data[1][nieuw_maximum]-threshold and not (laatste_maximum == nieuw_maximum):
            # we hebben een dal overgestoken en voegen dat toe aan de minima en beginnen aan het nieuwe dal
            minima.append(laatste_minima)
            laatste_maximum=nieuw_maximum
            laatste_minima=[i]
            huidig_minimum=i
           # print("Verwerk dal")

        elif data[1][i]<data[1][huidig_minimum]:
            # bij het dalen
            huidig_minimum=i
            nieuw_maximum=i
            laatste_minima=[i]
         #   print("nieuw minimum")

        elif data[1][i]== data[1][huidig_minimum]:
            laatste_minima.append(i)
          #  print("extra minimum")

        elif data[1][i]>data[1][huidig_minimum]+threshold and(data[1][i]>data[1][nieuw_maximum]):
            # hier beginnen we terug te stijgen
            nieuw_maximum=i
        #    print("nieuw maxima")

        else:
            #print("geen actie")
            pass

    minima.append(laatste_minima)
    m = []
    for i in minima:
        m.append(data[0][i])
    return m

def eerste_minimum(data, threshold=0.1):
    # retunt eerst punt dat boven threshold uitkomt
    indices=np.where(data[1]<0.1)
    i=0
    while i in list(indices[0]):
        i+=1
    return data[0][i]

def punt(data):
    return lambda index:  (data[0][index], data[1][index])

def gemiddelde(lijst, zonder=[]):
    #zonder zijn de indexen die niet worden meegeteld.
    lijst=list(lijst)
    for i in zonder[::-1]:
        lijst.pop(i)
    som=0
    for i in lijst:
        som+=i
    som=som/(len(lijst))
    return som

def std(lijst):
    s = np.std(lijst)
    if s ==0:
        s=1
    return s

def gemiddelden(lijst):
    gem=[]
    for l in lijst:
        gem.append(gemiddelde(l))
    return(gem)

def gemiddelden_std(lijst):
    l2 = []
    for l in lijst:
        l2.append(std(l))
    return l2

def verschillen(lijst):
    l1=lijst.copy()
    l1.pop(-1)
    l2 = lijst.copy()
    l2.pop(0)
    return np.array(l2)-np.array(l1)

def verschillend_std(lijst, sd):
    l=[]
    for i in range(len(lijst)-1):
        l.append( (sd[i]**2 + sd[i+1]**2) **0.5 )
    return l

def gewogen_gem(waarden, fouten):
    mg = np.sum( np.array(waarden) / np.array(fouten)**2 )
    fout = 1 /np.sum( 1 / np.array(fouten)**2)**.5
    mg = mg * fout**2
    return mg, fout
