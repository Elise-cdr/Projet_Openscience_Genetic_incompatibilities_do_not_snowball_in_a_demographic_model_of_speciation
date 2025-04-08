import numpy as np
import matplotlib.pyplot as plt
from Orr_sans_forcer import orr
from matplotlib.lines import Line2D


def graphe(Nb_sim,len_genome,nb_gen,sequence,inf):
    liste_moy=np.array([0 for i in range(nb_gen)])
    liste_moy_AD=np.array([0 for i in range(nb_gen)])
    liste_moy_DD=np.array([0 for i in range(nb_gen)])

    #liste qui à terme va contenir la moyenne des DMIs des simulations réalisées (indice correspond à une génération
    modèle='orr' if inf else 'orr-inf'

    Liste_tot=[[0]*Nb_sim for i in range(nb_gen)]

    for num_sim in range(Nb_sim):

        orr(len_genome,nb_gen,num_sim,sequence,inf)
        num = (3-len(str(num_sim)))*'0'+str(num_sim)
        nom_fichier_DMIs=f'Simulations/{num}_DMIs.csv'

        with open(nom_fichier_DMIs,'r') as f:
            lignes=f.readlines()[1:]


            for gen,DMIs in enumerate(lignes):#♥on ne veut pas prendre la première ligne

                DMIs = DMIs[:-1].split(',')
                liste_moy[gen]+=int(DMIs[1])
                liste_moy_AD[gen]+=int(DMIs[2])
                liste_moy_DD[gen]+=int(DMIs[3])

                Liste_tot[gen][num_sim]=int(DMIs[1])
        print(num_sim)


    with open(f'Simulations/{modèle}_timeComparison.csv','w') as file:
        file.write('generations,'+'avgDMIsperHyb,'+'stdError\n')
        for gen in range(nb_gen):
            mean=np.mean(Liste_tot[gen])
            se = np.std(Liste_tot[gen], ddof=1) / np.sqrt(len(Liste_tot[gen]))
            file.write(f'{gen},{mean},{se}\n')


    liste_moy=liste_moy/Nb_sim
    liste_moy_AD=liste_moy_AD/Nb_sim
    liste_moy_DD=liste_moy_DD/Nb_sim
    Generations=np.array([i for i in range(nb_gen)])

    bar_width=0.25
    #plt.bar(Generations,liste_moy, width=bar_width,label='DMIs',color='green',alpha=0.2)
    plt.bar(Generations,liste_moy_AD,width=bar_width,label='DMIs_AD',color='blue',alpha=0.2)
    plt.bar(Generations,liste_moy_DD,width=bar_width,label='DMIs_DD',color='red',alpha=0.2)


    plt.figtext(0.12, 0.87, f'Taille du génome: {len_genome}', fontsize=10)
    plt.figtext(0.12, 0.82, f'Modèle: {modèle}', fontsize=10)
    plt.figtext(0.12, 0.77, f'Nombre de simulations: {Nb_sim}', fontsize=10)

    plt.title(f'Evolution du nombre de DMIs au cours des générations')
    plt.axis([0, nb_gen+1, 0,liste_moy[-1]+1 ])
    plt.xlabel('Génération')
    plt.ylabel('Nombre de DMIs')
    plt.tight_layout()

    plt.legend()
    plt.show()

graphe(60,24,100,None,False)



















