# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 18:33:16 2025

@author: admin
"""

#%% Remarques
"""
    - Peut-être que parfois le code va tourner sans arrêt sans trouver
    de locus qui respecte toutes les conditions pour qu'il y ait mutation

    - Il faudra réécrire la partie input avec des argparse et tout ce bazar

    - J'ai fait un dico pour les DMIs finalement, il faudrait convertir les
    DMIs générés par leur code dans cette forme, pour faire des tests avec des
    génomes plus grands

"""
#%% Répertoire de travail
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir) # Définit le répertoire de travail comme celui où se trouve le script

#%% Importations

import random as rd
import argparse
from dmiGenerator import dmiGenerator
import json

#%% Fonctions du programme

def DMIs_list_to_dico(list_DMIs):
    """
    Transforme la liste de DMIs (obtenue avec dmiGenerator.py) en dico
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZΑΒΓΔΕΖΗΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩϴАБВГДЕЖЗИКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯԱԲԳԴԵԶԷԸԹԺԻԼԽԾԿՀՁՂՃՄՅՆՇՈՉՊՋՌՍՎՏՐՑՒՓՔՕՖႠႡႢႣႤႥႦჁႧႨႩႪႫႬჂႭႮႯႰႱႲჃႳჇႴႵႶႷႸႹႺႻႼႽႾჄႿჀჅჍⲀⲂⲄⲆⲈⲊⲌⲎⲐⲒⲔⲖⲘⲚⲜⲞⲠϤⲢⲤⲦⲨⲪⲬⲮⲰⳀϢϦϨϪϬϮ"

    DMIs_simple = {}
    for dmi in list_DMIs:
        DMIs_simple[(alphabet.find(dmi[0].upper()), alphabet.find(dmi[1].upper()))] = str(int(dmi[0].isupper()))+str(int(dmi[1].isupper()))

    # Puis on créé un dictionnaire complet pour ne pas avoir à chercher dans les 2 sens
    # donc dans notre exemple, il contient : 'AB', 'BA', 'aC', 'Ca', 'BC', 'CB'
    DMIs = {}
    for pos in DMIs_simple.keys():
        DMIs[pos] = DMIs_simple[pos]
        DMIs[pos[::-1]] = DMIs_simple[pos][::-1]
    #print(DMIs)
    return DMIs


def initialize_sim(len_genome, num_sim, nb_gen, DMIs, inf):
    """
    Créé les génomes de départs (ne contenant que des 0)
    Créé le fichier dans lequel seront stockés les données de la simulation
    Renvoie ces deux génome et le nom du fichier
    """

    if 'Simulations' not in os.listdir(script_dir):
        os.mkdir('Simulations')

    pop0 = [0]*len_genome
    pop1 = [0]*len_genome

    num = (3-len(str(num_sim)))*'0'+str(num_sim)
    # Ce truc chelou c'est pour écrire 002 quand num_sim = 2

    for i in range(2):
        file_name = f'Simulations/{num}_pop{i}.csv'
        with open(file_name, 'w') as file:
            file.write('0,'+',0'*(len_genome)+'\n')
            file.close()

    #Création du fichiers de stockage des DMIs
    with open(f'Simulations/{num}_DMIs.csv', 'w') as file:
        file.write('Generation,Nb_DMIs,Nb_AD,Nb_DD\n')
        file.write('0,0,0,0\n')

    # Il faut aussi stocker les paramètres dans un fichier json:
    params={'len_genome':len_genome, 'nb_gen':nb_gen,'num_sim':num_sim, 'Sites infinis':['Non', 'Oui'][inf], 'DMIs':str(DMIs)}

    nom_fichier= num + '_Paramètres' + '.json'
    with open(f'Simulations/{nom_fichier}','w') as f:
        json.dump(params, f, indent=4)#l'indentation permet de rendre le fichier plus lisible, json.dump permet de convertir le dico en json

    return pop0, pop1

def mutate(pop0, pop1, to_mute, len_genome, DMIs, inf):
    """
    Effectue une mutation sur la population n° to_mute (to_mute = 0 ou 1)
    Cette mutation respecte les règles suivantes :
        - Pas de mutation dans un locus déjà muté
        - Le locus ne doit pas non plus avoir été muté chez l'autre population
        - Pas de mutation créant une DMI dans cette population
    pop0 et pop1 sont ensuite renvoyées (bien dans cet ordre)
    """


    # Détermination du rôle de chaque pop en fonction de to_mute


    pop_to_mute = [pop0, pop1][to_mute]
    other_pop = [pop1, pop0][to_mute]

    # On créé ensuite une boucle while pour trouver un locus tel que
    # les conditions permettant la mutation soient vérifiées
    is_mutation = 0

    if inf :
        loci_not_mutable = set() # répertorie les locus où mutation = impossible

    while is_mutation == 0:


        if inf :
            if len(loci_not_mutable) == len_genome:
                return(pop0,pop1)
            else:# choix au hasard du locus parmi les locus pas déjà testés
                locus = rd.choice(list(set(range(len_genome)) - loci_not_mutable))

        else :
            locus = rd.choice(list(set(range(len_genome))))

        is_mutation = 1 # on considère a priori qu'on va muter
        result_mutation = [1,0][pop_to_mute[locus]]
        # result_mutation est un 1 s'il y avait un 0 à ce locus, et inversement

        # Mais on vérifie les différente conditions :
        # Gestion de pop déjà mutée là / autre pop mutée là (pour sites infinis):
        if inf and pop_to_mute[locus] == 1:
            is_mutation = 0
            # c'est optionnel mais pour aller plus vite on pourrait mettre
            # un pass parce que c'est inutile de faire la suite
        elif inf and other_pop[locus] == 1:
                is_mutation = 0

        # Gestion des DMIs :
        dmi_pbs = 0
        for (i,j) in DMIs.keys():
            if i == locus:
                if DMIs[(i,j)][0] == str(result_mutation):
                    if DMIs[(i,j)][1] == str(pop_to_mute[j]):
                        dmi_pbs += 1

        if dmi_pbs != 0:
            is_mutation = 0

        # Si toutes les conditions sont vérifiées, on mute :
        if is_mutation == 1 :
            pop_to_mute[locus] = result_mutation

        # Sinon, (et si sites infinis) on ajoute le locus dans la liste des loci où il ne peut pas y avoir de mutation
        elif inf:
            loci_not_mutable.add(locus)

    # Ce return dépendant de to_mute permet de garder l'info de qui est
    # pop0 et qui est pop1 : on renvoie forcément dans l'ordre pop0, pop1
    return [(pop_to_mute, other_pop), (other_pop, pop_to_mute)][to_mute]

def count_dmis(pop0, pop1, len_genome, DMIs):
    """
    Renvoie le nombre de DMIs entre pop0 et pop1
    """
    nb_dmis = 0
    nb_AD = 0

    # On parcourt toutes les paires
    # Il n'est nécessaire de le faire que dans 1 sens comme on a doublé le dico
    for loc0 in range(len_genome):
        for loc1 in range(len_genome):
            if (loc0,loc1) in DMIs.keys():
                if f'{pop0[loc0]}{pop1[loc1]}' == DMIs[(loc0,loc1)]:
                    nb_dmis += 1
                    if '0' in DMIs[(loc0,loc1)]:
                        nb_AD += 1

    return nb_dmis, nb_AD

def stock_infos(n, num_sim, nb_dmis, nb_AD, pop0, pop1):
    """
    Pour une génération, stock le génome de chacune des populations et le nombre de DMIs
    dans le fichier file, créé au début de l'expérience
    """
    num = (3-len(str(num_sim)))*'0'+str(num_sim)

    pops = [pop0, pop1]
    if type(nb_dmis)==str:
        nb_DD=''
    else:
        nb_DD = nb_dmis - nb_AD

    for i, pop in enumerate(pops) :
        pop_string = ','.join(map(str, pop))
        file_name = f'{num}_pop{i}.csv'

        with open(f'Simulations/{file_name}', 'a') as file:
            file.write(f'{n},,{pop_string}\n')
            file.close()


    with open(f'Simulations/{num}_DMIs.csv', 'a') as file:
        file.write(f'{n},{nb_dmis},{nb_AD},{nb_DD}\n')

def orr(len_genome, nb_gen, num_sim, sequence = None, inf = True):

    list_DMIs=dmiGenerator(False,False,sequence,False,False,len_genome)

    # Transformation de la liste de DMIs obtenu par dmiGnerator.py
    # en dico de la forme qu'on veut (avec des O et des 1)
    DMIs = DMIs_list_to_dico(list_DMIs)

    # Création des pops de départ et du fichier de stockage des infos
    pop0, pop1 = initialize_sim(len_genome, num_sim, nb_gen, DMIs, inf)

    for n in range(1, nb_gen):
        # Choix au hasard de la population à muter
        to_mute = rd.randint(0, 1)

        # Mutation de cette population
        (pop0, pop1) = mutate(pop0, pop1, to_mute, len_genome, DMIs, inf)
        # Comptage du nombre de DMIs entre les deux populations
        nb_dmis, nb_mi = count_dmis(pop0, pop1, len_genome, DMIs)

        # Stock des génomes des deux pops et du nb de DMIS dans un csv
        stock_infos(n, num_sim, nb_dmis, nb_mi, pop0, pop1)



#%% Programme principal

if __name__ == "__main__":

    # Entrée des arguments de  l'utilisateur
    parser = argparse.ArgumentParser(description="Nombre de DMIs au cours des générations")

    # Declaration des arguments
    parser.add_argument('len_genome', type=int, help='Taille du génome')
    parser.add_argument('nb_gen', type=int, help='Nombre de générations')
    parser.add_argument('num_sim', type=int, help="Le numéro de la simulation")
    parser.add_argument('--inf', type=int, default = 1, help='Hypothèse des sites infinis')
    parser.add_argument('--sequence', type=str, default=None, help="Séquence finale")

    # Analyse des arguments passés lors de l’appel du programme
    args = parser.parse_args()
    len_genome=args.len_genome
    nb_gen=args.nb_gen
    num_sim=args.num_sim
    inf=args.inf
    sequence=args.sequence

    orr(len_genome, nb_gen, num_sim, sequence, bool(inf))

#%%











































#orr(10, 6, 7, inf = False)