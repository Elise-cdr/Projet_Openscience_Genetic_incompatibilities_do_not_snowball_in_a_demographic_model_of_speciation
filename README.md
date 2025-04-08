# Projet_Openscience_Genetic_incompatibilities_do_not_snowball_in_a_demographic_model_of_speciation
## Modélisation de l'accumulation du nombre d'incompatibilités génétiques au cours des générations

### Présentation du script Orr
Le script Orr.py va simuler l'accumulation d'incompatibilités génétiques entre 2 populations évoluant indépendamment, il va permettre de produire 4 fichiers:
\-1 fichier json avec les paramètres entrés
\-1 fichier par population comportant le génome au cours des générations
\-1 fichier recensant le nombre de DMIs entre les populations à chaque génération, avec le détail des types de DMIs

### Comment l'exécuter
#### En ligne de commande
Pour l'exécuter, il faut avoir dans un répertoire principal les scripts Orr.py et dmiGenerator.py, ainsi qu'un dossier vide Simulations dans lequel seront stockées les résultats des manipulations.
Taper dans le répertoire la ligne de commande suivante  :

python3 Orr.py Elle est suivie des paramètres d'entrés: la longueur du génome, le nombre de génération, le numéro de la simulation (va servir pour nommer les fichiers), de manière facultative la séquence génomique qui a permis de construire l'ensemble de DMIs (sinon elle est créée par le script dmiGenerator par défaut), et également de manière facultative le fait de faire ou non l'hypothèse des sites infinis (l'hypothèse est faite par défaut)

Par exemple on entrera : python3 Orr.py 6 4 1 --inf False --sequence AbcDEf

pour faire une simulation sur un génome à 6 gènes, sur 4 générations, avec 1 comme numéro de simulation.

#### En appelant la fonction Orr

On peut aussi appeler la fonction Orr du script Orr pour la réutiliser dans un autre fichier en tapant dans le même répertoire: from Orr import Orr On rentre ensuite les mêmes arguments qu'en ligne de commande

### Sortie du script
On récupèrera:
- 001_Paramètres.json avec la longueur du génome, le nombre de génération, le numéro de la simulation, l'ensemble de DMIs sous forme de dictionnaire
- 001_pop0.csv avec 2 colonnes: 1 pour le numéro de la génération, 1 avec le génome à cette génération
- 001_pop1.csv même chose que le fichier précédent mais pour l'autre population
- 001_DMIs.csv avec 4 colonnes: le numéro de la génération, le nombre de DMIs à cette génération, le nombre de DMIs ancestrales-dérivées, puis dérivées-dérivées

### Visualisation graphique et fit:
On utilise le code Graphique_Orr.py pour tracer le nombre de DMIs de type ancestrales-dérivées, et celles dérivées-dérivées au cours des générations. La fonction graphe prend en paramètres le nombre de simulations du script précédent (utilisées pour faire une moyenne du nombre de DMIs), la taille du génome désirée, le nombre de générations, la séquence finale (permettant la synthèse de la liste de DMIs, sinon celle-ci est générée par défaut par dmiGenerator en rentrant None), le modèle désiré (True pour Orr, False pour Orr-Inf).

Cette fonction permet également de créer un fichier orr_timeComparison.csv avec en première colonne le numéro de la génération, en seconde la moyenne des DMIs sur les différentes simulations, en dernière colonne l'erreur standard. Ce fichier est ensuite utilisé pour trouver le meilleur fit avec le code de l'article AICw.R
